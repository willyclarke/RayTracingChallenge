# rtc-metal milestone 12: path tracing

Branch `metal/12-path-tracing`, merged into develop with `--no-ff` when done.
rtc-rust (f64, CPU) stays the reference.

## Goal

Port rtc-rust's chapter 17 path-traced indirect lighting
(`World::set_path_tracing(samples, depth)`) to the GPU, in both render modes:

- **deterministic** (no `--samples`): `samples` paths per pixel, an image
  with the same expected value as rtc-rust's;
- **sampled** (`--samples N`, `--window`): one path per pixel per frame,
  converging together with anti-aliasing, soft shadows and focal blur.

The model is rtc-rust's, not a physically based path tracer: every hit still
takes direct light from the scene's point or area light; diffuse surfaces
add `albedo · mean(incoming color)` gathered along cosine-weighted
hemisphere rays; the material ambient term is switched off while path
tracing is on. There are no emissive materials.

Success:
- a probe shows the GPU's `cosine_direction` agrees with rtc-rust's for the
  same `(r1, r2)`;
- `test_chap_17_path_tracing_gpu_render_matches_cpu` (the Cornell box,
  16 samples, depth 2) stays under calibrated error thresholds on both
  traversals;
- a sampled render of that scene converges towards the CPU image;
- every existing chapter compare and sampled test passes unchanged;
- a JSON scene with `settings.path_tracing` renders from the CLI and in
  `--window`.

## Out of scope

- **Emissive materials / light only from emitters.** No CPU reference; a
  possible later milestone with a matching rtc-rust change.
- **Motion blur** (still rejected), **Russian roulette**, importance
  sampling of the light beyond what rtc-rust does.
- **rtc-rust changes.** None expected: `World::gi_samples`, `gi_depth` and
  `cosine_direction` are public, and the scene loader already maps
  `settings.path_tracing` to `set_path_tracing`.

## What rtc-rust computes

`color_at_gi(ray, remaining, gi_samples, gi_depth)` at a hit
(`rtc-rust/src/world.rs`, `shade_hit_gi`, `indirect_color_gi`):

```
surface   = lighting(..., ambient_scale = 0 if world.gi_samples > 0 else 1)
indirect  = black if gi_samples == 0 or gi_depth <= 0 or remaining <= 0
                    or material.diffuse ≈ 0
            else albedo · (1/gi_samples) Σ_i color_at_gi(bounce_i, remaining-1, 1, gi_depth-1)
albedo    = pattern color at object_point (or material.color) · diffuse
bounce_i  = Ray(over_point, cosine_direction(normalv, gi_uniforms(over_point, i, gi_depth)))
reflected = color_at_gi(reflect ray, remaining-1, gi_samples, gi_depth) · reflective
refracted likewise; Schlick blends them as before; indirect is added unblended.
```

So `depth` counts diffuse bounces after the first hit, only the first
diffuse hit branches (`gi_samples` rays), a mirror or glass bounce keeps the
path-tracing state, and every bounce also uses up one of `remaining`.

`gi_uniforms` hashes the f64 bits of `over_point`; f32 points can't
reproduce it, so unlike milestone 11's lens samples **the GPU image cannot
match the CPU pixel for pixel**. Both are unbiased estimates of the same
value; tests compare them statistically.

## Decisions

- **Whole paths per pixel, never branching.** The GPU's `color_at` traces
  one path: at a diffuse hit with `gi_depth > 0` it pushes **one** bounce ray
  with weight `weight · albedo`. Averaging `gi_samples` such paths has the
  CPU's expected value: the first hit's direct light is the same in every
  path, so averaging it changes nothing, and the mean of single-sample
  bounces is the CPU's mean of `gi_samples` bounces. A mirror seen first
  then a diffuse hit also matches, because the CPU keeps branching
  `gi_samples` there too.
- **A deterministic path-traced render is `gi_samples` frames,** one path
  per pixel each, every frame its own command buffer; `Gpu::render` adds
  them and divides. (Amended during implementation: the first version
  looped all paths inside one dispatch, as for lens samples, and the
  chapter 11 glass scene's long kernel tripped macOS's GPU watchdog,
  "Impacting Interactivity", killing other tests' command buffers too.)
  Without path tracing it stays one frame, so every existing deterministic
  image is bit-identical. `color_at` keeps **one call site** (milestone
  11's open debt).
- **Sampled mode traces one path per frame**; frames do the averaging and
  `gi_samples` only switches the feature on.
- **Hemisphere numbers come from a new `random01` range,**
  `GI_DIM = 0x4000_0000`: the `b`-th diffuse bounce of a path draws
  `random01(pixel, path, GI_DIM + 2b)` and `GI_DIM + 2b + 1`, where
  `path = frame · n + i` for lens ray `i` of the frame's `n` (so each lens
  ray of a focal-blur frame bounces its own way). That range is
  far above the light sampler's `2 + 6·call + slot` and below `LENS_DIM`, so
  every existing sampled sequence is untouched.
- **`cosine_direction` is ported exactly** (same basis choice, same
  formula) so a probe can compare it with the CPU's.
- **The stack grows to `2·MAX_DEPTH + 1`.** A hit pushes at most three
  children (refracted, reflected, diffuse bounce), each with
  `remaining - 1`, so at most two siblings wait per level.

## Host side

`GpuSceneInfo` (`src/types.rs`) gains the two settings:

```rust
pub struct GpuSceneInfo {
    pub camera: GpuCamera,
    pub light: GpuLight,
    pub counts: GpuSceneCounts,
    pub frame: u32,
    pub sampling: u32,
    /// rtc-rust's `World::gi_samples`: paths per deterministic pixel;
    /// 0 = path tracing off (ambient on).
    pub gi_samples: u32,
    /// rtc-rust's `World::gi_depth`: diffuse bounces after the first hit.
    pub gi_depth: u32,
    pub _pad: [u32; 2],
}
```

Size goes 208 → 224 (`GpuMat4` is `align(16)`); the Rust `const` assert
and the MSL `static_assert` both move to 224. `GpuScene::from_world`
(`src/scene.rs`) drops the "path tracing is not supported so far" error and
fills the fields from `world.gi_samples` / `world.gi_depth` (a negative
depth becomes 0).

## Shader side (`shaders/rt.metal`, `shaders/render.metal`)

- `SceneInfo` mirrors the new fields.
- `cosine_direction(float4 normal, float r1, float r2)`: rtc-rust's basis
  (`axis = |n.x| > 0.9 ? y : x`, `tangent = normalize(cross(n, axis))`,
  `bitangent = cross(n, tangent)`) and `r = sqrt(r2)`, `phi = 2π·r1`.
- `lighting` takes an `ambient_scale` (1 everywhere except `color_at`
  with path tracing on, where it is 0); `probe_lighting` passes 1.
- `PendingRay` gains `int gi_depth` (diffuse bounces left) and
  `uint bounce` (diffuse bounces taken, for the random dimension).
  The camera ray starts with `gi_depth = scene.gi_samples > 0 ?
  scene.gi_depth : 0` and `bounce = 0`.
- `color_at` takes the path index. After the existing reflect/refract
  pushes, when `p.remaining > 0`, `p.gi_depth > 0` and `m.diffuse != 0`, it
  pushes `Ray{over_point, cosine_direction(normalv, u1, u2)}` with weight
  `p.weight · pattern_at(...) · m.diffuse`, `remaining - 1`,
  `gi_depth - 1`, `bounce + 1`. Reflected and refracted rays copy
  `gi_depth` and `bounce`.
- `render` passes `scene.frame · n + i` as the path index; frame 0
  overwrites `pixels`, later frames add to it, in both modes.

## Scenes

`scenes/chapter17_path_tracing.json`: rtc-rust's `cornell_box` with the
point light at (0, 1.8, -0.6), camera from (0, 1, -3.5) to (0, 1, 0), field
of view 39°, `"settings": { "path_tracing": { "samples": 16, "depth": 2 } }`
— rtc-rust's `test_chap_17_path_tracing_putting_it_together`. Planes and
cubes only; nothing the GPU lacks. If the blocks' bottom faces, coplanar
with the floor, speckle, lift the blocks 0.01 in the JSON (CPU and GPU
both read it).

## Testing

- **Probe** `probe_cosine_direction`: a handful of normals (axis-aligned,
  `|x| > 0.9`, oblique) × `(r1, r2)` pairs including `r2 = 0`; GPU vs
  `rtc_rust::world::cosine_direction` within 1e-5 per component.
- **Unchanged images**: all existing `*_gpu_render_matches_cpu` and sampled
  tests pass without threshold changes (path tracing off ⇒ same code path).
- **`test_chap_17_path_tracing_gpu_render_matches_cpu`**, both traversals,
  at the scene's 150×150: deterministic GPU vs `render_parallel`. Not
  `check_matches_cpu` (no pixel-exact budget is possible); instead two
  calibrated thresholds, measured once and written down like milestone 11's:
  - mean absolute error (catches noise-level problems; both images carry
    independent 16-sample noise, so it does not go to 0);
  - mean **signed** error per channel (catches bias — a missing bounce,
    ambient left on, a wrong weight — which mean noise would hide).
- **`test_chap_17_path_tracing_gpu_sampled_render_converges`**: like the
  focal-blur one: error after 256 frames ≤ half the error after 4 and
  under a calibrated threshold.
- **Bias sanity check while calibrating** (not a kept test): render the CPU
  at e.g. 256 samples once and confirm the GPU's signed error against it
  is near 0, so the thresholds aren't calibrated around a bug.
- **CLI / window**: render `scenes/chapter17_path_tracing.json` headless,
  with `--samples 64`, and in `--window`; check by eye.

## Risks

- **The METAL_RT miscompile** from milestone 11 came back with more
  `color_at` call sites; this design keeps one, but `color_at` grows. If the
  Metal traversal returns garbage, bisect the shader change before anything
  else and record the finding in the memory note.
- **Thread memory**: the stack doubles to 21 `PendingRay`s; if occupancy
  drops noticeably on existing scenes, size it from `gi_depth` at the
  host (function constant) instead — only if measured.
- **f32 self-intersection** on bounce rays from `over_point` (EPSILON 1e-5)
  near the Cornell box's wall/floor seams: would show as dark speckles in
  the GPU image and a signed-error bias; debug with a diff map.
