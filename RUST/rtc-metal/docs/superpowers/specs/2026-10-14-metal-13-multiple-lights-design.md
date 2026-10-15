# rtc-metal milestone 13: multiple lights

Branch `metal/13-multiple-lights`, merged into develop with `--no-ff` when
done. rtc-rust gets the feature first and stays the reference; rtc-metal
then matches it.

## Goal

Scenes with several lights, e.g. the Cornell box with a second light in
front of its opening. Both renderers sum each light's contribution;
the GPU matches the CPU pixel for pixel in deterministic renders and
statistically when sampled or path traced.

Success:
- rtc-rust: `World` holds a list of lights; scene files accept a
  `"lights": [...]` array; every existing rtc-rust test and example passes
  and every single-light image is bit-identical.
- rtc-metal: up to 8 lights; `test_bonus_multiple_lights_gpu_render_matches_cpu`
  passes on both traversals through `check_matches_cpu`; a sampled
  multi-light render converges; a path-traced two-light Cornell box
  matches statistically; every existing test passes unchanged.
- `scenes/cornell_two_lights.json` renders from the CLI and in `--window`.

## Out of scope

- **Spot lights on the GPU** (still rejected by `gpu_light`).
- **Light importance sampling** (one light picked per sample); every
  shading call sums all lights.
- **A light limit in rtc-rust**; only the GPU caps the count.
- **Shadow-ray sharing or culling** between lights.

## Decisions

- **Ambient per light, summed** — the book's own suggestion: `shade_hit`
  sums the full `lighting_with_ambient` (ambient + diffuse + specular)
  over every light. Two equal lights double the ambient; path tracing
  (ambient scale 0) is unaffected.
- **Scene format: an optional `"lights": [...]` next to `"light"`.** A
  scene has exactly one of the two; `"lights"` must be non-empty. Every
  existing scene file keeps working unchanged.
- **Single-light images stay bit-identical**: the sum starts from black,
  and `0 + x == x` in floating point, so one light gives today's value.
- **GPU lights live in `SceneInfo`** as a fixed array of
  `MAX_LIGHTS = 8` plus `light_count`, not in a new buffer: no new
  `[[buffer(n)]]` binding ripples through every kernel and probe, and
  `SceneInfo` (784 bytes) stays far below `setBytes`' 4 KB.
- **Light 0 keeps today's random dimensions** (`2 + 6·call + slot`), so
  every existing sampled image is unchanged. Light `k ≥ 1` uses
  `LIGHTS_DIM + (k − 1)·0x0100_0000 + 6·call + slot` with
  `LIGHTS_DIM = 0x2000_0000`: below `GI_DIM` (0x4000_0000) and
  `LENS_DIM`, and 7 lights × 16M dims never reach `GI_DIM`.

## rtc-rust

`World` (`src/world.rs`):

```rust
pub struct World {
    pub shapes: Vec<Box<dyn Shape>>,
    /// Every light in the scene; `shade_hit` sums their contributions.
    pub lights: Vec<Light>,
    ...
}
impl World {
    /// Replace all lights with `light`.
    pub fn set_light(&mut self, light: Light);
    /// Add another light.
    pub fn add_light(&mut self, light: Light);
}
```

- `shade_hit_gi` loops over `self.lights`:
  `surface = Σ light.lighting_with_ambient(…, light.intensity_at_time(comps.over_point, self, comps.time), ambient_scale)`,
  then adds indirect, reflected and refracted exactly as now. No lights →
  black (today's `None` branch).
- The ~50 uses of `world.light` migrate mechanically:
  `w.light = Some(x)` → `w.set_light(x)`;
  `w.light.as_ref().unwrap()` → `&w.lights[0]`; `w.light.is_some()` /
  `is_none()` → `!w.lights.is_empty()` / `is_empty()`; the three
  examples likewise. `Display` lists every light.

Scene loader (`src/scene.rs`):

```rust
pub struct SceneDescription {
    ...
    #[serde(default)]
    pub light: Option<LightDescription>,
    #[serde(default)]
    pub lights: Vec<LightDescription>,
    ...
}
```

`build` errors with `"a scene needs exactly one of \"light\" and \"lights\""`
when both or neither are given, and `"\"lights\" is empty"` for `[]`.
`help.rs` (`LIGHTS`, the overview) and `README.md` document the array.

## rtc-metal

Host (`src/types.rs`, `src/scene.rs`):

```rust
pub const MAX_LIGHTS: usize = 8;
pub struct GpuSceneInfo {
    pub camera: GpuCamera,
    pub lights: [GpuLight; MAX_LIGHTS],
    pub counts: GpuSceneCounts,
    pub frame: u32,
    pub sampling: u32,
    pub gi_samples: u32,
    pub gi_depth: u32,
    pub light_count: u32,
    pub _pad: u32,
}
```

Size: `96 + 8·80 + 24 + 6·4 = 784`. `GpuScene::from_world` errors with
`"scene has no light"` for none and `"more than 8 lights"` above the cap;
unused slots are zeroed. Tests reading `info.light` read `info.lights[0]`.

Shader (`shaders/rt.metal`):

- `SceneInfo` mirrors it (`Light lights[MAX_LIGHTS]; … uint light_count;
  uint pad;`, `static_assert` 784).
- `Sampler` gains `uint light` (aggregate initialisers that omit it get
  0, so probes are unchanged); `sample()` picks the dimension as in
  Decisions.
- `color_at` loops `for (uint k = 0; k < scene.light_count; k++)`,
  setting `sampler.light = k`, and adds
  `p.weight * lighting(m, scene.lights[k], …, intensity_at(scene.lights[k], …))`.
  `sampler.call++` stays once per shaded hit, after the loop.
- Probes keep taking one `Light`.

## Scenes

- `scenes/multiple_lights.json`: the chapter 7 scene (three spheres,
  floor and walls) lit by a red light (−10, 10, −10) and a blue light
  (10, 10, −10), each intensity 0.6 in its colour: coloured, overlapping
  shadows.
- `scenes/multiple_lights_area.json`: the soft-shadows scene's area light
  plus a point light from the other side.
- `scenes/cornell_two_lights.json`: `chapter17_path_tracing.json` (box
  walls, path tracing 16/2) with a second point light in front of the
  opening, (0.6, 1.5, −2.5), intensity 0.5.

## Testing

rtc-rust (`test_bonus_multiple_lights_*`):
- `shade_hit` with two lights equals the sum of the two one-light worlds'
  `shade_hit` (ambient included).
- A point shadowed from one light is still lit by the other.
- Loader: a `"lights"` array builds both lights; both keys, neither key,
  and `"lights": []` are errors with the messages above.
- Every existing rtc-rust test passes (they cover single-light
  bit-identity).

rtc-metal:
- `GpuSceneInfo` holds `light_count` and the lights in order; 9 lights is
  an error (unit test in `scene.rs`).
- `test_bonus_multiple_lights_gpu_render_matches_cpu`:
  `multiple_lights.json` through `check_matches_cpu` (both traversals).
- `test_bonus_multiple_lights_area_gpu_render_matches_cpu` (deterministic,
  `check_matches_cpu`) and `…_gpu_sampled_render_converges`
  (calibrated, like the soft-shadows test).
- `test_bonus_multiple_lights_path_tracing_gpu_render_matches_cpu`:
  `cornell_two_lights.json` through `check_converges_to_cpu`, calibrated
  1.5×, with a bias check against a high-sample CPU render while
  calibrating.
- All existing tests unchanged (light 0's dims and single-light sums are
  identical).
- CLI and `--window` on `cornell_two_lights.json`, checked by eye.

## Risks

- **Cost scales with lights**: every shaded hit traces one shadow ray per
  light (more for area lights). Path-traced multi-light scenes may push
  long frames towards the GPU watchdog (milestone 12); the test scenes are
  small, and deterministic path tracing already splits paths into frames.
- **The METAL_RT miscompile** (milestone 11) may return with the light
  loop inlining `intensity_at` more often; `render.metal` keeps one
  `color_at` call site, and every image test runs both traversals.
- **`SceneInfo` grows 224 → 784 bytes** per `setBytes`; negligible, but
  the Accumulator copies it every frame.
