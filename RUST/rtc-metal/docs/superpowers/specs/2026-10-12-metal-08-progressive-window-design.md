# rtc-metal milestone 8: progressive accumulation and a live window (bonus: soft shadows)

Branch `metal/08-progressive-window`, merged into develop with `--no-ff` when
done. rtc-rust (f64, CPU) stays the reference.

## Goal

Render on the GPU progressively: every frame adds one random sample per pixel
(jittered sub-pixel ray, jittered area-light sample) to an accumulation
buffer, and the running mean converges to the integral over the pixel area
and the continuous area light. rtc-rust's `usteps × vsteps` cell centres and
single pixel-centre ray approximate that integral, so the sampled image ends
close to the CPU image but not on it (soft_shadows, 10×10 light: mean
absolute error 4.6e-3 after 256 samples). Show it live in a
native macOS window with an orbit camera, or accumulate N samples headless and
write a PPM. Success: the deterministic path still matches every existing
compare test pixel for pixel; area lights match the CPU exactly in
deterministic mode; a sampled render converges to within a calibrated
mean-error threshold of the CPU image; the window opens, converges, pauses,
orbits and writes the PPM on close.

## Decisions

- **Sampled features:** antialiasing (uniform sub-pixel jitter) and area
  lights (the book's bonus chapter "Rendering soft shadows"). Sampling
  integrates over the continuous light and the pixel area; the CPU's cell
  centres are an approximation of that, not the target. Focal blur,
  motion blur, path tracing and spot lights stay out and are rejected with
  "not supported so far".
- **Window:** objc2 AppKit + MetalKit: an `MTKView` subclass and an
  `MTKViewDelegate`, defined with `define_class!`. No winit.
- **Behaviour:** accumulate to a sample target, then pause; mouse drag orbits
  the camera and restarts accumulation; scroll zooms; closing the window
  writes the PPM.
- **One sampling path:** headless `--samples N` and the window share the same
  `Accumulator`.

## Shader side (`shaders/rt.metal`, `shaders/render.metal`)

### Random numbers

A stateless PCG hash of `(pixel index, frame, dimension)` returns a float in
[0, 1). No RNG state lives on the GPU. Dimensions: 0–1 sub-pixel jitter, 2–3
light cell A and its jitter, 4–5 light cell B and its jitter, further
dimensions per bounce for reflected/refracted rays.

### Area lights

rtc-rust shades an area light as *(fraction of cells visible) × (mean of
diffuse + specular over the cell points)* — a product of two means. Picking
one cell per frame and using it for both factors would converge to the mean
of the product, a slightly different image. So each sampled frame draws two
**independent** cells: cell A for the shadow ray, cell B for diffuse and
specular. Independence makes the expected product equal the product of the
expectations, which is the CPU's formula. Both expectations are over the
continuous light (a uniform point anywhere on it, from the jittered cell), so
they equal the CPU's cell-centre means only in the limit of many cells.

`GpuLight` / MSL `Light` grows from 32 to 80 bytes: `corner`, `uvec`, `vvec`
(float4), `usteps`, `vsteps`, `samples` (uint) plus padding. A point light is
`usteps = vsteps = 1` with `corner` = its position and zero `uvec`/`vvec`, so
there is one lighting path. The `static_assert` sizes change with it.

### Deterministic vs sampled

`SceneInfo` gains `frame` (uint) and `sampling` (uint, 0 or 1).

- `sampling == 0`: the ray goes through the pixel centre and lighting loops
  over every cell centre (jitter 0.5) in one pass, including rtc-rust's
  adaptive corner rule for the shadow fraction (only when both steps ≥ 2).
  This is exactly rtc-rust
  with its default `jitter_by` of `[0.5]` and no antialias, so existing
  compare tests stay pixel-exact and area lights get an exact compare.
- `sampling == 1`: one random sample per pixel.

### Kernels

- `render` gains `device float4 *accum`. Deterministic: writes the pixel as
  today. Sampled: frame 0 writes the sample, later frames add to it (so
  `reset()` needs no clear).
- New `resolve`: `accum / (frame + 1)` into either a
  `texture2d<float, access::write>` (BGRA8Unorm drawable, clamped) or a
  float4 buffer.

### Known differences from the CPU (the convergence threshold absorbs them)

- The CPU supersamples only edge pixels on an n×n grid (`antialias.n`); the
  GPU jitters every pixel uniformly, so in-pixel texture detail averages
  slightly smoother. The GPU ignores `antialias.n`.
- The CPU's shadow test stops after the four corner cells when they agree,
  missing shadows that fall only on interior cells; the GPU's random cells do
  not.
- Light cells behind the surface's tangent plane (`dot(cell - point,
  normalv) < 0`) count as shadowed without a shadow ray. rtc-rust gets the
  same result geometrically: from its 1e-9 `over_point` the shadow ray
  re-enters the shape. From the GPU's 1e-5 `over_point` a ray dipping
  slightly below the tangent plane can miss the shape, so the GPU decides
  explicitly. The test uses the shading normal, so smooth triangles under an
  area light can differ from the CPU near terminators, where the shading and
  geometric normals disagree.

## CPU side

### `src/types.rs`, `src/scene.rs`

- `GpuLight` as above; `gpu_light()` converts area lights, still rejects spot
  lights.
- `GpuSceneInfo` gains `frame` and `sampling`.
- `GpuScene::from_world` rejects `focal_blur`, `motion_blur` and
  `path_tracing` with "not supported so far".
- New `GpuScene::set_camera(&Camera)` rebuilds only the `GpuCamera`; shape
  and BVH buffers are untouched.

### `src/gpu.rs`

- `render(&scene)` unchanged: one deterministic frame (probes, compares).
- New `Accumulator { accum, resolved, frame }`:
  - `Accumulator::new(&gpu, width, height)`
  - `add_frame(&gpu, &scene)`: encode `render` with `sampling = 1`, `frame`
  - `resolve_to_texture(&gpu, encoder, texture)` (window)
  - `resolve_to_buffer(&gpu) -> Vec<[f32; 4]>` (PPM)
  - `reset()`: `frame = 0`
- Errors stay `Result<_, String>`.

### `src/orbit.rs` (new)

CPU-side f64, reusing rtc-rust. Holds `to`, `up` and `from` as spherical
coordinates (radius, yaw, pitch) around `to`, initialised from the scene
camera. Drag changes yaw/pitch (pitch clamped to ±89° so `up` never lines up
with the view direction); dragging "grabs the scene" (right orbits the
camera left, down raises it); scroll scales the radius by exp(-0.1) per unit
(trackpad deltas × 0.1), clamped to [1e-3, 1e3 × the initial radius].
Zero-delta events do not restart accumulation. Produces
`view_transform(from, to, up)` for `set_camera`.

### `src/window.rs` (new)

New dependencies, same objc2 0.3 family as `objc2-metal` (exact versions
checked against objc2 0.6.5 when planning): `objc2-app-kit`,
`objc2-metal-kit`, `objc2-quartz-core`.

- `RtcView` (`MTKView` subclass): `framebufferOnly = NO`, `BGRA8Unorm`,
  drawable size = scene `width × height` pixels; the window is scaled down to
  fit the screen if needed, rendering stays full resolution. Overrides
  `mouseDown`, `mouseDragged`, `scrollWheel` → orbit → `set_camera` →
  `reset()` → unpause.
- `RtcDelegate` (`MTKViewDelegate`): `drawInMTKView` runs `add_frame`,
  `resolve_to_texture` into `currentDrawable.texture`, present, commit,
  frame += 1. At the sample target: `paused = YES`,
  `enableSetNeedsDisplay = YES`, so the GPU idles until the next input. No
  drawable → skip the frame, do not count it.
- App: `NSApplication`, activation policy Regular, an
  `NSApplicationDelegate` that terminates after the last window closes and
  in `applicationWillTerminate` reads back `resolve_to_buffer` and writes the
  PPM (whatever is accumulated, even short of the target). Window title:
  `<scene> — n/N samples, x ms/frame`.
- Everything on the main thread. Slow frames (large meshes at native
  resolution) just lower MTKView's frame rate.

### `src/main.rs`

```
rtc-metal scene.json [-o out.ppm] [--samples N] [--window]
```

- No `--samples`, no `--window`: unchanged (one deterministic frame).
- `--samples N`: N sampled frames headless, each committed and waited on,
  then `resolve_to_buffer` → PPM; progress about once a second on stderr and
  total time at the end.
- `--window`: progressive to N (default 256), PPM on close.
- The usage text notes that `antialias.n` is ignored on the GPU.

## Testing

Test names follow the book: GPU soft-shadow tests are
`test_bonus_soft_shadows_gpu_*`.

### Probes (`src/probes.rs`, `shaders/probes.metal`)

- RNG: 10k values all in [0, 1), mean 0.5 ± 0.01; neighbouring pixels and
  consecutive frames differ; same `(pixel, frame, dimension)` gives the same
  value.
- `point_on_light` with jitter 0.5 matches rtc-rust's cell points for every
  cell (ported from rtc-rust's bonus soft-shadow tests).
- Deterministic area-light `lighting` matches rtc-rust `lighting()` with a
  2×2 light.
- Sampled shading of one point in the default world's penumbra, under
  rtc-rust's 2×2 soft-shadow test light: the mean of 4096 sampled frames is
  within 0.01 of the deterministic shading of the same light subdivided
  64×64 (the continuous-light reference; checks the two-independent-cells
  argument).

### Whole-image compares (`src/compare.rs`)

- Every existing compare unchanged (deterministic path).
- New `soft_shadows_scene()` (the book's soft-shadow scene, shrunk):
  deterministic GPU vs CPU pixel-exact within the existing budget.
- Convergence: same scene, 256 sampled frames vs the CPU; mean absolute
  error below a threshold calibrated once and documented in the test (it
  does not go to 0: the CPU uses cell centres and pixel centres); and
  the 256-sample error is at most half the 4-sample error (a sampler that does
  not converge fails).
- `Accumulator`: after `reset()` frame 0 overwrites stale data; when every
  sample is the same colour (camera inside an ambient-only sphere), resolving
  after 4 frames gives exactly that colour (catches an off-by-one divisor).

### Orbit (`src/orbit.rs`)

Unit tests: spherical → `from` round-trips the scene camera's `from`; pitch
clamps at ±89°; zero drag gives the original view transform; dragging right
by a quarter turn orbits the camera left; positive dy raises the camera; one
zoom unit scales the radius by exp(-0.1); zooming far out stays finite and
within 1e3 × the initial radius.

### Window (manual)

Not automatable. Run a book scene and the chapter 15 teapot: opens,
converges, pauses at N, drag restarts accumulation and orbits the expected
way, scroll zooms, closing writes the PPM.

## Tasks

1. RNG + probe tests.
2. Area light in `GpuLight` / `lighting` (deterministic) + probes +
   soft-shadow compare scene.
3. Sampled `render`, `Accumulator`, `resolve`, `--samples`, convergence test.
4. `orbit.rs` + unit tests.
5. Window (`RtcView`, delegate, app lifecycle) + `--window`.
6. Usage/README, `scenes/soft_shadows.json`, manual check.

Each task leaves the build and tests green; tasks 1–3 already give headless
progressive rendering.

## Out of scope

Focal blur, motion blur, path tracing (a later milestone), spot lights,
camera controls beyond orbit and zoom, a worker thread for rendering,
resizing the window to change the render resolution, Metal's ray tracing API
(milestone 9).
