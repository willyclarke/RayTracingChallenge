# rtc-metal milestone 11: focal blur

Branch `metal/11-focal-blur`, merged into develop with `--no-ff` when done.
rtc-rust (f64, CPU) stays the reference. Milestone 12 will be path tracing.

## Goal

Render rtc-rust's chapter 17 focal blur (a thin-lens camera) on the GPU, in
both render modes:

- **deterministic** (no `--samples`): the image rtc-rust's
  `World::distributed` gives, up to f32 differences;
- **sampled** (`--samples N`, `--window`): random lens points that converge
  together with the existing anti-aliasing and area-light soft shadows.

Success:
- `test_chap_17_focal_blur_gpu_render_matches_cpu` passes on both
  traversals through `check_matches_cpu`;
- probes show the GPU's `lens_sample` and `ray_for_lens` agree with the CPU;
- a sampled render of a focal-blur scene converges towards the CPU image;
- a JSON scene with a `focal_blur` block renders from the CLI and in
  `--window`.

## Out of scope

- **Motion blur.** It needs a per-ray time and per-shape velocity in every
  intersection, on both traversals. `gpu_camera` keeps rejecting
  `motion_samples > 1` and shapes with `motion` stay rejected.
- **rtc-rust changes.** None expected; `Camera::lens_sample`,
  `ray_for_lens`, `aperture`, `focal_distance` and `dof_samples` are already
  public.

## Decisions

- **Port `Camera::lens_sample` exactly** rather than only sampling randomly,
  so the deterministic image keeps matching rtc-rust (the invariant every
  chapter compare relies on). Only the parts focal blur uses are ported:
  the hash, the stratified radius `sqrt((i + h(1)) / n)`, the angle
  `2π·h(2)` and the sub-pixel jitter `(h(3), h(4))`. The shutter time (and
  its gcd/stride permutation) is motion blur's and is left out.
- **The hash runs on `ulong`.** MSL's 64-bit integers wrap like Rust's
  `wrapping_mul`, and `(x & 0xFFFFFF) / 2^24` is exact in f32, so `h(k)` is
  bit-identical to the CPU; only `sqrt`, `cos` and `sin` round differently.
- **One dispatch per deterministic frame.** The kernel loops over the
  camera's `n` lens samples per pixel, like `distributed`. Fallback if a test
  scene hits the GPU watchdog: `n` dispatches summed through `Accumulator`.
- **Lens points in sampled mode come from two new `random01` dimensions**,
  `LENS_DIM` and `LENS_DIM + 1` with `LENS_DIM = 0x8000_0000`, far above
  the light sampler's `2 + 6·call + slot`. Existing sampled sequences are
  untouched. The disk point is `sqrt(u)` radius, `2π·v` angle, like the CPU.
- **Window zoom scales the focal distance** by the same factor as the orbit
  radius, so the subject stays in focus; orbiting leaves it alone.

## Host side

`GpuCamera` (`src/types.rs`) takes the three new fields from its padding
(84 bytes of fields today, padded to 96):

```rust
pub struct GpuCamera {
    pub inverse: GpuMat4,
    pub half_width: f32,
    pub half_height: f32,
    pub pixel_size: f32,
    pub width: u32,
    pub height: u32,
    /// Lens radius in world units; 0 = pinhole.
    pub aperture: f32,
    /// Distance at which the lens is in focus.
    pub focal_distance: f32,
    /// rtc-rust's `dof_samples`: lens samples per deterministic pixel.
    pub lens_samples: u32,
}
```

Size stays 96 (`64 + 8·4 = 96`); the MSL `static_assert` and the Rust
`const` assert keep guarding it. `gpu_camera` (`src/scene.rs`) fills the
fields and only rejects `motion_samples > 1`. With `aperture == 0` it sets
`lens_samples = 1`, matching the CPU's `is_distributed()` being false.

`window.rs`: after `Orbit::zoom`, multiply `camera.focal_distance` by the
ratio of new to old orbit radius before `set_camera`.

## Shader side (`shaders/rt.metal`, `shaders/render.metal`)

```metal
// rtc-rust's Camera::lens_sample, without the shutter time: lens point
// (x, y) in the unit disk and sub-pixel jitter (z, w).
inline float4 lens_sample(uint px, uint py, uint i, uint n);

// rtc-rust's Camera::ray_for_lens: from (lx, ly) * aperture on the lens
// (camera z = 0) through the pinhole ray's point at focal_distance.
inline Ray ray_for_lens(Camera c, uint px, uint py, float dx, float dy,
                        float lx, float ly);
```

`render`:
- `sampling == 0`, `aperture == 0`: unchanged (one ray through the centre).
- `sampling == 0`, `aperture > 0`: average `color_at` over
  `i in 0..lens_samples` of `ray_for_lens(c, x, y, s.z, s.w, s.x, s.y)`
  with `s = lens_sample(x, y, i, lens_samples)`.
- `sampling == 1`: jitter from dims 0–1 as today; lens point from
  `LENS_DIM`/`LENS_DIM + 1` when `aperture > 0`; one `ray_for_lens` per
  frame.

The other scene kernels (probes, trace) keep `ray_for_pixel`.

## Tests

TDD, CPU as reference, names per the book's chapter (17):

- `test_chap_17_gpu_lens_sample`: a probe kernel returns `lens_sample` for
  a spread of `(px, py, i, n)`; matches `Camera::lens_sample`'s
  `(lx, ly, dx, dy)` within f32 tolerance.
- `test_chap_17_gpu_ray_for_lens`: probe against `Camera::ray_for_lens`
  for a few cameras and lens offsets, including `aperture == 0` giving the
  pinhole ray.
- `test_chap_17_focal_blur_gpu_render_matches_cpu`: rtc-rust's
  `test_chap_17_focal_blur_putting_it_together` scene (five spheres on a
  checkered floor, aperture 0.25, 32 samples) at reduced size, through
  `check_matches_cpu` on both traversals. If the 1e-3 per-pixel budget is
  too tight for averaged lens rays, calibrate it from measurements and say
  so in the test.
- `test_chap_17_focal_blur_gpu_sampled_render_converges`: mean absolute
  error to the CPU image falls with more frames, in the style of
  `test_bonus_soft_shadows_gpu_sampled_render_converges`.
- `scene.rs`: a JSON `focal_blur` block loads into `GpuCamera`;
  `motion_blur` still errors.

## Risks

- **u64 hash mismatch** (MSL `ulong` multiply or shift behaving
  differently): the `lens_sample` probe catches it before any image test.
- **Long deterministic dispatches** with 32 lens samples on heavy scenes:
  keep compare scenes small; fall back to per-sample dispatches if needed.
- **f32 edge differences averaged over 32 rays** may exceed the compare
  budget at silhouettes; measure before loosening anything.
