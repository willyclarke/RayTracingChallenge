# rtc-metal milestone 11: focal blur Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** rtc-rust's chapter 17 thin-lens camera on the GPU: the
deterministic render reproduces `World::distributed`, sampled renders
converge with random lens points.

**Architecture:** `GpuCamera` carries `aperture`, `focal_distance` and
`lens_samples` in its existing padding. MSL gets exact ports of
`Camera::lens_sample` (minus shutter time) and `Camera::ray_for_lens`;
`render` loops over the lens samples when deterministic and draws one
random lens point per frame when sampling. The window scales the focal
distance with zoom.

**Tech Stack:** Rust 2024, objc2-metal, MSL compiled at runtime; rtc-rust
as the f64 reference (unchanged).

**Spec:** `RUST/rtc-metal/docs/superpowers/specs/2026-10-13-metal-11-focal-blur-design.md`

## Global Constraints

- Branch `metal/11-focal-blur` (exists, spec committed as `3c46e5d`);
  commit per task; never push.
- Commit messages `ENHANCEMENT:`/`TEST:`/`DOCS:` prefix, `(milestone 11)`
  suffix on the subject, ending with
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- Test names `test_chap_17_*` (focal blur is the book's chapter 17 bonus).
- Run tests from `RUST/rtc-metal` with `cargo test --release` (baseline:
  100 pass); clippy clean (`cargo clippy --release --all-targets`).
- Every Rust `#[repr(C)]` struct keeps its `static_assert(sizeof…)` twin in
  MSL in step. `GpuCamera` stays 96 bytes.
- Motion blur stays rejected: `gpu_camera` errors on `motion_samples > 1`,
  shapes with `motion` keep their existing errors.
- `LENS_DIM = 0x8000_0000` (sampled lens point dims `LENS_DIM`,
  `LENS_DIM + 1`).
- No rtc-rust changes.

## Review Focus

1. **`aperture == 0` with `focal_blur.samples > 1`** — the CPU's
   `is_distributed()` is false, so it renders one centre ray per pixel;
   the GPU must too (`lens_samples = 1`, pinhole). Test in Task 1.
2. **Negative aperture** — the CPU never calls `distributed` for it
   (`aperture > 0.0` is false) and renders a pinhole; `gpu_camera` must
   store `aperture = 0`, not a negative lens. Test in Task 1.
3. **`focal_blur.samples: 1` with `aperture > 0`** — one lens sample at
   `i = 0, n = 1` is *not* the lens centre (`r = sqrt(h(1))`); the GPU must
   follow the CPU, not special-case it. Case in the Task 2 probe.
4. **Focal blur on the Metal traversal** — lens rays start off the camera
   origin; both traversals must still match. Covered by Task 4's
   `check_matches_cpu` (runs both).
5. **Zoom at the radius clamp** — when `Orbit::zoom` clamps, the focal
   distance must scale by the *actual* radius ratio (1 when clamped), not
   `exp(-0.1·dy)`. Test in Task 6.

---

### Task 1: Focal blur fields in `GpuCamera`

**Files:**
- Modify: `src/types.rs` (`GpuCamera`, size assert comment)
- Modify: `shaders/rt.metal:568-576` (MSL `Camera`)
- Modify: `src/scene.rs:173-192` (`gpu_camera` + doc comment), tests in
  its `mod tests`

**Interfaces:**
- Produces: `GpuCamera { …, aperture: f32, focal_distance: f32,
  lens_samples: u32 }` (after `height`), mirrored in MSL `Camera` with
  the same names. `gpu_camera(&Camera) -> Result<GpuCamera, String>`
  fills them: `aperture = max(camera.aperture, 0)`, `focal_distance`,
  `lens_samples = camera.dof_samples` when `aperture > 0`, else `1`.

- [ ] **Step 1: Write the failing tests** in `src/scene.rs` tests:

```rust
#[test]
fn gpu_camera_takes_focal_blur() {
    let camera = Camera::new(4, 2, PI / 2.0).with_focal_blur(0.25, 8.5, 32);
    let c = gpu_camera(&camera).unwrap();
    assert_eq!((c.aperture, c.focal_distance, c.lens_samples), (0.25, 8.5, 32));
}

#[test]
fn gpu_camera_without_aperture_is_a_pinhole() {
    // Review focus 1 and 2: the CPU renders these with one centre ray.
    for aperture in [0.0, -0.1] {
        let camera = Camera::new(4, 2, PI / 2.0).with_focal_blur(aperture, 8.5, 32);
        let c = gpu_camera(&camera).unwrap();
        assert_eq!((c.aperture, c.lens_samples), (0.0, 1));
    }
}

#[test]
fn gpu_camera_rejects_motion_blur() {
    let camera = Camera::new(4, 2, PI / 2.0).with_motion_blur(4);
    assert!(gpu_camera(&camera).unwrap_err().contains("motion blur"));
}

#[test]
fn json_focal_blur_reaches_the_gpu_camera() {
    let json = r#"{
      "camera": { "width": 4, "height": 2, "field_of_view": 1.0,
                  "from": [0, 0, -5], "to": [0, 0, 0], "up": [0, 1, 0],
                  "focal_blur": { "aperture": 0.25, "focal_distance": 8.5, "samples": 16 } },
      "light": { "type": "point", "position": [-10, 10, -10], "intensity": [1, 1, 1] },
      "shapes": [ { "type": "sphere" } ]
    }"#;
    // Build as the existing JSON test at src/scene.rs:738 does.
    let (world, camera) = SceneDescription::from_json(json).unwrap()
        .build(std::path::Path::new(".")).unwrap();
    let c = GpuScene::from_world(&world, &camera).unwrap().info.camera;
    assert_eq!((c.aperture, c.focal_distance, c.lens_samples), (0.25, 8.5, 16));
}
```

- [ ] **Step 2: Run** — `cargo test --release gpu_camera json_focal_blur`
  → compile error (no such fields).
- [ ] **Step 3: Implement** the three fields in `GpuCamera` and MSL
  `Camera`, update the assert comment to `64 + 8 * 4 = 96 bytes`;
  `gpu_camera` rejects only `motion_samples > 1` with
  `"motion blur is not supported so far"` and fills the fields as in
  Interfaces. Update its doc comment: aperture/`dof_samples` are used,
  anti-aliasing samples are still ignored.
- [ ] **Step 4: Run** — `cargo test --release` → 104 pass (the GPU still
  ignores the fields, so every image test is unchanged); clippy clean.
- [ ] **Step 5: Commit** — `ENHANCEMENT: rtc-metal GpuCamera carries the
  focal blur lens (milestone 11)`.

### Task 2: `lens_sample` in MSL

**Files:**
- Modify: `shaders/rt.metal` (after `ray_for_pixel`, ~line 727)
- Modify: `shaders/probes.metal` (new probe after `probe_ray_for_pixel`)
- Test: `src/probes.rs` (after `test_chap_7_gpu_ray_for_pixel`)

**Interfaces:**
- Produces (MSL): `inline float4 lens_sample(uint px, uint py, uint i,
  uint n)` returning `(lx, ly, dx, dy)` = the first four values of
  rtc-rust's `Camera::lens_sample(px, py, i, n)`.

Port of the hash, all on `ulong` (Rust's `wrapping_mul` = MSL's wrapping
`ulong *`):

```metal
// h(k) of rtc-rust's Camera::lens_sample: exact, (x & 0xFFFFFF) / 2^24.
inline float lens_hash(uint px, uint py, uint i, ulong k) {
    ulong x = ulong(px) * 0x9E3779B9ul ^ ulong(py) * 0x85EBCA6Bul
            ^ ulong(i) * 0xC2B2AE35ul ^ k * 0x27D4EB2Ful;
    x ^= x >> 15; x *= 0x2C1B3C6Dul;
    x ^= x >> 12; x *= 0x297A2D39ul;
    x ^= x >> 15;
    return float(x & 0xFFFFFFul) * 0x1p-24f;
}
```

`lens_sample`: `r = sqrt((i + h(1)) / n)`, `theta = 2π·h(2)`, return
`(r cos θ, r sin θ, h(3), h(4))`. The time stratification is not ported.

- [ ] **Step 1: Write the failing test** `test_chap_17_gpu_lens_sample`:
  a `#[repr(C)] struct LensProbe { pixel: [u32; 2], i: u32, n: u32 }`
  (16 bytes); cases `(0,0,0,1)`, `(0,0,0,32)` (review focus 3 is the
  first), `(199,124,31,32)`, `(57,3,7,32)`, `(400,250,15,16)`,
  `(1,1,0,1)`; run `"probe_lens_sample"` → `Vec<Float4>`; for each, the
  four values from `Camera::lens_sample(px, py, i, n)` agree via
  `check_close` (existing helper) for `lx, ly`, and **exactly**
  (`got == want as f32`) for `dx, dy`, which proves the u64 hash is
  bit-identical.
- [ ] **Step 2: Run** — `cargo test --release test_chap_17_gpu_lens_sample`
  → FAIL (no kernel `probe_lens_sample`).
- [ ] **Step 3: Implement** `lens_hash`, `lens_sample`, and
  `kernel void probe_lens_sample(device const LensProbe *in [[buffer(0)]],
  device float4 *out [[buffer(1)]], uint i [[thread_position_in_grid]])`
  with `static_assert(sizeof(LensProbe) == 16, …)`.
- [ ] **Step 4: Run** — `cargo test --release` → 105 pass.
- [ ] **Step 5: Commit** — `ENHANCEMENT: rtc-metal ports rtc-rust's lens
  samples (milestone 11)`.

### Task 3: `ray_for_lens` in MSL

**Files:**
- Modify: `shaders/rt.metal` (after `lens_sample`)
- Modify: `shaders/probes.metal`
- Test: `src/probes.rs`

**Interfaces:**
- Consumes: Task 1's `Camera.aperture`, `Camera.focal_distance`.
- Produces (MSL): `inline Ray ray_for_lens(Camera c, uint px, uint py,
  float dx, float dy, float lx, float ly)` — rtc-rust's
  `Camera::ray_for_lens`: `pinhole = ray_for_pixel(c, px, py, dx, dy)`;
  if `c.aperture == 0` return it; else `focal = pinhole.origin +
  pinhole.direction * c.focal_distance`, `origin = c.inverse *
  float4(lx * a, ly * a, 0, 1)`, direction `normalize(focal - origin)`.

- [ ] **Step 1: Write the failing test** `test_chap_17_gpu_ray_for_lens`:
  `#[repr(C)] struct LensRayProbe { camera: GpuCamera, pixel: [u32; 2],
  offsets: [f32; 4] }` (`dx, dy, lx, ly`; 128 bytes with padding — assert
  `size_of == 128` in Rust). Cameras: rtc-rust's chapter 17 view
  (`Camera::new(400, 250, PI / 3.5)`, `view_transform((0, 2.2, -6),
  (-0.5, 1, 2), (0, 1, 0))`, `with_focal_blur(0.25, 8.616989, 32)`) and
  the same without focal blur. Cases on the focal camera: pixels (200,
  125), (0, 0), (399, 249) with offsets `(0.5, 0.5, 0, 0)`,
  `(0.3, 0.6, 0.7, -0.2)`, `(0.9, 0.1, -0.5, 0.5)`; on the pinhole camera
  `(0.3, 0.6, 0.7, -0.2)` must equal `ray_for_subpixel(.., 0.3, 0.6)`.
  Compare against `camera.ray_for_lens(...)` with `check_tuple` on origin
  and direction.
- [ ] **Step 2: Run** → FAIL (no kernel `probe_ray_for_lens`).
- [ ] **Step 3: Implement** `ray_for_lens` and `probe_ray_for_lens`
  (output `Ray`), with the `static_assert(sizeof(LensRayProbe) == 128)`.
- [ ] **Step 4: Run** — `cargo test --release` → 106 pass.
- [ ] **Step 5: Commit** — `ENHANCEMENT: rtc-metal casts rays from the
  lens (milestone 11)`.

### Task 4: Deterministic focal blur render

**Files:**
- Create: `scenes/chapter17_focal_blur.json`
- Modify: `shaders/render.metal` (`render`, `sampling == 0` branch; header
  comment)
- Test: `src/compare.rs` (after the chapter 16 compare)

**Interfaces:**
- Consumes: Task 2's `lens_sample`, Task 3's `ray_for_lens`.
- Produces: `scenes/chapter17_focal_blur.json`, used by Tasks 5 and 6.

The scene is rtc-rust's `test_chap_17_focal_blur_putting_it_together`
(`rtc-rust/src/world.rs:5203`) as JSON, with a `"comment"` saying so:
camera 400×250, `field_of_view` 0.8975979010 (π/3.5), from
`[0, 2.2, -6]`, to `[-0.5, 1, 2]`, up `[0, 1, 0]`, `"focal_blur":
{ "aperture": 0.25, "focal_distance": 8.616989, "samples": 32 }`
(distance from `from` to the second sphere at `(-0.75, 1, 2.5)`); point
light `(-8, 10, -10)` white; a plane with checkers `(0.35×3)` /
`(0.7×3)` and `specular 0`; five unit spheres translated to
`(-1.5 + 0.75 i, 1, 2.5 i)` for `i = 0..5`, colors `(0.9,0.2,0.2)`,
`(0.2,0.7,0.2)`, `(0.2,0.3,0.9)`, `(0.9,0.8,0.2)`, `(0.8,0.3,0.8)`, each
`specular 0.6`, `shininess 40`.

- [ ] **Step 1: Write the failing test**:

```rust
#[test]
fn test_chap_17_focal_blur_gpu_render_matches_cpu() -> Result<(), String> {
    // Every pixel averages 32 lens rays (rtc-rust's World::distributed).
    let (world, camera) = scene_file(include_str!("../scenes/chapter17_focal_blur.json"), 200, 125);
    check_matches_cpu(&world, camera)
}
```

- [ ] **Step 2: Run** — `cargo test --release test_chap_17_focal_blur` →
  FAIL: far over budget (the GPU renders the pinhole image).
- [ ] **Step 3: Implement** in `render`'s `s.on == 0` branch: if
  `c.aperture > 0`, sum `color_at(ray_for_lens(c, x, y, l.z, l.w, l.x,
  l.y), …)` over `i < c.lens_samples` with `l = lens_sample(x, y, i,
  c.lens_samples)`, divide by `lens_samples`; otherwise unchanged. Extend
  the kernel's header comment.
- [ ] **Step 4: Run** — `cargo test --release` → 107 pass. If the compare
  is over budget only at silhouettes (check the eprintln'd worst pixel),
  do not loosen `check_matches_cpu`: report the numbers and stop for a
  ruling (spec: "measure before loosening anything").
- [ ] **Step 5: Smoke-run the CLI** — `cargo run --release --
  scenes/chapter17_focal_blur.json -o chapter17_focal_blur.ppm` writes an
  image with the second sphere sharp and the others blurred (look at it).
- [ ] **Step 6: Commit** — `TEST: rtc-metal chapter 17 focal blur render
  matches the CPU (milestone 11)` (scene, shader, test).

### Task 5: Sampled focal blur

**Files:**
- Modify: `shaders/rt.metal` (`LENS_DIM` constant next to the `Sampler`
  dims comment, which gains a line about it)
- Modify: `shaders/render.metal` (`render`, sampled branch)
- Test: `src/compare.rs`

**Interfaces:**
- Consumes: Task 3's `ray_for_lens`, Task 4's scene file.
- Produces: `constant uint LENS_DIM = 0x80000000u;`

- [ ] **Step 1: Write the failing test**
  `test_chap_17_focal_blur_gpu_sampled_render_converges`, structured like
  `test_bonus_soft_shadows_gpu_sampled_render_converges`: the chapter 17
  scene at 200×125, CPU reference `render_parallel`, `err_4` after 4
  frames, `err_256` after 256; require `err_256 <= 0.5 * err_4` and
  `err_256 < THRESHOLD`. Start with `THRESHOLD = 1.0` and a comment
  `// Calibrated once: …` placeholder to fill in at Step 4.
- [ ] **Step 2: Run** → FAIL on the halving check (sampled frames are
  still pinhole, so the error to the blurred CPU image stays put).
- [ ] **Step 3: Implement**: in the sampled branch, when `c.aperture > 0`
  take `u = random01(pixel, s.frame, LENS_DIM)`, `v = random01(pixel,
  s.frame, LENS_DIM + 1)`, lens point `sqrt(u)·(cos 2πv, sin 2πv)`, and
  cast `ray_for_lens` with the dims 0–1 jitter; else `ray_for_pixel` as
  now.
- [ ] **Step 4: Calibrate** — run, read `err_4`/`err_256` from stderr,
  set `THRESHOLD` = 1.5 × `err_256` rounded up to 2 significant digits,
  write the measured values into the comment (as the soft-shadow test
  does: the error does not reach 0 because the CPU uses 32 stratified
  samples, the GPU the continuous lens).
- [ ] **Step 5: Run** — `cargo test --release` → 108 pass; clippy clean.
- [ ] **Step 6: Commit** — `ENHANCEMENT: rtc-metal samples the lens in
  progressive renders (milestone 11)`.

### Task 6: Window keeps the subject in focus

**Files:**
- Modify: `src/orbit.rs` (`radius` getter + test)
- Modify: `src/window.rs` (`RtcView::interact` or `Renderer::restart`)

**Interfaces:**
- Produces: `Orbit::radius(&self) -> f64`.

- [ ] **Step 1: Write the failing test** in `orbit.rs`:

```rust
#[test]
fn radius_reports_the_clamped_zoom() {
    // Review focus 5: the window scales the focal distance by the actual
    // radius ratio, which is 1 once zoom is clamped.
    let mut o = Orbit::new(p(0.0, 0.0, -5.0), p(0.0, 0.0, 0.0), up());
    assert!((o.radius() - 5.0).abs() < 1e-12);
    o.zoom(-1e6);
    let r = o.radius();
    o.zoom(-1.0);
    assert_eq!(o.radius(), r);
}
```

- [ ] **Step 2: Run** — `cargo test --release radius_reports` → compile
  error.
- [ ] **Step 3: Implement** `radius()`; in `RtcView::interact`, read
  `orbit.radius()` before and after `change`, and multiply
  `r.camera.focal_distance` by `after / before` before `restart()`
  (drags leave the ratio at 1).
- [ ] **Step 4: Run** — `cargo test --release` → 109 pass; clippy clean.
- [ ] **Step 5: Try it** — `cargo run --release --
  scenes/chapter17_focal_blur.json --window`: the second sphere is sharp
  and stays sharp while zooming; orbiting changes the view. (Manual; ask
  the user to look if you cannot.)
- [ ] **Step 6: Commit** — `ENHANCEMENT: rtc-metal window zoom keeps the
  focal plane on the subject (milestone 11)`.
