# rtc-metal milestone 12: path tracing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** rtc-rust's chapter 17 path-traced indirect lighting
(`World::set_path_tracing(samples, depth)`) on the GPU: deterministic
renders average `samples` paths per pixel, sampled renders one path per
frame, both converging to the CPU image's expected value.

**Architecture:** `GpuSceneInfo` carries `gi_samples`/`gi_depth`.
`color_at` traces unbranched paths: a diffuse hit with bounces left pushes
one cosine-weighted bounce ray weighted by its albedo, and ambient is off
while path tracing is on. `render` loops `lens_n · gi_n` paths per
deterministic pixel through its single `color_at` call site. Tests compare
statistically (mean absolute and mean signed error), since rtc-rust's
f64-bit hash of the hit point can't be reproduced in f32.

**Tech Stack:** Rust 2024, objc2-metal, MSL compiled at runtime; rtc-rust
as the f64 reference (unchanged).

**Spec:** `RUST/rtc-metal/docs/superpowers/specs/2026-10-14-metal-12-path-tracing-design.md`

## Global Constraints

- Branch `metal/12-path-tracing` (exists, spec committed as `fdd62bb`);
  commit per task; never push.
- Commit messages `ENHANCEMENT:`/`TEST:`/`DOCS:` prefix, `(milestone 12)`
  suffix on the subject, ending with
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- Test names `test_chap_17_*` (path tracing is the book's chapter 17 bonus).
- Run tests from `RUST/rtc-metal` with `cargo test --release` (baseline:
  109 pass); clippy clean (`cargo clippy --release --all-targets`).
- Every Rust `#[repr(C)]` struct keeps its `static_assert(sizeof…)` twin in
  MSL in step. `GpuSceneInfo` becomes 224 bytes.
- `GI_DIM = 0x4000_0000`; bounce `b` of a path draws dims `GI_DIM + 2b`
  and `GI_DIM + 2b + 1` with `frame = path index`.
- `render.metal` keeps exactly ONE `color_at` call site (milestone 11's
  METAL_RT miscompile).
- With path tracing off, no existing test threshold changes.
- No rtc-rust changes.

## Review Focus

1. **`samples > 0`, `depth: 0`** — rtc-rust turns ambient off but adds no
   indirect light; the image is deterministic, so the GPU must match it
   pixel for pixel through `check_matches_cpu`. Test in Task 3.
2. **Path tracing with mirrors and glass** — reflected/refracted rays keep
   `gi_depth`; a bug there (resetting or decrementing it) shows as a signed
   bias in reflections. Test in Task 3 on chapter 11's scene.
3. **Path tracing with focal blur** — `n = lens_n · gi_n`, lens sample
   `i / gi_n`; mixing up the index repeats or skips lens points. Test in
   Task 3 on chapter 17's focal blur scene.
4. **Negative `depth`** in JSON — CPU treats `gi_depth <= 0` as no
   bounces; `from_world` must store 0, not a wrapped `u32`. Test in Task 1.
5. **Metal traversal** — `color_at` grows; the METAL_RT pipeline once
   returned ~4e9 garbage. Every image test here runs both traversals.

---

### Task 1: Path-tracing settings in `GpuSceneInfo`

**Files:**
- Modify: `src/types.rs` (`GpuSceneInfo`, size assert 208 → 224)
- Modify: `src/scene.rs` (`GpuScene::from_world`: drop the "path tracing
  is not supported so far" error; fill the fields; tests module)
- Modify: `shaders/rt.metal` (`SceneInfo` mirror, `static_assert` 224)

**Interfaces:**
- Produces: `GpuSceneInfo { …, frame, sampling, gi_samples: u32,
  gi_depth: u32, _pad: [u32; 2] }`; MSL `SceneInfo` fields `gi_samples`,
  `gi_depth`, `pad[2]` (uint).

- [ ] **Step 1: Write the failing test** in `src/scene.rs` tests:

```rust
#[test]
fn test_chap_17_path_tracing_settings() {
    let scene = |settings: &str| {
        let json = format!(r#"{{ "camera": {{ "width": 4, "height": 4, "field_of_view": 1.0,
            "from": [0, 0, -5], "to": [0, 0, 0], "up": [0, 1, 0] }},
            "light": {{ "type": "point", "position": [0, 5, -5], "intensity": [1, 1, 1] }},
            "shapes": [ {{ "type": "sphere" }} ] {settings} }}"#);
        let (world, camera) = rtc_rust::scene::SceneDescription::from_json(&json)
            .unwrap().build(std::path::Path::new(".")).unwrap();
        let info = GpuScene::from_world(&world, &camera).unwrap().info;
        (info.gi_samples, info.gi_depth)
    };
    assert_eq!(scene(""), (0, 0));
    assert_eq!(scene(r#", "settings": { "path_tracing": { "samples": 16, "depth": 2 } }"#), (16, 2));
    assert_eq!(scene(r#", "settings": { "path_tracing": { "samples": 4, "depth": -3 } }"#), (4, 0));
}
```

(Adjust the JSON keys to whatever the existing scene tests in this module
use if any differ.)

- [ ] **Step 2: Run** — `cargo test --release path_tracing_settings` →
  compile error (no `gi_samples` field).
- [ ] **Step 3: Implement** the fields, `from_world` filling them
  (`world.gi_depth.max(0) as u32`), the MSL mirror and both size asserts.
- [ ] **Step 4: Run** — `cargo test --release` → 110 pass; clippy clean.
- [ ] **Step 5: Commit** — `ENHANCEMENT: rtc-metal GpuSceneInfo carries the path tracing settings (milestone 12)`

### Task 2: `cosine_direction` in MSL

**Files:**
- Modify: `shaders/rt.metal` (new `cosine_direction`, next to `lighting`)
- Modify: `shaders/probes.metal` (`CosineProbe`, `probe_cosine_direction`)
- Modify: `src/probes.rs` (`CosineProbe`, test)

**Interfaces:**
- Produces: `inline float4 cosine_direction(float4 normal, float r1, float r2)`
  — rtc-rust's `world::cosine_direction`: `axis = |n.x| > 0.9 ? (0,1,0) :
  (1,0,0)`, `tangent = normalize(cross(n, axis))`, `bitangent = cross(n,
  tangent)`, `phi = 2π·r1`, `r = sqrt(r2)`, result
  `normalize(tangent·r·cos phi + bitangent·r·sin phi + n·sqrt(1 − r2))`,
  `w = 0`.
- Produces: `struct CosineProbe { normal: Float4, r: [f32; 2], _pad: [u32; 2] }`
  (32 bytes) and kernel `probe_cosine_direction` → `float4` per input.

- [ ] **Step 1: Write the failing test** in `src/probes.rs`:
  `test_chap_17_gpu_cosine_direction` over normals `(0,1,0)`, `(1,0,0)`,
  `(-0.95, 0.312, 0)` normalised, `(0,0,-1)`, `(1,2,3)` normalised ×
  `(r1, r2)` in `{(0,0), (0.3,0.7), (0.99,0.01), (0.5,0.999), (0.123,0.456)}`;
  assert each component within `1e-5` of
  `rtc_rust::world::cosine_direction(normal, r1, r2)` and `dot ≥ -1e-6`.
- [ ] **Step 2: Run** — `cargo test --release gpu_cosine_direction` →
  fails (no kernel `probe_cosine_direction`).
- [ ] **Step 3: Implement** `cosine_direction` and the probe kernel.
- [ ] **Step 4: Run** — `cargo test --release` → 111 pass; clippy clean.
- [ ] **Step 5: Commit** — `ENHANCEMENT: rtc-metal ports cosine_direction (milestone 12)`

### Task 3: Path-traced `color_at` and the deterministic render

**Files:**
- Modify: `shaders/rt.metal` (`GI_DIM`, `lighting` gets `ambient_scale`,
  `PendingRay`, `MAX_DEPTH` stack, `color_at`)
- Modify: `shaders/render.metal` (path loop)
- Modify: `shaders/probes.metal` (two `lighting` calls pass `1.0`)
- Create: `scenes/chapter17_path_tracing.json`
- Modify: `src/compare.rs` (`check_converges_to_cpu`, four tests)

**Interfaces:**
- Consumes: Task 1's `SceneInfo.gi_samples/gi_depth`, Task 2's
  `cosine_direction`.
- Produces: `constant uint GI_DIM = 0x40000000u;`
  `lighting(…, float intensity, float ambient_scale, thread const Sampler &s,
  device const uint *tables)` (ambient term × `ambient_scale`);
  `PendingRay { Ray ray; float3 weight; int remaining; int gi_depth; uint bounce; }`;
  `color_at(…, thread Sampler &sampler, uint path)`.
- Produces: `fn check_converges_to_cpu(world: &World, camera: Camera,
  max_abs: f64, max_signed: f64) -> Result<(), String>` in `compare.rs`:
  deterministic `gpu.render` on both traversals vs `render_parallel`;
  prints and checks `mean_abs_error` and the largest per-channel
  `|mean(gpu − cpu)|` (new `mean_signed_error`, same shape as
  `mean_abs_error`).

Shader behaviour (from the spec):
- stack `PendingRay stack[2 * MAX_DEPTH + 1]`; camera ray starts with
  `gi_depth = scene.gi_samples > 0 ? int(scene.gi_depth) : 0`, `bounce = 0`;
- `ambient_scale = scene.gi_samples > 0 ? 0.0 : 1.0`;
- after the reflect/refract pushes (which copy `gi_depth` and `bounce`),
  if `p.gi_depth > 0 && m.diffuse != 0.0` push
  `Ray{over_point, cosine_direction(normalv, random01(sampler.pixel, path, GI_DIM + 2*p.bounce), random01(…, GI_DIM + 2*p.bounce + 1))}`
  with weight `p.weight * pattern_at(m, object_point, tables) * m.diffuse`,
  `remaining - 1`, `gi_depth - 1`, `bounce + 1` (the existing
  `p.remaining <= 0 → continue` already guards it);
- `render`: `lens_n = (s.on == 0 && c.aperture > 0) ? c.lens_samples : 1`,
  `gi_n = (s.on == 0 && scene.gi_samples > 0) ? scene.gi_samples : 1`,
  `n = lens_n * gi_n`; iteration `i` uses `lens_sample(…, i / gi_n, lens_n)`
  and path `s.on == 0 ? i : scene.frame`. Still one `color_at` call.

- [ ] **Step 1: Write `scenes/chapter17_path_tracing.json`** — rtc-rust's
  `cornell_box` (`rtc-rust/src/world.rs`, `fn cornell_box`): matte materials
  (specular 0) white 0.73, red (0.65, 0.05, 0.05), green (0.12, 0.45, 0.15);
  five planes and two cubes with the same transforms; point light
  (0, 1.8, -0.6) white; camera 150×150, field of view `39π/180`, from
  (0, 1, -3.5) to (0, 1, 0), up y; `"settings": { "path_tracing":
  { "samples": 16, "depth": 2 } }`; a `comment` naming
  `test_chap_17_path_tracing_putting_it_together`.
- [ ] **Step 2: Write the failing tests** in `src/compare.rs`:
  - `test_chap_17_path_tracing_without_bounces_gpu_render_matches_cpu`:
    the Cornell box with `world.set_path_tracing(4, 0)` →
    `check_matches_cpu` (ambient off, no indirect: deterministic).
  - `test_chap_17_path_tracing_gpu_render_matches_cpu`: the scene file at
    150×150 → `check_converges_to_cpu(…, MAX_ABS, MAX_SIGNED)`.
  - `test_chap_17_path_tracing_reflection_gpu_render_matches_cpu`:
    `scenes/chapter11.json` at 100×50 with `set_path_tracing(4, 1)`.
  - `test_chap_17_path_tracing_focal_blur_gpu_render_matches_cpu`:
    `scenes/chapter17_focal_blur.json` at 160×100 with
    `set_path_tracing(2, 1)`.
  Start each `check_converges_to_cpu` test with placeholder thresholds of
  `0.0` so it fails and prints the measured errors.
- [ ] **Step 3: Run** — `cargo test --release chap_17_path_tracing` → the
  `without_bounces` test fails (ambient still on / no
  `ambient_scale`), the others fail on thresholds.
- [ ] **Step 4: Implement** the shader changes above.
- [ ] **Step 5: Bias sanity check (not committed)** — temporarily raise the
  CPU's samples for the Cornell box to 256 (`world.set_path_tracing(256, 2)`
  for the CPU render only) and confirm the GPU's per-channel signed error
  against it is ≤ ~2e-3 on both traversals. If it isn't, it's a bug —
  diff-map before calibrating.
- [ ] **Step 6: Calibrate** — run each statistical test, set
  `MAX_ABS`/`MAX_SIGNED` to 1.5× the larger of the two traversals'
  measurements, rounded up to 2 significant digits, and write the measured
  values in a comment as in `test_chap_17_focal_blur_gpu_sampled_render_converges`.
  Sanity: with the bounce push disabled, the Cornell box test must fail
  its `MAX_SIGNED` (check once, then restore).
- [ ] **Step 7: Run** — `cargo test --release` → 115 pass, no existing
  threshold changed; clippy clean; `grep -c "color_at(" shaders/render.metal` → 1.
- [ ] **Step 8: Commit** — two commits:
  `ENHANCEMENT: rtc-metal traces diffuse bounces (milestone 12)` (shaders)
  and `TEST: rtc-metal chapter 17 path tracing matches the CPU (milestone 12)`
  (JSON, `compare.rs`).

### Task 4: Sampled path tracing converges

**Files:**
- Modify: `src/compare.rs` (one test)

**Interfaces:**
- Consumes: Task 3's shaders and scene; the existing `Accumulator`
  (`gpu.accumulator(&scene)`, `add_frame`, `resolve`) and `mean_abs_error`.

- [ ] **Step 1: Write the test**
  `test_chap_17_path_tracing_gpu_sampled_render_converges`, shaped like
  `test_chap_17_focal_blur_gpu_sampled_render_converges`: the Cornell box
  scene file at 150×150; CPU reference `render_parallel` with
  `world.set_path_tracing(64, 2)` (less reference noise); assert
  `err_256 <= 0.5 · err_4` and `err_256 < THRESHOLD`, with `THRESHOLD`
  starting at 0.0.
- [ ] **Step 2: Run** — fails on the threshold, printing `err_4`/`err_256`.
  If `err_256 > 0.5 · err_4`, stop: sampled paths aren't converging
  (check that `path = scene.frame` and that GI dims don't collide with
  `2 + 6·call + slot` or `LENS_DIM`).
- [ ] **Step 3: Calibrate** `THRESHOLD = 1.5 × err_256`, 2 significant
  digits, measurements in the comment.
- [ ] **Step 4: Run** — `cargo test --release` → 116 pass; clippy clean.
- [ ] **Step 5: Commit** — `TEST: rtc-metal sampled path tracing converges (milestone 12)`

### Task 5: CLI and window check

**Files:** none expected (any fix found lands in the owning file).

- [ ] **Step 1: Run** headless:
  `cargo run --release -- scenes/chapter17_path_tracing.json` and with `-o pt-sampled.ppm`
  `--samples 256`; open both PPMs and compare by eye with
  `../rtc-rust/test_chap_17_path_tracing_putting_it_together.ppm` (colour
  bleeding from the red/green walls onto the blocks, no ambient flatness,
  no speckles at the floor/block seams — if speckles, lift the blocks 0.01
  in the JSON and re-run Task 3's tests).
- [ ] **Step 2: Run** `--window` on the same scene; orbit and zoom, check
  the image restarts and converges. Report what you saw to the user.
- [ ] **Step 3: Commit** only if a fix was needed:
  `FIX: … (milestone 12)`.
