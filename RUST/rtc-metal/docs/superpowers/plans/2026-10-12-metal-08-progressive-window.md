# rtc-metal milestone 8: progressive accumulation and a live window — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Accumulate jittered samples (sub-pixel antialiasing, area-light soft shadows) on the GPU, headless (`--samples N`) or in a native macOS window with an orbit camera (`--window`), converging to the rtc-rust image.

**Architecture:** A stateless PCG hash gives each `(pixel, frame, dimension)` a random number. `render` either writes one deterministic frame (as today, now also handling area lights exactly like rtc-rust) or adds one sampled frame to an accumulation buffer. A new `Accumulator` owns that buffer, the scene buffers and the pipelines; `resolve` kernels divide by the frame count into a buffer (PPM) or the window's drawable. The window is an objc2 `MTKView` subclass plus delegates, driving the same `Accumulator`.

**Tech Stack:** Rust 2024, objc2 0.6.5, objc2-metal / objc2-foundation 0.3.2, new objc2-app-kit / objc2-metal-kit / objc2-quartz-core 0.3.2, MSL compiled at runtime, rtc-rust as path dependency and CPU reference.

**Spec:** `RUST/rtc-metal/docs/superpowers/specs/2026-10-12-metal-08-progressive-window-design.md`

## Global Constraints

- All commands run from `RUST/rtc-metal`. Tests: `cargo test --release`. Lint: `cargo clippy --release --all-targets` must report 0 warnings; run `cargo fmt` before each commit.
- Rust structs in `types.rs` mirror MSL structs byte for byte. Each size change updates both the Rust `const _: () = assert!(size_of::<…>() == N)` and the MSL `static_assert`.
- `GpuLight` / MSL `Light` becomes 80 bytes. `GpuSceneInfo` / `SceneInfo` becomes 208 bytes (96 camera + 80 light + 12 counts + `frame` + `sampling` = 196, padded to 16). `LightingProbe` grows from 256 to 304 bytes (144 material + 80 light + 4 × 16 + `intensity` + `pixel`, `frame`, `sampling` = 304).
- Deterministic mode (`sampling == 0`) must reproduce rtc-rust exactly: pixel-centre rays, cell centres (jitter 0.5), and rtc-rust's adaptive corner rule in `shadow_fraction` (only when `usteps >= 2 && vsteps >= 2`). Existing compare budget unchanged: at most `w * h / 2000` pixels off by more than 1e-3.
- Sampled mode: cell A (shadow ray) and cell B (diffuse + specular) are drawn independently; no corner shortcut.
- Random dimensions: 0–1 sub-pixel jitter; shading call k (the k-th ray `color_at` pops, k from 0) uses the six dimensions `2 + 6k … 7 + 6k` (table in Task 1 Interfaces).
- `antialias.n` is accepted and ignored on the GPU. `focal_blur`, `motion_blur`, path tracing and spot lights are rejected with "not supported so far".
- Window default sample target: 256. `--samples 0` is a usage error.
- Test names follow the book: `test_bonus_soft_shadows_gpu_…` for area-light and sampling tests.
- Branch `metal/08-progressive-window` off develop. Commit after each task with the repo's `ENHANCEMENT:`/`TEST:`/`DOCS:` prefix style and the `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>` trailer. Never push. Merge into develop only when the user asks. Do not stage `scenes/chapter10.json` (user's uncommitted local change) or the repo-root `.gitignore`.
- rtc-rust is not modified.

## Review Focus

1. **Area light with `usteps` or `vsteps` = 1** (rtc-rust skips the corner rule and counts every cell): deterministic GPU intensity equals the CPU's. Pinned in Task 2, step 1.
2. **Point-light scene in sampled mode** (most existing scenes): accumulation only smooths edges; the 64-sample image stays within a small mean error of the deterministic one, and shadows never go soft. Pinned in Task 3, step 1.
3. **Window closed (or Cmd-Q) before any frame finished**: no divide by zero, no garbage PPM — `resolve` returns `Err("no samples accumulated")` and nothing is written. Pinned in Task 3, step 1.
4. **Drag while frames are in flight**: the next frame after `reset()` overwrites `accum`, so no stale samples leak in. Pinned in Task 3, step 1 (reset-overwrites test).
5. **Scene larger than the screen** (the local 3456×2234 `chapter10.json`): window scaled to fit, render at full resolution, mouse deltas still feel right. Pinned in Task 5, step 4 (manual check).

---

### Task 1: Random numbers on the GPU

**Files:**
- Modify: `shaders/rt.metal` (new "Bonus: sampling" section near the top, after `Ray`)
- Modify: `shaders/probes.metal` (`RandomProbe`, `probe_random`)
- Test: `src/probes.rs`

**Interfaces:**
- Produces (MSL, used by Tasks 2–3):
  - `uint pcg_hash(uint v)`
  - `float random01(uint pixel, uint frame, uint dim)` → [0, 1)
  - `struct Sampler { uint pixel; uint frame; uint on; uint call; };` — `on` 0/1; `call` = index of the current shading call in `color_at`.
  - `float sample(thread const Sampler &s, uint slot)` = `random01(s.pixel, s.frame, 2 + 6 * s.call + slot)`. Slots: 0 cell A index, 1 cell A jitter u, 2 cell A jitter v, 3 cell B index, 4 cell B jitter u, 5 cell B jitter v. A and B never share a slot, which is what makes them independent.
  - Sub-pixel jitter uses dims 0 and 1.
- Rust probe struct `RandomProbe { pixel: u32, frame: u32, dim: u32 }` (12 bytes) → `f32`.

- [ ] **Step 0: Branch** — `git switch -c metal/08-progressive-window` from develop; commit the spec and this plan first as `DOCS: rtc-metal milestone 8 design and plan (progressive, window)`.

- [ ] **Step 1: Write the failing tests** in `src/probes.rs`
  - `test_bonus_soft_shadows_gpu_random_in_unit_interval`: 10 000 probes (pixel = i, frame = 7, dim = 3) → every value `>= 0.0 && < 1.0`; mean within 0.5 ± 0.01.
  - `test_bonus_soft_shadows_gpu_random_decorrelated`: for i in 0..1000, `(i, 0, 0)` vs `(i + 1, 0, 0)` and `(i, 0, 0)` vs `(i, 1, 0)` and `(i, 0, 0)` vs `(i, 0, 1)` differ in at least 990 of 1000 cases each.
  - `test_bonus_soft_shadows_gpu_random_deterministic`: the same 100 probes run twice give identical results.

- [ ] **Step 2: Run** `cargo test --release random` — FAIL (kernel `probe_random` not found).

- [ ] **Step 3: Implement** `pcg_hash` (PCG-RXS-M-XS 32), `random01`, `Sampler`, `sample` and `probe_random`:

```metal
inline uint pcg_hash(uint v) {
    uint s = v * 747796405u + 2891336453u;
    uint w = ((s >> ((s >> 28u) + 4u)) ^ s) * 277803737u;
    return (w >> 22u) ^ w;
}
// Top 24 bits → exactly representable floats in [0, 1).
inline float random01(uint pixel, uint frame, uint dim) {
    return float(pcg_hash(pixel ^ pcg_hash(frame ^ pcg_hash(dim))) >> 8) * 0x1p-24f;
}
```

- [ ] **Step 4: Run** `cargo test --release` — all pass.

- [ ] **Step 5: Commit** — `ENHANCEMENT: rtc-metal PCG random numbers for sampling (bonus: soft shadows)`

---

### Task 2: Area lights, deterministic

**Files:**
- Modify: `src/types.rs` (`GpuLight`, `GpuSceneInfo`, asserts)
- Modify: `src/scene.rs` (`gpu_light`, `gpu_camera`, `GpuScene::from_world`, new `GpuScene::set_camera`)
- Modify: `shaders/rt.metal` (`Light`, `SceneInfo`, `light_cell_point`, `intensity_at`, `lighting`, `color_at`)
- Modify: `shaders/render.metal` (pass a `Sampler` with `on = 0`)
- Modify: `shaders/probes.metal` (`LightingProbe`, `probe_lighting`, new `probe_point_on_light`, `probe_intensity_at`)
- Create: `scenes/soft_shadows.json`
- Test: `src/probes.rs`, `src/compare.rs`, `src/scene.rs`

**Interfaces:**
- Consumes: `Sampler`, `sample` (Task 1).
- Produces:
  - `GpuLight { corner: Float4, uvec: Float4, vvec: Float4, intensity: Float4, usteps: u32, vsteps: u32, samples: u32, _pad: u32 }` (80 bytes). Point light: `corner = position`, `uvec = vvec = 0`, steps 1, samples 1. Area light: rtc-rust's `corner`, `uvec`, `vvec` (already per-cell), `usteps`, `vsteps`, `samples`.
  - `GpuSceneInfo { camera, light, counts, frame: u32, sampling: u32 }` (208 bytes).
  - `pub fn gpu_camera(camera: &Camera) -> Result<GpuCamera, String>` no longer rejects `aa_samples > 1` (still rejects aperture and motion blur).
  - `impl GpuScene { pub fn set_camera(&mut self, camera: &Camera) -> Result<(), String> }` — replaces `info.camera` only.
  - MSL: `float4 light_cell_point(Light l, uint u, uint v, float ju, float jv)` = `corner + uvec * (u + ju) + vvec * (v + jv)`.
  - MSL: `float intensity_at(Light l, float4 point, device const Shape*, device const BvhNode*, SceneCounts, thread const Sampler &s)` — `s.on == 0`: rtc-rust's `shadow_fraction` with cell centres, including the corner rule when both steps ≥ 2; `s.on == 1`: one shadow ray to cell A (index `min(uint(sample(s,0) * samples), samples - 1)`, u = index % usteps, v = index / usteps, jitter `sample(s,1)`, `sample(s,2)`), returns 0 or 1. `is_shadowed` takes the light point instead of `light.position`.
  - MSL: `lighting(Material, Light, float4 point, float4 object_point, float4 eyev, float4 normalv, float intensity, thread const Sampler &s)` — ambient unchanged; `s.on == 0`: diffuse and specular averaged over all cell centres (rtc-rust `lighting_with_ambient`); `s.on == 1`: one cell B (slots 3, 4, 5).
  - `color_at(Ray, constant SceneInfo&, shapes, nodes, thread Sampler &s)` increments `s.call` after each shading.
  - `LightingProbe` gains `pixel: u32, frame: u32, sampling: u32` after `intensity` (Rust stays 304 bytes).
  - Probe structs: `PointOnLightProbe { light: GpuLight, u: u32, v: u32, ju: f32, jv: f32 }` → `Float4`; `IntensityProbe { light: GpuLight, point: Float4, pixel: u32, frame: u32, sampling: u32 }` → `f32` via `run_with_scene`.

- [ ] **Step 1: Write the failing tests**
  - `src/probes.rs`, `test_bonus_soft_shadows_gpu_point_on_light`: rtc-rust's `test_bonus_soft_shadows_4` light (4×2) and its five `(u, v)` cases with jitter 0.5, then `test_bonus_soft_shadows_7`'s cases with `ju = 0.3, jv = 0.7`; GPU point matches the expected CPU point (`check_tuple`).
  - `test_bonus_soft_shadows_gpu_intensity_at`: `World::default_world()` as a `GpuScene`, rtc-rust's `test_bonus_soft_shadows_5` light and points; GPU deterministic intensity equals `light.intensity_at(point, &world)`. Add a 4×1 light (`usteps = 4, vsteps = 1`, corner `(-1.5, -0.5, -5)`, uvec `(3,0,0)`, vvec `(0,1,0)`) over the same points (Review Focus 1).
  - `test_bonus_soft_shadows_gpu_lighting`: rtc-rust's `test_bonus_soft_shadows_9` cases (2×2 light, intensity 1.0); GPU `probe_lighting` (sampling 0) matches CPU `light.lighting(...)`.
  - `src/scene.rs` test `set_camera_replaces_only_the_camera`: after `set_camera` with a different transform, `info.camera` changed and `shapes`, `nodes`, `info.light`, `info.counts` are unchanged.
  - `src/compare.rs`, `test_bonus_soft_shadows_gpu_render_matches_cpu`: `scene_file(include_str!("../scenes/soft_shadows.json"), 200, 80)` through `check_matches_cpu`.
  - `scenes/soft_shadows.json` (the book's soft-shadow scene, no light cube — rtc-rust has no non-shadow-casting flag): camera 400×160, fov 0.7854, from (-3, 1, 2.5), to (0, 0.5, 0), up (0, 1, 0); area light corner (-1, 2, 4), uvec (2, 0, 0) usteps 10, vvec (0, 2, 0) vsteps 10, intensity (1.5, 1.5, 1.5); a plane (white, ambient 0.025, diffuse 0.67, specular 0); a sphere translated (0.5, 0.5, 0) scaled 0.5 (color (1, 0, 0), ambient 0.1, specular 0, diffuse 0.6, reflective 0.3); a sphere translated (-0.25, 0.33, 0) scaled 0.33 (color (0.5, 0.5, 1), ambient 0.1, specular 0, diffuse 0.6, reflective 0.3). Include a `comment` like the other scene files.

- [ ] **Step 2: Run** `cargo test --release soft_shadows` — FAIL (struct sizes / missing kernels).

- [ ] **Step 3: Implement** per Interfaces. `from_world` keeps rejecting path tracing and spot lights; the existing "point lights only" error goes away.

- [ ] **Step 4: Run** `cargo test --release` — every existing compare (chapters 7–15) and probe still passes, plus the new ones.

- [ ] **Step 5: Commit** — `ENHANCEMENT: rtc-metal renders area lights like rtc-rust (bonus: soft shadows)`

---

### Task 3: Sampled frames, `Accumulator`, `--samples`

**Files:**
- Modify: `shaders/render.metal` (`render` sampled branch, new `resolve_buffer`, `resolve_texture`)
- Modify: `shaders/probes.metal` (new `probe_shade_point`)
- Modify: `src/gpu.rs` (`Accumulator`, `Gpu::accumulator`, `Gpu::queue`)
- Modify: `src/main.rs` (`--samples`, PPM writing helper, usage text)
- Test: `src/probes.rs`, `src/compare.rs`, `src/gpu.rs`

**Interfaces:**
- Consumes: `GpuSceneInfo.frame/sampling`, `intensity_at`, `lighting`, `Sampler` (Tasks 1–2).
- Produces:
  - `render`: `Sampler s = {pixel = gid.y * width + gid.x, frame = scene.frame, on = scene.sampling, call = 0}`. `on == 1`: ray through `(px + random01(pixel, frame, 0), py + random01(pixel, frame, 1))` — `ray_for_pixel` gains `float jx, float jy` offsets (0.5 for deterministic); `accum[i] = frame == 0 ? c : accum[i] + c`. `on == 0`: writes `c` as today.
  - `kernel void resolve_buffer(device const float4 *accum [[buffer(0)]], device float4 *out [[buffer(1)]], constant uint3 &size_frames [[buffer(2)]], uint2 gid)` → `accum / frames`.
  - `kernel void resolve_texture(device const float4 *accum [[buffer(0)]], texture2d<float, access::write> out [[texture(0)]], constant uint3 &size_frames [[buffer(2)]], uint2 gid)` → `clamp(accum / frames, 0, 1)`, alpha 1.
  - Rust, `src/gpu.rs`:
    ```rust
    pub struct Accumulator { /* accum, shapes, nodes buffers; render + 2 resolve pipelines; width, height, frames: u32 */ }
    impl Gpu {
        pub fn accumulator(&self, scene: &GpuScene) -> Result<Accumulator, String>;
        pub fn queue(&self) -> &ProtocolObject<dyn MTLCommandQueue>;
    }
    impl Accumulator {
        pub fn frames(&self) -> u32;                 // samples accumulated so far
        pub fn reset(&mut self);                     // frames = 0; next frame overwrites
        pub fn encode_frame(&mut self, encoder: &Encoder, info: &GpuSceneInfo); // sampling = 1, frame = frames; frames += 1
        pub fn encode_resolve_texture(&self, encoder: &Encoder, texture: &ProtocolObject<dyn MTLTexture>);
        pub fn add_frame(&mut self, gpu: &Gpu, info: &GpuSceneInfo) -> Result<(), String>; // submit(encode_frame), waits
        pub fn resolve(&self, gpu: &Gpu) -> Result<Vec<[f32; 4]>, String>; // Err("no samples accumulated") when frames == 0
    }
    ```
    Pipelines are created once in `accumulator()` (the per-call `pipeline()` builds a new pipeline state). Make `Encoder` `pub(crate)`. The spec's `Accumulator::new(&gpu, w, h)` becomes `gpu.accumulator(&scene)` so the scene buffers are uploaded once and kept.
  - Probe `ShadeProbe { material: GpuMaterial, light: GpuLight, point: Float4, eyev: Float4, normalv: Float4, pixel: u32, frame: u32, sampling: u32 }` → `Float4`: `intensity_at` then `lighting` with the same `Sampler` (`call = 0`), via `run_with_scene`.
  - `src/main.rs`: `fn write_ppm(pixels: Vec<[f32; 4]>, width: usize, height: usize, out: &str)` (the current canvas loop, moved). Usage: `rtc-metal [scene.json] [-o out.ppm] [--samples N] [--window]`, with a line that `antialias.n` is ignored on the GPU.

- [ ] **Step 1: Write the failing tests**
  - `src/probes.rs`, `test_bonus_soft_shadows_gpu_sampled_shading_converges`: `World::default_world()` scene; rtc-rust's 2×2 light from `test_bonus_soft_shadows_5`; point `(1.0, -1.0, 2.0)` (penumbra, CPU intensity 0.25), eyev `(0, 0, -1)`, normalv `(0, 0, -1)`, default material; 4096 probes with frames 0..4096, pixel 0, sampling 1; mean of the outputs within 0.01 (each channel) of the single deterministic probe (sampling 0).
  - `src/gpu.rs`, `accumulator_reset_overwrites`: chapter 7 scene (`crate::compare::chapter_7_scene`, make it `pub(crate)`) at 20×10; 3 frames, `reset()`, 1 frame; `resolve` equals a fresh accumulator's 1-frame `resolve` exactly (Review Focus 4).
  - `src/gpu.rs`, `accumulator_without_frames_errors`: `resolve` on a fresh accumulator is `Err("no samples accumulated")` (Review Focus 3).
  - `src/compare.rs`, `test_bonus_soft_shadows_gpu_sampled_render_converges`: `soft_shadows.json` at 200×80; mean absolute error per channel vs CPU `render_parallel` after 4 and after 256 sampled frames; assert `err_256 <= 0.5 * err_4` and `err_256 < THRESHOLD`. Calibrate `THRESHOLD` once: run, take the measured `err_256`, set `THRESHOLD` to it ×1.5 rounded up to 2 significant digits, and document the measured value in a comment.
  - `src/compare.rs`, `test_chap_7_gpu_sampled_render_stays_close`: chapter 7 scene (point light) at 160×80, 64 sampled frames vs the deterministic GPU frame; mean absolute error < 0.01 (Review Focus 2).

- [ ] **Step 2: Run** `cargo test --release` — new tests FAIL.

- [ ] **Step 3: Implement** kernels, `Accumulator`, `ShadeProbe`, then `main.rs`: with `--samples N` (N ≥ 1, else usage), build `gpu.accumulator(&scene)`, loop `add_frame` N times printing `\r{n}/{N} samples` to stderr at most once a second, then `resolve` → `write_ppm`, and print total ms. Without `--samples`/`--window`, keep the current deterministic path unchanged.

- [ ] **Step 4: Run** `cargo test --release` and `cargo clippy --release --all-targets` — all pass, 0 warnings. Then `cargo run --release -- ../rtc-rust/scenes/dice-light-area.json --samples 64 -o <scratchpad>/dice.ppm`, convert with `sips -s format png`, look at it: soft shadows, smooth edges.

- [ ] **Step 5: Commit** — `ENHANCEMENT: rtc-metal accumulates sampled frames headless (--samples)`

---

### Task 4: Orbit camera maths

**Files:**
- Create: `src/orbit.rs`
- Modify: `src/lib.rs` (`pub mod orbit;`)

**Interfaces:**
- Produces:
  ```rust
  pub struct Orbit { to: Tuple, up: Tuple, radius: f64, yaw: f64, pitch: f64 }
  impl Orbit {
      pub fn new(from: Tuple, to: Tuple, up: Tuple) -> Self;
      pub fn from(&self) -> Tuple;
      pub fn drag(&mut self, dx: f64, dy: f64);   // points; 0.01 rad per point; pitch clamped to ±89°
      pub fn zoom(&mut self, dy: f64);            // radius *= exp(-0.1 * dy), kept ≥ 1e-3
      pub fn transform(&self) -> Matrix4;         // view_transform(self.from(), to, up)
  }
  ```
  Yaw rotates around world +y, pitch is the elevation above the xz-plane through `to` (every scene in the repo has `up = (0, 1, 0)`; document that assumption on the struct).

- [ ] **Step 1: Write the failing tests** (in `src/orbit.rs`)
  - `new_round_trips_from`: `Orbit::new((-3, 1, 2.5), (0, 0.5, 0), (0, 1, 0)).from()` approx-equals `(-3, 1, 2.5)`; and `transform()` approx-equals `view_transform` of the original vectors.
  - `drag_zero_is_identity`: `drag(0, 0)` leaves `from()` unchanged.
  - `pitch_is_clamped`: `drag(0, 1e6)` then `drag(0, -1e6)` keep pitch within ±89° (check the angle of `from() - to` against the xz-plane).
  - `zoom_keeps_positive_radius`: `zoom(1e6)` leaves `(from() - to).magnitude() >= 1e-3`.
  - `yaw_quarter_turn`: from `(0, 0, -5)` to the origin, `drag(PI / 2 / 0.01, 0)` puts `from()` at `(±5, 0, 0)` (sign chosen so dragging right moves the camera right around the scene; state it in a comment).

- [ ] **Step 2: Run** `cargo test --release orbit` — FAIL.
- [ ] **Step 3: Implement** `Orbit`.
- [ ] **Step 4: Run** `cargo test --release` — pass.
- [ ] **Step 5: Commit** — `ENHANCEMENT: rtc-metal orbit camera maths`

---

### Task 5: Live window (`--window`)

**Files:**
- Modify: `Cargo.toml` — add `objc2-app-kit = "0.3.2"`, `objc2-metal-kit = { version = "0.3.2", features = ["MTKView"] }` (plus any feature the compiler asks for), `objc2-quartz-core = "0.3.2"`; enable the extra `objc2-metal` features the window code needs (`MTLTexture`, `MTLDrawable`).
- Create: `src/window.rs`
- Modify: `src/lib.rs` (`pub mod window;`), `src/main.rs` (`--window`)

**Interfaces:**
- Consumes: `Gpu`, `Gpu::queue`, `Accumulator` (Task 3), `GpuScene::set_camera` (Task 2), `Orbit` (Task 4), `write_ppm` (Task 3 — move it to `lib` as `rtc_metal::ppm::write_ppm` if `window.rs` needs it; otherwise pass a `Box<dyn FnOnce(Vec<[f32;4]>)>` "on close" callback from `main.rs`).
- Produces: `pub fn run(gpu: Gpu, scene: GpuScene, camera: Camera, orbit: Orbit, samples: u32, title: String, on_close: Box<dyn FnOnce(Vec<[f32; 4]>)>) -> !`
- Shared state: one `Rc<RefCell<Renderer>>` with `gpu, scene, acc, camera, orbit, target, last_frame: Instant`, held in the ivars of both classes below (main thread only).
- `RtcView` (`define_class!`, superclass `MTKView`): overrides `mouseDown:` (remember location), `mouseDragged:` (`orbit.drag(dx, dy)` with the event's delta in points), `scrollWheel:` (`orbit.zoom(scrollingDeltaY)`); each then `camera.set_transform(&orbit.transform())`, `scene.set_camera(&camera)`, `acc.reset()`, `setPaused(false)`, `setEnableSetNeedsDisplay(false)`. `acceptsFirstResponder` → YES.
- `RtcDelegate` (`define_class!`, `NSObject` + `MTKViewDelegate` + `NSApplicationDelegate`):
  - `drawInMTKView:` — if `acc.frames() >= target`: `setPaused(true)`, `setEnableSetNeedsDisplay(true)`, return. Get `currentDrawable` (none → return without counting). One command buffer: compute encoder → `encode_frame(&scene.info)` → `encode_resolve_texture(drawable.texture())` → end → `presentDrawable` → `commit` (no wait). Update the window title `"{title} — {n}/{target} samples, {ms:.1} ms/frame"`.
  - `drawableSizeWillChange` — nothing (drawable size is fixed).
  - `applicationShouldTerminateAfterLastWindowClosed:` → YES.
  - `applicationWillTerminate:` — `acc.resolve(&gpu)`: `Ok` → `on_close(pixels)`; `Err` → print it to stderr (Review Focus 3).
- View setup: `MTKView` with `gpu`'s device, `colorPixelFormat = BGRA8Unorm`, `framebufferOnly = false`, `autoResizeDrawable = false`, `drawableSize = (width, height)`, `preferredFramesPerSecond = 60`. Window content size in points: `(width, height) / backingScaleFactor`, scaled down uniformly to fit `NSScreen::mainScreen().visibleFrame()` if larger. `NSApplication` activation policy `Regular`, `activateIgnoringOtherApps(true)`, `makeKeyAndOrderFront`, `run()`.
- `main.rs`: `--window` loads the scene with `SceneDescription::from_path` + `build(dir)` (keeps `camera.from/to/up` for `Orbit::new`), target = `--samples` or 256, `on_close` writes the PPM to `-o` / the default name.

- [ ] **Step 1: Write** `src/window.rs` per Interfaces. (No automated test: AppKit needs a running app on the main thread. The logic it drives is covered by Tasks 3–4.)
- [ ] **Step 2: Build** `cargo build --release` and `cargo clippy --release --all-targets` — 0 warnings; `cargo test --release` still green.
- [ ] **Step 3: Manual check — book scene.** `cargo run --release -- scenes/soft_shadows.json --window -o <scratchpad>/soft.ppm`: window opens with the image refining; title counts to 256 and stops (Activity Monitor GPU idle); dragging right orbits the camera right and restarts accumulation; scrolling zooms; closing writes `soft.ppm` (check with `sips`).
- [ ] **Step 4: Manual check — large scenes.** `scenes/chapter15.json --window` (teapots: slower frames, still responsive) and the local 3456×2234 `scenes/chapter10.json --window` (window fits the screen, Review Focus 5). Cmd-Q before the first frame prints "no samples accumulated" and writes nothing.
- [ ] **Step 5: Commit** — `ENHANCEMENT: rtc-metal live progressive window with orbit camera (--window)`
