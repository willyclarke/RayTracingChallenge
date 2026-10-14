# rtc-metal milestone 13: multiple lights Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Scenes with several lights in rtc-rust (the reference) and
rtc-metal, each light's full Phong contribution (ambient included) summed
at every hit.

**Architecture:** rtc-rust's `World.light: Option<Light>` becomes
`lights: Vec<Light>`, `shade_hit` sums over it, and scene files gain an
optional `"lights": [...]`. rtc-metal's `SceneInfo` carries a fixed
`lights[8]` plus `light_count`; `color_at` loops over them, light `k ≥ 1`
drawing sampled random numbers from its own dimension range.

**Tech Stack:** Rust 2024 (rtc-rust: rayon, serde), objc2-metal, MSL
compiled at runtime.

**Spec:** `RUST/rtc-metal/docs/superpowers/specs/2026-10-14-metal-13-multiple-lights-design.md`

## Global Constraints

- Branch `metal/13-multiple-lights` (exists, spec committed as `c74d475`);
  commit per task; never push.
- Commit messages `ENHANCEMENT:`/`TEST:`/`DOCS:` prefix; rtc-metal subjects
  end `(milestone 13)`; rtc-rust subjects start `rtc-rust`; all end with
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- New test names `test_bonus_multiple_lights_*` in both crates.
- **Both crates build and pass at every commit** (rtc-metal path-depends on
  rtc-rust). Baselines: rtc-rust `cargo test --release --lib` 344 pass
  (7 ignored), `cargo test --release --doc` 37 pass, examples build;
  rtc-metal `cargo test --release --lib` 116 pass. Clippy clean in both
  (`cargo clippy --release --all-targets`); `cargo fmt`.
- Single-light images stay bit-identical: sums start from black.
- `MAX_LIGHTS = 8`; `LIGHTS_DIM = 0x2000_0000`; light `k ≥ 1` slot dims
  `LIGHTS_DIM + (k − 1)·0x0100_0000 + 6·call + slot`; light 0 keeps
  `2 + 6·call + slot`.
- Error texts, verbatim: `a scene needs exactly one of "light" and "lights"`,
  `"lights" is empty`, `scene has no light`, `more than 8 lights`.
- `render.metal` keeps exactly one `color_at` call site.

## Review Focus

1. **Two area lights, sampled** — light 1 must not reuse light 0's random
   dims (correlated shadows: more noise, no bias, invisible to mean-error
   tests). Pinned by a probe of `sample()`'s dimension in Task 4.
2. **Both `"light"` and `"lights"`, or `"lights": []`** — a clear load
   error, not a silent pick. Test in Task 2.
3. **A 9th light on the GPU** — error, not truncation. Test in Task 3.
4. **Shadowed from one light, lit by the other** — the second light's
   diffuse must survive the first's shadow. Test in Task 1 (CPU) and by
   the pixel-exact compare in Task 4 (GPU).
5. **Multi-light on the Metal traversal** — the light loop inlines
   `intensity_at` more; every image test runs both traversals.

---

### Task 1: rtc-rust `World` holds a list of lights

**Files:**
- Modify: `rtc-rust/src/world.rs` (field, `set_light`, `add_light`,
  `shade_hit_gi`, `Display`, tests), `rtc-rust/src/light.rs` (tests using
  `w.light`), `rtc-rust/src/scene.rs` (`build` sets the light; test at
  ~788), `rtc-rust/examples/{wireframe,bvh_stress,hexagons}.rs`,
  `rtc-rust/CLAUDE.md` (the `world` line: "the shapes and lights")
- Modify: every rtc-metal `world.light` use (`src/scene.rs`, `accel.rs`,
  `compare.rs`, `gpu.rs`, `probes.rs`) — mechanical; `GpuScene::from_world`
  takes `world.lights.first()` until Task 3

**Interfaces:**
- Produces: `pub lights: Vec<Light>` on `World`;
  `pub fn set_light(&mut self, light: Light)` (replaces all);
  `pub fn add_light(&mut self, light: Light)` (appends).

- [ ] **Step 1: Write the failing tests** in `world.rs`'s test module:
  - `test_bonus_multiple_lights_shade_hit_sums_lights`: the default world
    (`World::default_world()` or whatever the chapter 7 tests use) with
    light A `point_light((-10,10,-10), WHITE)` and B
    `point_light((10,10,-10), color(0.5,0.25,0.1))`; for the first hit of
    ray `(0,0,-5)→(0,0,1)`, `shade_hit` with `[A, B]` (via `set_light(A)`,
    `add_light(B)`) `approx_eq`s `shade_hit` with only A plus `shade_hit`
    with only B.
  - `test_bonus_multiple_lights_shadow_from_one_light`: a sphere at the
    origin, a plane `y = -1`; light A straight above at `(0,10,0)`, light
    B at `(10,10,0)`. The floor point `(0,-1,0)` is shadowed from A
    (`is_shadowed(A.position, p)`) but its `shade_hit` (ray from
    `(0,5,-5)` aimed at it) is brighter than the ambient-only colour,
    and equals A-alone + B-alone.
- [ ] **Step 2: Run** — `cargo test --release --lib multiple_lights` in
  rtc-rust → compile error (`add_light` / `lights` missing).
- [ ] **Step 3: Implement** the field, both methods, the summing
  `shade_hit_gi` (start from `BLACK`, add each
  `light.lighting_with_ambient(…, light.intensity_at_time(comps.over_point, self, comps.time), ambient_scale)`),
  `Display` listing every light (`light: none` when empty), and migrate
  every `.light` use in both crates as the spec lists.
- [ ] **Step 4: Run** — rtc-rust: `cargo test --release --lib` → 346 pass,
  `--doc` 37, `cargo build --release --examples` ok; rtc-metal:
  `cargo test --release --lib` → 116 pass; clippy clean in both.
- [ ] **Step 5: Commit** — `ENHANCEMENT: rtc-rust World sums several lights`

### Task 2: rtc-rust scene files accept `"lights"`

**Files:**
- Modify: `rtc-rust/src/scene.rs` (`SceneDescription.light` →
  `Option<LightDescription>`, new `lights: Vec<LightDescription>`, `build`,
  tests), `rtc-rust/src/help.rs` (`LIGHTS` and overview text: "one light,
  or a "lights" list"), `rtc-rust/README.md` (scene format section)

**Interfaces:**
- Consumes: Task 1's `set_light`/`add_light`.
- Produces: JSON `"lights": [<light>, …]` building `World.lights` in order.

- [ ] **Step 1: Write the failing tests** in `scene.rs`:
  - `test_bonus_multiple_lights_scene_lights_array`: a scene with
    `"lights": [point (-10,10,-10) white, area light]` builds
    `world.lights.len() == 2`, `lights[0].position ≈ (-10,10,-10)`, and the
    second is the area light (`usteps`/`vsteps` as given).
  - `test_bonus_multiple_lights_scene_light_key_errors`: both keys,
    neither key, and `"lights": []` each make `build` (or `from_json`)
    return an error whose text contains, respectively,
    `exactly one of "light" and "lights"` (both cases) and
    `"lights" is empty`.
- [ ] **Step 2: Run** — `cargo test --release --lib multiple_lights_scene` →
  fail (unknown field `lights`).
- [ ] **Step 3: Implement** the two fields (`#[serde(default)]`) and the
  checks in `build`; update `help.rs` and `README.md` with a `"lights"`
  example.
- [ ] **Step 4: Run** — rtc-rust lib → 348 pass, doc 37; `cargo run
  --release -- help lights` shows the array; rtc-metal lib → 116 pass;
  clippy clean.
- [ ] **Step 5: Commit** — `ENHANCEMENT: rtc-rust scene files take a "lights" list`

### Task 3: `GpuSceneInfo` carries up to 8 lights

**Files:**
- Modify: `rtc-metal/src/types.rs` (`MAX_LIGHTS`, `GpuSceneInfo`, size
  assert 784), `src/scene.rs` (`from_world`, tests; `info.light` →
  `info.lights[0]`), `shaders/rt.metal` (`SceneInfo` mirror,
  `static_assert` 784; `color_at` uses `scene.lights[0]` for now)

**Interfaces:**
- Produces: `pub const MAX_LIGHTS: usize = 8;`
  `GpuSceneInfo { camera, lights: [GpuLight; MAX_LIGHTS], counts, frame,
  sampling, gi_samples, gi_depth, light_count: u32, _pad: u32 }`;
  MSL `constant uint MAX_LIGHTS = 8;` and `SceneInfo.lights`,
  `SceneInfo.light_count`.

- [ ] **Step 1: Write the failing test** in `scene.rs`:
  `test_bonus_multiple_lights_reach_scene_info` — a world with 3 point
  lights at x = 1, 2, 3 → `info.light_count == 3`,
  `info.lights[k] == gpu_light(&world.lights[k])` for k < 3, slots 3..8
  zeroed; the same world with 9 lights → `from_world` error containing
  `more than 8 lights`; no lights → `scene has no light`.
- [ ] **Step 2: Run** — `cargo test --release multiple_lights_reach` →
  compile error (no `lights` field).
- [ ] **Step 3: Implement** the layout on both sides and `from_world`.
- [ ] **Step 4: Run** — rtc-metal lib → 117 pass; clippy clean.
- [ ] **Step 5: Commit** — `ENHANCEMENT: rtc-metal GpuSceneInfo carries up to 8 lights (milestone 13)`

### Task 4: The GPU sums the lights

**Files:**
- Modify: `shaders/rt.metal` (`LIGHTS_DIM`, `Sampler.light`, `sample()`,
  `color_at` light loop), `shaders/probes.metal` (`probe_sample`),
  `src/probes.rs` (sample-dim probe test), `src/compare.rs` (three tests)
- Create: `scenes/multiple_lights.json`, `scenes/multiple_lights_area.json`

**Interfaces:**
- Consumes: Task 3's `SceneInfo.lights`/`light_count`; the existing
  `probe_random(pixel, frame, dim)`; `check_matches_cpu`, `scene_file`,
  `mean_abs_error`, `Accumulator`.
- Produces: `constant uint LIGHTS_DIM = 0x20000000u;`
  `struct Sampler { uint pixel; uint frame; uint on; uint call; uint light; }`;
  kernel `probe_sample(SampleProbe{pixel, frame, call, light, slot})` →
  `float`.

- [ ] **Step 1: Write the scenes.**
  - `multiple_lights.json`: copy `scenes/chapter7.json`'s camera and
    shapes; `"lights"`: red `(−10,10,−10)` intensity `[0.6,0.15,0.15]`,
    blue `(10,10,−10)` intensity `[0.15,0.15,0.6]`.
  - `multiple_lights_area.json`: copy `scenes/soft_shadows.json`; keep its
    area light first in `"lights"`, add a point light mirrored across x
    from it at intensity `[0.4,0.4,0.4]`.
- [ ] **Step 2: Write the failing tests**:
  - `probes.rs` `test_bonus_multiple_lights_gpu_sample_dims`: for
    `(pixel, frame, call, light, slot)` in a handful of cases with
    `light ∈ {0, 1, 7}`, `call ∈ {0, 3}`, `slot ∈ {0, 5}`, the
    `probe_sample` value equals `probe_random` at dim
    `2 + 6·call + slot` (light 0) or
    `0x2000_0000 + (light−1)·0x0100_0000 + 6·call + slot`.
  - `compare.rs`: `test_bonus_multiple_lights_gpu_render_matches_cpu`
    (`multiple_lights.json` at 200×100, `check_matches_cpu`);
    `test_bonus_multiple_lights_area_gpu_render_matches_cpu`
    (`multiple_lights_area.json` at 200×80, `check_matches_cpu`);
    `test_bonus_multiple_lights_area_gpu_sampled_render_converges`, shaped
    like `test_bonus_soft_shadows_gpu_sampled_render_converges`,
    `THRESHOLD` starting at 0.0.
- [ ] **Step 3: Run** — `cargo test --release multiple_lights` → probe
  fails (no kernel), compares fail (GPU renders only light 0).
- [ ] **Step 4: Implement** the sampler and the loop: per shaded hit,
  `for k < light_count { sampler.light = k; color += p.weight *
  lighting(m, scene.lights[k], …, intensity_at(scene.lights[k], …), …); }`
  then one `sampler.call++`.
- [ ] **Step 5: Calibrate** the sampled test's `THRESHOLD` = 1.5 ×
  `err_256`, 2 significant digits, measurements in the comment.
- [ ] **Step 6: Run** — rtc-metal lib → 121 pass in **several** parallel
  runs (watch for GPU watchdog kills); clippy clean;
  `grep -c "color_at(" shaders/render.metal` → 1.
- [ ] **Step 7: Commit** — `ENHANCEMENT: rtc-metal sums every light (milestone 13)`
  (shaders, probe) and `TEST: rtc-metal multiple lights match the CPU (milestone 13)`
  (scenes, compare tests).

### Task 5: Two-light Cornell box, path traced

**Files:**
- Create: `scenes/cornell_two_lights.json`
- Modify: `src/compare.rs` (one test)

**Interfaces:**
- Consumes: `check_converges_to_cpu(world, reference, camera, max_abs, max_signed)`.

- [ ] **Step 1: Write the scene** — `chapter17_path_tracing.json`'s box
  walls, blocks and settings at camera 150×150 (not the 3456×2234 working
  size); `"lights"`: the inside light `(0,1.8,−0.6)` white, plus
  `(0.6,1.5,−2.5)` intensity `[0.5,0.5,0.5]`; comment explains the second
  light shines in through the opening.
- [ ] **Step 2: Write the test**
  `test_bonus_multiple_lights_path_tracing_gpu_render_matches_cpu`: the
  scene at 150×149 (diagonal-edge ties, milestone 12), reference world
  with `set_path_tracing(256, 2)`, thresholds 0.0.
- [ ] **Step 3: Run** — fails on thresholds, printing the errors. Bias
  check: signed error against the 256-path reference must be ≲ 2e-4
  (milestone 12's level); if not, diff-map before calibrating.
- [ ] **Step 4: Calibrate** 1.5×, 2 significant digits, measurements in
  the comment.
- [ ] **Step 5: Run** — rtc-metal lib → 122 pass in several parallel
  runs; clippy clean.
- [ ] **Step 6: Check by eye** — `cargo run --release --
  scenes/cornell_two_lights.json --samples 256 -o …` (second light's
  pool on the floor near the opening, second set of shadows) and
  `--window` launches cleanly.
- [ ] **Step 7: Commit** — `TEST: rtc-metal two-light Cornell box matches the CPU (milestone 13)`
