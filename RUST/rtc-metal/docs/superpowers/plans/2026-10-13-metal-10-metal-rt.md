# rtc-metal milestone 10: Metal ray tracing API Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A second GPU traversal, Metal's `intersection_query` over one
primitive acceleration structure, selectable per render next to the BVH
walk, with both matching the CPU.

**Architecture:** `GpuScene` gains world-space triangle vertices, box
bounds and a primitive-id → shape-index map; `accel.rs` builds the
structure; `hit_shapes`/`is_shadowed` branch on `counts.traversal`.
Planes stay in the unbounded loop and `refractive_indices` stays on the
BVH walk.

**Tech Stack:** Rust 2024, objc2-metal 0.3.2 (default features already
include the acceleration-structure bindings), MSL compiled at runtime,
`<metal_raytracing>`.

**Spec:** `RUST/rtc-metal/docs/superpowers/specs/2026-10-13-metal-10-metal-rt-design.md`

## Global Constraints

- Branch `metal/10-metal-rt`; commit per task; never push.
- Commit messages `ENHANCEMENT:`/`TEST:`/`DOCS:` prefix, ending with
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- Test names `test_chap_N_*` after the chapter of the scene they use.
- `TRAVERSAL_BVH = 0`, `TRAVERSAL_METAL = 1`; default traversal stays BVH.
- Run tests from `RUST/rtc-metal` with `cargo test --release`; clippy clean
  (`cargo clippy --release --all-targets`).
- Every Rust `#[repr(C)]` struct keeps its `static_assert(sizeof…)` twin in
  MSL in step.

**Plan-level refinement of the spec (ruled while planning):** the spec puts
`traversal` in `SceneInfo`, but probes receive only `SceneCounts`, which
`hit_shapes` and `is_shadowed` already take. `traversal` therefore lives in
`SceneCounts` next to `triangle_count`, reachable from renders as
`scene.counts.traversal`. `SceneCounts` grows 12 → 20 bytes; `SceneInfo`
stays 208 (the counts sit in its tail padding: 176 + 20 + 8 = 204 → 208).

## Review Focus

1. A ray that starts inside a sphere/CSG box (refraction, shadow rays from
   inside a die's pip): the query must still report the box candidate and
   our code must commit the exit root. Pinned by Task 4's
   `test_chap_11_metal_ray_from_inside_sphere`; shadow rays leaving a
   die's surface inside its CSG box are covered by the chapter 16
   agreement test.
2. A box candidate whose first root lies beyond the committed distance must
   not commit (committing past it is undefined). Pinned by Task 4's
   "nearer triangle in front of a sphere" probe.
3. A scene with planes only, or with nothing at all, on `METAL`: the
   placeholder structure is bound and the query skipped. Pinned by Task 2's
   planes-only build and Task 4's planes-only trace.
4. The window's `Accumulator` must bind the structure too, or `--window
   --traversal metal` reads an unbound slot. Pinned by Task 5's sampled
   test on `METAL`.
5. Triangle hits exactly on shared edges may land on the neighbour triangle;
   budgeted (1 in 1000 rays) in Task 4, never loosened elsewhere.

---

### Task 1: Acceleration-structure inputs in `GpuScene`

**Files:**
- Modify: `src/bvh.rs` (keep each bounded entry's box)
- Modify: `src/scene.rs` (`Leaf` world vertices, `GpuScene` fields)
- Test: `src/scene.rs` (new `#[cfg(test)] mod tests` at the end)

**Interfaces:**
- Produces on `GpuScene`:
  `pub rt_vertices: Vec<[f32; 3]>` (3 per triangle, world space, in
  `rt_ids` triangle order),
  `pub rt_boxes: Vec<([f32; 3], [f32; 3])>` (in `rt_ids` box order),
  `pub rt_ids: Vec<u32>` (shape indices of every bounded triangle in
  `shapes` order, then of every other bounded entry in `shapes` order),
  `pub rt_triangle_count: u32`.
- Produces on `Bvh`: `pub bounds: Vec<([f32; 3], [f32; 3])>`, one per
  bounded entry, i.e. `bounds[k]` is the box of
  `shapes[unbounded_count + k]`.
- Produces on `bvh::Leaf`: `pub vertices: Option<[[f32; 3]; 3]>`, world
  triangle vertices; `Some` exactly for triangle leaves.

- [ ] **Step 1: Write the failing test `test_chap_15_rt_inputs`**

World: a plane, a sphere at `translation(3, 0, 0)`, a group with transform
`translation(0, 1, 0) * scaling(2, 2, 2)` holding one `Triangle` with
points (0,1,0), (-1,0,0), (1,0,0), and one CSG (union of two spheres);
`world.build_bounds()`. Assert:
- `rt_triangle_count == 1`, `rt_ids.len() == 3` (triangle, sphere, CSG);
- `shapes[rt_ids[0]].kind == SHAPE_TRIANGLE`, the next two are
  `SHAPE_SPHERE` and `SHAPE_CSG` in `shapes` order;
- `rt_vertices == [[0,3,0], [-2,1,0], [2,1,0]]` exactly;
- `rt_boxes[i]` equals `bvh.bounds` for that entry (check the sphere's
  box contains (2,-1,-1)..(4,1,1));
- no `rt_ids` entry is `< info.counts.unbounded_count` (the plane) or
  `>= info.counts.shape_count` (CSG leaves).

- [ ] **Step 2: Run it** — `cargo test --release test_chap_15_rt_inputs`
  Expected: compile error, no field `rt_vertices`.

- [ ] **Step 3: Implement**
  - `flatten`'s leaf arm fills `Leaf.vertices` for `Triangle` and
    `TriangleUV` (downcast as `gpu_shape` does): `full * p1/p2/p3` in f64,
    then cast. Every other `Leaf` construction gets `vertices: None`.
  - `bvh::build` records `bounds` and keeps `vertices` alongside (a
    `pub vertices: Vec<Option<[[f32; 3]; 3]>>` parallel to `bounds` is
    fine).
  - `from_world` splits the bounded range `unbounded_count..shape_count`
    into triangles then the rest to fill the four fields.

- [ ] **Step 4: Run** — `cargo test --release` → all pass (83 + 1).

- [ ] **Step 5: Commit** — `ENHANCEMENT: rtc-metal collects
  acceleration-structure inputs (milestone 10)`

### Task 2: `accel::build`

**Files:**
- Create: `src/accel.rs`; Modify: `src/lib.rs` (`pub mod accel;`)
- Test: `src/accel.rs` (`#[cfg(test)] mod tests`)

**Interfaces:**
- Consumes: Task 1's `GpuScene` fields.
- Produces:
  ```rust
  pub struct Accel {
      pub structure: Retained<ProtocolObject<dyn MTLAccelerationStructure>>,
      pub rt_ids: Retained<ProtocolObject<dyn MTLBuffer>>,
  }
  pub fn build(gpu: &Gpu, scene: &GpuScene) -> Result<Accel, String>;
  ```
  `Gpu` gains `pub(crate) fn input_buffer`/`output_buffer` visibility if
  needed (they exist, private).

- [ ] **Step 1: Write failing tests** `test_chap_15_accel_builds_triangles`,
  `test_chap_12_accel_builds_boxes`, `test_chap_16_accel_builds_mixed`
  (chapter16.json at 40×25, which has a plane, CSGs; add one triangle) and
  `test_chap_9_accel_builds_planes_only` (one plane). Each asserts
  `build` is `Ok` and `structure.size() > 0`.

- [ ] **Step 2: Run** — Expected: compile error, no module `accel`.

- [ ] **Step 3: Implement `build`**
  - Triangle geometry only when `rt_triangle_count > 0`: vertex buffer from
    `rt_vertices` (stride 12, `MTLAttributeFormat::Float3`), no index
    buffer, `triangleCount = rt_triangle_count`, `opaque = true`.
  - Box geometry only when there are boxes: buffer of
    `MTLAxisAlignedBoundingBox` from `rt_boxes`, stride 24,
    `boundingBoxCount`, `opaque = false`.
  - Neither: one box geometry with a single zero-size box at the origin.
  - Sizes → `newAccelerationStructureWithSize`, private scratch buffer of
    `buildScratchBufferSize`, one acceleration-structure encoder build,
    commit, wait, check the command buffer status like `Gpu::submit`.
  - `rt_ids` buffer: the `rt_ids` Vec, or one zero word when empty.

- [ ] **Step 4: Run** — `cargo test --release` → all pass.

- [ ] **Step 5: Commit** — `ENHANCEMENT: rtc-metal builds a primitive
  acceleration structure (milestone 10)`

### Task 3: Bind the structure everywhere; shaders accept it

**Files:**
- Modify: `src/types.rs` (`GpuSceneCounts`), `src/scene.rs` (fill
  counts), `src/gpu.rs` (bindings), `shaders/rt.metal`,
  `shaders/render.metal`, `shaders/probes.metal`, every `GpuSceneCounts {…}`
  literal in `src/probes.rs`.

**Interfaces:**
- Produces: `GpuSceneCounts { shape_count, unbounded_count, node_count,
  triangle_count: u32, traversal: u32 }` and MSL twin (20 bytes);
  `pub const TRAVERSAL_BVH: u32 = 0; pub const TRAVERSAL_METAL: u32 = 1;`
  in `types.rs`, same names as MSL `constant uint`s.
- MSL: `hit_shapes(Ray, shapes, nodes, tables, rt_ids, accel, counts,
  thread Hit &)` and `is_shadowed(light_position, point, shapes, nodes,
  tables, rt_ids, accel, counts)`, with `device const uint *rt_ids` and
  `primitive_acceleration_structure accel`, threaded through
  `light_reaches`, `cell_lit`, `intensity_at`, `color_at`, `trace`.
- Bindings: `render` kernel `rt_ids [[buffer(5)]]`, `accel [[buffer(6)]]`;
  probe kernels that reach `hit_shapes`/`is_shadowed` take
  `rt_ids [[buffer(6)]]`, `accel [[buffer(7)]]`. `run_with_scene` binds
  6 and 7 for every scene probe (unused bindings are harmless).
- `scene_buffers` returns a struct `SceneBuffers { shapes, nodes, tables,
  accel: Accel }` instead of the tuple; `Accumulator` holds one.

- [ ] **Step 1: Write the failing test `test_chap_7_traversal_flag_reaches_shader`**
  in `probes.rs`: a probe kernel `probe_counts` (out = `counts.traversal`,
  `counts.triangle_count`) run through `run_with_scene` on chapter 7 with
  `scene.info.counts.traversal = TRAVERSAL_METAL` returns `[1, 0]`.
- [ ] **Step 2: Run** — Expected: compile error (no field `traversal`).
- [ ] **Step 3: Implement** the struct change, the argument threading and
  the bindings (`setAccelerationStructure_atBufferIndex`). The shader
  ignores `traversal` so far.
- [ ] **Step 4: Run** — `cargo test --release` → all pass; clippy clean.
- [ ] **Step 5: Commit** — `ENHANCEMENT: rtc-metal binds the acceleration
  structure to every scene kernel (milestone 10)`

### Task 4: `TRAVERSAL_METAL` in `hit_shapes` and `is_shadowed`

**Files:**
- Modify: `shaders/rt.metal` (`#include <metal_raytracing>`, both
  functions), `src/compare.rs` (`scene_file` → `pub(crate)`)
- Test: `src/probes.rs`

**Interfaces:**
- Consumes: Task 3's signatures and `counts.traversal`.
- Produces: helper `fn traces_agree(world: &World, camera: &Camera,
  light: Tuple) -> Result<(), String>` in `probes.rs`: one `probe_trace`
  per camera pixel on `TRAVERSAL_BVH` and on `TRAVERSAL_METAL` from the
  same `GpuScene`; requires identical hit/miss and `shadowed` everywhere,
  identical `index` and `root` except where both hits are triangles
  (counted; Err if count > rays / 1000), `|t_m - t_b| <= 1e-4 * max(1, t_b)`.

- [ ] **Step 1: Write failing tests**
  - `test_chap_15_metal_traversal_matches_bvh`: `scenes/chapter15.json`
    at 160×80 via `traces_agree`.
  - `test_chap_16_metal_traversal_matches_bvh`: `scenes/chapter16.json`
    at 160×100.
  - `test_chap_12_metal_traversal_matches_bvh` and
    `test_chap_13_metal_traversal_matches_bvh`: their JSON scenes at
    160×80 (glass: rays starting inside boxes).
  - `test_chap_9_metal_traversal_planes_only`: one plane, `traces_agree`.
  - `test_chap_8_metal_shadow_occluders`: light at (0,10,0), point
    (0,0,0) via `probe_trace` rays straight down onto a floor plane:
    (a) a sphere at (0,5,0) → shadowed; (b) a triangle spanning
    y = 5 above the point → shadowed; (c) a sphere at (0,12,0), beyond
    the light → not shadowed. All on `TRAVERSAL_METAL`.
  - `test_chap_15_metal_nearer_triangle_wins`: a triangle at z = -1 in
    front of a sphere at the origin, ray from (0,0,-5) toward +z on
    `METAL`: `index` is the triangle's, `t == 4` within 1e-5.
  - `test_chap_11_metal_ray_from_inside_sphere`: unit sphere, ray from the
    origin toward +z on `METAL`: hit, `t == 1` within 1e-5, `root == 1`.
- [ ] **Step 2: Watch them fail.** Until Step 3 the shader ignores
  `traversal`, so the agreement tests would pass trivially. Add a stub
  first: `TRAVERSAL_METAL` in `hit_shapes` returns `false` (a miss) and in
  `is_shadowed` returns `false`. Run `cargo test --release metal_`.
  Expected: every Step 1 test FAILS (the stub skips the plane loop too, so
  the planes-only test fails as well). Replace the stub in Step 3.
- [ ] **Step 3: Implement** as the spec's "hit_shapes, TRAVERSAL_METAL"
  and "is_shadowed, TRAVERSAL_METAL" sections: unbounded loop with a
  `SceneCounts` copy whose `node_count = 0`, then the query when
  `counts.node_count > 0`. Box candidates: `i = rt_ids[triangle_count +
  q.get_candidate_primitive_id()]`, CSG via `intersect_csg`. Commit only
  `0 <= t < q.get_committed_distance()`.
- [ ] **Step 4: Run** — `cargo test --release` → all pass. If a chapter
  15 agreement exceeds the edge budget, debug (diff map of disagreeing
  rays) before touching the budget.
- [ ] **Step 5: Commit** — `ENHANCEMENT: rtc-metal traverses with Metal
  intersection queries (milestone 10)`

### Task 5: Whole images on both traversals, CLI, timing

**Files:**
- Modify: `src/compare.rs`, `src/main.rs`, `src/window.rs` (if it builds
  its own `GpuScene` info), spec's Timing table.

**Interfaces:**
- Consumes: everything above.
- Produces: CLI `--traversal bvh|metal` (sets
  `scene.info.counts.traversal` before any render or accumulator);
  `GPU render: {ms} ms ({bvh|metal})`.

- [ ] **Step 1: Write the failing change** — `check_matches_cpu` renders on
  `TRAVERSAL_BVH` and `TRAVERSAL_METAL` and checks each against the CPU
  with the existing budget, printing one line per traversal; add
  `test_chap_7_gpu_sampled_metal_matches_bvh`: 16 sampled frames of chapter
  7 at 160×80 through `Accumulator` on each traversal; the resolved images
  differ by at most 1e-3 in 99.9% of pixels.
- [ ] **Step 2: Run** — `cargo test --release` → Expected: PASS if Task 4
  is right (these tests are coverage, not new behaviour; say so in the
  ledger). Any failure is a Task 4 bug: debug it there.
- [ ] **Step 3: Implement** the `--traversal` flag (usage text included;
  an unknown value prints usage).
- [ ] **Step 4: Verify manually**
  - `cargo run --release -- ../rtc-rust/scenes/dice-light-area.json
    --traversal metal -o /tmp/…/dice-metal.ppm` renders; compare by eye
    with the bvh render.
  - Fill the spec's Timing table: median of 3 `GPU render` lines per cell
    for `scenes/chapter15.json`, `dice-light-area.json` and `cover.json`
    (write "not supported: <error>" if a scene is rejected).
  - The live `--window --traversal metal` check is the user's.
- [ ] **Step 5: Run** `cargo test --release` and clippy → green.
- [ ] **Step 6: Commit** — `TEST: rtc-metal renders every chapter on both
  traversals (milestone 10)`
