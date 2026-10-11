# rtc-metal milestone 7: triangles, OBJ and BVH — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Render the book's chapter 15 triangle and OBJ scenes on the GPU, through a BVH, matching rtc-rust within the existing pixel budget.

**Architecture:** Triangles become a new `SHAPE_TRIANGLE` kind in the existing tagged `Shape` array. `GpuScene` flattening also computes world-space boxes, and a new `bvh.rs` builds a median-split BVH over the bounded leaves, storing the shapes reordered as `[unbounded…, bounded…]`. The shader's three linear loops become *linear over the unbounded shapes, then a stack walk of the BVH*.

**Tech Stack:** Rust 2024 (rustc 1.97), objc2-metal, MSL compiled at runtime, rtc-rust as path dependency and CPU reference.

**Spec:** `RUST/rtc-metal/docs/superpowers/specs/2026-10-11-metal-07-triangles-bvh-design.md`

## Global Constraints

- All commands run from `RUST/rtc-metal`. Tests: `cargo test --release`. Lint: `cargo clippy --release --all-targets` must report 0 warnings; run `cargo fmt` before each commit.
- Rust structs in `types.rs` mirror MSL structs byte for byte. Each size change updates both the Rust `const _: () = assert!(size_of::<…>() == N)` and the MSL `static_assert`.
- `GpuShape` / MSL `Shape` become 320 bytes. `GpuBvhNode` / MSL `BvhNode` is 48 bytes. `GpuSceneInfo` stays 144 bytes.
- BVH: median split on the longest axis of the centroids' box, at most 4 shapes per leaf, depth at most 32 (the shader's traversal stack size, `BVH_STACK = 32`).
- GPU triangle test: reject only `det == 0.0`. The u/v range checks handle the rest.
- Compare budget unchanged: at most `w * h / 2000` pixels off by more than 1e-3.
- Test names follow book chapters: `test_chap_15_…` for triangle/OBJ/BVH tests.
- Branch `metal/07-triangles-obj-bvh`. Commit after each task with the repo's `ENHANCEMENT:`/`TEST:` prefix style and the `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>` trailer. Never push. Merge into develop only when the user asks.
- rtc-rust is not modified.

## Review Focus

1. **Degenerate triangles** (zero area; OBJ models have them at poles): GPU and CPU both report a miss, never NaN. Pinned in Task 1, step 1.
2. **Rays through a shared triangle edge** (f32 cracks): a ray aimed exactly at the shared edge of a two-triangle quad hits one of the two triangles. Pinned in Task 4, step 1.
3. **Scenes with no bounded shapes** (only planes) and empty scenes: render without error, with `node_count == 0`. Pinned in Task 3, step 1.
4. **Transforms that rotate an infinite box** (a rotated plane or uncapped cylinder gives NaN corners): the shape is classed unbounded, not given a NaN box. Pinned in Task 3, step 1.
5. **Many leaves with identical boxes** (all centroids equal): the builder terminates with leaves of at most 4 shapes and depth within 32. Pinned in Task 2, step 1.

---

### Task 1: Triangles on the GPU (still linear loops)

**Files:**
- Modify: `src/types.rs` (`GpuShape`, `SHAPE_TRIANGLE`, size asserts)
- Modify: `src/scene.rs` (`gpu_leaf` accepts `Triangle`, `TriangleUV`)
- Modify: `shaders/rt.metal` (`Shape`, `intersect_triangle`, `triangle_uv`, `local_intersect`, `local_normal_at`, `normal_at`, `color_at`)
- Modify: `shaders/probes.metal` (callers of `normal_at`; new `probe_triangle`)
- Test: `src/probes.rs`

**Interfaces:**
- Produces (Rust): `pub const SHAPE_TRIANGLE: u32 = 5;` and `GpuShape` gaining, after `closed`, `pub p1: Float4, pub e1: Float4, pub e2: Float4, pub n1: Float4, pub n2: Float4, pub n3: Float4` (object space). Non-triangles leave them zeroed. A flat rtc-rust `Triangle` sets `n1 = n2 = n3 = normal`. `TriangleUV` copies `n1..n3`.
- Produces (MSL):
  - `inline uint intersect_triangle(Shape s, Ray r, thread float2 &ts)`: returns 0 or 1 roots.
  - `inline float2 triangle_uv(Shape s, Ray r)`: the Möller–Trumbore (u, v) of the object-space ray `r`.
  - `inline float4 local_normal_at(Shape s, float4 object_point, float2 uv)` and `inline float4 normal_at(Shape s, float4 world_point, float2 uv)`: `uv` is used only by `SHAPE_TRIANGLE`, as `n2*u + n3*v + n1*(1-u-v)`, not normalized before the world transform (as on the CPU).
- In `color_at` and `probe_trace`, `uv = s.kind == SHAPE_TRIANGLE ? triangle_uv(s, ray_transform(s.inverse, ray)) : float2(0)`.

- [ ] **Step 1: Write the failing tests** in `src/probes.rs`

`test_chap_15_gpu_intersect_triangle` uses the existing `check_local_intersect` on `Triangle::default()` (book points (0,1,0), (−1,0,0), (1,0,0)) with the book's rays:
- parallel: origin (0,−1,−2), direction (0,1,0)
- misses past the p1–p3 edge: origin (1,1,−2)
- misses past the p1–p2 edge: origin (−1,1,−2)
- misses past the p2–p3 edge: origin (0,−1,−2)
- hits at t = 2: origin (0,0.5,−2)

All except the parallel one use direction (0,0,1). Add a degenerate `Triangle::new` with p1, p2, p3 collinear at (0,0,0), (1,0,0), (2,0,0), and assert that GPU and CPU both report 0 hits for rays at it (Review Focus 1).

`test_chap_15_gpu_triangle_uv_and_normal` uses a new probe:
- Rust `TriangleProbe { shape: GpuShape, ray: GpuRay }` and `TriangleResult { normal: Float4, uv: [f32; 2], t: f32, count: u32 }`.
- MSL `probe_triangle` sets `out = {normal_at(s, point, uv) flipped toward -ray.direction, uv, t, count}` for a world-space ray.

Cases:
- The book's smooth triangle: p1 (0,1,0), p2 (−1,0,0), p3 (1,0,0); n1 (0,1,0), n2 (−1,0,0), n3 (1,0,0). Ray origin (−0.2,0.3,−2), direction (0,0,1). Expect u = 0.45, v = 0.25 and normal ≈ (−0.5547, 0.83205, 0) after normalization.
- A 7×7 grid of rays at that triangle and at a flat `Triangle`, each under the transform `translation(1,−2,0.5) * rotation_y(0.7) * scaling(1,2,3)`. Compare each with a one-shape rtc-rust `World`:
  - `count`, and `t` with `check_close`, against `world.intersect(&ray)`.
  - `uv` against `hit.u`, `hit.v`.
  - `normal` with `check_tuple` against `prepare_computations(hit, …).normalv`.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cargo test --release --lib chap_15`
Expected: FAIL. `gpu_shape` errors with "only spheres, planes, cubes, cylinders and cones", and `probe_triangle` doesn't exist yet.

- [ ] **Step 3: Implement**
  - Update `types.rs` and the MSL `Shape` (320 bytes).
  - Add the `gpu_leaf` triangle branches. The error message becomes "only spheres, planes, cubes, cylinders, cones and triangles are supported so far". Update the `contains("only spheres, planes")` assert in `scene.rs`, which still matches.
  - Add the MSL functions listed under Interfaces. `local_intersect` gets `case SHAPE_TRIANGLE`.
  - Add the `uv` argument at every `normal_at` call: `color_at`, `probe_normal_at` (pass `float2(0)`), `probe_trace`.

- [ ] **Step 4: Run all tests**

Run: `cargo test --release`
Expected: PASS. That includes the chapter 7–14 compares, now with 320-byte shapes.

- [ ] **Step 5: Commit** — `ENHANCEMENT: rtc-metal intersects flat and smooth triangles (book chapter 15)`

---

### Task 2: BVH builder (CPU only)

**Files:**
- Create: `src/bvh.rs` (and `pub mod bvh;` in `src/lib.rs`)
- Modify: `src/types.rs` (`GpuBvhNode`)
- Test: unit tests inside `src/bvh.rs`

**Interfaces:**
- Produces in `types.rs` (layout from the spec, 48 bytes):
  ```rust
  #[repr(C, align(16))] pub struct GpuBvhNode {
      pub min: Float4, pub max: Float4,
      pub start: u32, pub count: u32, pub right: u32, pub _pad: u32 }
  ```
- Produces in `bvh.rs`:
  ```rust
  pub const MAX_DEPTH: usize = 32;
  pub const LEAF_SIZE: usize = 4;
  /// A flattened leaf. `bounds` is its f32 world box, None when unbounded.
  pub struct Leaf { pub shape: GpuShape, pub id: usize, pub bounds: Option<([f32; 3], [f32; 3])> }
  pub struct Bvh {
      pub shapes: Vec<GpuShape>,   // [unbounded…, bounded…], leaves' ranges contiguous
      pub ids: Vec<usize>,         // rtc-rust shape id per entry of `shapes`
      pub nodes: Vec<GpuBvhNode>,  // depth-first; empty when nothing is bounded
      pub unbounded_count: u32,
  }
  pub fn build(leaves: Vec<Leaf>) -> Result<Bvh, String>; // Err if depth > MAX_DEPTH
  ```
- Node semantics: `count > 0` is a leaf covering `shapes[start..start + count]`. Otherwise it's an inner node with left child `index + 1` and right child `right`. Node min/max carry `w = 0`.

- [ ] **Step 1: Write the failing tests** in `src/bvh.rs`. Leaves use a dummy `GpuShape` (all zero bytes via `unsafe { std::mem::zeroed() }` is fine in tests) with `id = i`.

  - `test_chap_15_bvh_covers_every_leaf_once`: 2 unbounded leaves plus 37 unit boxes at random-ish integer positions (deterministic, e.g. `((i * 7) % 11, (i * 3) % 5, i % 13)`):
    - `ids[..2]` are the unbounded ids and `unbounded_count == 2`.
    - Every bounded id appears in exactly one leaf range.
    - Leaf ranges are disjoint and together cover `shapes[2..]`.
  - `test_chap_15_bvh_boxes_nest`: on the same input, each inner node's box contains both children's boxes, and each leaf's box contains its shapes' boxes. Leaf `count <= LEAF_SIZE`.
  - `test_chap_15_bvh_identical_boxes_terminate` (Review Focus 5): 100 leaves with the same box. Building succeeds, every leaf has `count <= 4`, and depth ≤ `MAX_DEPTH`.
  - `test_chap_15_bvh_empty_and_unbounded_only`: no leaves gives empty `nodes` and `shapes`. Three unbounded leaves give empty `nodes` and `unbounded_count == 3`.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cargo test --release --lib bvh`
Expected: FAIL to compile (`bvh` module missing).

- [ ] **Step 3: Implement `build`.** Partition unbounded leaves first, keeping their order. Then recurse over the bounded slice:

  ```
  node(range, depth):
    if depth > MAX_DEPTH: Err
    box = union of leaf boxes in range
    if range.len() <= LEAF_SIZE: push leaf {box, start, count}; return
    axis = longest axis of the box of the leaves' box centroids
    sort range by centroid[axis] (f32::total_cmp); mid = range.len() / 2
    push inner {box, count: 0}; recurse left half; set right = nodes.len(); recurse right half
  ```

  Splitting by count, not position, guarantees termination for identical centroids.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cargo test --release --lib bvh`
Expected: PASS.

- [ ] **Step 5: Commit** — `ENHANCEMENT: rtc-metal median-split BVH builder (book chapter 15)`

---

### Task 3: GpuScene carries the BVH (shader loops still linear)

**Files:**
- Modify: `src/scene.rs` (`flatten` builds forward transform and bounds; `GpuScene { info, shapes, nodes, ids }`)
- Modify: `src/types.rs` (`GpuSceneCounts`; `GpuSceneInfo.counts` replaces `shape_count`)
- Modify: `src/gpu.rs` (`render` binds nodes at buffer 3; `run_with_shapes` becomes `run_with_scene`)
- Modify: `shaders/rt.metal`, `shaders/render.metal`, `shaders/probes.metal` (`SceneCounts`, new kernel arguments)
- Modify: `src/probes.rs` (refractive and trace probes go through `GpuScene`)
- Test: `src/scene.rs` tests, `src/probes.rs`

**Interfaces:**
- Consumes: `bvh::{Leaf, Bvh, build}` from Task 2.
- Produces (Rust):
  ```rust
  #[repr(C)] pub struct GpuSceneCounts { pub shape_count: u32, pub unbounded_count: u32, pub node_count: u32 }
  pub struct GpuSceneInfo { pub camera: GpuCamera, pub light: GpuLight, pub counts: GpuSceneCounts } // still 144
  pub struct GpuScene { pub info: GpuSceneInfo, pub shapes: Vec<GpuShape>, pub nodes: Vec<GpuBvhNode>, pub ids: Vec<usize> }
  impl Gpu { pub fn run_with_scene<I: Copy, O: Copy>(&self, kernel: &str, input: &[I], scene: &GpuScene) -> Result<Vec<O>, String> }
  ```
  `run_with_scene` binds shapes at 2, nodes at 3 and `scene.info.counts` via `setBytes` at 4. An empty `shapes` or `nodes` vector gets a one-element placeholder buffer, as `render` already does for shapes.
- Produces (MSL): `struct SceneCounts { uint shape_count; uint unbounded_count; uint node_count; };` and `struct BvhNode` (48 bytes, mirroring `GpuBvhNode`). `SceneInfo.counts` replaces `shape_count`. `render` takes `device const BvhNode *nodes [[buffer(3)]]`. Probe kernels that read the scene take `device const Shape *shapes [[buffer(2)]], device const BvhNode *nodes [[buffer(3)]], constant SceneCounts &counts [[buffer(4)]]`. `RefractiveProbe` and `TraceProbe` drop their `shape_count` field, so their size asserts change: RefractiveProbe 48, TraceProbe 48.
- Bounds: `flatten` carries `parent_transform: &Matrix4` alongside `parent_inverse`, with `full = parent_transform * *shape.transform()`. The leaf box is `shape.bounds()`. If any of its object-space min/max coordinates is non-finite, the leaf is unbounded. Otherwise transform it with `full` and cast outward (min: `f32::next_down` when the cast rounded up; max: `f32::next_up` when it rounded down). If any result is non-finite, the leaf is also unbounded.

- [ ] **Step 1: Write the failing tests** in `src/scene.rs`

  - `test_chap_15_scene_bounds_and_order`. The world has:
    - a plane;
    - a plane under `rotation_x(0.3)`, which has NaN corners if transformed naively (Review Focus 4);
    - an uncapped `Cylinder` under `rotation_z(0.5)`;
    - a sphere under `translation(5,0,0)`;
    - a group under `scaling(2,2,2)` holding a sphere.

    Assert:
    - `info.counts.unbounded_count == 3` and `info.counts.shape_count == 5`.
    - `scene.ids` lists the three unbounded ids first.
    - The single node, a leaf of 2, has a box containing [4,6]×[−1,1]×[−1,1] and [−2,2]³.
    - The f32 box contains the f64 box (compare as f64).
  - `test_chap_15_scene_without_bounded_shapes_renders` (Review Focus 3): a world with one plane and a light, and a 4×2 camera. `GpuScene::from_world` gives `nodes.is_empty()` and `counts.node_count == 0`, and `Gpu::new()?.render(&scene)` returns 8 pixels. Do the same for a world with no shapes.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cargo test --release --lib scene`
Expected: FAIL to compile (`counts`, `nodes`, `ids` missing).

- [ ] **Step 3: Implement**
  - `flatten` produces `Vec<Leaf>`, and `from_world` calls `bvh::build` and fills `info.counts`.
  - Update `render` and `run_with_scene`.
  - In MSL, the three loops still run `for (i = 0; i < counts.shape_count; i++)` over the reordered array. Their signatures change from `(…, device const Shape *shapes, uint count, …)` to `(…, device const Shape *shapes, device const BvhNode *nodes, SceneCounts counts, …)`, so Task 4 changes only their bodies.
  - In the probes, `test_chap_11_gpu_refractive_indices` and `test_chap_13_gpu_cone_far_origin_no_acne` build `GpuScene::from_world` and map shape ids to GPU indices with `scene.ids` (`hit_index = scene.ids.iter().position(|&id| id == object_id)`). The refractive world gets a point light, since `from_world` requires one.

- [ ] **Step 4: Run all tests**

Run: `cargo test --release`
Expected: PASS (every existing probe and compare unchanged in outcome).

- [ ] **Step 5: Commit** — `ENHANCEMENT: rtc-metal GpuScene builds a BVH and uploads nodes (book chapter 15)`

---

### Task 4: Shader BVH traversal

**Files:**
- Modify: `shaders/rt.metal` (`BVH_STACK`, `intersect_box`, bodies of `hit_shapes`, `is_shadowed`, `refractive_indices`)
- Test: `src/probes.rs`

**Interfaces:**
- Consumes: `SceneCounts`, `BvhNode` and the loop signatures from Task 3.
- Produces (MSL): `constant uint BVH_STACK = 32;` and `inline bool intersect_box(float4 mn, float4 mx, Ray r, float tmin, float tmax)`, a slab test that is true when the ray's interval overlap with the box is non-empty. Use `1.0 / direction` with IEEE infinities, and `min`/`max` per axis.
- Each loop first tests `shapes[0 .. counts.unbounded_count]`. Then, when `counts.node_count > 0`, it walks: start at node 0; when a node box fails, pop. A leaf tests its shapes and then pops. An inner node pushes `right` and continues at `index + 1`. Interval per loop:
  - `hit_shapes`: `[0, hit.t]`, re-read after every leaf.
  - `is_shadowed`: `[0, distance]`; return `true` on the first blocking root.
  - `refractive_indices`: `[-INFINITY, INFINITY]`, visiting every candidate. Its per-shape body stays as it is.

  Share the stack-step logic through a small inline helper of your choosing (MSL has no closures). The body that differs per loop stays in each loop.

- [ ] **Step 1: Write the failing tests** in `src/probes.rs`

  - `test_chap_15_gpu_trace_through_bvh`: load `../rtc-rust/models/icosahedron.obj` with `World::load_obj_file` under `scaling(1.5,1.5,1.5)`, with a floor plane under `translation(0,−2,0)` and a point light at (−10,10,−10). Use `Camera::new(20, 20, PI / 3.0)` with `view_transform` from (0,1.5,−5) to (0,0,0), and fire `ray_for_pixel` for all 400 pixels. Per ray, `probe_trace` via `run_with_scene` must match the CPU: hit or miss, which shape (via `scene.ids`), `t`, the normal, and shadowed. Also assert that at least 100 rays hit the icosahedron. This passes on the linear loops already, and it's the guard the traversal must keep passing.
  - `test_chap_15_gpu_bvh_shared_edge_has_no_crack` (Review Focus 2): a quad from two `TriangleUV::flat` triangles, (−1,0,−1), (1,0,−1), (1,0,1) and (−1,0,−1), (1,0,1), (−1,0,1), plus 12 small spheres elsewhere so a real tree is built. 50 rays straight down (direction (0,−1,0)) at points along the diagonal from (−1,0,−1) to (1,0,1) must each hit one of the two triangles at t ≈ 1 from origin y = 1.
  - `test_chap_15_gpu_bvh_matches_linear`: a GPU-only regression for the traversal. Add a `probe_trace_linear` kernel that sets `counts.unbounded_count = counts.shape_count` and `counts.node_count = 0` before tracing (so everything is tested linearly) and compare `probe_trace` against it on the chapter 14 hexagon world from `compare.rs` (make `chapter_14_scene` `pub(crate)`): identical index and t for 1000 camera rays.

- [ ] **Step 2: Run the tests**

Run: `cargo test --release --lib chap_15_gpu`
Expected: the trace and shared-edge tests pass on linear loops. `test_chap_15_gpu_bvh_matches_linear` passes trivially too, since both paths are still linear. Record that, then implement.

- [ ] **Step 3: Implement the traversal** as described under Interfaces.

- [ ] **Step 4: Run all tests**

Run: `cargo test --release`
Expected: PASS, including every chapter 7–14 compare (now through the BVH) and the three tests above.

- [ ] **Step 5: Commit** — `ENHANCEMENT: rtc-metal traverses the BVH for hits, shadows and refraction (book chapter 15)`

---

### Task 5: Chapter 15 scene and whole-image compare

**Files:**
- Create: `scenes/chapter15.json`
- Modify: `src/compare.rs` (`test_chap_15_gpu_render_matches_cpu`, timing output)

**Interfaces:**
- Consumes: everything above. Scene JSON format as rtc-rust's `help.rs` documents (`"type": "obj", "file": …`, `"settings": { "bvh_threshold": N }`).

- [ ] **Step 1: Write the scene and the test**
  - `scenes/chapter15.json` contains:
    - a floor plane with `reflective` 0.1;
    - `"../rtc-rust/models/teapot.obj"` (flat), material color (0.9, 0.6, 0.3);
    - `"../rtc-rust/models/teapot_smooth.obj"` (smooth), color (0.3, 0.6, 0.9), `reflective` 0.2;
    - both teapots scaled about 0.1, rotated `rotate_x: -1.5707963268` if the models turn out to be z-up (check in step 3), and placed side by side, lifted 0.01 off the floor;
    - the light at (−10,10,−10);
    - `"settings": { "bvh_threshold": 8 }`.

    If the CPU reference takes over ~20 s at 400×200 in release, swap the smooth teapot for `torus_smooth.obj`.
  - `test_chap_15_gpu_render_matches_cpu` renders through `scene_file(include_str!("../scenes/chapter15.json"), 400, 200)` and calls `check_matches_cpu`. In `check_matches_cpu`, time the GPU and CPU renders with `std::time::Instant` and add both to the existing `eprintln!` line, e.g. `gpu 0.12 s, cpu 4.3 s`.

- [ ] **Step 2: Run it**

Run: `cargo test --release --lib chap_15_gpu_render -- --nocapture`
Expected: PASS within budget, with the timing line printed. If it is over budget, debug with a diff map, as in earlier milestones (the f32 edge-crack risk in the spec), before changing any tolerance.

- [ ] **Step 3: Look at it.** Render with `cargo run --release -- scenes/chapter15.json -o <scratchpad>/chapter15.ppm`, convert with `sips -s format png`, and check that both teapots appear, the smooth one without facets.

- [ ] **Step 4: Run all tests and lint**

Run: `cargo fmt && cargo clippy --release --all-targets && cargo test --release`
Expected: 0 clippy warnings, all tests PASS.

- [ ] **Step 5: Commit** — `ENHANCEMENT: rtc-metal renders OBJ teapots (book chapter 15)`
