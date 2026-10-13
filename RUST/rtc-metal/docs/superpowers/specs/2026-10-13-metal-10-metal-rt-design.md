# rtc-metal milestone 10: Metal's ray tracing API

Branch `metal/10-metal-rt`, merged into develop with `--no-ff` when done.
rtc-rust (f64, CPU) stays the reference.

## Goal

Learn Metal's ray tracing API by adding a second traversal next to the
hand-written BVH walk, chosen per render. Both paths must give the image
the CPU gives, and a timing comparison shows what the API buys on this
machine.

The development machine is an M2 Max, which has no ray tracing hardware
(that arrived with M3). There Metal runs the same API on Apple's tuned
software traversal, so a speed-up is not a success criterion; the
measurement is.

Success:
- every whole-image chapter compare (7, 9–16) passes on both traversals;
- per-ray traces of the triangle, CSG and analytic scenes agree between
  traversals;
- `--traversal metal` renders and converges in `--window`;
- the spec's timing table is filled in for both traversals.

## Decisions

- **Keep both traversals.** `SceneInfo.traversal`: `TRAVERSAL_BVH` (0) is
  today's walk, `TRAVERSAL_METAL` (1) the new one. The host always builds
  and binds the acceleration structure, so switching is one field and every
  pixel takes the same branch. Default stays `bvh`.
- **Inline `intersection_query`**, not `intersector` with intersection
  function tables: the kernel loops over candidates itself and calls the
  existing `intersect`/`intersect_csg`, so no linked functions, function
  tables or pipeline changes.
- **One primitive acceleration structure per scene, two geometries:**
  - triangle geometry: every triangle's world-space vertices, opaque;
  - bounding-box geometry: every other bounded `shapes` entry (sphere,
    cube, cylinder, cone, each CSG root) with its world box from the BVH
    leaf, not opaque.
  No instancing: groups are already flattened and motion is rejected.
- **Planes stay in the unbounded linear loop** on both traversals. An
  infinite box cannot go in an acceleration structure.
- **`refractive_indices` keeps the BVH walk** on both traversals. It needs
  roots behind the ray origin, which a query will not report, and runs only
  for transparent hits.
- **Triangle hits come from Metal** (f32, world space) rather than our
  object-space Möller–Trumbore. Their `uv` for smooth normals is still
  computed by `triangle_uv`, so shading code is unchanged. Expect small
  edge differences against the BVH path; the per-ray test budgets them.

## Shader side (`shaders/rt.metal`)

```metal
#include <metal_raytracing>
using namespace metal::raytracing;

constant uint TRAVERSAL_BVH = 0;
constant uint TRAVERSAL_METAL = 1;
```

`SceneCounts` gains `triangle_count`: the first `triangle_count` words of
`rt_ids` map triangle `primitive_id`s to `shapes` indices, the rest map box
`primitive_id`s. `SceneInfo` gains `traversal`.

`hit_shapes` and `is_shadowed` take two more arguments,
`device const uint *rt_ids` and `primitive_acceleration_structure accel`,
threaded through `intensity_at`, `color_at` and the probes like `tables`.

### `hit_shapes`, `TRAVERSAL_METAL`

1. Run the unbounded loop (planes) as today; it leaves `hit.t`.
2. If `counts.node_count > 0`, run
   `intersection_query<triangle_data> q(ray, accel, params)` over
   `[0, hit.t]` (`ray.min_distance = 0`, `max_distance = hit.t`, infinite
   when nothing hit yet). For each candidate (only boxes are candidates;
   opaque triangles commit inside the query):
   - `i = rt_ids[counts.triangle_count + q.get_candidate_primitive_id()]`;
   - CSG: `intersect_csg`, else `intersect`; take the first root with
     `0 ≤ t < q.get_committed_distance()`; if any, record
     `Hit{t, index, root}` (for CSG the leaf and its root) and
     `q.commit_bounding_box_intersection(t)`.
3. After the loop, if `q.get_committed_intersection_type()` is
   `triangle`, set `hit = Hit{q.get_committed_distance(),
   rt_ids[q.get_committed_primitive_id()], 0}`. A committed box already set
   `hit`.

### `is_shadowed`, `TRAVERSAL_METAL`

Unbounded loop as today, then the same query over `[0, distance)` with
`params.accept_any_intersection(true)`. A box candidate commits if any root
lies in `[0, distance)`; a committed triangle or box means shadowed.

### Unchanged

`refractive_indices`, `color_at`'s iteration, lighting, patterns, noise,
CSG membership, sampling, `render`/`resolve` kernels apart from the two new
arguments.

## Host side (`src/`)

### `GpuScene` (`scene.rs`, `bvh.rs`)

`GpuScene` gains the acceleration-structure inputs, built during
`from_world` after `bvh::build`:

```rust
pub rt_vertices: Vec<[f32; 3]>,           // 3 per triangle, world space
pub rt_boxes: Vec<([f32; 3], [f32; 3])>,  // per non-triangle bounded entry
pub rt_ids: Vec<u32>,                     // triangle shape indices, then box ones
```

`Bvh` keeps each bounded entry's world box so `rt_boxes` reuses the BVH's
outward-rounded f32 boxes. Triangle vertices are rtc-rust's `p1, p2, p3`
transformed by the triangle's full world transform (parents included) in
f64, then rounded to f32. `info.counts.triangle_count = triangles`.

### Acceleration structure (`src/accel.rs`, new)

```rust
pub struct Accel {
    pub structure: Retained<ProtocolObject<dyn MTLAccelerationStructure>>,
    pub rt_ids: Buffer,
}
pub fn build(gpu: &Gpu, scene: &GpuScene) -> Result<Accel, String>;
```

- `MTLPrimitiveAccelerationStructureDescriptor` with a
  `MTLAccelerationStructureTriangleGeometryDescriptor` (vertex buffer of
  packed `float3`, stride 12, no index buffer, `opaque = true`) when there
  are triangles, and a
  `MTLAccelerationStructureBoundingBoxGeometryDescriptor`
  (`MTLAxisAlignedBoundingBox`, stride 24, `opaque = false`) when there are
  boxes. Triangles are geometry 0 only when present; the shader never reads
  `geometry_id`, it relies on Metal numbering `primitive_id` per geometry.
- Sizes from `accelerationStructureSizesWithDescriptor`, a private scratch
  buffer, one `MTLAccelerationStructureCommandEncoder` build, waited on.
  No compaction or refit (static scenes).
- A scene with nothing bounded gets a structure over one zero-size box at
  the origin, so the argument is always bound; the shader skips the query
  when `counts.node_count == 0`.

### Binding (`gpu.rs`)

`Gpu::render`, `Accumulator` and `run_with_scene` build an `Accel` with
the other scene buffers and bind `rt_ids` and the structure
(`setAccelerationStructure_atBufferIndex`): buffers 5 and 6 for `render`,
6 and 7 for probes. A primitive structure bound directly needs no
`useResource`.

### CLI (`main.rs`)

`--traversal bvh|metal` sets `info.traversal` for the one-shot render, the
sampled render and `--window`. The existing `GPU render: … ms` line prints
the traversal name.

## Testing

Test names follow the chapter of the scene they use.

1. **Accel build:** `accel::build` succeeds for a triangles-only, a
   boxes-only, a mixed and a planes-only scene (unit tests in `accel.rs`).
2. **Per-ray traversal agreement:** `probe_trace` on the chapter 15
   teapot, chapter 16 dice and chapter 12/13 shapes with
   `traversal = METAL` against `traversal = BVH`: same miss/hit and shadow
   on every ray, same `index` except on triangle edges (both hits on
   triangles of the same mesh), `t` within 1e-4 relative. Edge
   disagreements budgeted at 1 in 1000 rays.
3. **Shadow occluders:** on `METAL`, a point whose only occluder is a
   sphere (a box primitive) and one whose only occluder is a triangle are
   both shadowed, and an occluder beyond the light is not.
4. **Whole image:** `check_matches_cpu` renders on both traversals and
   checks each against the CPU with the existing budget; every
   `test_chap_N_gpu_render_matches_cpu` therefore covers both.
5. **Manual:** `--window --traversal metal` on `dice-light-area.json` and
   the chapter 15 scene.

## Timing (filled in at the end)

Full-size headless renders, release build, median of 3, `GPU render` line:

| scene | size | bvh | metal |
|---|---|---|---|
| `scenes/chapter15.json` | | | |
| `../rtc-rust/scenes/dice-light-area.json` | 3456×2234 | | |
| `../rtc-rust/scenes/cover.json` | | | |

## Order of work

1. `rt_vertices`/`rt_boxes`/`rt_ids` in `GpuScene`, `accel::build`, unit
   tests; nothing bound yet.
2. Bind `rt_ids` and the structure everywhere, add `traversal` and
   `triangle_count`; the shader accepts them and ignores them. Tests green.
3. `hit_shapes` and `is_shadowed` on `TRAVERSAL_METAL`; per-ray agreement
   and shadow probes.
4. `check_matches_cpu` on both traversals.
5. `--traversal` CLI, window check, timing table.

Each step leaves the build and tests green.

## Out of scope

Intersection function tables and `intersector<>`, instance acceleration
structures, compaction and refit, motion, moving `refractive_indices` onto
the query, changing the default traversal, and the milestone 9 deferred
minors.
