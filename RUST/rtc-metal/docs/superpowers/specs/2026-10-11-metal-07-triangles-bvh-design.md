# rtc-metal milestone 7: triangles, OBJ models and a BVH (book chapter 15)

Branch `metal/07-triangles-obj-bvh`, merged into develop with `--no-ff` when
done. rtc-rust (f64, CPU) stays the reference.

## Goal

Render the book's chapter 15 scenes on the GPU: flat and smooth triangles,
including OBJ models (which rtc-rust loads as `TriangleUV`), fast enough for
meshes of 1k–10k triangles. Success: GPU images match the CPU within the
existing budget (at most 0.05% of pixels off by more than 1e-3), backed by
probe tests, and every existing compare test (chapters 7–14) still passes.

## Decisions

- **BVH:** built by rtc-metal over the flattened leaves' world-space boxes,
  not taken from rtc-rust's `divide()` tree.
- **Split rule:** median split on the longest centroid axis, at most 4 shapes
  per leaf.
- **Triangles:** a new `SHAPE_TRIANGLE` kind in the existing tagged `Shape`
  array. `Shape` grows from 224 to 320 bytes for every kind.

## CPU side (`src/scene.rs`, new `src/bvh.rs`, `src/types.rs`)

### Triangles

`gpu_leaf` accepts rtc-rust `Triangle` and `TriangleUV` as `SHAPE_TRIANGLE`.
`GpuShape` gains six `Float4` fields after the existing ones: `p1`, `e1`,
`e2`, `n1`, `n2`, `n3` (object space, as in rtc-rust). A flat `Triangle`
stores its face `normal` as `n1 = n2 = n3`. The MSL `Shape` and the size
asserts on both sides change with it (320 bytes).

### World-space boxes

`flatten` carries the forward transform (parent × own) alongside the inverse.
A leaf's box is `shape.bounds().transform(full)` in f64, then cast to f32
rounded outward (min rounded down, max rounded up), so the f32 box always
contains the f64 box. A box with any infinite side (planes, uncapped
cylinders and cones) marks the shape **unbounded**.

### Build

`GpuScene` reorders `shapes` as `[unbounded…, bounded…]`. The builder runs
over the bounded part:

1. Take the box of the shapes' box centroids and pick its longest axis.
2. Sort the shapes by centroid on that axis and split at the middle.
3. Recurse until a node holds at most 4 shapes.

Nodes are stored in a flat array in depth-first order:

```rust
#[repr(C, align(16))]
pub struct GpuBvhNode {
    pub min: Float4,
    pub max: Float4,
    pub start: u32, // leaf: first shape index into `shapes`
    pub count: u32, // > 0: leaf of shapes[start..start + count]; 0: inner
    pub right: u32, // inner: index of the right child; left is this + 1
    pub _pad: u32,
} // 48 bytes
```

`GpuSceneInfo` gains `unbounded_count` and `node_count`. The render kernel
takes the node array as buffer 3. The builder returns an error if tree depth
exceeds the shader's traversal stack (32 entries). Median split on 10k
triangles gives a depth of about 12. A scene with no bounded shapes has
`node_count = 0`.

### Probes

The trace probes that took a raw shape list build through the same code
path, so they exercise the traversal as well.

## Shader side (`shaders/rt.metal`)

### Triangle intersection

`local_intersect` gets a `SHAPE_TRIANGLE` case: Möller–Trumbore in object
space, as in rtc-rust's `moller_trumbore`, returning at most one root. The
CPU rejects rays with `|det| ≈ 0` (CPU EPSILON 1e-9). The GPU rejects only
`det == 0` and lets the u/v range checks reject near-parallel rays. Probes
confirm this matches the CPU.

### Smooth normals

`Hit` is unchanged. For a triangle hit, `color_at` reruns Möller–Trumbore on
the hit shape with the same ray to get u, v. `normal_at` then interpolates
`n2·u + n3·v + n1·(1 − u − v)` in object space and maps it to world space as
today (transpose(inverse), normalize once), matching rtc-rust's
`local_normal_at` plus `normal_to_world`.

### Traversal

All three loops (`hit_shapes`, `is_shadowed`, `refractive_indices`) test the
unbounded shapes linearly and then walk the BVH with a 32-entry stack and a
slab test on each node box. They differ only in the ray interval and the
stopping rule:

| Loop                 | Interval                | Stopping rule                     |
|----------------------|-------------------------|-----------------------------------|
| `hit_shapes`         | [0, nearest t so far]   | interval shrinks as hits are found |
| `is_shadowed`        | [0, distance to light]  | returns on the first hit           |
| `refractive_indices` | (−∞, ∞)                 | visits every shape the line crosses |

`refractive_indices` needs the whole line because it counts roots behind the
ray origin too. MSL has no closures, so the three loops share a small
`bvh_step`-style routine instead of a callback.

### Behaviour kept on purpose

In the refraction container logic each triangle counts as its own solid,
exactly as on the CPU, which tracks containers by object id. A glass mesh
looks equally odd on both, so the compare tests stay meaningful. Fixing that
would be a CPU change and is out of scope.

### Known risk

f32 cracks along shared triangle edges: a ray that just misses both
neighbours shows as a speckle. If the compare test shows these, debug with
diff maps (for example, widen the u/v range by a hair) rather than guess in
advance.

## Testing

- **BVH builder unit tests (CPU only, no GPU):**
  - Every bounded shape lands in exactly one leaf, and the unbounded ones
    come first.
  - Each node's box contains its children's boxes or its shapes' boxes.
  - Leaves hold at most 4 shapes, and depth stays within 32.
  - The f32 box contains the f64 box.
- **Triangle probes (GPU against CPU):**
  - The book's chapter 15 cases: a parallel ray, misses past each edge, a
    hit at the right t.
  - The smooth-triangle cases: u/v at a hit, the interpolated normal.
  - A grid of rays against transformed flat and smooth triangles, comparing
    hit/miss, t and the normal with rtc-rust's `local_intersect` and
    `normal_to_world`.
- **Traversal probe:** `probe_trace` through the BVH on a few hundred rays
  aimed at a mesh. Hit, normal and shadow must match the CPU per ray.
- **Regression:** the chapter 7–14 compare tests pass unchanged, now through
  the BVH. Chapter 14 has around 100 bounded leaves, so it is a real tree.
- **New compare, `test_chap_15_gpu_render_matches_cpu`:** `scenes/chapter15.json`
  with the flat `teapot.obj` and the smooth `teapot_smooth.obj` (or
  `torus_smooth.obj` if the teapot makes the CPU reference too slow) on a
  floor, one slightly reflective. The models load via the scene file's `obj`
  type, with `bvh_threshold` set so the CPU reference renders in seconds. The
  pixel budget is unchanged.
- **Timing (printed, not asserted):** the chapter 15 test prints the GPU and
  CPU render times side by side.

## Out of scope

SAH splitting, slimmer shape records (shared transform and material tables),
a correct container model for glass meshes, and texture/UV patterns.
