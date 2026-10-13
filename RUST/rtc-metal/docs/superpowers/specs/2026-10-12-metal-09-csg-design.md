# rtc-metal milestone 9: CSG (book chapter 16) and the perturbed pattern

Branch `metal/09-csg`, merged into develop with `--no-ff` when done.
rtc-rust (f64, CPU) stays the reference. Metal's ray tracing API, once
"optional milestone 9", moves to milestone 10.

## Goal

Render rtc-rust's dice scenes (`scenes/dice-light-area.json`,
`dice-sentence.json`) on the GPU. That takes two features:

- **Chapter 16 CSG:** union, intersection and difference of primitives and
  of nested CSG trees. Each die is one CSG tree of 23 leaves,
  `(cube ∩ sphere) − pip − pip − …` with 21 pip spheres, all opaque, no
  groups inside the CSG.
- **The perturbed pattern** (rtc-rust's `PerturbedPattern` over Perlin
  noise, tests `test_chap_17_*`): the dice bodies are perturbed stripes.

Success:
- a GPU probe of the CSG rule agrees with `intersection_allowed` on all 24
  rows of the book's table;
- GPU Perlin noise and perturbed colours match the CPU at sample points;
- per-pixel hit/normal/shadow traces of CSG scenes match the CPU;
- a whole-image compare of a reduced dice scene matches the CPU like the
  other chapters;
- both dice scenes render in `--window`.

## Decisions

- **One BVH entry per CSG root.** The CPU flattens a CSG tree's leaves into
  the shapes array behind the BVH-indexed range and puts a single
  `SHAPE_CSG` entry, with the tree's world bounds, into the BVH. The CSG
  node's own bounds test on the CPU becomes the BVH box test here.
- **Membership-change filter instead of the book's recursive filter.** For
  well-formed operands the book's rule table says exactly one thing: keep a
  hit if and only if it changes whether the ray is inside the combined
  solid. Union keeps `lhit && !inr`, which is "union membership flipped";
  intersection and difference work out the same way. A nested node's
  surviving hits are exactly the hits where the inner solid's membership
  flips, so the outer node sees the same `inl`/`inr` toggles the book's
  recursion produces. The GPU therefore keeps one inside bit per leaf, and
  for each hit, in t order, evaluates the tree before and after flipping
  that leaf's bit. Leaf bits are parity (toggle on every root), as on the
  CPU, so open shapes and planes behave as in `filter_intersections`.
- **All roots, sorted, negative t included.** The CPU's
  `filter_intersections` walks the full ordered list starting with
  `inl = inr = false`, roots behind the ray origin included. The GPU does
  the same, then picks the first surviving root with t ≥ 0.
- **One new buffer, `tables`** (`u32` words), holds every CSG program and
  every distinct Perlin permutation. `GpuShape` and `GpuMaterial` refer to
  their entries by word offset.
- **Perturbed pattern:** exactly rtc-rust's `PerturbedPattern` over one of
  the four base patterns (stripe, gradient, ring, checkers). The base
  pattern keeps its own transform. Nested perturbed, blended or other
  patterns stay rejected.
- **rtc-rust gets one accessor,** `Perlin::permutation(&self) -> &[u8]`,
  because the table is private and the seed isn't stored. Nothing else in
  rtc-rust changes.
- **Rejected at upload with "not supported yet":** transparent leaves in a
  CSG (the dice are opaque, and the CPU's own container walk is unbalanced
  for CSG); groups inside a CSG (the CPU's parity-over-subtree semantics
  aren't union, and the dice don't need it); motion on a CSG node; trees
  over the limits below.
- **Limits:** at most 32 leaves per tree (one `uint` of inside bits) and at
  most 64 roots per tree, counted from each leaf kind's maximum (sphere,
  cube 2; plane, triangle 1; cylinder, cone 4). A die is 23 leaves and at
  most 46 roots.

## CPU side (`src/scene.rs`, `src/types.rs`)

### Tables

```rust
#[derive(Default)]
pub struct Tables {
    pub words: Vec<u32>,
    // permutation -> its offset in `words`, so dice sharing a seed share one
    noise: Vec<(Vec<u8>, u32)>,
}
```

`GpuScene` gains `tables: Vec<u32>`. `gpu_material` and `gpu_shape` take
`&mut Tables`, so a perturbed material can add its permutation. Probes pass
`&mut Tables::default()` when they don't upload one.

A permutation is 256 words, one per entry of `perm[0..256]`. The shader
indexes `perm[i & 255]`, which equals rtc-rust's 512-entry `perm[i]`.

### CSG programs

Each CSG root compiles to a postfix program in `tables.words`:

```
words[p]     = n            // instruction count
words[p + 1] = leaf_base    // first leaf, relative to counts.shape_count
words[p + 2 .. p + 2 + n]   // instructions
```

An instruction is `CSG_LEAF` (0), `CSG_UNION` (1), `CSG_INTERSECTION` (2)
or `CSG_DIFFERENCE` (3), as new constants in `types.rs` and `rt.metal`.
Leaves appear in program order, so the k-th `CSG_LEAF` is inside bit k and
shape `counts.shape_count + leaf_base + k`; it needs no operand. A die
compiles to `L L ∩ L − L − … L −`.

`GpuShape` has no spare field, so a `SHAPE_CSG` entry stores its program
offset `p` in `closed`, documented on both the Rust and the MSL struct. Its
`inverse` and `material` are unused: the leaves carry the full
world-to-object transform, parents' and CSG nodes' transforms included,
exactly as `flatten` already does for groups.

### Flattening

`flatten` gets a CSG arm. It compiles the tree (`compile_csg`, recursive:
left, right, then the operation) into `tables.words` and pushes the leaves
to `CsgLeaves { shapes: Vec<GpuShape>, ids: Vec<usize> }`. It pushes one
`Leaf` of kind `SHAPE_CSG` with `world_bounds(csg, full)`, using the CSG
bounds `World::build_bounds` stores. After `bvh::build`, `from_world`
appends the CSG leaves to `shapes` and their ids to `ids`.
`counts.shape_count` stays the BVH-indexed count, so the existing loops
never visit a CSG leaf directly.

### Material

`GpuMaterial` / MSL `Material` grow from 144 to 224 bytes:

| field | meaning |
|---|---|
| `pattern` | base pattern kind (unchanged) |
| `pattern_inverse` | object → pattern space of the outermost pattern (the perturbed one when perturbed) |
| `base_inverse` (new, mat4) | perturbed only: the base pattern's own inverse; identity otherwise |
| `noise_table` (new, u32) | word offset of the permutation, `NOISE_NONE` = `0xffffffff` when not perturbed |
| `noise_amplitude`, `noise_frequency` (new, f32) | `PerturbedPattern::amplitude`, `frequency` |

`GpuShape` grows from 320 to 400 bytes, and every probe struct embedding a
`Material` or `Shape` grows by 80.

`gpu.rs` uploads `tables` (one zero word when empty, like `nodes`): at
`[[buffer(4)]]` for `render` and the `Accumulator`, and at
`[[buffer(5)]]` for `run_with_scene` probes.

## Shader side (`shaders/rt.metal`)

### Noise and the perturbed pattern

```metal
inline float perlin_noise(device const uint *tables, uint table, float3 p);
inline float3 perlin_vector(device const uint *tables, uint table, float3 p);
```

These are line-for-line ports of `Perlin::noise` (`fade`, `lerp`, `grad`,
`int(floor(x)) & 255`) and `Perlin::vector`, with offsets (31.4, 15.9, 26.5)
and (−27.1, 82.8, −18.2). `pattern_at(Material m, float4 object_point,
device const uint *tables)` applies `pattern_inverse`. When
`noise_table != NOISE_NONE`, it then perturbs
`p += perlin_vector(tables, noise_table, p.xyz * frequency) * amplitude`
and applies `base_inverse`, then evaluates the base pattern as today.
`tables` is threaded through `lighting`, `color_at` and the probes; probes
without a scene pass `nullptr` (their materials are never perturbed).

### CSG

```metal
struct CsgHit { float t; uint leaf; uint root; };
constant uint CSG_MAX_LEAVES = 32;
constant uint CSG_MAX_HITS = 64;

// Is the ray inside the tree's solid, given each leaf's inside bit?
inline bool csg_inside(device const uint *tables, uint p, uint bits);

// Sorted surviving roots of the CSG entry `s` (all t, negative included).
inline uint intersect_csg(Shape s, Ray r, device const Shape *shapes,
                          device const uint *tables, SceneCounts counts,
                          thread CsgHit *out);
```

- `csg_inside` runs the postfix program with a bit stack in a `uint`
  (depth ≤ 32): `CSG_LEAF` pushes bit k, an operation pops two and pushes
  `a | b`, `a & b` or `a & ~b`.
- `intersect_csg` tests every leaf with the existing `intersect` and
  insertion-sorts `(t, leaf, root)` into a thread-local array of
  `CSG_MAX_HITS`. It walks the array, compares `csg_inside` before and after
  flipping the leaf's bit, and compacts the survivors in place.
- `hit_shapes`: a `SHAPE_CSG` entry contributes its first surviving root
  with t ≥ 0. `Hit.index` is the **leaf's** shape index, so `normal_at`,
  materials and patterns work unchanged.
- `is_shadowed`: a `SHAPE_CSG` entry occludes if any surviving root lies in
  [0, distance).
- `refractive_indices` skips `SHAPE_CSG` entries. It only runs for
  transparent hits, and CSG leaves are opaque, so the only divergence from
  the CPU is glass inside a die's volume, which no scene has.

`intersect(Shape, Ray, ts)` stays as it is for every other kind. The CSG
branch sits in its three callers instead, because it needs `shapes`,
`tables` and `counts`.

## Testing

CSG tests are named `test_chap_16_*`; noise and perturbed tests
`test_chap_17_*`, as in rtc-rust.

1. **Noise probe:** `probe_noise` vs `Perlin::noise` and `Perlin::vector`
   at a few hundred points over several seeds, including negative
   coordinates; f32 tolerance.
2. **Perturbed pattern probe:** `pattern_at` of a perturbed-stripe material
   vs rtc-rust's colour at object points away from stripe edges.
3. **Program compile:** unit tests on `compile_csg`: a single node, a nested
   tree and a die give the expected words, leaf order and leaf transforms;
   transparent leaves, groups inside a CSG and the 32/64 limits are
   rejected.
4. **Rule probe:** `probe_csg_inside` evaluates the programs of a union, an
   intersection and a difference for each bit pattern. Comparing before and
   after the hit side's flip matches `intersection_allowed` on all 24 rows.
5. **Intersections:** the book's "a ray hits/misses a CSG object" scenarios
   plus a nested tree and a die, via `probe_trace` against rtc-rust's
   `World::intersect`/hit/`is_shadowed`, and the BVH and linear traces
   agree.
6. **Whole image:** `scenes/chapter16.json`, a checkered floor and two
   marble dice from `dicelib`'s builders, one point light, small
   resolution, compared with the CPU through `check_matches_cpu`.
7. **Manual:** `rtc-metal --window ../rtc-rust/scenes/dice-light-area.json`
   and `dice-sentence.json` render and converge.

## Order of work

1. `tables` buffer plumbing (empty), all tests green.
2. Perturbed pattern: `Perlin::permutation`, `Tables::noise`, material
   layout, noise and pattern probes.
3. `compile_csg` + `flatten` CSG arm + unit tests.
4. `csg_inside` + rule probe.
5. `intersect_csg` in `hit_shapes`, `is_shadowed`, `refractive_indices` +
   trace probes.
6. `scenes/chapter16.json` + whole-image compare; manual window check.

Each step leaves the build and tests green.

## Out of scope

Metal's ray tracing API (milestone 10), a per-CSG BVH to cull pips,
transparent CSG leaves, groups inside a CSG, bump maps, blended, nested and
texture-map patterns, focal blur and path tracing.
