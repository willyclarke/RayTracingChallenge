
// Probe kernels: one thread per test case, each calling a single rt.metal
// function so Rust tests (src/probes.rs) can compare it with rtc-rust.
// Gpu::run dispatches exactly one thread per input, so no bounds checks.

struct TransformProbe {
    float4x4 m;
    Ray ray;
};

kernel void probe_ray_transform(device const TransformProbe *in [[buffer(0)]],
                                device Ray *out                 [[buffer(1)]],
                                uint i                          [[thread_position_in_grid]])
{
    out[i] = ray_transform(in[i].m, in[i].ray);
}

// out[i] = (t0, t1, hit count, unused)
kernel void probe_intersect_sphere(device const Ray *in [[buffer(0)]],
                                   device float4 *out   [[buffer(1)]],
                                   uint i               [[thread_position_in_grid]])
{
    float2 ts = float2(0.0);
    uint count = intersect_sphere(in[i], ts);
    out[i] = float4(ts, float(count), 0.0);
}

struct NormalProbe {
    Shape shape;
    float4 point;
};

kernel void probe_normal_at(device const NormalProbe *in [[buffer(0)]],
                            device float4 *out           [[buffer(1)]],
                            uint i                       [[thread_position_in_grid]])
{
    out[i] = normal_at(in[i].shape, in[i].point, float2(0.0));
}

// The ray is in object space; shape.inverse is unused.
struct IntersectProbe {
    Shape shape;
    Ray ray;
};
static_assert(sizeof(IntersectProbe) == 432, "IntersectProbe must match probes.rs");

struct IntersectResult {
    float4 ts;  // ascending; ts[count..] unspecified
    uint count;
};
static_assert(sizeof(IntersectResult) == 32, "IntersectResult must match probes.rs");

kernel void probe_local_intersect(device const IntersectProbe *in [[buffer(0)]],
                                  device IntersectResult *out     [[buffer(1)]],
                                  uint i                          [[thread_position_in_grid]])
{
    float4 ts = float4(0.0);
    uint count = local_intersect(in[i].shape, in[i].ray, ts);
    out[i] = IntersectResult{ts, count};
}

struct LightingProbe {
    Material material;
    Light light;
    float4 point;
    float4 object_point;
    float4 eyev;
    float4 normalv;
    float intensity;
    uint pixel;  // the Sampler; sampling 0 is rtc-rust's deterministic grid
    uint frame;
    uint sampling;
};
static_assert(sizeof(LightingProbe) == 384, "LightingProbe must match probes.rs");

kernel void probe_lighting(device const LightingProbe *in [[buffer(0)]],
                           device float4 *out             [[buffer(1)]],
                           uint i                         [[thread_position_in_grid]])
{
    LightingProbe p = in[i];
    Sampler s = Sampler{p.pixel, p.frame, p.sampling, 0};
    out[i] = float4(lighting(p.material, p.light, p.point, p.object_point, p.eyev, p.normalv,
                             p.intensity, s, nullptr),
                    0.0);
}

struct PointOnLightProbe {
    Light light;
    uint u;
    uint v;
    float ju;
    float jv;
};
static_assert(sizeof(PointOnLightProbe) == 96, "PointOnLightProbe must match probes.rs");

kernel void probe_point_on_light(device const PointOnLightProbe *in [[buffer(0)]],
                                 device float4 *out                 [[buffer(1)]],
                                 uint i                             [[thread_position_in_grid]])
{
    PointOnLightProbe p = in[i];
    out[i] = light_cell_point(p.light, p.u, p.v, p.ju, p.jv);
}

struct IntensityProbe {
    Light light;
    float4 point;
    float4 normalv;
    uint pixel;  // the Sampler
    uint frame;
    uint sampling;
};
static_assert(sizeof(IntensityProbe) == 128, "IntensityProbe must match probes.rs");

kernel void probe_intensity_at(device const IntensityProbe *in [[buffer(0)]],
                               device float *out               [[buffer(1)]],
                               device const Shape *shapes      [[buffer(2)]],
                               device const BvhNode *nodes     [[buffer(3)]],
                               constant SceneCounts &counts    [[buffer(4)]],
                               device const uint *tables       [[buffer(5)]],
                               device const uint *rt_ids [[buffer(6), function_constant(METAL_RT)]],
                               raytracing::primitive_acceleration_structure triangles [[buffer(7), function_constant(METAL_RT)]],
                               raytracing::primitive_acceleration_structure boxes [[buffer(8), function_constant(METAL_RT)]],
                               uint i                          [[thread_position_in_grid]])
{
    RtScene rt;
    if (METAL_RT) {
        rt = RtScene{rt_ids, triangles, boxes};
    }
    IntensityProbe p = in[i];
    Sampler s = Sampler{p.pixel, p.frame, p.sampling, 0};
    out[i] = intensity_at(p.light, p.point, p.normalv, shapes, nodes, tables, rt, counts, s);
}

struct ShadeProbe {
    Material material;
    Light light;
    float4 point;
    float4 eyev;
    float4 normalv;
    uint pixel;  // the Sampler
    uint frame;
    uint sampling;
};
static_assert(sizeof(ShadeProbe) == 368, "ShadeProbe must match probes.rs");

// One shading call as color_at makes it: intensity_at and lighting share the
// Sampler, which keeps cell A (shadow) and cell B (diffuse) independent.
kernel void probe_shade_point(device const ShadeProbe *in    [[buffer(0)]],
                              device float4 *out             [[buffer(1)]],
                              device const Shape *shapes     [[buffer(2)]],
                              device const BvhNode *nodes    [[buffer(3)]],
                              constant SceneCounts &counts   [[buffer(4)]],
                              device const uint *tables      [[buffer(5)]],
                              device const uint *rt_ids [[buffer(6), function_constant(METAL_RT)]],
                              raytracing::primitive_acceleration_structure triangles [[buffer(7), function_constant(METAL_RT)]],
                              raytracing::primitive_acceleration_structure boxes [[buffer(8), function_constant(METAL_RT)]],
                              uint i                         [[thread_position_in_grid]])
{
    RtScene rt;
    if (METAL_RT) {
        rt = RtScene{rt_ids, triangles, boxes};
    }
    ShadeProbe p = in[i];
    Sampler s = Sampler{p.pixel, p.frame, p.sampling, 0};
    float intensity = intensity_at(p.light, p.point, p.normalv, shapes, nodes, tables, rt, counts, s);
    out[i] = float4(lighting(p.material, p.light, p.point, p.point, p.eyev, p.normalv,
                             intensity, s, tables),
                    0.0);
}

struct CameraProbe {
    Camera camera;
    uint2 pixel;
};
static_assert(sizeof(CameraProbe) == 112, "CameraProbe must match probes.rs");

kernel void probe_ray_for_pixel(device const CameraProbe *in [[buffer(0)]],
                                device Ray *out              [[buffer(1)]],
                                uint i                       [[thread_position_in_grid]])
{
    out[i] = ray_for_pixel(in[i].camera, in[i].pixel.x, in[i].pixel.y, 0.5, 0.5);
}

struct LensRayProbe {
    Camera camera;
    uint2 pixel;
    float4 offsets;  // dx, dy, lx, ly
};
static_assert(sizeof(LensRayProbe) == 128, "LensRayProbe must match probes.rs");

kernel void probe_ray_for_lens(device const LensRayProbe *in [[buffer(0)]],
                               device Ray *out               [[buffer(1)]],
                               uint i                        [[thread_position_in_grid]])
{
    float4 o = in[i].offsets;
    out[i] = ray_for_lens(in[i].camera, in[i].pixel.x, in[i].pixel.y, o.x, o.y, o.z, o.w);
}

struct LensProbe {
    uint2 pixel;
    uint i;
    uint n;
};
static_assert(sizeof(LensProbe) == 16, "LensProbe must match probes.rs");

kernel void probe_lens_sample(device const LensProbe *in [[buffer(0)]],
                              device float4 *out         [[buffer(1)]],
                              uint i                     [[thread_position_in_grid]])
{
    out[i] = lens_sample(in[i].pixel.x, in[i].pixel.y, in[i].i, in[i].n);
}

struct RefractiveProbe {
    Ray ray;
    uint hit_index;
    uint hit_root;
};
static_assert(sizeof(RefractiveProbe) == 48, "RefractiveProbe must match probes.rs");

// out[i] = (n1, n2) at root hit_root of shapes[hit_index] along ray
kernel void probe_refractive_indices(device const RefractiveProbe *in [[buffer(0)]],
                                     device float2 *out               [[buffer(1)]],
                                     device const Shape *shapes       [[buffer(2)]],
                                     device const BvhNode *nodes      [[buffer(3)]],
                                     constant SceneCounts &counts     [[buffer(4)]],
                                     device const uint *tables        [[buffer(5)]],
                                     uint i                           [[thread_position_in_grid]])
{
    RefractiveProbe p = in[i];
    Shape s = shapes[p.hit_index];
    float4 ts;
    intersect(s, p.ray, ts);
    Hit hit = Hit{ts[p.hit_root], p.hit_index, p.hit_root};
    out[i] = refractive_indices(p.ray, hit, shapes, nodes, tables, counts);
}

// in[i] = (cos_i, n1, n2, unused)
kernel void probe_schlick(device const float4 *in [[buffer(0)]],
                          device float *out       [[buffer(1)]],
                          uint i                  [[thread_position_in_grid]])
{
    out[i] = schlick(in[i].x, in[i].y, in[i].z);
}

struct PatternProbe {
    Shape shape;
    float4 point;
};

// out[i] = the pattern color of in[i].shape at world point in[i].point
kernel void probe_pattern_at(device const PatternProbe *in [[buffer(0)]],
                             device float4 *out            [[buffer(1)]],
                             uint i                        [[thread_position_in_grid]])
{
    Shape s = in[i].shape;
    out[i] = float4(pattern_at(s.material, s.inverse * in[i].point, nullptr), 0.0);
}

// probe_pattern_at for perturbed materials, which read the scene's tables.
kernel void probe_pattern_at_scene(device const PatternProbe *in [[buffer(0)]],
                                   device float4 *out            [[buffer(1)]],
                                   device const Shape *shapes    [[buffer(2)]],
                                   device const BvhNode *nodes   [[buffer(3)]],
                                   constant SceneCounts &counts  [[buffer(4)]],
                                   device const uint *tables     [[buffer(5)]],
                                   uint i                        [[thread_position_in_grid]])
{
    Shape s = in[i].shape;
    out[i] = float4(pattern_at(s.material, s.inverse * in[i].point, tables), 0.0);
}

struct NoiseProbe {
    float4 point;
    uint table;
    uint pad[3];
};
static_assert(sizeof(NoiseProbe) == 32, "NoiseProbe must match probes.rs");

// out[i] = (Perlin::noise, Perlin::vector) at in[i].point, with the
// permutation at tables[in[i].table].
kernel void probe_noise(device const NoiseProbe *in  [[buffer(0)]],
                        device float4 *out           [[buffer(1)]],
                        device const Shape *shapes   [[buffer(2)]],
                        device const BvhNode *nodes  [[buffer(3)]],
                        constant SceneCounts &counts [[buffer(4)]],
                        device const uint *tables    [[buffer(5)]],
                        uint i                       [[thread_position_in_grid]])
{
    float3 p = in[i].point.xyz;
    out[i] = float4(perlin_noise(tables, in[i].table, p), perlin_vector(tables, in[i].table, p));
}

struct CsgInsideProbe {
    uint program;
    uint bits;
};

// Chapter 16: out[i] = csg_inside of the program at in[i].program.
kernel void probe_csg_inside(device const CsgInsideProbe *in [[buffer(0)]],
                             device uint *out                [[buffer(1)]],
                             device const Shape *shapes      [[buffer(2)]],
                             device const BvhNode *nodes     [[buffer(3)]],
                             constant SceneCounts &counts    [[buffer(4)]],
                             device const uint *tables       [[buffer(5)]],
                             uint i                          [[thread_position_in_grid]])
{
    out[i] = csg_inside(tables, in[i].program, in[i].bits) ? 1 : 0;
}

struct TraceProbe {
    Ray ray;
    float4 light_position;
};
static_assert(sizeof(TraceProbe) == 48, "TraceProbe must match probes.rs");

struct TraceResult {
    float4 normal;  // world-space, flipped toward the eye as in color_at
    float t;
    uint index;
    uint root;
    uint shadowed;
};
static_assert(sizeof(TraceResult) == 32, "TraceResult must match probes.rs");

// The first steps of color_at for one ray: nearest hit, normal, shadow test.
// index = 0xffffffff on a miss (counts.shape_count is the first CSG leaf).
inline TraceResult trace(TraceProbe p, device const Shape *shapes, device const BvhNode *nodes,
                         device const uint *tables, RtScene rt, SceneCounts counts) {
    Hit hit;
    if (!hit_shapes(p.ray, shapes, nodes, tables, rt, counts, hit)) {
        return TraceResult{float4(0.0), INFINITY, 0xffffffff, 0, 0};
    }
    float4 point = p.ray.origin + p.ray.direction * hit.t;
    Shape s = shapes[hit.index];
    float2 uv = s.kind == SHAPE_TRIANGLE ? triangle_uv(s, ray_transform(s.inverse, p.ray))
                                         : float2(0.0);
    float4 normalv = normal_at(s, point, uv);
    if (dot(normalv, -p.ray.direction) < 0.0) {
        normalv = -normalv;
    }
    bool shadowed =
        is_shadowed(p.light_position, point + normalv * EPSILON, shapes, nodes, tables, rt, counts);
    return TraceResult{normalv, hit.t, hit.index, hit.root, shadowed ? 1u : 0u};
}

kernel void probe_trace(device const TraceProbe *in  [[buffer(0)]],
                        device TraceResult *out      [[buffer(1)]],
                        device const Shape *shapes   [[buffer(2)]],
                        device const BvhNode *nodes  [[buffer(3)]],
                        constant SceneCounts &counts [[buffer(4)]],
                        device const uint *tables    [[buffer(5)]],
                        device const uint *rt_ids [[buffer(6), function_constant(METAL_RT)]],
                        raytracing::primitive_acceleration_structure triangles [[buffer(7), function_constant(METAL_RT)]],
                        raytracing::primitive_acceleration_structure boxes [[buffer(8), function_constant(METAL_RT)]],
                        uint i                       [[thread_position_in_grid]])
{
    RtScene rt;
    if (METAL_RT) {
        rt = RtScene{rt_ids, triangles, boxes};
    }
    out[i] = trace(in[i], shapes, nodes, tables, rt, counts);
}

// probe_trace with the BVH switched off: every shape counts as unbounded,
// so the loops test them all linearly. The reference for the traversal.
kernel void probe_trace_linear(device const TraceProbe *in  [[buffer(0)]],
                               device TraceResult *out      [[buffer(1)]],
                               device const Shape *shapes   [[buffer(2)]],
                               device const BvhNode *nodes  [[buffer(3)]],
                               constant SceneCounts &counts [[buffer(4)]],
                               device const uint *tables    [[buffer(5)]],
                               device const uint *rt_ids [[buffer(6), function_constant(METAL_RT)]],
                               raytracing::primitive_acceleration_structure triangles [[buffer(7), function_constant(METAL_RT)]],
                               raytracing::primitive_acceleration_structure boxes [[buffer(8), function_constant(METAL_RT)]],
                               uint i                       [[thread_position_in_grid]])
{
    RtScene rt;
    if (METAL_RT) {
        rt = RtScene{rt_ids, triangles, boxes};
    }
    SceneCounts linear = counts;
    linear.unbounded_count = counts.shape_count;
    linear.node_count = 0;
    linear.traversal = TRAVERSAL_BVH;
    out[i] = trace(in[i], shapes, nodes, tables, rt, linear);
}

struct TriangleProbe {
    Shape shape;
    Ray ray;  // world space
};
static_assert(sizeof(TriangleProbe) == 432, "TriangleProbe must match probes.rs");

struct TriangleResult {
    float4 normal;  // world-space, flipped toward the eye as in color_at
    float2 uv;
    float t;
    uint count;
};
static_assert(sizeof(TriangleResult) == 32, "TriangleResult must match probes.rs");

// Chapter 15: hit count, t, (u, v) and the smooth normal of a world-space
// ray at one triangle. A miss leaves every field zero.
kernel void probe_triangle(device const TriangleProbe *in [[buffer(0)]],
                           device TriangleResult *out     [[buffer(1)]],
                           uint i                         [[thread_position_in_grid]])
{
    TriangleProbe p = in[i];
    float4 ts = float4(0.0);
    uint count = intersect(p.shape, p.ray, ts);
    if (count == 0) {
        out[i] = TriangleResult{float4(0.0), float2(0.0), 0.0, 0};
        return;
    }
    float2 uv = triangle_uv(p.shape, ray_transform(p.shape.inverse, p.ray));
    float4 point = p.ray.origin + p.ray.direction * ts.x;
    float4 normalv = normal_at(p.shape, point, uv);
    if (dot(normalv, -p.ray.direction) < 0.0) {
        normalv = -normalv;
    }
    out[i] = TriangleResult{normalv, uv, ts.x, count};
}

struct BoxProbe {
    float4 min;
    float4 max;
    Ray ray;
    float tmin;
    float tmax;
};
static_assert(sizeof(BoxProbe) == 80, "BoxProbe must match probes.rs");

// out[i] = 1 when the ray crosses the box within [tmin, tmax], else 0
kernel void probe_intersect_box(device const BoxProbe *in [[buffer(0)]],
                                device uint *out          [[buffer(1)]],
                                uint i                    [[thread_position_in_grid]])
{
    BoxProbe p = in[i];
    out[i] = intersect_box(p.min, p.max, p.ray, p.tmin, p.tmax) ? 1 : 0;
}

struct RandomProbe {
    uint pixel;
    uint frame;
    uint dim;
};
static_assert(sizeof(RandomProbe) == 12, "RandomProbe must match probes.rs");

kernel void probe_random(device const RandomProbe *in [[buffer(0)]],
                         device float *out            [[buffer(1)]],
                         uint i                       [[thread_position_in_grid]])
{
    out[i] = random01(in[i].pixel, in[i].frame, in[i].dim);
}

// Milestone 10: the scene counts as a scene kernel sees them.
kernel void probe_counts(device const uint *in         [[buffer(0)]],
                         device uint4 *out             [[buffer(1)]],
                         constant SceneCounts &counts  [[buffer(4)]],
                         uint i                        [[thread_position_in_grid]])
{
    out[i] = uint4(counts.traversal, counts.triangle_count, counts.box_count, 0);
}
