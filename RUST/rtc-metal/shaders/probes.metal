
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
static_assert(sizeof(IntersectProbe) == 352, "IntersectProbe must match probes.rs");

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
};
static_assert(sizeof(LightingProbe) == 256, "LightingProbe must match probes.rs");

kernel void probe_lighting(device const LightingProbe *in [[buffer(0)]],
                           device float4 *out             [[buffer(1)]],
                           uint i                         [[thread_position_in_grid]])
{
    LightingProbe p = in[i];
    out[i] = float4(lighting(p.material, p.light, p.point, p.object_point, p.eyev, p.normalv, p.intensity), 0.0);
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
    out[i] = ray_for_pixel(in[i].camera, in[i].pixel.x, in[i].pixel.y);
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
                                     uint i                           [[thread_position_in_grid]])
{
    RefractiveProbe p = in[i];
    Shape s = shapes[p.hit_index];
    float4 ts;
    intersect(s, p.ray, ts);
    Hit hit = Hit{ts[p.hit_root], p.hit_index, p.hit_root};
    out[i] = refractive_indices(p.ray, hit, shapes, nodes, counts);
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
    out[i] = float4(pattern_at(s.material, s.inverse * in[i].point), 0.0);
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
// index = counts.shape_count on a miss.
kernel void probe_trace(device const TraceProbe *in  [[buffer(0)]],
                        device TraceResult *out      [[buffer(1)]],
                        device const Shape *shapes   [[buffer(2)]],
                        device const BvhNode *nodes  [[buffer(3)]],
                        constant SceneCounts &counts [[buffer(4)]],
                        uint i                       [[thread_position_in_grid]])
{
    TraceProbe p = in[i];
    Hit hit;
    if (!hit_shapes(p.ray, shapes, nodes, counts, hit)) {
        out[i] = TraceResult{float4(0.0), INFINITY, counts.shape_count, 0, 0};
        return;
    }
    float4 point = p.ray.origin + p.ray.direction * hit.t;
    Shape s = shapes[hit.index];
    float2 uv = s.kind == SHAPE_TRIANGLE ? triangle_uv(s, ray_transform(s.inverse, p.ray))
                                         : float2(0.0);
    float4 normalv = normal_at(s, point, uv);
    if (dot(normalv, -p.ray.direction) < 0.0) {
        normalv = -normalv;
    }
    bool shadowed = is_shadowed(p.light_position, point + normalv * EPSILON, shapes, nodes, counts);
    out[i] = TraceResult{normalv, hit.t, hit.index, hit.root, shadowed ? 1u : 0u};
}

struct TriangleProbe {
    Shape shape;
    Ray ray;  // world space
};
static_assert(sizeof(TriangleProbe) == 352, "TriangleProbe must match probes.rs");

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
