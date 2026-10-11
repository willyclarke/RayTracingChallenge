
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
    out[i] = normal_at(in[i].shape, in[i].point);
}

struct IntersectProbe {
    Ray ray;
    uint kind;
};

// out[i] = (t0, t1, hit count, unused)
kernel void probe_local_intersect(device const IntersectProbe *in [[buffer(0)]],
                                  device float4 *out              [[buffer(1)]],
                                  uint i                          [[thread_position_in_grid]])
{
    float2 ts = float2(0.0);
    uint count = local_intersect(in[i].kind, in[i].ray, ts);
    out[i] = float4(ts, float(count), 0.0);
}

struct LightingProbe {
    Material material;
    Light light;
    float4 point;
    float4 eyev;
    float4 normalv;
    float intensity;
};
static_assert(sizeof(LightingProbe) == 144, "LightingProbe must match probes.rs");

kernel void probe_lighting(device const LightingProbe *in [[buffer(0)]],
                           device float4 *out             [[buffer(1)]],
                           uint i                         [[thread_position_in_grid]])
{
    LightingProbe p = in[i];
    out[i] = float4(lighting(p.material, p.light, p.point, p.eyev, p.normalv, p.intensity), 0.0);
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
    uint shape_count;
};
static_assert(sizeof(RefractiveProbe) == 48, "RefractiveProbe must match probes.rs");

// out[i] = (n1, n2) at root hit_root of shapes[hit_index] along ray
kernel void probe_refractive_indices(device const RefractiveProbe *in [[buffer(0)]],
                                     device float2 *out               [[buffer(1)]],
                                     device const Shape *shapes       [[buffer(2)]],
                                     uint i                           [[thread_position_in_grid]])
{
    RefractiveProbe p = in[i];
    Shape s = shapes[p.hit_index];
    float2 ts;
    local_intersect(s.kind, ray_transform(s.inverse, p.ray), ts);
    Hit hit = Hit{ts[p.hit_root], p.hit_index, p.hit_root};
    out[i] = refractive_indices(p.ray, hit, shapes, p.shape_count);
}

// in[i] = (cos_i, n1, n2, unused)
kernel void probe_schlick(device const float4 *in [[buffer(0)]],
                          device float *out       [[buffer(1)]],
                          uint i                  [[thread_position_in_grid]])
{
    out[i] = schlick(in[i].x, in[i].y, in[i].z);
}
