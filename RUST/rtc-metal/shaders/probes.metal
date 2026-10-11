
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
static_assert(sizeof(LightingProbe) == 128, "LightingProbe must match probes.rs");

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
