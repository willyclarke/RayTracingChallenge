#include <metal_stdlib>
using namespace metal;

// Ray-tracer core shared by every kernel. Structs mirror src/types.rs.
// Points have w = 1 and vectors w = 0, as in rtc-rust's Tuple, so one
// float4x4 * float4 transforms both correctly.

struct Ray {
    float4 origin;
    float4 direction;
};
static_assert(sizeof(Ray) == 32, "Ray must match GpuRay in types.rs");

inline Ray ray_transform(float4x4 m, Ray r) {
    return Ray{m * r.origin, m * r.direction};
}

// Unit sphere at the origin, ray in object space (book chapter 5). Returns
// the hit count, 0 or 2; on a hit `ts` holds both t values, ascending.
inline uint intersect_sphere(Ray r, thread float2 &ts) {
    float3 sphere_to_ray = r.origin.xyz;
    float a = dot(r.direction.xyz, r.direction.xyz);
    float b = 2.0 * dot(r.direction.xyz, sphere_to_ray);
    float c = dot(sphere_to_ray, sphere_to_ray) - 1.0;
    float discriminant = b * b - 4.0 * a * c;
    if (discriminant < 0.0) {
        return 0;
    }
    float s = sqrt(discriminant);
    ts = float2((-b - s) / (2.0 * a), (-b + s) / (2.0 * a));
    return 2;
}
