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

// Shape kinds; must match the SHAPE_* constants in src/types.rs.
constant uint SHAPE_SPHERE = 0;
constant uint SHAPE_PLANE = 1;

inline Ray ray_transform(float4x4 m, Ray r) {
    return Ray{m * r.origin, m * r.direction};
}

// Unit sphere at the origin, ray in object space (book chapter 5). Returns
// the hit count, 0 or 2; on a hit `ts` holds both t values, ascending.
//
// Not the book's b^2 - 4ac: a flattened sphere (scaled 0.01) puts the
// object-space origin hundreds of units out, where b^2 and 4ac are huge and
// nearly equal and f32 loses most digits of their difference. Instead
// (Ray Tracing Gems ch. 7, "Precision Improvements for Ray/Sphere
// Intersection"):
//   - the discriminant comes from the closest-approach vector f, computed
//     directly: b'^2 - ac = a * (1 - |f|^2), with b' = b/2;
//   - the roots come from q = -(b' + sign(b') sqrt(disc)) as c/q and q/a,
//     so neither subtracts two nearly equal numbers.
inline uint intersect_sphere(Ray r, thread float2 &ts) {
    float3 o = r.origin.xyz;
    float3 d = r.direction.xyz;
    float a = dot(d, d);
    float b_half = dot(o, d);
    float3 f = o - (b_half / a) * d;
    float discriminant = a * (1.0 - dot(f, f));
    if (discriminant < 0.0) {
        return 0;
    }
    float c = dot(o, o) - 1.0;
    float q = -(b_half + copysign(sqrt(discriminant), b_half));
    float t0 = c / q;
    float t1 = q / a;
    ts = t0 < t1 ? float2(t0, t1) : float2(t1, t0);
    return 2;
}

// The xz plane, ray in object space (book chapter 9). Returns the hit count,
// 0 or 1; on a hit both `ts` slots hold the one t, so callers can treat
// every primitive as "ascending pair, take the first non-negative".
inline uint intersect_plane(Ray r, thread float2 &ts) {
    // rtc-rust's parallel test, same threshold. f32 represents 1e-9 fine;
    // it only has to keep t = -y / 0 out.
    if (abs(r.direction.y) < 1e-9) {
        return 0;
    }
    ts = float2(-r.origin.y / r.direction.y);
    return 1;
}

// Dispatch on the primitive kind; hits in object space, ascending.
inline uint local_intersect(uint kind, Ray r, thread float2 &ts) {
    switch (kind) {
    case SHAPE_PLANE:
        return intersect_plane(r, ts);
    default: // SHAPE_SPHERE
        return intersect_sphere(r, ts);
    }
}

// ---- Chapter 6: normals and Phong lighting --------------------------------

struct Material {
    float4 color;
    float ambient;
    float diffuse;
    float specular;
    float shininess;
};
static_assert(sizeof(Material) == 32, "Material must match GpuMaterial");

// A primitive of `kind` seen through its world-to-object transform. One
// tagged array holds every kind; functions switch on `kind`.
struct Shape {
    float4x4 inverse;
    Material material;
    uint kind;
};
static_assert(sizeof(Shape) == 112, "Shape must match GpuShape");

struct Light {
    float4 position;
    float4 intensity;
};
static_assert(sizeof(Light) == 32, "Light must match GpuLight");

// Object-space normal of a primitive at `object_point`.
inline float4 local_normal_at(uint kind, float4 object_point) {
    switch (kind) {
    case SHAPE_PLANE:
        return float4(0.0, 1.0, 0.0, 0.0); // constant (book chapter 9)
    default: // SHAPE_SPHERE
        return object_point - float4(0.0, 0.0, 0.0, 1.0);
    }
}

// World-space normal: to object space, local normal, back with the inverse
// transpose (book chapter 6).
inline float4 normal_at(Shape s, float4 world_point) {
    float4 object_normal = local_normal_at(s.kind, s.inverse * world_point);
    float4 world_normal = transpose(s.inverse) * object_normal;
    world_normal.w = 0.0;
    return normalize(world_normal);
}

// Mirrors rtc-rust's Light::lighting for a point light (one sample point).
// `intensity` is the fraction of the light reaching `point`: 0 in shadow.
inline float3 lighting(Material m, Light light, float4 point, float4 eyev,
                       float4 normalv, float intensity) {
    float3 effective_color = m.color.rgb * light.intensity.rgb;
    float3 ambient = effective_color * m.ambient;
    if (intensity <= 0.0) {
        return ambient;
    }
    float4 to_light = light.position - point;
    float light_dot_normal = dot(to_light, normalv) / length(to_light);
    float diffuse = 0.0;
    float specular = 0.0;
    if (light_dot_normal >= 0.0) {
        diffuse = light_dot_normal;
        if (m.specular > 0.0) {
            float4 reflectv = -reflect(normalize(to_light), normalv);
            float reflect_dot_eye = dot(reflectv, eyev);
            if (reflect_dot_eye > 0.0) {
                specular = pow(reflect_dot_eye, m.shininess);
            }
        }
    }
    return ambient + (effective_color * m.diffuse * diffuse
                      + light.intensity.rgb * m.specular * specular) * intensity;
}

// ---- Chapter 7: camera and world ------------------------------------------

struct Camera {
    float4x4 inverse;
    float half_width;
    float half_height;
    float pixel_size;
    uint width;
    uint height;
};
static_assert(sizeof(Camera) == 96, "Camera must match GpuCamera");

struct SceneInfo {
    Camera camera;
    Light light;
    uint shape_count;
};
static_assert(sizeof(SceneInfo) == 144, "SceneInfo must match GpuSceneInfo");

// Ray through the centre of pixel (px, py). The camera looks toward -z with
// the canvas at z = -1, so +x is to the left.
inline Ray ray_for_pixel(Camera c, uint px, uint py) {
    float world_x = c.half_width - (float(px) + 0.5) * c.pixel_size;
    float world_y = c.half_height - (float(py) + 0.5) * c.pixel_size;
    float4 pixel = c.inverse * float4(world_x, world_y, -1.0, 1.0);
    float4 origin = c.inverse * float4(0.0, 0.0, 0.0, 1.0);
    return Ray{origin, normalize(pixel - origin)};
}

// Nearest non-negative hit over all shapes. Returns false on a miss.
inline bool hit_shapes(Ray r, device const Shape *shapes, uint count,
                       thread float &t_hit, thread uint &hit_index) {
    t_hit = INFINITY;
    for (uint i = 0; i < count; i++) {
        float2 ts;
        if (local_intersect(shapes[i].kind, ray_transform(shapes[i].inverse, r), ts) == 0) {
            continue;
        }
        // ts is ascending, so the first non-negative one is this sphere's hit
        float t = ts.x >= 0.0 ? ts.x : ts.y;
        if (t >= 0.0 && t < t_hit) {
            t_hit = t;
            hit_index = i;
        }
    }
    return t_hit < INFINITY;
}

// ---- Chapter 8: shadows ---------------------------------------------------

// rtc-rust uses 1e-9 in f64. In f32 a point ~10 units from the origin is
// only resolved to ~1e-6, so the shadow-acne offset has to be much larger.
constant float EPSILON = 1e-4;

// Any-hit test: is anything between `point` and the light?
inline bool is_shadowed(float4 light_position, float4 point,
                        device const Shape *shapes, uint count) {
    float4 v = light_position - point;
    float distance = length(v);
    Ray r = Ray{point, v / distance};
    for (uint i = 0; i < count; i++) {
        float2 ts;
        if (local_intersect(shapes[i].kind, ray_transform(shapes[i].inverse, r), ts) != 0
            && ((ts.x >= 0.0 && ts.x < distance) || (ts.y >= 0.0 && ts.y < distance))) {
            return true;
        }
    }
    return false;
}

inline float3 color_at(Ray r, constant SceneInfo &scene, device const Shape *shapes) {
    float t;
    uint index;
    if (!hit_shapes(r, shapes, scene.shape_count, t, index)) {
        return float3(0.0);
    }
    Shape s = shapes[index];
    float4 point = r.origin + r.direction * t;
    float4 eyev = -r.direction;
    float4 normalv = normal_at(s, point);
    if (dot(normalv, eyev) < 0.0) {
        normalv = -normalv; // hit from inside
    }
    float4 over_point = point + normalv * EPSILON;
    float intensity =
        is_shadowed(scene.light.position, over_point, shapes, scene.shape_count) ? 0.0 : 1.0;
    return lighting(s.material, scene.light, over_point, eyev, normalv, intensity);
}
