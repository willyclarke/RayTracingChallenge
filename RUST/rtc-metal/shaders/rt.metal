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
    float reflective;       // chapter 11
    float transparency;
    float refractive_index;
};
static_assert(sizeof(Material) == 48, "Material must match GpuMaterial");

// A primitive of `kind` seen through its world-to-object transform. One
// tagged array holds every kind; functions switch on `kind`.
struct Shape {
    float4x4 inverse;
    Material material;
    uint kind;
};
static_assert(sizeof(Shape) == 128, "Shape must match GpuShape");

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

// The nearest non-negative hit: root `root` (0 or 1, of the ascending pair
// local_intersect returns) of shapes[index], at distance t.
struct Hit {
    float t;
    uint index;
    uint root;
};

// Nearest non-negative hit over all shapes. Returns false on a miss.
inline bool hit_shapes(Ray r, device const Shape *shapes, uint count, thread Hit &hit) {
    hit.t = INFINITY;
    for (uint i = 0; i < count; i++) {
        float2 ts;
        if (local_intersect(shapes[i].kind, ray_transform(shapes[i].inverse, r), ts) == 0) {
            continue;
        }
        // ts is ascending, so the first non-negative one is this shape's hit
        uint root = ts.x >= 0.0 ? 0 : 1;
        float t = ts[root];
        if (t >= 0.0 && t < hit.t) {
            hit = Hit{t, i, root};
        }
    }
    return hit.t < INFINITY;
}

// ---- Chapter 8: shadows ---------------------------------------------------

// rtc-rust uses 1e-9 in f64. In f32 a point ~10 units from the origin is
// only resolved to ~1e-6, so the shadow-acne offset has to be much larger.
// But not too large: refraction magnifies the offset of each under_point
// into visible edge shifts. At 1e-4, the chapter 11 test scene's glass ball
// put 0.23% of pixels off the CPU image; at 1e-5, 0.02%. Scenes much larger
// than ~10 units may need a larger value or an offset that scales with |point|.
constant float EPSILON = 1e-5;

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

// ---- Chapter 11: reflection and refraction -------------------------------

// Schlick's approximation to the Fresnel reflectance: the fraction of light
// reflected at a boundary from index n1 into n2, at cos(angle) = eyev . normalv.
inline float schlick(float cos_i, float n1, float n2) {
    float cos = cos_i;
    // total internal reflection can only occur if n1 > n2
    if (n1 > n2) {
        float n = n1 / n2;
        float sin2_t = n * n * (1.0 - cos * cos);
        if (sin2_t > 1.0) {
            return 1.0;
        }
        cos = sqrt(1.0 - sin2_t);
    }
    float r0 = (n1 - n2) / (n1 + n2);
    r0 *= r0;
    return r0 + (1.0 - r0) * pow(1.0 - cos, 5.0);
}

// (n1, n2): refractive indices on the near and far side of `hit`.
//
// rtc-rust walks the t-sorted intersection list keeping a list of the shapes
// the ray is inside, in entry order. With no sorted list here, the same state
// is rebuilt per shape: the ray is inside a shape just before the hit when an
// odd number of its roots (negative ones too) come first, and among those
// the innermost - the end of rtc-rust's list - is the one entered last.
inline float2 refractive_indices(Ray r, Hit hit, device const Shape *shapes, uint count) {
    float inner_t = -INFINITY;  // entry t of the innermost other medium
    float inner_n = 1.0;        // its refractive index; vacuum when none
    float hit_entry_t = -INFINITY;
    for (uint i = 0; i < count; i++) {
        float2 ts;
        uint n = local_intersect(shapes[i].kind, ray_transform(shapes[i].inverse, r), ts);
        if (i == hit.index) {
            // the hit's own roots before it: compare by index, not t, so the
            // hit root itself never counts
            if (hit.root == 1) {
                hit_entry_t = ts.x;
            }
            continue;
        }
        uint before = 0;
        for (uint k = 0; k < n; k++) {
            if (ts[k] < hit.t) {
                before++;
            }
        }
        // ts is ascending, so the last root before the hit is ts[before - 1]
        if (before % 2 == 1 && ts[before - 1] > inner_t) {
            inner_t = ts[before - 1];
            inner_n = shapes[i].material.refractive_index;
        }
    }
    float hit_n = shapes[hit.index].material.refractive_index;
    if (hit.root == 1) {
        // leaving the hit shape: it is n1 only if it is the innermost
        return float2(hit_entry_t > inner_t ? hit_n : inner_n, inner_n);
    }
    return float2(inner_n, hit_n); // entering it
}

// rtc-rust's render_parallel gives every camera ray remaining = 10.
constant int MAX_DEPTH = 10;

// A ray still to be traced, and how much its color counts in the pixel.
struct PendingRay {
    Ray ray;
    float3 weight;
    int remaining;
};

// rtc-rust's recursive color_at / shade_hit / reflected_color /
// refracted_color, as a loop: MSL has no recursion. The pixel is a sum over
// the tree of rays, each hit adding weight * its surface color, so pending
// rays go on a stack with their weights. A ray pushes at most two children
// and only while remaining > 0, so the stack never holds more than
// MAX_DEPTH + 1 rays (one pending sibling per level, plus the current one).
inline float3 color_at(Ray camera_ray, constant SceneInfo &scene, device const Shape *shapes) {
    PendingRay stack[MAX_DEPTH + 1];
    int top = 0;
    stack[top++] = PendingRay{camera_ray, float3(1.0), MAX_DEPTH};
    float3 color = float3(0.0);
    while (top > 0) {
        PendingRay p = stack[--top];
        Hit hit;
        if (!hit_shapes(p.ray, shapes, scene.shape_count, hit)) {
            continue;
        }
        Shape s = shapes[hit.index];
        Material m = s.material;
        float4 point = p.ray.origin + p.ray.direction * hit.t;
        float4 eyev = -p.ray.direction;
        float4 normalv = normal_at(s, point);
        if (dot(normalv, eyev) < 0.0) {
            normalv = -normalv; // hit from inside
        }
        float4 over_point = point + normalv * EPSILON;
        float intensity =
            is_shadowed(scene.light.position, over_point, shapes, scene.shape_count) ? 0.0 : 1.0;
        color += p.weight * lighting(m, scene.light, over_point, eyev, normalv, intensity);

        if (p.remaining <= 0) {
            continue;
        }
        float reflected = m.reflective;
        if (m.transparency > 0.0) {
            float2 n = refractive_indices(p.ray, hit, shapes, scene.shape_count);
            float cos_i = dot(eyev, normalv);
            float refracted = m.transparency;
            if (m.reflective > 0.0) {
                float reflectance = schlick(cos_i, n.x, n.y);
                reflected *= reflectance;
                refracted *= 1.0 - reflectance;
            }
            // Snell's law; nothing refracts on total internal reflection
            float n_ratio = n.x / n.y;
            float sin2_t = n_ratio * n_ratio * (1.0 - cos_i * cos_i);
            if (sin2_t <= 1.0) {
                float cos_t = sqrt(1.0 - sin2_t);
                float4 direction = normalv * (n_ratio * cos_i - cos_t) - eyev * n_ratio;
                float4 under_point = point - normalv * EPSILON;
                stack[top++] = PendingRay{Ray{under_point, direction}, p.weight * refracted,
                                          p.remaining - 1};
            }
        }
        if (m.reflective > 0.0) {
            Ray reflect_ray = Ray{over_point, reflect(p.ray.direction, normalv)};
            stack[top++] = PendingRay{reflect_ray, p.weight * reflected, p.remaining - 1};
        }
    }
    return color;
}
