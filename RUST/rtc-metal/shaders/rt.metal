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
constant uint SHAPE_CUBE = 2;
constant uint SHAPE_CYLINDER = 3;
constant uint SHAPE_CONE = 4;
constant uint SHAPE_TRIANGLE = 5;

// rtc-rust uses 1e-9 in f64. In f32 a point ~10 units from the origin is
// only resolved to ~1e-6, so the shadow-acne offset has to be much larger.
// But not too large: refraction magnifies the offset of each under_point
// into visible edge shifts. At 1e-4, the chapter 11 test scene's glass ball
// put 0.23% of pixels off the CPU image; at 1e-5, 0.02%. Scenes much larger
// than ~10 units may need a larger value or an offset that scales with |point|.
// The cylinder and cone normals use it too, to tell a cap point from a wall
// point (rtc-rust's 1e-9 is below f32 resolution there as well).
constant float EPSILON = 1e-5;

// Pattern kinds (book chapter 10); must match the PATTERN_* constants in
// src/types.rs. PATTERN_NONE uses the material's plain color.
constant uint PATTERN_NONE = 0;
constant uint PATTERN_STRIPE = 1;
constant uint PATTERN_GRADIENT = 2;
constant uint PATTERN_RING = 3;
constant uint PATTERN_CHECKERS = 4;

struct Material {
    float4 color;
    float ambient;
    float diffuse;
    float specular;
    float shininess;
    float reflective;       // chapter 11
    float transparency;
    float refractive_index;
    uint pattern;           // chapter 10: PATTERN_*
    float4x4 pattern_inverse; // object space to pattern space
    float4 pattern_a;       // the pattern's two colors
    float4 pattern_b;
};
static_assert(sizeof(Material) == 144, "Material must match GpuMaterial");

// A primitive of `kind` seen through its world-to-object transform. One
// tagged array holds every kind; functions switch on `kind`.
struct Shape {
    float4x4 inverse;
    Material material;
    uint kind;
    float minimum;  // cylinder and cone: y extent (book chapter 13)
    float maximum;
    uint closed;    // cylinder and cone: nonzero for end caps
    // Triangle (book chapter 15), object space: vertex 1, the edges to
    // vertices 2 and 3, and the vertex normals.
    float4 p1;
    float4 e1;
    float4 e2;
    float4 n1;
    float4 n2;
    float4 n3;
};
static_assert(sizeof(Shape) == 320, "Shape must match GpuShape");

struct Light {
    float4 position;
    float4 intensity;
};
static_assert(sizeof(Light) == 32, "Light must match GpuLight");

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

// Insert t into the ascending ts[0..n). For primitives whose roots come
// out of order (cylinder and cone walls and caps); at most 4 roots.
inline void push_sorted(thread float4 &ts, thread uint &n, float t) {
    uint k = n++;
    for (; k > 0 && ts[k - 1] > t; k--) {
        ts[k] = ts[k - 1];
    }
    ts[k] = t;
}

// The unit cube, ray in object space (book chapter 12): the slab method,
// intersecting the ray's t intervals between x, y and z = -1 and 1. Returns
// 0 or 2 hits, ascending.
inline uint intersect_cube(Ray r, thread float2 &ts) {
    float3 t_lo = (-1.0 - r.origin.xyz) / r.direction.xyz;
    float3 t_hi = (1.0 - r.origin.xyz) / r.direction.xyz;
    float3 t_near = min(t_lo, t_hi);
    float3 t_far = max(t_lo, t_hi);
    float tmin = max3(t_near.x, t_near.y, t_near.z);
    float tmax = min3(t_far.x, t_far.y, t_far.z);
    if (tmin > tmax) {
        return 0;
    }
    ts = float2(tmin, tmax);
    return 2;
}

// Caps of a closed cylinder (radius 1) or cone (radius |y|) at y = minimum
// and y = maximum: discs on those planes (book chapter 13).
inline void intersect_caps(Shape s, Ray r, bool cone, thread float4 &ts, thread uint &n) {
    // rtc-rust's approx_eq(direction.y, 0) threshold
    if (s.closed == 0 || abs(r.direction.y) < 1e-9) {
        return;
    }
    float caps[2] = {s.minimum, s.maximum};
    for (float y : caps) {
        float t = (y - r.origin.y) / r.direction.y;
        float x = r.origin.x + t * r.direction.x;
        float z = r.origin.z + t * r.direction.z;
        if (x * x + z * z <= (cone ? y * y : 1.0)) {
            push_sorted(ts, n, t);
        }
    }
}

// Wall of the cylinder x^2 + z^2 = 1 or the double cone x^2 + z^2 = y^2
// (book chapter 13), clipped to minimum < y < maximum. Both are the
// quadratic a t^2 + 2 b' t + c = 0 with
//   cylinder: a = dx^2 + dz^2,         b' = ox dx + oz dz,         c = ox^2 + oz^2 - 1
//   cone:     a = dx^2 - dy^2 + dz^2,  b' = ox dx - oy dy + oz dz, c = ox^2 - oy^2 + oz^2
//
// As for the sphere (Ray Tracing Gems ch. 7), the textbook b'^2 - ac fails
// in f32 when the object-space origin is far out: both terms are huge and
// nearly equal. A cone scaled 0.5 had t off by ~3e-5 at ~17 units, more
// than EPSILON, so over_point fell inside and the cone shadowed itself.
// Instead, by Lagrange's identity, with w = o x d (small for a ray aimed at
// the shape, however far out o is):
//   cylinder: b'^2 - ac = a - w.y^2
//   cone:     b'^2 - ac = w.x^2 - w.y^2 + w.z^2
// The roots then come from q = -(b' + sign(b') sqrt(disc)) as q/a and c/q,
// so neither subtracts nearly equal numbers; that also keeps the near root
// of a ray almost parallel to a cone's side, where a ~ 0.
//
// The cone's dz^2 - dy^2 is written (dz - dy)(dz + dy): exactly 0 when
// |dy| = |dz|, and accurate near it. As a sum of squares the compiler fuses
// it into an fma that leaves a ~ 3e-8 for the book's parallel ray.
inline void intersect_wall(Shape s, Ray r, bool cone, thread float4 &ts, thread uint &n) {
    float3 o = r.origin.xyz;
    float3 d = r.direction.xyz;
    float3 w = cross(o, d);
    float a, b_half, c, discriminant;
    if (cone) {
        a = d.x * d.x + (d.z - d.y) * (d.z + d.y);
        b_half = o.x * d.x - o.y * d.y + o.z * d.z;
        c = o.x * o.x + (o.z - o.y) * (o.z + o.y);
        discriminant = w.x * w.x + (w.z - w.y) * (w.z + w.y);
    } else {
        a = d.x * d.x + d.z * d.z;
        b_half = o.x * d.x + o.z * d.z;
        c = o.x * o.x + o.z * o.z - 1.0;
        discriminant = a - w.y * w.y;
    }
    float t[2];
    uint roots;
    // rtc-rust's approx_eq(a, 0) threshold
    if (abs(a) < 1e-9) {
        // cylinder: parallel to the y axis, so no wall hits. cone: parallel
        // to one half, so a single root. rtc-rust copies the book's
        // t = -c / (2b) with b = 2b'; the root of 2b' t + c = 0 is -c / 2b'
        // (a known erratum of the book), but the CPU is the reference here.
        float b = 2.0 * b_half;
        if (!cone || abs(b) < 1e-9) {
            return;
        }
        t[0] = -c / (2.0 * b);
        roots = 1;
    } else {
        if (discriminant < 0.0) {
            return;
        }
        float q = -(b_half + copysign(sqrt(discriminant), b_half));
        t[0] = q / a;
        t[1] = c / q;
        roots = 2;
    }
    for (uint k = 0; k < roots; k++) {
        float y = o.y + t[k] * d.y;
        if (s.minimum < y && y < s.maximum) {
            push_sorted(ts, n, t[k]);
        }
    }
}

// Triangle, ray in object space (book chapter 15): Moller-Trumbore. Returns
// false on a miss; on a hit `tuv` is (t, u, v), where u and v are the
// barycentric weights of vertices 2 and 3.
inline bool moller_trumbore(Shape s, Ray r, thread float3 &tuv) {
    float3 dir_cross_e2 = cross(r.direction.xyz, s.e2.xyz);
    float det = dot(s.e1.xyz, dir_cross_e2);
    // Only an exactly parallel ray (or collinear vertices) divides by zero.
    // The u/v range checks reject near-parallel rays, and a fixed threshold
    // would drop hits on tiny mesh triangles, whose det scales with |e1||e2|.
    if (det == 0.0) {
        return false;
    }
    float f = 1.0 / det;
    float3 p1_to_origin = r.origin.xyz - s.p1.xyz;
    float u = f * dot(p1_to_origin, dir_cross_e2);
    if (u < 0.0 || u > 1.0) {
        return false;
    }
    float3 origin_cross_e1 = cross(p1_to_origin, s.e1.xyz);
    float v = f * dot(r.direction.xyz, origin_cross_e1);
    if (v < 0.0 || u + v > 1.0) {
        return false;
    }
    tuv = float3(f * dot(s.e2.xyz, origin_cross_e1), u, v);
    return true;
}

// 0 or 1 roots; the hit's t goes in ts.x.
inline uint intersect_triangle(Shape s, Ray r, thread float2 &ts) {
    float3 tuv;
    if (!moller_trumbore(s, r, tuv)) {
        return 0;
    }
    ts = float2(tuv.x, 0.0);
    return 1;
}

// The (u, v) of the object-space ray's hit, for the smooth normal; zero
// on a miss.
inline float2 triangle_uv(Shape s, Ray r) {
    float3 tuv = float3(0.0);
    moller_trumbore(s, r, tuv);
    return tuv.yz;
}

// Dispatch on the primitive kind: up to 4 hits in object space, ascending
// in ts[0..count). Only a closed cone has more than 2 (both halves plus
// both caps).
inline uint local_intersect(Shape s, Ray r, thread float4 &ts) {
    float2 pair;
    uint n = 0;
    switch (s.kind) {
    case SHAPE_PLANE:
        n = intersect_plane(r, pair);
        break;
    case SHAPE_CUBE:
        n = intersect_cube(r, pair);
        break;
    case SHAPE_TRIANGLE:
        n = intersect_triangle(s, r, pair);
        break;
    case SHAPE_CYLINDER:
    case SHAPE_CONE:
        intersect_caps(s, r, s.kind == SHAPE_CONE, ts, n);
        intersect_wall(s, r, s.kind == SHAPE_CONE, ts, n);
        return n;
    default: // SHAPE_SPHERE
        n = intersect_sphere(r, pair);
        break;
    }
    ts.xy = pair;
    return n;
}

// Intersect s with a world-space ray.
inline uint intersect(Shape s, Ray r, thread float4 &ts) {
    return local_intersect(s, ray_transform(s.inverse, r), ts);
}

// ---- Chapter 6: normals and Phong lighting --------------------------------


// Object-space normal of a primitive at `object_point`. `uv` is the hit's
// barycentric pair, used only by triangles.
inline float4 local_normal_at(Shape s, float4 object_point, float2 uv) {
    float4 p = object_point;
    switch (s.kind) {
    case SHAPE_PLANE:
        return float4(0.0, 1.0, 0.0, 0.0); // constant (book chapter 9)
    case SHAPE_TRIANGLE:
        // interpolated vertex normals (flat triangles have three equal ones)
        return s.n2 * uv.x + s.n3 * uv.y + s.n1 * (1.0 - uv.x - uv.y);
    case SHAPE_CUBE: {
        // the face of the largest coordinate (book chapter 12)
        float3 a = abs(p.xyz);
        float m = max3(a.x, a.y, a.z);
        if (m == a.x) {
            return float4(p.x, 0.0, 0.0, 0.0);
        }
        if (m == a.y) {
            return float4(0.0, p.y, 0.0, 0.0);
        }
        return float4(0.0, 0.0, p.z, 0.0);
    }
    case SHAPE_CYLINDER:
    case SHAPE_CONE: {
        // caps, else the wall (book chapter 13)
        bool cone = s.kind == SHAPE_CONE;
        float dist = p.x * p.x + p.z * p.z;
        if (dist < (cone ? s.maximum * s.maximum : 1.0) && p.y >= s.maximum - EPSILON) {
            return float4(0.0, 1.0, 0.0, 0.0);
        }
        if (dist < (cone ? s.minimum * s.minimum : 1.0) && p.y <= s.minimum + EPSILON) {
            return float4(0.0, -1.0, 0.0, 0.0);
        }
        float y = 0.0;
        if (cone) {
            y = p.y > 0.0 ? -sqrt(dist) : sqrt(dist);
        }
        return float4(p.x, y, p.z, 0.0);
    }
    default: // SHAPE_SPHERE
        return p - float4(0.0, 0.0, 0.0, 1.0);
    }
}

// World-space normal: to object space, local normal, back with the inverse
// transpose (book chapter 6).
inline float4 normal_at(Shape s, float4 world_point, float2 uv) {
    float4 object_normal = local_normal_at(s, s.inverse * world_point, uv);
    float4 world_normal = transpose(s.inverse) * object_normal;
    world_normal.w = 0.0;
    return normalize(world_normal);
}

// ---- Chapter 10: patterns -------------------------------------------------

// Parity of floor(x): rtc-rust's `floor() as i32 % 2 == 0`, where a negative
// odd value leaves -1, not 0, so it is odd too.
inline bool floor_is_even(float x) {
    return (int(floor(x)) & 1) == 0;
}

// The pattern's color at a point in the shape's object space; rtc-rust's
// Pattern::color_at_local.
inline float3 pattern_at(Material m, float4 object_point) {
    float4 p = m.pattern_inverse * object_point;
    float3 a = m.pattern_a.rgb;
    float3 b = m.pattern_b.rgb;
    switch (m.pattern) {
    case PATTERN_STRIPE:
        return floor_is_even(p.x) ? a : b;
    case PATTERN_GRADIENT:
        return a + (b - a) * (p.x - floor(p.x));
    case PATTERN_RING:
        return floor_is_even(sqrt(p.x * p.x + p.z * p.z)) ? a : b;
    case PATTERN_CHECKERS:
        return floor_is_even(floor(p.x) + floor(p.y) + floor(p.z)) ? a : b;
    default: // PATTERN_NONE
        return m.color.rgb;
    }
}

// Mirrors rtc-rust's Light::lighting for a point light (one sample point).
// `object_point` is `point` in the shape's object space, for the pattern.
// `intensity` is the fraction of the light reaching `point`: 0 in shadow.
inline float3 lighting(Material m, Light light, float4 point, float4 object_point,
                       float4 eyev, float4 normalv, float intensity) {
    float3 effective_color = pattern_at(m, object_point) * light.intensity.rgb;
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

// The nearest non-negative hit: root `root` (0 to 3, of the ascending roots
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
        float4 ts;
        uint n = intersect(shapes[i], r, ts);
        // ts is ascending, so the first non-negative one is this shape's hit
        uint root = 0;
        while (root < n && ts[root] < 0.0) {
            root++;
        }
        if (root < n && ts[root] < hit.t) {
            hit = Hit{ts[root], i, root};
        }
    }
    return hit.t < INFINITY;
}

// ---- Chapter 8: shadows ---------------------------------------------------

// Any-hit test: is anything between `point` and the light?
inline bool is_shadowed(float4 light_position, float4 point,
                        device const Shape *shapes, uint count) {
    float4 v = light_position - point;
    float distance = length(v);
    Ray r = Ray{point, v / distance};
    for (uint i = 0; i < count; i++) {
        float4 ts;
        uint n = intersect(shapes[i], r, ts);
        for (uint k = 0; k < n; k++) {
            if (ts[k] >= 0.0 && ts[k] < distance) {
                return true;
            }
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
        float4 ts;
        uint n = intersect(shapes[i], r, ts);
        if (i == hit.index) {
            // the hit's own roots before it: compare by index, not t, so the
            // hit root itself never counts. An odd root leaves the shape,
            // and the root before it entered.
            if (hit.root % 2 == 1) {
                hit_entry_t = ts[hit.root - 1];
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
    if (hit.root % 2 == 1) {
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
        float2 uv = s.kind == SHAPE_TRIANGLE
                        ? triangle_uv(s, ray_transform(s.inverse, p.ray))
                        : float2(0.0);
        float4 normalv = normal_at(s, point, uv);
        if (dot(normalv, eyev) < 0.0) {
            normalv = -normalv; // hit from inside
        }
        float4 over_point = point + normalv * EPSILON;
        float intensity =
            is_shadowed(scene.light.position, over_point, shapes, scene.shape_count) ? 0.0 : 1.0;
        // rtc-rust samples the pattern at over_point too: on a plane the raw
        // hit's y is ~0 with either sign, and checkers would speckle
        float4 object_point = s.inverse * over_point;
        color += p.weight
                 * lighting(m, scene.light, over_point, object_point, eyev, normalv, intensity);

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
