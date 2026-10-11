//! Probe tests: run single `rt.metal` functions on the GPU, one thread per
//! case, and compare with the rtc-rust CPU version of the same function.

use crate::gpu::Gpu;
use crate::scene::{gpu_camera, gpu_light, gpu_material};
use crate::types::{
    Float4, GpuCamera, GpuLight, GpuMat4, GpuMaterial, GpuRay, GpuShape, SHAPE_PLANE, SHAPE_SPHERE,
    float4,
};
use rtc_rust::camera::Camera;
use rtc_rust::light::Light;
use rtc_rust::material::Material;
use rtc_rust::matrix::Matrix4;
use rtc_rust::ray::Ray;
use rtc_rust::shape::Shape;
use rtc_rust::shapes::plane::Plane;
use rtc_rust::shapes::sphere::Sphere;
use rtc_rust::tuple::Tuple;
use rtc_rust::world::{World, prepare_computations, schlick};

/// Relative tolerance for f32 GPU results against f64 CPU results.
const EPSILON: f64 = 1e-5;

fn check_close(gpu: f32, cpu: f64, what: &str) -> Result<(), String> {
    if (gpu as f64 - cpu).abs() <= EPSILON * cpu.abs().max(1.0) {
        Ok(())
    } else {
        Err(format!("{what}: gpu {gpu} vs cpu {cpu}"))
    }
}

fn check_tuple(gpu: Float4, cpu: Tuple, what: &str) -> Result<(), String> {
    for (g, c) in gpu.iter().zip([cpu.x, cpu.y, cpu.z, cpu.w]) {
        check_close(*g, c, what)?;
    }
    Ok(())
}

/// Input to `probe_ray_transform`.
#[repr(C)]
#[derive(Clone, Copy)]
struct TransformProbe {
    m: GpuMat4,
    ray: GpuRay,
}

#[test]
fn test_chap_5_gpu_ray_transform() -> Result<(), String> {
    let ray = Ray::new(Tuple::point(1.0, 2.0, 3.0), Tuple::vector(0.0, 1.0, 0.0));
    let matrices = [
        Matrix4::translation(3.0, 4.0, 5.0), // book: translating a ray
        Matrix4::scaling(2.0, 3.0, 4.0),     // book: scaling a ray
        // every element non-trivial, so a transposed upload cannot pass
        Matrix4::translation(1.0, -2.0, 0.5)
            * Matrix4::rotation_y(0.7)
            * Matrix4::shearing(1.0, 0.0, 0.5, 0.0, 0.0, 0.25),
    ];
    let cases: Vec<TransformProbe> = matrices
        .iter()
        .map(|m| TransformProbe {
            m: m.into(),
            ray: (&ray).into(),
        })
        .collect();

    let out: Vec<GpuRay> = Gpu::new()?.run("probe_ray_transform", &cases)?;

    for (i, (m, got)) in matrices.iter().zip(&out).enumerate() {
        let want = ray.transform(*m);
        check_tuple(got.origin, want.origin, &format!("case {i} origin"))?;
        check_tuple(
            got.direction,
            want.direction,
            &format!("case {i} direction"),
        )?;
    }
    Ok(())
}

#[test]
fn test_chap_5_gpu_intersect_sphere() -> Result<(), String> {
    // The book's five cases: two hits, tangent, miss, inside, behind.
    let mut rays: Vec<Ray> = [
        (0.0, 0.0, -5.0),
        (0.0, 1.0, -5.0),
        (0.0, 2.0, -5.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 5.0),
    ]
    .iter()
    .map(|&(x, y, z)| Ray::new(Tuple::point(x, y, z), Tuple::vector(0.0, 0.0, 1.0)))
    .collect();
    // A 16x16 fan of unnormalized rays aimed off-centre: hits and misses,
    // with direction lengths != 1 so the `a` term matters.
    for i in 0..256 {
        let origin = Tuple::point(
            -3.0 + 6.0 * (i % 16) as f64 / 15.0,
            -3.0 + 6.0 * (i / 16) as f64 / 15.0,
            -5.0,
        );
        rays.push(Ray::new(origin, Tuple::point(0.3, -0.2, 0.0) - origin));
    }
    let cases: Vec<GpuRay> = rays.iter().map(GpuRay::from).collect();

    // Each result is (t0, t1, hit count, unused).
    let out: Vec<Float4> = Gpu::new()?.run("probe_intersect_sphere", &cases)?;

    let sphere = Sphere::new();
    for (i, (ray, got)) in rays.iter().zip(&out).enumerate() {
        let xs = sphere.local_intersect(ray);
        if got[2] as usize != xs.count() {
            return Err(format!(
                "ray {i}: gpu {} hits vs cpu {}",
                got[2],
                xs.count()
            ));
        }
        if !xs.is_empty() {
            check_close(got[0], xs[0].t, &format!("ray {i} t0"))?;
            check_close(got[1], xs[1].t, &format!("ray {i} t1"))?;
        }
    }
    Ok(())
}

/// Input to `probe_normal_at`.
#[repr(C)]
#[derive(Clone, Copy)]
struct NormalProbe {
    shape: GpuShape,
    point: Float4,
}

/// A `GpuShape` of `kind` with transform `m` (material irrelevant here).
fn probe_shape(kind: u32, m: &Matrix4) -> GpuShape {
    GpuShape {
        inverse: (&m.inverse().unwrap()).into(),
        material: gpu_material(&Material::new()).unwrap(),
        kind,
    }
}

#[test]
fn test_chap_6_gpu_sphere_normal_at() -> Result<(), String> {
    use std::f64::consts::{FRAC_1_SQRT_2, PI};
    let s3 = 3f64.sqrt() / 3.0;
    let cases = [
        (Matrix4::identity(), Tuple::point(1.0, 0.0, 0.0)),
        (Matrix4::identity(), Tuple::point(0.0, 1.0, 0.0)),
        (Matrix4::identity(), Tuple::point(0.0, 0.0, 1.0)),
        (Matrix4::identity(), Tuple::point(s3, s3, s3)),
        (
            Matrix4::translation(0.0, 1.0, 0.0),
            Tuple::point(0.0, 1.0 + FRAC_1_SQRT_2, -FRAC_1_SQRT_2),
        ),
        (
            Matrix4::scaling(1.0, 0.5, 1.0) * Matrix4::rotation_z(PI / 5.0),
            Tuple::point(0.0, FRAC_1_SQRT_2, -FRAC_1_SQRT_2),
        ),
    ];
    let input: Vec<NormalProbe> = cases
        .iter()
        .map(|(m, p)| NormalProbe {
            shape: probe_shape(SHAPE_SPHERE, m),
            point: float4(*p),
        })
        .collect();

    let out: Vec<Float4> = Gpu::new()?.run("probe_normal_at", &input)?;

    for (i, ((m, p), got)) in cases.iter().zip(&out).enumerate() {
        let mut sphere = Sphere::new();
        sphere.set_transform(*m);
        check_tuple(
            *got,
            sphere.normal_at_no_intersect(*p),
            &format!("case {i}"),
        )?;
    }
    Ok(())
}

/// Input to `probe_lighting`.
#[repr(C)]
#[derive(Clone, Copy)]
struct LightingProbe {
    material: GpuMaterial,
    light: GpuLight,
    point: Float4,
    eyev: Float4,
    normalv: Float4,
    intensity: f32,
}

#[test]
fn test_chap_6_gpu_lighting() -> Result<(), String> {
    use std::f64::consts::FRAC_1_SQRT_2 as H;
    let material = Material::new();
    let point = Tuple::point(0.0, 0.0, 0.0);
    let normalv = Tuple::vector(0.0, 0.0, -1.0);
    // (eye, light position, intensity): the book's chapter 6 cases, then
    // chapter 8's "lighting with the surface in shadow".
    let cases = [
        (
            Tuple::vector(0.0, 0.0, -1.0),
            Tuple::point(0.0, 0.0, -10.0),
            1.0,
        ),
        (
            Tuple::vector(0.0, H, -H),
            Tuple::point(0.0, 0.0, -10.0),
            1.0,
        ),
        (
            Tuple::vector(0.0, 0.0, -1.0),
            Tuple::point(0.0, 10.0, -10.0),
            1.0,
        ),
        (
            Tuple::vector(0.0, -H, -H),
            Tuple::point(0.0, 10.0, -10.0),
            1.0,
        ),
        (
            Tuple::vector(0.0, 0.0, -1.0),
            Tuple::point(0.0, 0.0, 10.0),
            1.0,
        ),
        (
            Tuple::vector(0.0, 0.0, -1.0),
            Tuple::point(0.0, 0.0, -10.0),
            0.0,
        ),
    ];
    let lights: Vec<Light> = cases
        .iter()
        .map(|(_, pos, _)| Light::point_light(*pos, Tuple::color(1.0, 1.0, 1.0)))
        .collect();
    let input: Vec<LightingProbe> = cases
        .iter()
        .zip(&lights)
        .map(|((eyev, _, intensity), light)| LightingProbe {
            material: gpu_material(&material).unwrap(),
            light: gpu_light(light).unwrap(),
            point: float4(point),
            eyev: float4(*eyev),
            normalv: float4(normalv),
            intensity: *intensity as f32,
        })
        .collect();

    let out: Vec<Float4> = Gpu::new()?.run("probe_lighting", &input)?;

    let sphere = Sphere::new(); // carries the default material
    for (i, (((eyev, _, intensity), light), got)) in cases.iter().zip(&lights).zip(&out).enumerate()
    {
        let want = light.lighting(&sphere, point, point, *eyev, normalv, *intensity);
        for (c, (g, w)) in got.iter().zip([want.x, want.y, want.z]).enumerate() {
            check_close(*g, w, &format!("case {i} channel {c}"))?;
        }
    }
    Ok(())
}

/// Input to `probe_ray_for_pixel`.
#[repr(C)]
#[derive(Clone, Copy)]
struct CameraProbe {
    camera: GpuCamera,
    pixel: [u32; 2],
}

#[test]
fn test_chap_7_gpu_ray_for_pixel() -> Result<(), String> {
    use std::f64::consts::PI;
    let plain = Camera::new(201, 101, PI / 2.0);
    let moved = Camera::new(201, 101, PI / 2.0)
        .with_transform(Matrix4::rotation_y(PI / 4.0) * Matrix4::translation(0.0, -2.0, 5.0));
    let cases = [(&plain, 100, 50), (&plain, 0, 0), (&moved, 100, 50)];
    let input: Vec<CameraProbe> = cases
        .iter()
        .map(|(camera, x, y)| CameraProbe {
            camera: gpu_camera(camera).unwrap(),
            pixel: [*x as u32, *y as u32],
        })
        .collect();

    let out: Vec<GpuRay> = Gpu::new()?.run("probe_ray_for_pixel", &input)?;

    for (i, ((camera, x, y), got)) in cases.iter().zip(&out).enumerate() {
        let want = camera.ray_for_pixel(*x, *y);
        check_tuple(got.origin, want.origin, &format!("case {i} origin"))?;
        check_tuple(
            got.direction,
            want.direction,
            &format!("case {i} direction"),
        )?;
    }
    Ok(())
}

#[test]
fn test_chap_5_gpu_intersect_sphere_far_origin() -> Result<(), String> {
    use std::f64::consts::{FRAC_1_SQRT_2, SQRT_2};
    // A flattened sphere (a book chapter 7 wall, scaled 0.01) puts the
    // object-space ray origin hundreds of units out. The textbook quadratic
    // then subtracts nearly equal huge numbers and loses most f32 digits.
    // The first case is the chapter 7 wall ray that showed shadow acne.
    let rays = [
        Ray::new(
            Tuple::point(FRAC_1_SQRT_2, -500.0 * SQRT_2, -0.15),
            Tuple::vector(-0.03, 83.6, 0.02),
        ),
        Ray::new(
            Tuple::point(0.3, 0.2, -1000.0),
            Tuple::vector(0.0, 0.0, 1.0),
        ),
        Ray::new(
            Tuple::point(500.0, -400.0, 300.0),
            Tuple::vector(-5.0, 4.0, -3.0),
        ),
    ];
    let cases: Vec<GpuRay> = rays.iter().map(GpuRay::from).collect();

    let out: Vec<Float4> = Gpu::new()?.run("probe_intersect_sphere", &cases)?;

    let sphere = Sphere::new();
    for (i, (ray, got)) in rays.iter().zip(&out).enumerate() {
        let xs = sphere.local_intersect(ray);
        if got[2] as usize != xs.count() || xs.is_empty() {
            return Err(format!(
                "ray {i}: gpu {} hits vs cpu {}",
                got[2],
                xs.count()
            ));
        }
        check_close(got[0], xs[0].t, &format!("ray {i} t0"))?;
        check_close(got[1], xs[1].t, &format!("ray {i} t1"))?;
    }
    Ok(())
}

#[test]
fn test_chap_9_gpu_plane_normal_at() -> Result<(), String> {
    // Constant in object space; the second case checks the world transform.
    let cases = [
        (Matrix4::identity(), Tuple::point(0.0, 0.0, 0.0)),
        (Matrix4::identity(), Tuple::point(10.0, 0.0, -10.0)),
        (Matrix4::identity(), Tuple::point(-5.0, 0.0, 150.0)),
        (
            Matrix4::translation(0.0, 0.0, 5.0) * Matrix4::rotation_x(1.0),
            Tuple::point(1.0, 2.0, 3.0),
        ),
    ];
    let input: Vec<NormalProbe> = cases
        .iter()
        .map(|(m, p)| NormalProbe {
            shape: probe_shape(SHAPE_PLANE, m),
            point: float4(*p),
        })
        .collect();

    let out: Vec<Float4> = Gpu::new()?.run("probe_normal_at", &input)?;

    for (i, ((m, p), got)) in cases.iter().zip(&out).enumerate() {
        let mut plane = Plane::new();
        plane.set_transform(*m);
        check_tuple(*got, plane.normal_at_no_intersect(*p), &format!("case {i}"))?;
    }
    Ok(())
}

/// Input to `probe_local_intersect`.
#[repr(C, align(16))]
#[derive(Clone, Copy)]
struct IntersectProbe {
    ray: GpuRay,
    kind: u32,
}

#[test]
fn test_chap_9_gpu_intersect_plane() -> Result<(), String> {
    // The book's cases: parallel, coplanar, from above, from below.
    let rays = [
        Ray::new(Tuple::point(0.0, 10.0, 0.0), Tuple::vector(0.0, 0.0, 1.0)),
        Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0)),
        Ray::new(Tuple::point(0.0, 1.0, 0.0), Tuple::vector(0.0, -1.0, 0.0)),
        Ray::new(Tuple::point(0.0, -1.0, 0.0), Tuple::vector(0.0, 1.0, 0.0)),
    ];
    let input: Vec<IntersectProbe> = rays
        .iter()
        .map(|r| IntersectProbe {
            ray: r.into(),
            kind: SHAPE_PLANE,
        })
        .collect();

    // Each result is (t0, t1, hit count, unused).
    let out: Vec<Float4> = Gpu::new()?.run("probe_local_intersect", &input)?;

    let plane = Plane::new();
    for (i, (ray, got)) in rays.iter().zip(&out).enumerate() {
        let xs = plane.local_intersect(ray);
        if got[2] as usize != xs.count() {
            return Err(format!(
                "ray {i}: gpu {} hits vs cpu {}",
                got[2],
                xs.count()
            ));
        }
        if !xs.is_empty() {
            check_close(got[0], xs[0].t, &format!("ray {i} t"))?;
        }
    }
    Ok(())
}

/// Input to `probe_refractive_indices`: the hit is root `hit_root` (0 or 1,
/// ascending) of shape `hit_index` along `ray`.
#[repr(C)]
#[derive(Clone, Copy)]
struct RefractiveProbe {
    ray: GpuRay,
    hit_index: u32,
    hit_root: u32,
    shape_count: u32,
}

#[test]
fn test_chap_11_gpu_refractive_indices() -> Result<(), String> {
    // The book's n1/n2 table: glass spheres A (scaled 2), B and C (offset
    // along z) overlapping, so the ray crosses six boundaries. The real
    // geometry puts them at the book's t = 2, 2.75, 3.25, 4.75, 5.25, 6.
    let mut world = World::new();
    for (transform, refractive_index) in [
        (Matrix4::scaling(2.0, 2.0, 2.0), 1.5),
        (Matrix4::translation(0.0, 0.0, -0.25), 2.0),
        (Matrix4::translation(0.0, 0.0, 0.25), 2.5),
    ] {
        let mut s = Sphere::glass();
        s.set_transform(transform);
        s.data.material.refractive_index = refractive_index;
        world.add_shape(Box::new(s));
    }
    let shapes: Vec<GpuShape> = world
        .shapes
        .iter()
        .map(|s| GpuShape {
            material: gpu_material(s.material()).unwrap(),
            ..probe_shape(SHAPE_SPHERE, s.transform())
        })
        .collect();
    let ray = Ray::new(Tuple::point(0.0, 0.0, -4.0), Tuple::vector(0.0, 0.0, 1.0));
    let xs = world.intersect(&ray);
    if xs.count() != 6 {
        return Err(format!("cpu: {} intersections, expected 6", xs.count()));
    }
    let index_of = |id| world.shapes.iter().position(|s| s.id() == id).unwrap();
    let input: Vec<RefractiveProbe> = (0..xs.count())
        .map(|k| RefractiveProbe {
            ray: (&ray).into(),
            hit_index: index_of(xs[k].object_id) as u32,
            // how many of this shape's hits come before this one
            hit_root: (0..k)
                .filter(|&j| xs[j].object_id == xs[k].object_id)
                .count() as u32,
            shape_count: shapes.len() as u32,
        })
        .collect();

    // Each result is (n1, n2).
    let out: Vec<[f32; 2]> =
        Gpu::new()?.run_with_shapes("probe_refractive_indices", &input, &shapes)?;

    for (k, got) in out.iter().enumerate() {
        let comps = prepare_computations(xs[k], &ray, &world.shapes, &xs);
        check_close(got[0], comps.n1, &format!("hit {k} n1"))?;
        check_close(got[1], comps.n2, &format!("hit {k} n2"))?;
    }
    Ok(())
}

#[test]
fn test_chap_11_gpu_schlick() -> Result<(), String> {
    use std::f64::consts::FRAC_1_SQRT_2 as H;
    // The book's cases, all on a glass sphere: total internal reflection,
    // a perpendicular view, and a small angle with n2 > n1. Each case is
    // (ray, which of the t-sorted intersections is the hit).
    let cases = [
        (
            Ray::new(Tuple::point(0.0, 0.0, H), Tuple::vector(0.0, 1.0, 0.0)),
            1,
        ),
        (
            Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 1.0, 0.0)),
            1,
        ),
        (
            Ray::new(Tuple::point(0.0, 0.99, -2.0), Tuple::vector(0.0, 0.0, 1.0)),
            0,
        ),
    ];
    let mut world = World::new();
    world.add_shape(Box::new(Sphere::glass()));
    let comps: Vec<_> = cases
        .iter()
        .map(|(ray, k)| {
            let xs = world.intersect(ray);
            prepare_computations(xs[*k], ray, &world.shapes, &xs)
        })
        .collect();
    // Each input is (cos_i, n1, n2, unused).
    let input: Vec<Float4> = comps
        .iter()
        .map(|c| [c.eyev.dot(c.normalv) as f32, c.n1 as f32, c.n2 as f32, 0.0])
        .collect();

    let out: Vec<f32> = Gpu::new()?.run("probe_schlick", &input)?;

    for (i, (c, got)) in comps.iter().zip(&out).enumerate() {
        check_close(*got, schlick(c), &format!("case {i}"))?;
    }
    Ok(())
}
