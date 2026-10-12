//! Probe tests: run single `rt.metal` functions on the GPU, one thread per
//! case, and compare with the rtc-rust CPU version of the same function.

use crate::gpu::Gpu;
use crate::scene::{GpuScene, gpu_camera, gpu_light, gpu_material, gpu_shape};
use crate::types::{
    Float4, GpuCamera, GpuLight, GpuMat4, GpuMaterial, GpuRay, GpuShape, SHAPE_PLANE, SHAPE_SPHERE,
    float4,
};
use rtc_rust::camera::Camera;
use rtc_rust::light::Light;
use rtc_rust::material::Material;
use rtc_rust::matrix::Matrix4;
use rtc_rust::pattern::Pattern;
use rtc_rust::ray::Ray;
use rtc_rust::shape::Shape;
use rtc_rust::shapes::cone::Cone;
use rtc_rust::shapes::cube::Cube;
use rtc_rust::shapes::cylinder::Cylinder;
use rtc_rust::shapes::plane::Plane;
use rtc_rust::shapes::sphere::Sphere;
use rtc_rust::shapes::triangle::Triangle;
use rtc_rust::shapes::triangleuv::TriangleUV;
use rtc_rust::tuple::Tuple;
use rtc_rust::world::{World, prepare_computations, schlick, view_transform};
use std::f64::consts::PI;

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
        minimum: 0.0,
        maximum: 0.0,
        closed: 0,
        p1: [0.0; 4],
        e1: [0.0; 4],
        e2: [0.0; 4],
        n1: [0.0; 4],
        n2: [0.0; 4],
        n3: [0.0; 4],
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
    object_point: Float4,
    eyev: Float4,
    normalv: Float4,
    intensity: f32,
    /// The `Sampler`: `sampling` 0 is rtc-rust's deterministic grid.
    pixel: u32,
    frame: u32,
    sampling: u32,
}

const _: () = assert!(size_of::<LightingProbe>() == 304);

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
            object_point: float4(point),
            eyev: float4(*eyev),
            normalv: float4(normalv),
            intensity: *intensity as f32,
            pixel: 0,
            frame: 0,
            sampling: 0,
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

/// Input to `probe_local_intersect`: `ray` is in object space, so only
/// the shape's kind and extent matter.
#[repr(C)]
#[derive(Clone, Copy)]
struct IntersectProbe {
    shape: GpuShape,
    ray: GpuRay,
}

/// Output of `probe_local_intersect`.
#[repr(C, align(16))]
#[derive(Clone, Copy)]
struct IntersectResult {
    ts: Float4, // ascending; ts[count..] unspecified
    count: u32,
}

/// Intersect `rays` (object space) with `shape` on the GPU and compare every
/// t with rtc-rust's `local_intersect`.
fn check_local_intersect(shape: &dyn Shape, rays: &[Ray]) -> Result<(), String> {
    let gpu_shape = gpu_shape(shape)?;
    let input: Vec<IntersectProbe> = rays
        .iter()
        .map(|r| IntersectProbe {
            shape: gpu_shape,
            ray: r.into(),
        })
        .collect();

    let out: Vec<IntersectResult> = Gpu::new()?.run("probe_local_intersect", &input)?;

    for (i, (ray, got)) in rays.iter().zip(&out).enumerate() {
        let xs = shape.local_intersect(ray);
        if got.count as usize != xs.count() {
            return Err(format!(
                "ray {i}: gpu {} hits vs cpu {}",
                got.count,
                xs.count()
            ));
        }
        for k in 0..xs.count() {
            check_close(got.ts[k], xs[k].t, &format!("ray {i} t{k}"))?;
        }
    }
    Ok(())
}

/// Object-space normals of `shape` (identity transform) at `points`, GPU
/// against rtc-rust.
fn check_local_normals(shape: &dyn Shape, points: &[Tuple]) -> Result<(), String> {
    let gpu_shape = gpu_shape(shape)?;
    let input: Vec<NormalProbe> = points
        .iter()
        .map(|p| NormalProbe {
            shape: gpu_shape,
            point: float4(*p),
        })
        .collect();

    let out: Vec<Float4> = Gpu::new()?.run("probe_normal_at", &input)?;

    for (i, (p, got)) in points.iter().zip(&out).enumerate() {
        check_tuple(
            *got,
            shape.normal_at_no_intersect(*p),
            &format!("point {i}"),
        )?;
    }
    Ok(())
}

/// A 16x16 fan of unnormalized rays from z = -5 (plus `y_offset`) toward
/// (0.3, -0.2, 0): hits and misses, direction lengths != 1.
fn ray_fan(y_offset: f64) -> Vec<Ray> {
    (0..256)
        .map(|i| {
            let origin = Tuple::point(
                -3.0 + 6.0 * (i % 16) as f64 / 15.0,
                -3.0 + 6.0 * (i / 16) as f64 / 15.0 + y_offset,
                -5.0,
            );
            Ray::new(origin, Tuple::point(0.3, -0.2, 0.0) - origin)
        })
        .collect()
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
    check_local_intersect(&Plane::new(), &rays)
}

fn ray(origin: (f64, f64, f64), direction: (f64, f64, f64)) -> Ray {
    Ray::new(
        Tuple::point(origin.0, origin.1, origin.2),
        Tuple::vector(direction.0, direction.1, direction.2),
    )
}

#[test]
fn test_chap_12_gpu_intersect_cube() -> Result<(), String> {
    // The book's hits (each face, and from inside) and misses, then the fan.
    let mut rays = vec![
        ray((5.0, 0.5, 0.0), (-1.0, 0.0, 0.0)),
        ray((-5.0, 0.5, 0.0), (1.0, 0.0, 0.0)),
        ray((0.5, 5.0, 0.0), (0.0, -1.0, 0.0)),
        ray((0.5, -5.0, 0.0), (0.0, 1.0, 0.0)),
        ray((0.5, 0.0, 5.0), (0.0, 0.0, -1.0)),
        ray((0.5, 0.0, -5.0), (0.0, 0.0, 1.0)),
        ray((0.0, 0.5, 0.0), (0.0, 0.0, 1.0)),
        ray((-2.0, 0.0, 0.0), (0.2673, 0.5345, 0.8018)),
        ray((0.0, -2.0, 0.0), (0.8018, 0.2673, 0.5345)),
        ray((0.0, 0.0, -2.0), (0.5345, 0.8018, 0.2673)),
        ray((2.0, 0.0, 2.0), (0.0, 0.0, -1.0)),
        ray((0.0, 2.0, 2.0), (0.0, -1.0, 0.0)),
        ray((2.0, 2.0, 0.0), (-1.0, 0.0, 0.0)),
    ];
    rays.extend(ray_fan(0.0));
    check_local_intersect(&Cube::new(), &rays)
}

#[test]
fn test_chap_12_gpu_cube_normal_at() -> Result<(), String> {
    // The book's table: each face, and two corners.
    let points = [
        (1.0, 0.5, -0.8),
        (-1.0, -0.2, 0.9),
        (-0.4, 1.0, -0.1),
        (0.3, -1.0, -0.7),
        (-0.6, 0.3, 1.0),
        (0.4, 0.4, -1.0),
        (1.0, 1.0, 1.0),
        (-1.0, -1.0, -1.0),
    ]
    .map(|(x, y, z)| Tuple::point(x, y, z));
    check_local_normals(&Cube::new(), &points)
}

/// A cylinder or cone truncated to `minimum..maximum`.
fn truncated<S: Shape + Default>(set: impl Fn(&mut S)) -> S {
    let mut s = S::default();
    set(&mut s);
    s
}

#[test]
fn test_chap_13_gpu_intersect_cylinder() -> Result<(), String> {
    // The book's misses and hits on the infinite cylinder, then truncated
    // and capped cases, then the fan against all three.
    let book = [
        ray((1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
        ray((0.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
        ray((0.0, 0.0, -5.0), (1.0, 1.0, 1.0)),
        ray((1.0, 0.0, -5.0), (0.0, 0.0, 1.0)),
        ray((0.0, 0.0, -5.0), (0.0, 0.0, 1.0)),
        ray((0.5, 0.0, -5.0), (0.1, 1.0, 1.0)),
    ];
    let truncated_cases = [
        ray((0.0, 1.5, 0.0), (0.1, 1.0, 0.0)),
        ray((0.0, 3.0, -5.0), (0.0, 0.0, 1.0)),
        ray((0.0, 0.0, -5.0), (0.0, 0.0, 1.0)),
        ray((0.0, 2.0, -5.0), (0.0, 0.0, 1.0)),
        ray((0.0, 1.0, -5.0), (0.0, 0.0, 1.0)),
        ray((0.0, 1.5, -2.0), (0.0, 0.0, 1.0)),
    ];
    let capped_cases = [
        ray((0.0, 3.0, 0.0), (0.0, -1.0, 0.0)),
        ray((0.0, 3.0, -2.0), (0.0, -1.0, 2.0)),
        ray((0.0, 4.0, -2.0), (0.0, -1.0, 1.0)),
        ray((0.0, 0.0, -2.0), (0.0, 1.0, 2.0)),
        ray((0.0, -1.0, -2.0), (0.0, 1.0, 1.0)),
    ];
    let infinite = Cylinder::new();
    let open = truncated(|c: &mut Cylinder| (c.minimum, c.maximum) = (1.0, 2.0));
    let closed = Cylinder {
        closed: true,
        ..open.clone()
    };
    check_local_intersect(&infinite, &[&book[..], &ray_fan(0.0)].concat())?;
    check_local_intersect(&open, &[&truncated_cases[..], &ray_fan(1.5)].concat())?;
    check_local_intersect(&closed, &[&capped_cases[..], &ray_fan(1.5)].concat())
}

#[test]
fn test_chap_13_gpu_cylinder_normal_at() -> Result<(), String> {
    // The book's wall and cap tables.
    let wall = [
        (1.0, 0.0, 0.0),
        (0.0, 5.0, -1.0),
        (0.0, -2.0, 1.0),
        (-1.0, 1.0, 0.0),
    ];
    let caps = [
        (0.0, 1.0, 0.0),
        (0.5, 1.0, 0.0),
        (0.0, 1.0, 0.5),
        (0.0, 2.0, 0.0),
        (0.5, 2.0, 0.0),
        (0.0, 2.0, 0.5),
    ];
    let point = |&(x, y, z): &(f64, f64, f64)| Tuple::point(x, y, z);
    check_local_normals(
        &Cylinder::new(),
        &wall.iter().map(point).collect::<Vec<_>>(),
    )?;
    let closed = truncated(|c: &mut Cylinder| {
        (c.minimum, c.maximum, c.closed) = (1.0, 2.0, true);
    });
    check_local_normals(&closed, &caps.iter().map(point).collect::<Vec<_>>())
}

#[test]
fn test_chap_13_gpu_intersect_cone() -> Result<(), String> {
    // The book's wall hits, the ray parallel to one half, and the caps;
    // then the fan against a closed double cone, where a ray can cross
    // both halves and both caps: 4 roots.
    let mut book: Vec<Ray> = [
        ray((0.0, 0.0, -5.0), (0.0, 0.0, 1.0)),
        ray((1.0, 1.0, -5.0), (-0.5, -1.0, 1.0)),
        ray((0.0, 0.0, -1.0), (0.0, 1.0, 1.0)),
    ]
    .iter()
    .map(|r| Ray::new(r.origin, r.direction.normalize()))
    .collect();
    // The book's tangent ray (one double root), left unnormalized: the
    // normalized f32 direction makes the discriminant round to a tiny
    // negative and the GPU misses. A ray exactly grazing a cone is a yes/no
    // decision f32 cannot make reliably, as on shadow edges.
    book.push(ray((0.0, 0.0, -5.0), (1.0, 1.0, 1.0)));
    let caps = [
        ray((0.0, 0.0, -5.0), (0.0, 1.0, 0.0)),
        ray((0.0, 0.0, -0.25), (0.0, 1.0, 1.0)),
        ray((0.0, 0.0, -0.25), (0.0, 1.0, 0.0)),
    ];
    let caps: Vec<Ray> = caps
        .iter()
        .map(|r| Ray::new(r.origin, r.direction.normalize()))
        .collect();
    let half = truncated(|c: &mut Cone| {
        (c.minimum, c.maximum, c.closed) = (-0.5, 0.5, true);
    });
    let double = truncated(|c: &mut Cone| {
        (c.minimum, c.maximum, c.closed) = (-1.5, 1.0, true);
    });
    check_local_intersect(&Cone::new(), &[&book[..], &ray_fan(0.0)].concat())?;
    check_local_intersect(&half, &caps)?;
    // through the apex region: two wall roots, two cap roots
    let four = ray((-0.5, -2.0, -0.1), (0.3, 3.0, 0.05));
    if double.local_intersect(&four).count() != 4 {
        return Err("cpu: expected 4 roots on the double cone".into());
    }
    check_local_intersect(&double, &[&[four][..], &ray_fan(0.0)].concat())
}

#[test]
fn test_chap_13_gpu_cone_normal_at() -> Result<(), String> {
    // The book's wall table, then caps. Without the book's apex point: its
    // local normal is (0, 0, 0), which normal_at normalizes to NaN on both.
    let wall = [(1.0, 1.0, 1.0), (-1.0, -1.0, 0.0)];
    let caps = [(0.3, 1.0, 0.2), (-0.2, -1.5, 0.4), (0.5, 0.2, -0.1)];
    let point = |&(x, y, z): &(f64, f64, f64)| Tuple::point(x, y, z);
    check_local_normals(&Cone::new(), &wall.iter().map(point).collect::<Vec<_>>())?;
    let closed = truncated(|c: &mut Cone| {
        (c.minimum, c.maximum, c.closed) = (-1.5, 1.0, true);
    });
    check_local_normals(&closed, &caps.iter().map(point).collect::<Vec<_>>())
}

/// Input to `probe_refractive_indices`: the hit is root `hit_root` (0 or 1,
/// ascending) of `shapes[hit_index]` along `ray`.
#[repr(C)]
#[derive(Clone, Copy)]
struct RefractiveProbe {
    ray: GpuRay,
    hit_index: u32,
    hit_root: u32,
}

#[test]
fn test_chap_11_gpu_refractive_indices() -> Result<(), String> {
    // The book's n1/n2 table: glass spheres A (scaled 2), B and C (offset
    // along z) overlapping, so the ray crosses six boundaries. The real
    // geometry puts them at the book's t = 2, 2.75, 3.25, 4.75, 5.25, 6.
    let mut world = World::new();
    world.light = Some(Light::point_light(
        Tuple::point(-10.0, 10.0, -10.0),
        Tuple::color(1.0, 1.0, 1.0),
    ));
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
    let scene = GpuScene::from_world(&world, &Camera::default())?;
    let ray = Ray::new(Tuple::point(0.0, 0.0, -4.0), Tuple::vector(0.0, 0.0, 1.0));
    let xs = world.intersect(&ray);
    if xs.count() != 6 {
        return Err(format!("cpu: {} intersections, expected 6", xs.count()));
    }
    // The BVH reorders the shapes; scene.ids maps them back.
    let index_of = |id| scene.ids.iter().position(|&i| i == id).unwrap();
    let input: Vec<RefractiveProbe> = (0..xs.count())
        .map(|k| RefractiveProbe {
            ray: (&ray).into(),
            hit_index: index_of(xs[k].object_id) as u32,
            // how many of this shape's hits come before this one
            hit_root: (0..k)
                .filter(|&j| xs[j].object_id == xs[k].object_id)
                .count() as u32,
        })
        .collect();

    // Each result is (n1, n2).
    let out: Vec<[f32; 2]> =
        Gpu::new()?.run_with_scene("probe_refractive_indices", &input, &scene)?;

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

/// Input to `probe_trace`.
#[repr(C, align(16))]
#[derive(Clone, Copy)]
struct TraceProbe {
    ray: GpuRay,
    light_position: Float4,
}

/// Output of `probe_trace`.
#[repr(C, align(16))]
#[derive(Clone, Copy)]
struct TraceResult {
    normal: Float4,
    t: f32,
    index: u32,
    root: u32,
    shadowed: u32,
}

#[test]
fn test_chap_13_gpu_cone_far_origin_no_acne() -> Result<(), String> {
    // The chapter 13 scene's cone on its base is scaled 0.5 in x and z, so
    // camera rays reach it from ~17 units out in object space. There the
    // textbook discriminant b^2 - 4ac cancelled and put t off by ~3e-5,
    // more than EPSILON: over_point landed inside the cone, which then
    // shadowed itself (speckled acne). Every hit and shadow decision on a
    // patch of its wall must match the CPU.
    let mut desc =
        rtc_rust::scene::SceneDescription::from_json(include_str!("../scenes/chapter13.json"))
            .unwrap();
    (desc.camera.width, desc.camera.height) = (400, 200);
    let (world, camera) = desc.build(std::path::Path::new(".")).unwrap();
    let scene = GpuScene::from_world(&world, &camera)?;
    let light = world.light.as_ref().unwrap().position;
    let pixels: Vec<(usize, usize)> = (190..205)
        .flat_map(|x| (78..90).map(move |y| (x, y)))
        .collect();
    let rays: Vec<Ray> = pixels
        .iter()
        .map(|&(x, y)| camera.ray_for_pixel(x, y))
        .collect();
    let input: Vec<TraceProbe> = rays
        .iter()
        .map(|r| TraceProbe {
            ray: r.into(),
            light_position: float4(light),
        })
        .collect();

    let out: Vec<TraceResult> = Gpu::new()?.run_with_scene("probe_trace", &input, &scene)?;

    let cone = world.shapes[4].id(); // the cone on its base in chapter13.json
    let mut on_cone = 0;
    for (((x, y), ray), got) in pixels.iter().zip(&rays).zip(&out) {
        let xs = world.intersect(ray);
        let hit = xs.hit().ok_or(format!("({x}, {y}): cpu miss"))?;
        let comps = prepare_computations(hit, ray, &world.shapes, &xs);
        // got.index is into the BVH-ordered shapes; a miss is past the end.
        let got_id = scene.ids.get(got.index as usize);
        if got_id != Some(&hit.object_id) || hit.object_id != cone {
            continue; // only the cone's wall is under test
        }
        on_cone += 1;
        check_close(got.t, hit.t, &format!("({x}, {y}) t"))?;
        check_tuple(got.normal, comps.normalv, &format!("({x}, {y}) normal"))?;
        let shadowed = world.is_shadowed(comps.over_point, light);
        if (got.shadowed != 0) != shadowed {
            return Err(format!(
                "({x}, {y}): gpu shadowed {} vs cpu {shadowed}",
                got.shadowed != 0
            ));
        }
    }
    if on_cone < 50 {
        return Err(format!("only {on_cone} pixels on the cone; patch moved?"));
    }
    Ok(())
}

/// Input to `probe_pattern_at`.
#[repr(C)]
#[derive(Clone, Copy)]
struct PatternProbe {
    shape: GpuShape,
    point: Float4,
}

/// Pattern colors of `shape` at world `points`, GPU against rtc-rust's
/// `color_at_shape` (shape transform, then pattern transform).
fn check_pattern(shape: &dyn Shape, points: &[Tuple]) -> Result<(), String> {
    let gpu_shape = gpu_shape(shape)?;
    let input: Vec<PatternProbe> = points
        .iter()
        .map(|p| PatternProbe {
            shape: gpu_shape,
            point: float4(*p),
        })
        .collect();

    let out: Vec<Float4> = Gpu::new()?.run("probe_pattern_at", &input)?;

    let pattern = shape
        .material()
        .pattern
        .as_ref()
        .ok_or("shape has no pattern")?;
    for (i, (p, got)) in points.iter().zip(&out).enumerate() {
        let want = pattern.color_at_shape(shape, *p);
        for (c, (g, w)) in got.iter().zip([want.x, want.y, want.z]).enumerate() {
            check_close(*g, w, &format!("point {i} channel {c}"))?;
        }
    }
    Ok(())
}

/// A sphere with `pattern`, the book's object and pattern transforms
/// optional: (object transform, pattern transform).
fn patterned_sphere(mut pattern: Box<dyn Pattern>, transforms: (Matrix4, Matrix4)) -> Sphere {
    pattern.set_transform(transforms.1);
    let mut s = Sphere::new();
    s.set_transform(transforms.0);
    s.data.material.pattern = Some(pattern);
    s
}

/// A grid of points off the integer lattice, so no point sits on a pattern
/// boundary where f32 and f64 may floor differently, plus `extra`.
fn pattern_points(extra: &[(f64, f64, f64)]) -> Vec<Tuple> {
    let mut points: Vec<Tuple> = extra
        .iter()
        .map(|&(x, y, z)| Tuple::point(x, y, z))
        .collect();
    for i in 0..343 {
        let c = |k: usize| -2.95 + 0.83 * k as f64;
        points.push(Tuple::point(c(i % 7), c(i / 7 % 7), c(i / 49)));
    }
    points
}

#[test]
fn test_chap_10_gpu_pattern_at() -> Result<(), String> {
    use rtc_rust::patterns::checkerspattern::CheckersPattern;
    use rtc_rust::patterns::gradientpattern::GradientPattern;
    use rtc_rust::patterns::ringpattern::RingPattern;
    use rtc_rust::patterns::stripepattern::StripePattern;
    let (white, black) = (Tuple::color(1.0, 1.0, 1.0), Tuple::color(0.0, 0.0, 0.0));
    // The book's tables: stripes constant in y and z and alternating in x,
    // gradient, rings, checkers repeating in x, y and z.
    let book = pattern_points(&[
        (0.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 2.0, 0.0),
        (0.0, 0.0, 1.0),
        (0.0, 0.0, 2.0),
        (0.9, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (-0.1, 0.0, 0.0),
        (-1.0, 0.0, 0.0),
        (-1.1, 0.0, 0.0),
        (0.25, 0.0, 0.0),
        (0.5, 0.0, 0.0),
        (0.75, 0.0, 0.0),
        (0.708, 0.0, 0.708),
        (0.99, 0.0, 0.0),
        (1.01, 0.0, 0.0),
        (0.0, 0.99, 0.0),
        (0.0, 1.01, 0.0),
        (0.0, 0.0, 0.99),
        (0.0, 0.0, 1.01),
    ]);
    let identity = (Matrix4::identity(), Matrix4::identity());
    type NewPattern = fn(Tuple, Tuple) -> Box<dyn Pattern>;
    let patterns: [NewPattern; 4] = [
        |a, b| Box::new(StripePattern::new(a, b)),
        |a, b| Box::new(GradientPattern::new(a, b)),
        |a, b| Box::new(RingPattern::new(a, b)),
        |a, b| Box::new(CheckersPattern::new(a, b)),
    ];
    for (k, new) in patterns.iter().enumerate() {
        check_pattern(&patterned_sphere(new(white, black), identity), &book)
            .map_err(|e| format!("pattern {k}: {e}"))?;
    }
    // The book's object and pattern transformation cases, then both at once
    // (non-uniform and rotated, so mixing up the order cannot pass).
    let points = pattern_points(&[(1.5, 0.0, 0.0), (2.5, 0.0, 0.0), (3.5, 0.5, 1.0)]);
    let transforms = [
        (Matrix4::scaling(2.0, 2.0, 2.0), Matrix4::identity()),
        (Matrix4::identity(), Matrix4::scaling(2.0, 2.0, 2.0)),
        (
            Matrix4::scaling(2.0, 2.0, 2.0),
            Matrix4::translation(0.5, 0.0, 0.0),
        ),
        (
            Matrix4::rotation_y(0.4) * Matrix4::scaling(1.5, 0.7, 1.0),
            Matrix4::translation(0.3, -0.2, 0.1) * Matrix4::rotation_z(0.9),
        ),
    ];
    let (red, blue) = (Tuple::color(0.9, 0.1, 0.2), Tuple::color(0.1, 0.3, 0.8));
    for (k, new) in patterns.iter().enumerate() {
        for (j, t) in transforms.iter().enumerate() {
            check_pattern(&patterned_sphere(new(red, blue), *t), &points)
                .map_err(|e| format!("pattern {k}, transforms {j}: {e}"))?;
        }
    }
    Ok(())
}

#[test]
fn test_chap_10_gpu_lighting_with_pattern() -> Result<(), String> {
    use rtc_rust::patterns::stripepattern::StripePattern;
    // The book's case: only ambient, so the result is the stripe color.
    let mut material = Material::new();
    material.pattern = Some(Box::new(StripePattern::new(
        Tuple::color(1.0, 1.0, 1.0),
        Tuple::color(0.0, 0.0, 0.0),
    )));
    (material.ambient, material.diffuse, material.specular) = (1.0, 0.0, 0.0);
    let light = Light::point_light(Tuple::point(0.0, 0.0, -10.0), Tuple::color(1.0, 1.0, 1.0));
    let points = [Tuple::point(0.9, 0.0, 0.0), Tuple::point(1.1, 0.0, 0.0)];
    let (eyev, normalv) = (Tuple::vector(0.0, 0.0, -1.0), Tuple::vector(0.0, 0.0, -1.0));
    let input: Vec<LightingProbe> = points
        .iter()
        .map(|p| LightingProbe {
            material: gpu_material(&material).unwrap(),
            light: gpu_light(&light).unwrap(),
            point: float4(*p),
            object_point: float4(*p),
            eyev: float4(eyev),
            normalv: float4(normalv),
            intensity: 1.0,
            pixel: 0,
            frame: 0,
            sampling: 0,
        })
        .collect();

    let out: Vec<Float4> = Gpu::new()?.run("probe_lighting", &input)?;

    let mut sphere = Sphere::new();
    sphere.set_material(material);
    for (i, (p, got)) in points.iter().zip(&out).enumerate() {
        let want = light.lighting(&sphere, *p, *p, eyev, normalv, 1.0);
        for (c, (g, w)) in got.iter().zip([want.x, want.y, want.z]).enumerate() {
            check_close(*g, w, &format!("point {i} channel {c}"))?;
        }
    }
    Ok(())
}

#[test]
fn test_chap_15_gpu_intersect_triangle() -> Result<(), String> {
    // The book's rays at its triangle: parallel, three edge misses, a hit.
    let rays = [
        ray((0.0, -1.0, -2.0), (0.0, 1.0, 0.0)),
        ray((1.0, 1.0, -2.0), (0.0, 0.0, 1.0)),
        ray((-1.0, 1.0, -2.0), (0.0, 0.0, 1.0)),
        ray((0.0, -1.0, -2.0), (0.0, 0.0, 1.0)),
        ray((0.0, 0.5, -2.0), (0.0, 0.0, 1.0)),
    ];
    check_local_intersect(&Triangle::default(), &rays)?;

    // Collinear vertices: zero determinant, so no hit (Review Focus 1).
    let degenerate = Triangle::new(
        Tuple::point(0.0, 0.0, 0.0),
        Tuple::point(1.0, 0.0, 0.0),
        Tuple::point(2.0, 0.0, 0.0),
    );
    let rays = [
        ray((1.0, 0.0, -2.0), (0.0, 0.0, 1.0)),
        ray((1.0, 1.0, -2.0), (0.0, -0.5, 1.0)),
        ray((-1.0, 0.0, 0.0), (1.0, 0.0, 0.0)),
    ];
    check_local_intersect(&degenerate, &rays)?;
    let gpu = gpu_shape(&degenerate)?;
    let input: Vec<IntersectProbe> = rays
        .iter()
        .map(|r| IntersectProbe {
            shape: gpu,
            ray: r.into(),
        })
        .collect();
    let out: Vec<IntersectResult> = Gpu::new()?.run("probe_local_intersect", &input)?;
    for (i, (r, got)) in rays.iter().zip(&out).enumerate() {
        let cpu = degenerate.local_intersect(r).count();
        if got.count != 0 || cpu != 0 {
            return Err(format!("degenerate ray {i}: gpu {} cpu {cpu}", got.count));
        }
    }
    Ok(())
}

/// Input to `probe_triangle`: `ray` is in world space.
#[repr(C)]
#[derive(Clone, Copy)]
struct TriangleProbe {
    shape: GpuShape,
    ray: GpuRay,
}

/// Output of `probe_triangle`.
#[repr(C, align(16))]
#[derive(Clone, Copy)]
struct TriangleResult {
    normal: Float4, // flipped toward the eye
    uv: [f32; 2],
    t: f32,
    count: u32,
}

#[test]
fn test_chap_15_gpu_triangle_uv_and_normal() -> Result<(), String> {
    let (p1, p2, p3) = (
        Tuple::point(0.0, 1.0, 0.0),
        Tuple::point(-1.0, 0.0, 0.0),
        Tuple::point(1.0, 0.0, 0.0),
    );
    let (n1, n2, n3) = (
        Tuple::vector(0.0, 1.0, 0.0),
        Tuple::vector(-1.0, 0.0, 0.0),
        Tuple::vector(1.0, 0.0, 0.0),
    );

    // The book's smooth triangle: u = 0.45, v = 0.25.
    let smooth = TriangleUV::new(p1, p2, p3, n1, n2, n3);
    let book_ray = ray((-0.2, 0.3, -2.0), (0.0, 0.0, 1.0));
    let input = [TriangleProbe {
        shape: gpu_shape(&smooth)?,
        ray: (&book_ray).into(),
    }];
    let out: Vec<TriangleResult> = Gpu::new()?.run("probe_triangle", &input)?;
    let got = out[0];
    check_close(got.uv[0], 0.45, "book u")?;
    check_close(got.uv[1], 0.25, "book v")?;
    // The book's normal, up to the flip toward the eye the CPU also applies.
    let mut world = World::new();
    world.add_shape(Box::new(smooth));
    let xs = world.intersect(&book_ray);
    let hit = xs.hit().ok_or("cpu miss on the book ray")?;
    let comps = prepare_computations(hit, &book_ray, &world.shapes, &xs);
    let expected = Tuple::vector(-0.5547, 0.83205, 0.0);
    let flipped = Tuple::vector(0.5547, -0.83205, 0.0);
    let book = if (comps.normalv - expected).magnitude() < 1e-3 {
        expected
    } else {
        flipped
    };
    for (g, c) in got.normal.iter().zip([book.x, book.y, book.z, 0.0]) {
        if (*g as f64 - c).abs() > 1e-4 {
            return Err(format!("book normal {:?} vs {c}", got.normal));
        }
    }
    check_tuple(got.normal, comps.normalv, "book normal vs cpu")?;

    // A grid of rays at smooth and flat triangles under one transform.
    let m = Matrix4::translation(1.0, -2.0, 0.5)
        * Matrix4::rotation_y(0.7)
        * Matrix4::scaling(1.0, 2.0, 3.0);
    let mut smooth = TriangleUV::new(p1, p2, p3, n1, n2, n3);
    smooth.set_transform(m);
    let mut flat = Triangle::new(p1, p2, p3);
    flat.set_transform(m);
    for shape in [Box::new(smooth) as Box<dyn Shape>, Box::new(flat)] {
        let gpu = gpu_shape(shape.as_ref())?;
        let mut world = World::new();
        world.add_shape(shape);
        // Aim a fan from a fixed eye at points around the triangle's
        // centroid, so some rays hit and some miss.
        let eye = Tuple::point(0.5, -1.0, -6.0);
        let target = m * Tuple::point(0.0, 0.33, 0.0);
        let mut rays = Vec::new();
        for i in 0..7 {
            for j in 0..7 {
                let aim =
                    target + Tuple::vector((i as f64 - 3.0) * 0.3, (j as f64 - 3.0) * 0.5, 0.0);
                rays.push(Ray::new(eye, (aim - eye).normalize()));
            }
        }
        let input: Vec<TriangleProbe> = rays
            .iter()
            .map(|r| TriangleProbe {
                shape: gpu,
                ray: r.into(),
            })
            .collect();
        let out: Vec<TriangleResult> = Gpu::new()?.run("probe_triangle", &input)?;
        let mut hits = 0;
        for (i, (r, got)) in rays.iter().zip(&out).enumerate() {
            let xs = world.intersect(r);
            if got.count as usize != xs.count() {
                return Err(format!(
                    "ray {i}: gpu {} hits vs cpu {}",
                    got.count,
                    xs.count()
                ));
            }
            let Some(hit) = xs.hit() else { continue };
            hits += 1;
            check_close(got.t, hit.t, &format!("ray {i} t"))?;
            check_close(got.uv[0], hit.u, &format!("ray {i} u"))?;
            check_close(got.uv[1], hit.v, &format!("ray {i} v"))?;
            let comps = prepare_computations(hit, r, &world.shapes, &xs);
            check_tuple(got.normal, comps.normalv, &format!("ray {i} normal"))?;
        }
        if !(5..49).contains(&hits) {
            return Err(format!(
                "{hits} of 49 rays hit; grid should mix hits and misses"
            ));
        }
    }
    Ok(())
}

/// `probe_trace` inputs for `rays` toward `light`.
fn trace_probes(rays: &[Ray], light: Tuple) -> Vec<TraceProbe> {
    rays.iter()
        .map(|r| TraceProbe {
            ray: r.into(),
            light_position: float4(light),
        })
        .collect()
}

#[test]
fn test_chap_15_gpu_trace_through_bvh() -> Result<(), String> {
    // An OBJ mesh on a floor: every camera ray's hit or miss, shape, t,
    // normal and shadow must match the CPU, so the BVH never drops or
    // reorders a hit (book chapter 15).
    let mut world = World::new();
    let light = Tuple::point(-10.0, 10.0, -10.0);
    world.light = Some(Light::point_light(light, Tuple::color(1.0, 1.0, 1.0)));
    let mut floor = Plane::new();
    floor.set_transform(Matrix4::translation(0.0, -2.0, 0.0));
    let floor = world.add_shape(Box::new(floor));
    world
        .load_obj_file(
            "../rtc-rust/models/icosahedron.obj",
            Matrix4::scaling(1.5, 1.5, 1.5),
        )
        .map_err(|e| e.to_string())?;
    world.build_bounds();
    let camera = Camera::new(20, 20, PI / 3.0).with_transform(view_transform(
        Tuple::point(0.0, 1.5, -5.0),
        Tuple::point(0.0, 0.0, 0.0),
        Tuple::vector(0.0, 1.0, 0.0),
    ));
    let scene = GpuScene::from_world(&world, &camera)?;
    let rays: Vec<Ray> = (0..20)
        .flat_map(|y| (0..20).map(move |x| (x, y)))
        .map(|(x, y)| camera.ray_for_pixel(x, y))
        .collect();

    let out: Vec<TraceResult> =
        Gpu::new()?.run_with_scene("probe_trace", &trace_probes(&rays, light), &scene)?;

    let mut on_mesh = 0;
    for (i, (ray, got)) in rays.iter().zip(&out).enumerate() {
        let xs = world.intersect(ray);
        let got_id = scene.ids.get(got.index as usize);
        let Some(hit) = xs.hit() else {
            if got_id.is_some() {
                return Err(format!("ray {i}: gpu hit {got_id:?}, cpu miss"));
            }
            continue;
        };
        if got_id != Some(&hit.object_id) {
            return Err(format!(
                "ray {i}: gpu shape {got_id:?} vs cpu {}",
                hit.object_id
            ));
        }
        if hit.object_id != floor {
            on_mesh += 1;
        }
        let comps = prepare_computations(hit, ray, &world.shapes, &xs);
        check_close(got.t, hit.t, &format!("ray {i} t"))?;
        check_tuple(got.normal, comps.normalv, &format!("ray {i} normal"))?;
        let shadowed = world.is_shadowed(comps.over_point, light);
        if (got.shadowed != 0) != shadowed {
            return Err(format!(
                "ray {i}: gpu shadowed {} vs cpu {shadowed}",
                got.shadowed != 0
            ));
        }
    }
    // The unit icosahedron scaled 1.5 fills ~15% of the 60 degree view
    // from ~5.2 units: 56 of the 400 rays.
    if on_mesh < 50 {
        return Err(format!("only {on_mesh} rays hit the icosahedron"));
    }
    Ok(())
}

#[test]
fn test_chap_15_gpu_bvh_shared_edge_has_no_crack() -> Result<(), String> {
    // A quad split along its diagonal into two triangles, in a real tree
    // (12 spheres off to the side). Rays straight down the shared edge must
    // hit one of the two: the edge is on both triangles' u + v = 1 border,
    // and the end rays may start on or near the leaf box's faces. The
    // quad's leaf likely holds a sphere too, so its box is larger; on-face
    // cases are covered by test_chap_15_gpu_intersect_box.
    let mut world = World::new();
    world.light = Some(Light::point_light(
        Tuple::point(-10.0, 10.0, -10.0),
        Tuple::color(1.0, 1.0, 1.0),
    ));
    let a = world.add_shape(Box::new(TriangleUV::flat(
        Tuple::point(-1.0, 0.0, -1.0),
        Tuple::point(1.0, 0.0, -1.0),
        Tuple::point(1.0, 0.0, 1.0),
    )));
    let b = world.add_shape(Box::new(TriangleUV::flat(
        Tuple::point(-1.0, 0.0, -1.0),
        Tuple::point(1.0, 0.0, 1.0),
        Tuple::point(-1.0, 0.0, 1.0),
    )));
    for k in 0..12 {
        let mut s = Sphere::new();
        s.set_transform(
            Matrix4::translation(3.0 + (k % 4) as f64, 0.5, 3.0 + (k / 4) as f64)
                * Matrix4::scaling(0.3, 0.3, 0.3),
        );
        world.add_shape(Box::new(s));
    }
    let camera = Camera::new(1, 1, PI / 3.0);
    let scene = GpuScene::from_world(&world, &camera)?;
    if scene.nodes.len() < 3 {
        return Err(format!("only {} BVH nodes; no tree", scene.nodes.len()));
    }
    let rays: Vec<Ray> = (0..50)
        .map(|k| {
            let s = -1.0 + 2.0 * k as f64 / 49.0;
            ray((s, 1.0, s), (0.0, -1.0, 0.0))
        })
        .collect();
    let light = world.light.as_ref().unwrap().position;

    let out: Vec<TraceResult> =
        Gpu::new()?.run_with_scene("probe_trace", &trace_probes(&rays, light), &scene)?;

    for (i, got) in out.iter().enumerate() {
        let got_id = scene.ids.get(got.index as usize);
        if got_id != Some(&a) && got_id != Some(&b) {
            return Err(format!("ray {i}: missed the quad (gpu shape {got_id:?})"));
        }
        check_close(got.t, 1.0, &format!("ray {i} t"))?;
    }
    Ok(())
}

#[test]
fn test_chap_15_gpu_bvh_matches_linear() -> Result<(), String> {
    // GPU against GPU: the BVH traversal must find exactly the hit, t and
    // shadow that testing every shape linearly finds.
    let (world, camera) = crate::compare::chapter_14_scene(40, 25);
    let scene = GpuScene::from_world(&world, &camera)?;
    let light = world.light.as_ref().unwrap().position;
    let rays: Vec<Ray> = (0..25)
        .flat_map(|y| (0..40).map(move |x| (x, y)))
        .map(|(x, y)| camera.ray_for_pixel(x, y))
        .collect();
    let input = trace_probes(&rays, light);
    let gpu = Gpu::new()?;

    let bvh: Vec<TraceResult> = gpu.run_with_scene("probe_trace", &input, &scene)?;
    let linear: Vec<TraceResult> = gpu.run_with_scene("probe_trace_linear", &input, &scene)?;

    for (i, (b, l)) in bvh.iter().zip(&linear).enumerate() {
        if (b.index, b.t, b.shadowed) != (l.index, l.t, l.shadowed) {
            return Err(format!(
                "ray {i}: bvh (index {}, t {}, shadowed {}) vs linear ({}, {}, {})",
                b.index, b.t, b.shadowed, l.index, l.t, l.shadowed
            ));
        }
    }
    Ok(())
}

/// Input to `probe_intersect_box`.
#[repr(C, align(16))]
#[derive(Clone, Copy)]
struct BoxProbe {
    min: Float4,
    max: Float4,
    ray: GpuRay,
    tmin: f32,
    tmax: f32,
}

#[test]
fn test_chap_15_gpu_intersect_box() -> Result<(), String> {
    // The BVH's slab test on the box [-1, 1]^3. Axis-parallel rays have a
    // zero direction component, so 1 / direction is infinite; with the
    // origin on a face, the slab products would be 0 * inf = NaN.
    let cases = [
        // (origin, direction, tmin, tmax, crosses)
        ((0.0, 5.0, 0.0), (0.0, -1.0, 0.0), 0.0, f32::INFINITY, true),
        ((0.5, 5.0, -0.5), (0.0, -1.0, 0.0), 0.0, f32::INFINITY, true),
        ((1.0, 5.0, 0.0), (0.0, -1.0, 0.0), 0.0, f32::INFINITY, true), // on x = 1
        (
            (-1.0, 5.0, -1.0),
            (0.0, -1.0, 0.0),
            0.0,
            f32::INFINITY,
            true,
        ), // on an edge
        (
            (-1.0, 5.0, 1.0),
            (0.0, 1.0, 0.0),
            f32::NEG_INFINITY,
            f32::INFINITY,
            true,
        ),
        ((1.5, 5.0, 0.0), (0.0, -1.0, 0.0), 0.0, f32::INFINITY, false), // outside x
        (
            (0.0, 5.0, -2.0),
            (0.0, -1.0, 0.0),
            0.0,
            f32::INFINITY,
            false,
        ), // outside z
        ((0.0, 5.0, 0.0), (0.0, 1.0, 0.0), 0.0, f32::INFINITY, false),  // behind
        (
            (0.0, 5.0, 0.0),
            (0.0, 1.0, 0.0),
            f32::NEG_INFINITY,
            f32::INFINITY,
            true,
        ),
        ((0.0, 5.0, 0.0), (0.0, -1.0, 0.0), 0.0, 3.0, false), // ends before
        ((0.0, 5.0, 0.0), (0.0, -1.0, 0.0), 0.0, 4.0, true),  // ends on y = 1
        ((-5.0, 0.0, 0.0), (1.0, 1.0, 0.0), 0.0, f32::INFINITY, false), // passes by
        ((-5.0, 0.0, 0.0), (1.0, 0.2, 0.1), 0.0, f32::INFINITY, true),
    ];
    let input: Vec<BoxProbe> = cases
        .iter()
        .map(|&(o, d, tmin, tmax, _)| BoxProbe {
            min: [-1.0, -1.0, -1.0, 0.0],
            max: [1.0, 1.0, 1.0, 0.0],
            ray: (&ray(o, d)).into(),
            tmin,
            tmax,
        })
        .collect();

    let out: Vec<u32> = Gpu::new()?.run("probe_intersect_box", &input)?;

    for (i, (case, got)) in cases.iter().zip(&out).enumerate() {
        if (*got != 0) != case.4 {
            return Err(format!("case {i} {case:?}: gpu {}", *got != 0));
        }
    }
    Ok(())
}

/// Input to `probe_random`.
#[repr(C)]
#[derive(Clone, Copy)]
struct RandomProbe {
    pixel: u32,
    frame: u32,
    dim: u32,
}

fn random_on_gpu(probes: &[RandomProbe]) -> Result<Vec<f32>, String> {
    Gpu::new()?.run("probe_random", probes)
}

#[test]
fn test_bonus_soft_shadows_gpu_random_in_unit_interval() -> Result<(), String> {
    let input: Vec<RandomProbe> = (0..10_000)
        .map(|pixel| RandomProbe {
            pixel,
            frame: 7,
            dim: 3,
        })
        .collect();
    let out = random_on_gpu(&input)?;
    if let Some(bad) = out.iter().find(|v| !(0.0..1.0).contains(*v)) {
        return Err(format!("{bad} is outside [0, 1)"));
    }
    let mean = out.iter().sum::<f32>() / out.len() as f32;
    if (mean - 0.5).abs() > 0.01 {
        return Err(format!("mean {mean} is not 0.5 +- 0.01"));
    }
    Ok(())
}

#[test]
fn test_bonus_soft_shadows_gpu_random_decorrelated() -> Result<(), String> {
    let probe = |pixel, frame, dim| RandomProbe { pixel, frame, dim };
    let base: Vec<RandomProbe> = (0..1000).map(|i| probe(i, 0, 0)).collect();
    let neighbours = [
        (
            "pixel",
            (0..1000).map(|i| probe(i + 1, 0, 0)).collect::<Vec<_>>(),
        ),
        ("frame", (0..1000).map(|i| probe(i, 1, 0)).collect()),
        ("dim", (0..1000).map(|i| probe(i, 0, 1)).collect()),
    ];
    let base_out = random_on_gpu(&base)?;
    for (axis, input) in &neighbours {
        let out = random_on_gpu(input)?;
        let differing = base_out.iter().zip(&out).filter(|(a, b)| a != b).count();
        if differing < 990 {
            return Err(format!("{axis}: only {differing}/1000 values differ"));
        }
    }
    Ok(())
}

#[test]
fn test_bonus_soft_shadows_gpu_random_deterministic() -> Result<(), String> {
    let input: Vec<RandomProbe> = (0..100)
        .map(|i| RandomProbe {
            pixel: i,
            frame: i * 3,
            dim: i % 8,
        })
        .collect();
    if random_on_gpu(&input)? != random_on_gpu(&input)? {
        return Err("two runs differ".to_string());
    }
    Ok(())
}

/// Input to `probe_point_on_light`.
#[repr(C)]
#[derive(Clone, Copy)]
struct PointOnLightProbe {
    light: GpuLight,
    u: u32,
    v: u32,
    ju: f32,
    jv: f32,
}

const _: () = assert!(size_of::<PointOnLightProbe>() == 96);

#[test]
fn test_bonus_soft_shadows_gpu_point_on_light() -> Result<(), String> {
    use rtc_rust::light::Sequence;
    // rtc-rust's test_bonus_soft_shadows_4 (cell centres) and _7 (jitter
    // 0.3, 0.7) light and cells.
    let mut light = Light::area_light(
        Tuple::point(0.0, 0.0, 0.0),
        Tuple::vector(2.0, 0.0, 0.0),
        4,
        Tuple::vector(0.0, 0.0, 1.0),
        2,
        Tuple::color(1.0, 1.0, 1.0),
    );
    let cells = [(0, 0), (1, 0), (0, 1), (2, 0), (3, 1)];
    for (ju, jv) in [(0.5, 0.5), (0.3, 0.7)] {
        light.jitter_by = Sequence::new(vec![ju, jv]);
        let input: Vec<PointOnLightProbe> = cells
            .iter()
            .map(|&(u, v)| PointOnLightProbe {
                light: gpu_light(&light).unwrap(),
                u,
                v,
                ju: ju as f32,
                jv: jv as f32,
            })
            .collect();
        let out: Vec<Float4> = Gpu::new()?.run("probe_point_on_light", &input)?;
        for (&(u, v), got) in cells.iter().zip(&out) {
            let want = light.point_on_light(u as usize, v as usize);
            check_tuple(*got, want, &format!("jitter ({ju}, {jv}) cell ({u}, {v})"))?;
        }
    }
    Ok(())
}

/// Input to `probe_intensity_at`.
#[repr(C)]
#[derive(Clone, Copy)]
struct IntensityProbe {
    light: GpuLight,
    point: Float4,
    normalv: Float4,
    pixel: u32,
    frame: u32,
    sampling: u32,
}

const _: () = assert!(size_of::<IntensityProbe>() == 128);

#[test]
fn test_bonus_soft_shadows_gpu_intensity_at() -> Result<(), String> {
    let world = World::default_world();
    let scene = GpuScene::from_world(&world, &Camera::default())?;
    let white = Tuple::color(1.0, 1.0, 1.0);
    let lights = [
        // rtc-rust's test_bonus_soft_shadows_5: 2x2, with the corner rule.
        Light::area_light(
            Tuple::point(-0.5, -0.5, -5.0),
            Tuple::vector(1.0, 0.0, 0.0),
            2,
            Tuple::vector(0.0, 1.0, 0.0),
            2,
            white,
        ),
        // 4x1: every cell is tested, no corner rule.
        Light::area_light(
            Tuple::point(-1.5, -0.5, -5.0),
            Tuple::vector(3.0, 0.0, 0.0),
            4,
            Tuple::vector(0.0, 1.0, 0.0),
            1,
            white,
        ),
    ];
    let points = [
        Tuple::point(0.0, 0.0, 2.0),
        Tuple::point(1.0, -1.0, 2.0),
        Tuple::point(1.5, 0.0, 2.0),
        Tuple::point(1.25, 1.25, 3.0),
        Tuple::point(0.0, 0.0, -2.0),
    ];
    for (k, light) in lights.iter().enumerate() {
        let input: Vec<IntensityProbe> = points
            .iter()
            .map(|p| IntensityProbe {
                light: gpu_light(light).unwrap(),
                point: float4(*p),
                // Facing the light, as rtc-rust's intensity_at has no normal.
                normalv: [0.0, 0.0, -1.0, 0.0],
                pixel: 0,
                frame: 0,
                sampling: 0,
            })
            .collect();
        let out: Vec<f32> = Gpu::new()?.run_with_scene("probe_intensity_at", &input, &scene)?;
        for (p, got) in points.iter().zip(&out) {
            let want = light.intensity_at(*p, &world);
            check_close(*got, want, &format!("light {k} point {p}"))?;
        }
    }
    Ok(())
}

#[test]
fn test_bonus_soft_shadows_gpu_lighting() -> Result<(), String> {
    use std::f64::consts::FRAC_1_SQRT_2 as H;
    // rtc-rust's test_bonus_soft_shadows_9, then the same with the default
    // (specular) material so the specular average is covered too.
    let light = Light::area_light(
        Tuple::point(-0.5, -0.5, -5.0),
        Tuple::vector(1.0, 0.0, 0.0),
        2,
        Tuple::vector(0.0, 1.0, 0.0),
        2,
        Tuple::color(1.0, 1.0, 1.0),
    );
    let mut matte = Material::new();
    (matte.ambient, matte.diffuse, matte.specular) = (0.1, 0.9, 0.0);
    let eye = Tuple::point(0.0, 0.0, -5.0);
    let points = [Tuple::point(0.0, 0.0, -1.0), Tuple::point(0.0, H, -H)];
    for (k, material) in [matte, Material::new()].into_iter().enumerate() {
        let input: Vec<LightingProbe> = points
            .iter()
            .map(|p| LightingProbe {
                material: gpu_material(&material).unwrap(),
                light: gpu_light(&light).unwrap(),
                point: float4(*p),
                object_point: float4(*p),
                eyev: float4((eye - *p).normalize()),
                normalv: float4(Tuple::vector(p.x, p.y, p.z)),
                intensity: 1.0,
                pixel: 0,
                frame: 0,
                sampling: 0,
            })
            .collect();
        let out: Vec<Float4> = Gpu::new()?.run("probe_lighting", &input)?;
        let mut sphere = Sphere::new();
        sphere.set_material(material);
        for (p, got) in points.iter().zip(&out) {
            let (eyev, normalv) = ((eye - *p).normalize(), Tuple::vector(p.x, p.y, p.z));
            let want = light.lighting(&sphere, *p, *p, eyev, normalv, 1.0);
            for (c, (g, w)) in got.iter().zip([want.x, want.y, want.z]).enumerate() {
                check_close(*g, w, &format!("material {k} point {p} channel {c}"))?;
            }
        }
    }
    Ok(())
}

#[test]
fn test_bonus_soft_shadows_gpu_intensity_at_terminator() -> Result<(), String> {
    // A light 1e-3 rad below the tangent plane at the front of the default
    // world's unit sphere. rtc-rust's shadow ray from its over_point, 1e-9
    // off the surface, re-enters the sphere: shadowed. A shadow ray from the
    // GPU's over_point, 1e-5 off, would miss the sphere below
    // sqrt(2e-5) = 4.5e-3 rad; the GPU must still call it shadowed.
    let world = World::default_world();
    let scene = GpuScene::from_world(&world, &Camera::default())?;
    let (point, normalv) = (Tuple::point(0.0, 0.0, -1.0), Tuple::vector(0.0, 0.0, -1.0));
    let light = Light::point_light(Tuple::point(5.0, 0.0, -0.995), Tuple::color(1.0, 1.0, 1.0));
    let input = [IntensityProbe {
        light: gpu_light(&light).unwrap(),
        point: float4(point + normalv * EPSILON),
        normalv: float4(normalv),
        pixel: 0,
        frame: 0,
        sampling: 0,
    }];
    let out: Vec<f32> = Gpu::new()?.run_with_scene("probe_intensity_at", &input, &scene)?;
    let cpu_over_point = point + normalv * rtc_rust::math::EPSILON;
    check_close(
        out[0],
        light.intensity_at(cpu_over_point, &world),
        "terminator",
    )
}

/// Input to `probe_shade_point`: `intensity_at` then `lighting`, one Sampler.
#[repr(C)]
#[derive(Clone, Copy)]
struct ShadeProbe {
    material: GpuMaterial,
    light: GpuLight,
    point: Float4,
    eyev: Float4,
    normalv: Float4,
    pixel: u32,
    frame: u32,
    sampling: u32,
}

const _: () = assert!(size_of::<ShadeProbe>() == 288);

#[test]
fn test_bonus_soft_shadows_gpu_sampled_shading_converges() -> Result<(), String> {
    // rtc-rust's test_bonus_soft_shadows_5 light; (1, -1, 2) is in the
    // default world's penumbra. Sampling picks a uniform point on the whole
    // light, so the mean of many frames approaches the integral over the
    // light, not the 2x2 cell-centre estimate (intensity 0.25 there, while
    // about 49% of the light's area is visible). The reference is the
    // deterministic shading of the same light cut into 64x64 cells.
    let world = World::default_world();
    let scene = GpuScene::from_world(&world, &Camera::default())?;
    let light = |steps| {
        Light::area_light(
            Tuple::point(-0.5, -0.5, -5.0),
            Tuple::vector(1.0, 0.0, 0.0),
            steps,
            Tuple::vector(0.0, 1.0, 0.0),
            steps,
            Tuple::color(1.0, 1.0, 1.0),
        )
    };
    let probe = |light: &Light, frame, sampling| ShadeProbe {
        material: gpu_material(&Material::new()).unwrap(),
        light: gpu_light(light).unwrap(),
        point: [1.0, -1.0, 2.0, 1.0],
        eyev: [0.0, 0.0, -1.0, 0.0],
        normalv: [0.0, 0.0, -1.0, 0.0],
        pixel: 0,
        frame,
        sampling,
    };
    let gpu = Gpu::new()?;
    let want: Vec<Float4> =
        gpu.run_with_scene("probe_shade_point", &[probe(&light(64), 0, 0)], &scene)?;
    let coarse = light(2);
    let input: Vec<ShadeProbe> = (0..4096).map(|f| probe(&coarse, f, 1)).collect();
    let out: Vec<Float4> = gpu.run_with_scene("probe_shade_point", &input, &scene)?;
    for c in 0..3 {
        let mean = out.iter().map(|o| o[c] as f64).sum::<f64>() / out.len() as f64;
        let want = want[0][c] as f64;
        eprintln!("channel {c}: sampled mean {mean:.4}, 64x64 deterministic {want:.4}");
        if (mean - want).abs() > 0.01 {
            return Err(format!(
                "channel {c}: sampled mean {mean} vs deterministic {want}"
            ));
        }
    }
    Ok(())
}
