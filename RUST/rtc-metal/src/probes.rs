//! Probe tests: run single `rt.metal` functions on the GPU, one thread per
//! case, and compare with the rtc-rust CPU version of the same function.

use crate::gpu::Gpu;
use crate::types::{Float4, GpuMat4, GpuRay};
use rtc_rust::matrix::Matrix4;
use rtc_rust::ray::Ray;
use rtc_rust::shape::Shape;
use rtc_rust::shapes::sphere::Sphere;
use rtc_rust::tuple::Tuple;

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
