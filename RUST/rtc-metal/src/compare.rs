//! Whole-image tests: render a scene on the GPU and on the CPU (rtc-rust's
//! `render_parallel`) and compare the pixels.

use crate::gpu::Gpu;
use crate::scene::GpuScene;
use rtc_rust::camera::Camera;
use rtc_rust::light::Light;
use rtc_rust::material::Material;
use rtc_rust::matrix::Matrix4;
use rtc_rust::shape::Shape;
use rtc_rust::shapes::sphere::Sphere;
use rtc_rust::tuple::Tuple;
use rtc_rust::world::{World, view_transform};
use std::f64::consts::PI;

/// The book's chapter 7 "putting it together" scene: three spheres on a
/// floor and two walls that are themselves flattened spheres.
fn chapter_7_scene(width: usize, height: usize) -> (World, Camera) {
    let mut world = World::new();
    world.light = Some(Light::point_light(
        Tuple::point(-10.0, 10.0, -10.0),
        Tuple::color(1.0, 1.0, 1.0),
    ));
    let mut wall = Material::new();
    wall.color = Tuple::color(1.0, 0.9, 0.9);
    wall.specular = 0.0;
    let flat = Matrix4::scaling(10.0, 0.01, 10.0);
    let spheres = [
        (flat, wall.clone()),
        (
            Matrix4::translation(0.0, 0.0, 5.0)
                * Matrix4::rotation_y(-PI / 4.0)
                * Matrix4::rotation_x(PI / 2.0)
                * flat,
            wall.clone(),
        ),
        (
            Matrix4::translation(0.0, 0.0, 5.0)
                * Matrix4::rotation_y(PI / 4.0)
                * Matrix4::rotation_x(PI / 2.0)
                * flat,
            wall,
        ),
        (
            Matrix4::translation(-0.5, 1.0, 0.5),
            ball(Tuple::color(0.1, 1.0, 0.5)),
        ),
        (
            Matrix4::translation(1.5, 0.5, -0.5) * Matrix4::scaling(0.5, 0.5, 0.5),
            ball(Tuple::color(0.5, 1.0, 0.1)),
        ),
        (
            Matrix4::translation(-1.5, 0.33, -0.75) * Matrix4::scaling(0.33, 0.33, 0.33),
            ball(Tuple::color(1.0, 0.8, 0.1)),
        ),
    ];
    for (transform, material) in spheres {
        let mut s = Sphere::new();
        s.set_transform(transform);
        s.set_material(material);
        world.add_shape(Box::new(s));
    }
    let camera = Camera::new(width, height, PI / 3.0).with_transform(view_transform(
        Tuple::point(0.0, 1.5, -5.0),
        Tuple::point(0.0, 1.0, 0.0),
        Tuple::vector(0.0, 1.0, 0.0),
    ));
    (world, camera)
}

fn ball(color: Tuple) -> Material {
    let mut m = Material::new();
    m.color = color;
    m.diffuse = 0.7;
    m.specular = 0.3;
    m
}

#[test]
fn test_chap_7_gpu_render_matches_cpu() -> Result<(), String> {
    let (world, camera) = chapter_7_scene(200, 100);
    let scene = GpuScene::from_world(&world, &camera)?;
    let gpu = Gpu::new()?.render(&scene)?;
    let cpu = world.render_parallel(camera);

    let (w, h) = (scene.width(), scene.height());
    let mut worst = (0.0f64, 0, 0);
    let mut off = 0; // pixels more than one 8-bit step apart
    for y in 0..h {
        for x in 0..w {
            let c = cpu.pixel_at(x, y);
            let g = gpu[y * w + x];
            let diff = [c.x, c.y, c.z]
                .iter()
                .zip(g)
                .map(|(c, g)| (c - g as f64).abs())
                .fold(0.0, f64::max);
            if diff > worst.0 {
                worst = (diff, x, y);
            }
            if diff > 1.0 / 255.0 {
                off += 1;
            }
        }
    }
    eprintln!(
        "worst pixel diff {:.2e} at ({}, {}); {off} of {} pixels off by more than 1/255",
        worst.0,
        worst.1,
        worst.2,
        w * h
    );
    if worst.0 > 1e-3 {
        return Err(format!(
            "worst diff {:.2e} at ({}, {})",
            worst.0, worst.1, worst.2
        ));
    }
    Ok(())
}
