//! Whole-image tests: render a scene on the GPU and on the CPU (rtc-rust's
//! `render_parallel`) and compare the pixels.

use crate::gpu::Gpu;
use crate::scene::GpuScene;
use rtc_rust::camera::Camera;
use rtc_rust::light::Light;
use rtc_rust::material::Material;
use rtc_rust::matrix::Matrix4;
use rtc_rust::shape::Shape;
use rtc_rust::shapes::plane::Plane;
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

/// The book's chapter 9 "putting it together" scene: chapter 7's spheres,
/// now on a real plane floor in front of a plane wall.
fn chapter_9_scene(width: usize, height: usize) -> (World, Camera) {
    let (chapter_7, camera) = chapter_7_scene(width, height);
    let mut world = World::new();
    world.light = chapter_7.light.clone();
    let mut wall = Material::new();
    wall.color = Tuple::color(1.0, 0.9, 0.9);
    wall.specular = 0.0;
    let mut floor = Plane::new();
    floor.set_material(wall.clone());
    world.add_shape(Box::new(floor));
    let mut back = Plane::new();
    back.set_transform(Matrix4::translation(0.0, 0.0, 5.0) * Matrix4::rotation_x(PI / 2.0));
    back.set_material(wall);
    world.add_shape(Box::new(back));
    // chapter 7's last three shapes are the balls
    for ball in &chapter_7.shapes[3..] {
        let mut s = Sphere::new();
        s.set_transform(*ball.transform());
        s.set_material(ball.material().clone());
        world.add_shape(Box::new(s));
    }
    (world, camera)
}

/// Render on both GPU and CPU and compare. A pixel is "off" when a channel
/// differs by more than `1e-3`; at most 0.05% of pixels may be off.
///
/// Why any at all: a pixel on a hard shadow edge is a yes/no decision on a
/// shadow ray grazing an occluder (discriminant ~ 0), and f32 rounds some of
/// those the other way than rtc-rust's f64. At 1600x800 the chapter 9 scene
/// has 11 such pixels, all on shadow edges. Real bugs are far over budget
/// (the chapter 7 wall acne was 14% of pixels).
fn check_matches_cpu(world: &World, camera: Camera) -> Result<(), String> {
    let scene = GpuScene::from_world(world, &camera)?;
    let gpu = Gpu::new()?.render(&scene)?;
    let cpu = world.render_parallel(camera);

    let (w, h) = (scene.width(), scene.height());
    let mut worst = (0.0f64, 0, 0);
    let mut off = 0;
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
            if diff > 1e-3 {
                off += 1;
            }
        }
    }
    let budget = w * h / 2000;
    eprintln!(
        "worst pixel diff {:.2e} at ({}, {}); {off} of {} pixels off by > 1e-3 (budget {budget})",
        worst.0,
        worst.1,
        worst.2,
        w * h
    );
    if off > budget {
        return Err(format!(
            "{off} pixels off by > 1e-3 (budget {budget}); worst {:.2e} at ({}, {})",
            worst.0, worst.1, worst.2
        ));
    }
    Ok(())
}

#[test]
fn test_chap_7_gpu_render_matches_cpu() -> Result<(), String> {
    let (world, camera) = chapter_7_scene(200, 100);
    check_matches_cpu(&world, camera)
}

#[test]
fn test_chap_9_gpu_render_matches_cpu() -> Result<(), String> {
    // Higher resolution than chapter 7: more shadow-edge pixels to exercise
    // the budget, and the plane runs out to the horizon.
    let (world, camera) = chapter_9_scene(800, 400);
    check_matches_cpu(&world, camera)
}
