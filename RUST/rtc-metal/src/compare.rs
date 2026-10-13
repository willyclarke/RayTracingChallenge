//! Whole-image tests: render a scene on the GPU and on the CPU (rtc-rust's
//! `render_parallel`) and compare the pixels.

use crate::gpu::Gpu;
use crate::scene::GpuScene;
use rtc_rust::camera::Camera;
use rtc_rust::light::Light;
use rtc_rust::material::Material;
use rtc_rust::matrix::Matrix4;
use rtc_rust::scene::SceneDescription;
use rtc_rust::shape::Shape;
use rtc_rust::shapes::group::Group;
use rtc_rust::shapes::plane::Plane;
use rtc_rust::shapes::sphere::Sphere;
use rtc_rust::tuple::Tuple;
use rtc_rust::world::{World, hexagon, view_transform};
use std::f64::consts::PI;

/// The book's chapter 7 "putting it together" scene: three spheres on a
/// floor and two walls that are themselves flattened spheres.
pub(crate) fn chapter_7_scene(width: usize, height: usize) -> (World, Camera) {
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

/// Chapter 9's scene with chapter 11 materials: a dim mirror floor, a
/// mirror ball, and a hollow glass ball (an air bubble inside glass, so
/// refraction crosses nested media) in front of it.
fn chapter_11_scene(width: usize, height: usize) -> (World, Camera) {
    let (chapter_9, camera) = chapter_9_scene(width, height);
    let mut world = World::new();
    world.light = chapter_9.light.clone();
    let mut materials: Vec<Material> = chapter_9
        .shapes
        .iter()
        .map(|s| s.material().clone())
        .collect();
    materials[0].reflective = 0.3; // floor
    let mirror = &mut materials[2]; // middle ball
    mirror.color = Tuple::color(0.1, 0.1, 0.1);
    mirror.reflective = 0.9;
    let glass = &mut materials[3]; // right ball
    glass.color = Tuple::color(0.0, 0.0, 0.0);
    glass.diffuse = 0.1;
    glass.specular = 1.0;
    glass.shininess = 300.0;
    glass.reflective = 0.9;
    glass.transparency = 0.9;
    glass.refractive_index = 1.5;
    // chapter 9's first two shapes are the planes
    for (i, (shape, material)) in chapter_9.shapes.iter().zip(materials).enumerate() {
        let mut s: Box<dyn Shape> = if i < 2 {
            Box::new(Plane::new())
        } else {
            Box::new(Sphere::new())
        };
        s.set_transform(*shape.transform());
        s.set_material(material);
        world.add_shape(s);
    }
    let mut bubble = Sphere::new();
    bubble.set_transform(*chapter_9.shapes[3].transform() * Matrix4::scaling(0.5, 0.5, 0.5));
    let air = &mut bubble.data.material;
    air.color = Tuple::color(0.0, 0.0, 0.0);
    air.ambient = 0.0;
    air.diffuse = 0.0;
    air.specular = 1.0;
    air.shininess = 300.0;
    air.reflective = 0.9;
    air.transparency = 0.9;
    air.refractive_index = 1.0000034;
    world.add_shape(Box::new(bubble));
    (world, camera)
}

/// The book's chapter 14 hexagon (spheres for corners, open cylinders for
/// edges, three levels of groups) as rtc-rust's own chapter 14 scene builds
/// it, several times over inside a tilted, scaled top group, on a mirror
/// floor. `divide` then sorts the top group's children into BVH subgroups,
/// so the flattening also walks the identity-transform groups it inserts.
pub(crate) fn chapter_14_scene(width: usize, height: usize) -> (World, Camera) {
    let mut world = World::new();
    world.light = Some(Light::point_light(
        Tuple::point(-10.0, 10.0, -10.0),
        Tuple::color(1.0, 1.0, 1.0),
    ));
    let mut floor = Plane::new();
    let mut m = Material::new();
    m.color = Tuple::color(1.0, 0.9, 0.9);
    m.specular = 0.0;
    m.reflective = 0.2;
    floor.set_material(m);
    world.add_shape(Box::new(floor));

    let mut top = Group::new();
    top.set_transform(Matrix4::translation(0.0, 1.2, 0.0) * Matrix4::rotation_x(-PI / 6.0));
    let top = world.add_shape(Box::new(top));
    for ix in 0..4 {
        let mut m = Material::new();
        m.color = Tuple::color(ix as f64 / 4.0, 0.5, 0.7);
        m.specular = 0.6;
        for iz in 0..2 {
            let scale = 0.4 + 0.1 * iz as f64;
            hexagon(
                &mut world,
                top,
                Matrix4::translation(-2.4 + 1.6 * ix as f64, 0.0, 1.2 * iz as f64)
                    * Matrix4::rotation_y(ix as f64 * PI / 12.0)
                    * Matrix4::scaling(scale, scale, scale),
                m.clone(),
            );
        }
    }
    world.divide(top, 4);
    world.build_bounds();
    let camera = Camera::new(width, height, PI / 3.0).with_transform(view_transform(
        Tuple::point(0.0, 2.5, -6.0),
        Tuple::point(0.0, 0.8, 0.0),
        Tuple::vector(0.0, 1.0, 0.0),
    ));
    (world, camera)
}

/// A scene from `scenes/`, at `width` x `height` instead of its own size.
/// Model paths resolve relative to `scenes/`, as the CLI does for a scene file.
pub(crate) fn scene_file(json: &str, width: usize, height: usize) -> (World, Camera) {
    let mut desc = SceneDescription::from_json(json).unwrap();
    (desc.camera.width, desc.camera.height) = (width, height);
    desc.build(std::path::Path::new("scenes")).unwrap()
}

/// Render on both GPU and CPU and compare. A pixel is "off" when a channel
/// differs by more than `1e-3`; at most 0.05% of pixels may be off.
///
/// Why any at all: a pixel on a hard shadow edge is a yes/no decision on a
/// shadow ray grazing an occluder (discriminant ~ 0), and f32 rounds some of
/// those the other way than rtc-rust's f64. Refraction adds edges seen
/// through glass, which f32 shifts by a fraction of a pixel. At 1600x800
/// the chapter 11 scene has 290 such pixels, nearly all where the CPU image
/// itself changes sharply between neighbours. Real bugs are far over budget
/// (the chapter 7 wall acne was 14% of pixels).
fn check_matches_cpu(world: &World, camera: Camera) -> Result<(), String> {
    let scene = GpuScene::from_world(world, &camera)?;
    // Compiling the MSL library is setup, not render time.
    let gpu = Gpu::new()?;
    let start = std::time::Instant::now();
    let gpu = gpu.render(&scene)?;
    let gpu_secs = start.elapsed().as_secs_f64();
    let start = std::time::Instant::now();
    let cpu = world.render_parallel(camera);
    let cpu_secs = start.elapsed().as_secs_f64();

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
        "worst pixel diff {:.2e} at ({}, {}); {off} of {} pixels off by > 1e-3 (budget {budget}); gpu {gpu_secs:.2} s, cpu {cpu_secs:.2} s",
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

#[test]
fn test_chap_10_gpu_render_matches_cpu() -> Result<(), String> {
    // Chapter 9's layout with a pattern on every surface.
    let (world, camera) = scene_file(include_str!("../scenes/chapter10.json"), 400, 200);
    check_matches_cpu(&world, camera)
}

#[test]
fn test_chap_11_gpu_render_matches_cpu() -> Result<(), String> {
    let (world, camera) = chapter_11_scene(400, 200);
    check_matches_cpu(&world, camera)
}

#[test]
fn test_chap_12_gpu_render_matches_cpu() -> Result<(), String> {
    // Cubes: a room seen from inside, a flattened table top, glass and mirror.
    let (world, camera) = scene_file(include_str!("../scenes/chapter12.json"), 400, 200);
    check_matches_cpu(&world, camera)
}

#[test]
fn test_chap_13_gpu_render_matches_cpu() -> Result<(), String> {
    // Cylinders and cones, including a glass closed double cone (4 roots).
    let (world, camera) = scene_file(include_str!("../scenes/chapter13.json"), 400, 200);
    check_matches_cpu(&world, camera)
}

#[test]
fn test_chap_14_gpu_render_matches_cpu() -> Result<(), String> {
    // Hexagons three groups deep under a tilted top group plus BVH subgroups.
    let (world, camera) = chapter_14_scene(400, 200);
    check_matches_cpu(&world, camera)
}

#[test]
fn test_chap_15_gpu_render_matches_cpu() -> Result<(), String> {
    // Flat and smooth OBJ teapots, each a BVH of triangles. About 12 pixels
    // differ: the CPU (EPSILON 1e-9) and GPU (1e-5) over_point offsets decide
    // shadows differently on the smooth teapot's shadow terminator, where
    // interpolated normals put the light behind the flat facet and the CPU's
    // shadow ray re-hits the same triangle, plus silhouette edges. Investigated:
    // no triangle cracks.
    let (world, camera) = scene_file(include_str!("../scenes/chapter15.json"), 400, 200);
    check_matches_cpu(&world, camera)
}

#[test]
fn test_chap_16_gpu_render_matches_cpu() -> Result<(), String> {
    // Two CSG dice (23 leaves each) with perturbed marble bodies.
    let (world, camera) = scene_file(include_str!("../scenes/chapter16.json"), 400, 250);
    check_matches_cpu(&world, camera)
}

#[test]
fn test_bonus_soft_shadows_gpu_render_matches_cpu() -> Result<(), String> {
    // A 10x10 area light: penumbrae under both balls and in their reflections.
    let (world, camera) = scene_file(include_str!("../scenes/soft_shadows.json"), 200, 80);
    check_matches_cpu(&world, camera)
}

/// Mean absolute difference per channel (RGB) between two row-major images;
/// the worst channel.
fn mean_abs_error(a: &[[f32; 4]], b: impl Fn(usize) -> [f64; 3]) -> f64 {
    let mut sums = [0.0f64; 3];
    for (i, p) in a.iter().enumerate() {
        let q = b(i);
        for c in 0..3 {
            sums[c] += (p[c] as f64 - q[c]).abs();
        }
    }
    sums.iter()
        .fold(0.0, |m, s| f64::max(m, s / a.len() as f64))
}

#[test]
fn test_bonus_soft_shadows_gpu_sampled_render_converges() -> Result<(), String> {
    // Calibrated once: measured err_4 = 2.604e-2, err_256 = 4.586e-3;
    // THRESHOLD = 1.5x err_256, rounded up to 2 significant digits. The error
    // does not go to 0: the CPU reference uses cell centres and one ray per
    // pixel centre, the sampled frames integrate over the whole light and
    // pixel.
    const THRESHOLD: f64 = 6.9e-3;
    let (world, camera) = scene_file(include_str!("../scenes/soft_shadows.json"), 200, 80);
    let scene = GpuScene::from_world(&world, &camera)?;
    let cpu = world.render_parallel(camera);
    let w = scene.width();
    let cpu_at = |i: usize| {
        let c = cpu.pixel_at(i % w, i / w);
        [c.x, c.y, c.z]
    };
    let gpu = Gpu::new()?;
    let mut acc = gpu.accumulator(&scene)?;
    for _ in 0..4 {
        acc.add_frame(&gpu, &scene.info)?;
    }
    let err_4 = mean_abs_error(&acc.resolve(&gpu)?, cpu_at);
    for _ in 4..256 {
        acc.add_frame(&gpu, &scene.info)?;
    }
    let err_256 = mean_abs_error(&acc.resolve(&gpu)?, cpu_at);
    eprintln!("err_4 {err_4:.3e}, err_256 {err_256:.3e}");
    if err_256 > 0.5 * err_4 {
        return Err(format!("err_256 {err_256} not <= half of err_4 {err_4}"));
    }
    if err_256 >= THRESHOLD {
        return Err(format!("err_256 {err_256} >= {THRESHOLD}"));
    }
    Ok(())
}

#[test]
fn test_chap_7_gpu_sampled_render_stays_close() -> Result<(), String> {
    // A point light: sampling only jitters the ray within the pixel, so 64
    // frames are the deterministic image with anti-aliased edges.
    let (world, camera) = chapter_7_scene(160, 80);
    let scene = GpuScene::from_world(&world, &camera)?;
    let gpu = Gpu::new()?;
    let deterministic = gpu.render(&scene)?;
    let mut acc = gpu.accumulator(&scene)?;
    for _ in 0..64 {
        acc.add_frame(&gpu, &scene.info)?;
    }
    let err = mean_abs_error(&acc.resolve(&gpu)?, |i| {
        let [r, g, b, _] = deterministic[i];
        [r as f64, g as f64, b as f64]
    });
    eprintln!("chapter 7 sampled vs deterministic: {err:.3e}");
    if err >= 0.01 {
        return Err(format!("mean absolute error {err} >= 0.01"));
    }
    Ok(())
}
