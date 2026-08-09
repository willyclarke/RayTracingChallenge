// Measures what divide() (BVH) adds on top of the flat AABB cull.
//
// Scene: one FLAT group holding a grid of spheres (the case BVH is built for —
// a single group with many children the walker would otherwise scan linearly).
// Renders twice: AABB-only (build_bounds) vs. AABB + divide(). Reports
// prim-tests/ray and node-visits for each, and checks the images match.
//
// Run: cargo run --release --features stats --example bvh_stress

use rtc_rust::camera::Camera;
use rtc_rust::light::Light;
use rtc_rust::material::Material;
use rtc_rust::matrix::Matrix4;
use rtc_rust::shape::Shape;
use rtc_rust::shapes::group::Group;
use rtc_rust::shapes::sphere::Sphere;
use rtc_rust::tuple::Tuple;
use rtc_rust::world::{World, read_stats, reset_stats, view_transform};

const GRID: i32 = 11; // GRID x GRID spheres in one flat group

fn build_scene() -> (World, usize) {
    let mut world = World::new();
    world.light = Some(Light::point_light(
        Tuple::point(-10.0, 10.0, -10.0),
        Tuple::color(1.0, 1.0, 1.0),
    ));

    let g_id = world.add_shape(Box::new(Group::new()));
    let mut m = Material::new();
    m.color = Tuple::color(0.3, 0.6, 0.9);

    for i in 0..GRID {
        for j in 0..GRID {
            let x = (i - GRID / 2) as f64 * 1.5;
            let y = (j - GRID / 2) as f64 * 1.5;
            let mut s = Sphere::new();
            s.set_transform(Matrix4::translation(x, y, 0.0) * Matrix4::scaling(0.5, 0.5, 0.5));
            s.set_material(m.clone());
            world.add_child(g_id, Box::new(s));
        }
    }
    (world, g_id)
}

fn render_and_count(label: &str, world: &World, camera: Camera) {
    reset_stats();
    let _img = world.render_single(camera);
    let (nodes, prims) = read_stats();
    let rays = (camera.hsize * camera.vsize) as f64;
    println!(
        "{label:14}  nodes/ray: {:8.2}   prim-tests/ray: {:6.2}",
        nodes as f64 / rays,
        prims as f64 / rays
    );
}

fn main() {
    let from = Tuple::point(0.0, 0.0, -20.0);
    let to = Tuple::point(0.0, 0.0, 0.0);
    let up = Tuple::vector(0.0, 1.0, 0.0);
    let camera = Camera::new(400, 400, std::f64::consts::PI / 3.0)
        .with_transform(view_transform(from, to, up));

    println!("scene: 1 flat group of {} spheres\n", GRID * GRID);

    // AABB only
    let (mut a, _) = build_scene();
    a.build_bounds();
    render_and_count("AABB only", &a, camera);
    let img_a = a.render_single(camera);

    // AABB + divide
    let (mut b, g_id) = build_scene();
    b.divide(g_id, 4);
    b.build_bounds();
    render_and_count("AABB + divide", &b, camera);
    let img_b = b.render_single(camera);

    // correctness: divide must not change the image beyond float rounding.
    // NOTE: exact equality is the wrong bar here — BVH inserts identity-transform
    // sub-groups, and normal_to_world normalizes once per ancestor, so the normal
    // (hence color) drifts by ~1e-14 per added level. Compare against epsilon.
    const EPS: f64 = 1e-9;
    let (mut n_diff, mut max_d) = (0usize, 0.0_f64);
    for y in 0..camera.vsize {
        for x in 0..camera.hsize {
            let p = img_a.pixel_at(x, y);
            let q = img_b.pixel_at(x, y);
            let d = (p.x - q.x)
                .abs()
                .max((p.y - q.y).abs())
                .max((p.z - q.z).abs());
            max_d = max_d.max(d);
            if d > EPS {
                n_diff += 1;
            }
        }
    }
    println!(
        "\npixels differing > {EPS:e}: {n_diff}   (max channel diff: {max_d:e})  {}",
        if n_diff == 0 {
            "PASS - divide preserved the image"
        } else {
            "FAIL"
        }
    );
    // let _ = img_b.write_ppm("checkdivide.ppm"); // optional; comment out to keep I/O off the graph
}
