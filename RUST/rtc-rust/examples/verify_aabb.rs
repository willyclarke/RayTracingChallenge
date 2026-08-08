// Correctness check for the AABB group cull.
//
// Renders the SAME scene twice:
//   * reference: without build_bounds() -> groups keep infinite boxes -> NO cull
//   * culled:    with build_bounds()    -> AABB cull active
// If the cull is correct these must be pixel-identical. Any nonzero diff means
// a bounding box is culling a ray that has a real hit.

use rtc_rust::camera::Camera;
use rtc_rust::light::Light;
use rtc_rust::material::Material;
use rtc_rust::matrix::Matrix4;
use rtc_rust::tuple::Tuple;
use rtc_rust::world::{World, hexagon, view_transform};

fn build_scene() -> World {
    let mut world = World::new();
    world.light = Some(Light::point_light(
        Tuple::point(-10.0, 10.0, -10.0),
        Tuple::color(1.0, 1.0, 1.0),
    ));
    let mut m = Material::new();
    m.color = Tuple::color(0.3, 0.6, 0.9);
    hexagon(
        &mut world,
        0,
        Matrix4::rotation_x(-std::f64::consts::PI / 6.0),
        m,
    );
    world
}

fn main() {
    let from = Tuple::point(0.0, 2.5, -5.0);
    let to = Tuple::point(0.0, 0.0, 0.0);
    let up = Tuple::vector(0.0, 1.0, 0.0);
    let camera =
        Camera::new(400, 400, std::f64::consts::PI / 3.0).with_transform(view_transform(from, to, up));

    // reference: no build_bounds -> infinite group boxes -> cull is a no-op
    let reference = build_scene().render_single(camera);

    // culled: AABB active
    let mut culled_world = build_scene();
    culled_world.build_bounds();
    let culled = culled_world.render_single(camera);

    let mut max_diff = 0.0_f64;
    let mut n_diff = 0usize;
    for y in 0..camera.vsize {
        for x in 0..camera.hsize {
            let a = reference.pixel_at(x, y);
            let b = culled.pixel_at(x, y);
            let d = (a.x - b.x).abs().max((a.y - b.y).abs()).max((a.z - b.z).abs());
            if d > 0.0 {
                n_diff += 1;
                max_diff = max_diff.max(d);
            }
        }
    }

    println!("pixels differing: {n_diff}   max channel diff: {max_diff:e}");
    if n_diff == 0 {
        println!("PASS: AABB cull is pixel-identical to the reference render.");
    } else {
        println!("FAIL: cull changed the image -> a box is culling real hits.");
        std::process::exit(1);
    }
}
