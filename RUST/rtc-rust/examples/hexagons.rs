use rtc_rust::camera::Camera;
use rtc_rust::light::Light;
use rtc_rust::material::Material;
use rtc_rust::matrix::Matrix4;
use rtc_rust::tuple::Tuple;
use rtc_rust::world::{World, hexagon, view_transform};
use rtc_rust::world::{read_stats, reset_stats};

fn main() {
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

    let from = Tuple::point(0.0, 2.5, -5.0);
    let to = Tuple::point(0.0, 0.0, 0.0);
    let up = Tuple::vector(0.0, 1.0, 0.0);
    let camera = Camera::new(1000, 1000, std::f64::consts::PI / 3.0)
        .with_transform(view_transform(from, to, up));

    let start = std::time::Instant::now();

    world.build_bounds();
    for _ in 0..100 {
        let _ = world.render_single(camera); // Camera is Copy, so the loop is fine
    }
    eprintln!("render took {:?}", start.elapsed()); // your baseline number

    reset_stats();
    let _image = world.render_single(camera); // single-threaded: NO contention
    let (nodes, prims) = read_stats();
    let rays = (camera.hsize * camera.vsize) as f64;
    eprintln!(
        "nodes visited: {nodes}  prim tests: {prims}  ({:.2} prim-tests/ray)",
        prims as f64 / rays
    );
    // let _ = _image.write_ppm("hexagons.ppm");      // optional; comment out to keep I/O off the graph
}
