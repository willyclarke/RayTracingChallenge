// Renders a scene and overlays the AABB / BVH group boxes as a wireframe.
//
// The boxes are a 2D debug overlay: each group's 8 world-space corners are
// projected to screen and its 12 edges drawn on top of the render, so nested
// and divide()-created boxes are all visible regardless of occlusion.
//
// Run: cargo run --release --example wireframe   (writes wireframe.ppm)

use rtc_rust::camera::Camera;
use rtc_rust::canvas::Canvas;
use rtc_rust::light::Light;
use rtc_rust::material::Material;
use rtc_rust::matrix::Matrix4;
use rtc_rust::shape::Shape;
use rtc_rust::shapes::group::Group;
use rtc_rust::shapes::sphere::Sphere;
use rtc_rust::tuple::Tuple;
use rtc_rust::world::{World, view_transform};

/// Toggle the BVH subdivision so you can see divide() carve up the boxes.
const SHOW_BVH: bool = true;

fn build() -> (World, usize) {
    let mut w = World::new();
    w.light = Some(Light::point_light(
        Tuple::point(-10.0, 10.0, -10.0),
        Tuple::color(1.0, 1.0, 1.0),
    ));
    let g_id = w.add_shape(Box::new(Group::new()));
    let mut m = Material::new();
    m.color = Tuple::color(0.3, 0.6, 0.9);
    let n = 5;
    for i in 0..n {
        for j in 0..n {
            let x = (i - n / 2) as f64 * 1.6;
            let y = (j - n / 2) as f64 * 1.6;
            let mut s = Sphere::new();
            s.set_transform(Matrix4::translation(x, y, 0.0) * Matrix4::scaling(0.4, 0.4, 0.4));
            s.set_material(m.clone());
            w.add_child(g_id, Box::new(s));
        }
    }
    (w, g_id)
}

/// Project a world-space point to floating pixel coords (None if behind camera).
fn project(camera: &Camera, p: Tuple) -> Option<(i32, i32)> {
    let pc = camera.transform() * p; // world -> camera space
    if pc.z >= 0.0 {
        return None; // behind the eye (camera looks down -z)
    }
    let proj_x = -pc.x / pc.z; // onto the z = -1 image plane
    let proj_y = -pc.y / pc.z;
    let px = (camera.half_width - proj_x) / camera.pixel_size - 0.5;
    let py = (camera.half_height - proj_y) / camera.pixel_size - 0.5;
    Some((px.round() as i32, py.round() as i32))
}

/// Bresenham line, clipped to the canvas.
fn draw_line(c: &mut Canvas, mut x0: i32, mut y0: i32, x1: i32, y1: i32, color: Tuple) {
    let dx = (x1 - x0).abs();
    let dy = -(y1 - y0).abs();
    let sx = if x0 < x1 { 1 } else { -1 };
    let sy = if y0 < y1 { 1 } else { -1 };
    let mut err = dx + dy;
    loop {
        if x0 >= 0 && y0 >= 0 && (x0 as usize) < c.width() && (y0 as usize) < c.height() {
            c.write_pixel(x0 as usize, y0 as usize, color);
        }
        if x0 == x1 && y0 == y1 {
            break;
        }
        let e2 = 2 * err;
        if e2 >= dy {
            err += dy;
            x0 += sx;
        }
        if e2 <= dx {
            err += dx;
            y0 += sy;
        }
    }
}

fn main() {
    let (mut world, g_id) = build();
    if SHOW_BVH {
        world.divide(g_id, 4);
    }
    world.build_bounds();

    let camera = Camera::new(600, 600, std::f64::consts::PI / 3.0).with_transform(view_transform(
        Tuple::point(0.0, 0.0, -18.0),
        Tuple::point(0.0, 0.0, 0.0),
        Tuple::vector(0.0, 1.0, 0.0),
    ));

    let mut canvas = world.render_single(camera);

    // overlay the group boxes
    let wire = Tuple::color(1.0, 1.0, 0.0); // yellow
    const EDGES: [(usize, usize); 12] = [
        (0, 1), (2, 3), (4, 5), (6, 7), // x edges
        (0, 2), (1, 3), (4, 6), (5, 7), // y edges
        (0, 4), (1, 5), (2, 6), (3, 7), // z edges
    ];
    let boxes = world.group_world_boxes();
    for corners in &boxes {
        let scr: Vec<Option<(i32, i32)>> = corners.iter().map(|&c| project(&camera, c)).collect();
        for &(a, b) in &EDGES {
            if let (Some(pa), Some(pb)) = (scr[a], scr[b]) {
                draw_line(&mut canvas, pa.0, pa.1, pb.0, pb.1, wire);
            }
        }
    }

    let _ = canvas.write_ppm("wireframe.ppm");
    println!(
        "wrote wireframe.ppm  ({} group boxes drawn, BVH {})",
        boxes.len(),
        if SHOW_BVH { "on" } else { "off" }
    );
}
