//! World
//!
//! This module defines what rays can hit.
//!

use rayon::prelude::*;

use crate::camera::Camera;
use crate::canvas::Canvas;
use crate::intersection::{Intersection, Intersections};
use crate::light::Light;
use crate::log::*;
use crate::material::Material;
use crate::matrix::Matrix4;
use crate::ray::Ray;
use crate::shape::Shape;
use crate::shapes::sphere::Sphere;
use crate::tuple::Tuple;
use std::fmt;
use std::sync::atomic::{AtomicUsize, Ordering};

pub struct Computations<'a> {
    pub t: f64,
    pub object: &'a dyn Shape,
    pub point: Tuple,
    pub over_point: Tuple,
    pub eyev: Tuple,
    pub normalv: Tuple,
    pub inside: bool,
}

/// Encapsulating some precomputed information relating to the intersection.
///
/// Precomputes the following:
/// * the point (in world space) where the intersection occurred
/// * the eye vector (pointing back toward the eye or camera)
/// * the normal vector
///
pub fn prepare_computations<'a>(
    intersection: Intersection,
    ray: &Ray,
    shape: &'a dyn Shape,
) -> Computations<'a> {
    let t = intersection.t;
    let point = ray.position(t);
    let eyev = -ray.direction;
    let mut normalv = shape.normal_at(point);
    let inside = if normalv.dot(eyev) < 0.0 {
        normalv = -normalv;
        true
    } else {
        false
    };
    let over_point = point + normalv * crate::math::EPSILON;

    Computations {
        t,
        object: shape,
        point,
        over_point,
        eyev,
        normalv,
        inside,
    }
}

pub struct World {
    pub shapes: Vec<Box<dyn Shape>>,
    pub light: Option<Light>,
    next_id: AtomicUsize,
}

impl World {
    pub fn add_shape(&mut self, mut shape: Box<dyn Shape>) {
        shape.set_id(self.next_id.fetch_add(1, Ordering::Relaxed));
        self.shapes.push(shape);
    }

    pub fn set_light(&mut self, light: Light) {
        self.light = Some(light)
    }

    pub fn default_world() -> Self {
        let mut w = Self::new();

        let light = Light::point_light(
            Tuple::point(-10.0, 10.0, -10.0),
            Tuple::color(1.0, 1.0, 1.0),
        );

        w.light = Some(light);

        let mut m = Material::new();
        m.color = Tuple::color(0.8, 1.0, 0.6);
        m.diffuse = 0.7;
        m.specular = 0.2;
        let mut s1 = Sphere::new();
        s1.set_material(m);
        w.add_shape(Box::new(s1));

        let mut s2 = Sphere::new();
        s2.set_transform(Matrix4::scaling(0.5, 0.5, 0.5));
        w.add_shape(Box::new(s2));

        w
    }

    pub fn color_at(&self, ray: &Ray) -> Tuple {
        let xs = self.intersect(ray);
        match xs.hit() {
            None => Tuple::color(0.0, 0.0, 0.0),
            Some(hit) => match self.shapes.iter().find(|s| s.id() == hit.object_id) {
                None => Tuple::color(0.0, 0.0, 0.0),
                Some(shape) => {
                    let comps = prepare_computations(hit, ray, shape.as_ref());
                    self.shade_hit(&comps)
                }
            },
        }
    }

    pub fn intersect(&self, ray: &Ray) -> Intersections {
        let mut xs = Intersections::new();
        for shape in &self.shapes {
            for i in shape.intersect(ray).iter() {
                xs.push(i);
            }
        }
        xs
    }

    pub fn is_shadowed(&self, point: Tuple) -> bool {
        let light = match self.light {
            Some(l) => l,
            None => return false,
        };
        let v = light.position - point;
        let distance = Tuple::magnitude(v);
        let direction = Tuple::normalize(v);
        let r = Ray::new(point, direction);
        let intersections = self.intersect(&r);
        let hit = intersections.hit();
        let h = match hit {
            Some(x) => x,
            None => return false,
        };

        h.t < distance
    }

    pub fn new() -> Self {
        Self {
            shapes: Vec::new(),
            light: None,
            next_id: AtomicUsize::new(1),
        }
    }

    pub fn shade_hit(&self, comps: &Computations) -> Tuple {
        let shadowed = self.is_shadowed(comps.over_point);
        match self.light {
            Some(light) => light.lighting(
                *comps.object.material(),
                comps.over_point,
                comps.eyev,
                comps.normalv,
                shadowed,
            ),
            None => Tuple::color(0.0, 0.0, 0.0),
        }
    }

    pub fn render_single(&self, camera: Camera) -> Canvas {
        let mut image = Canvas::new(camera.hsize, camera.vsize);

        for y in 0..camera.vsize {
            for x in 0..camera.hsize {
                let ray = camera.ray_for_pixel(x, y);
                let color = self.color_at(&ray);
                image.write_pixel(x, y, color);
            }
        }

        image
    }

    pub fn render_parallel(&self, camera: Camera) -> Canvas {
        let width = camera.hsize;
        let height = camera.vsize;

        let mut pixels = vec![Tuple::color(0.0, 0.0, 0.0); width * height];
        pixels.par_iter_mut().enumerate().for_each(|(i, pixel)| {
            let x = i % width;
            let y = i / width;
            *pixel = self.color_at(&camera.ray_for_pixel(x, y));
        });

        let mut image = Canvas::new(width, height);
        for (i, color) in pixels.into_iter().enumerate() {
            image.write_pixel(i % width, i / width, color);
        }
        image
    }

    pub fn render(&self, camera: Camera) -> Canvas {
        // self.render_single(camera) // change this one line to switch
        self.render_parallel(camera) // change this one line to switch
    }
}

impl Default for World {
    fn default() -> Self {
        Self::new()
    }
}

impl fmt::Display for World {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.light {
            Some(l) => writeln!(
                f,
                "{}World{} light: {}{}{}",
                Color::Yellow,
                Color::Reset,
                Color::Green,
                l,
                Color::Reset
            )?,
            None => writeln!(
                f,
                "{}World{} light: {}none{}",
                Color::Yellow,
                Color::Reset,
                Color::Red,
                Color::Reset
            )?,
        }
        writeln!(
            f,
            "  shapes: {}{}{}",
            Color::Green,
            self.shapes.len(),
            Color::Reset
        )?;
        Ok(())
    }
}

/// A transformation matrix—like scaling, rotation, and translation—that orients the world relative
/// to your eye, thus allowing you to line everything up and get exactly the shot that you need
pub fn view_transform(from: Tuple, to: Tuple, up: Tuple) -> Matrix4 {
    let forward = (to - from).normalize();
    let upn = up.normalize();
    let left = forward.cross(upn);
    let true_up = left.cross(forward);

    let orientation = Matrix4::new([
        [left.x, left.y, left.z, 0.0],
        [true_up.x, true_up.y, true_up.z, 0.0],
        [-forward.x, -forward.y, -forward.z, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]);

    orientation * Matrix4::translation(-from.x, -from.y, -from.z)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::intersection::Intersection;
    use crate::math::approx_eq;
    use crate::shape::Shape;
    use crate::{loge, logi, tuple::Tuple};

    /// Chap 7 - Creating a world
    #[test]
    fn test_chap_7_1() -> Result<(), String> {
        let w = World::new();
        let chk = w.shapes.is_empty() && w.light.is_none();
        if chk {
            Ok(())
        } else {
            Err("Creating a world".into())
        }
    }

    /// Chap 7 - The default world
    #[test]
    fn test_chap_7_2() -> Result<(), String> {
        let w = World::default_world();

        let chk = !w.shapes.is_empty();
        let chk = chk
            && w.light.unwrap().approx_eq(Light::point_light(
                Tuple::point(-10.0, 10.0, -10.0),
                Tuple::color(1.0, 1.0, 1.0),
            ));

        let chk = chk && w.shapes.len() == 2;
        let chk = chk
            && w.shapes[0]
                .material()
                .color
                .approx_eq(Tuple::color(0.8, 1.0, 0.6));

        // The transform should be a scaling of 0.5 in all directions.
        let xform = w.shapes[1].transform();
        let chk = chk && approx_eq(xform[(0, 0)], 0.5);
        let chk = chk && approx_eq(xform[(1, 1)], 0.5);
        let chk = chk && approx_eq(xform[(2, 2)], 0.5);

        // The id's should be sphere 1 and sphere 2.
        let chk = chk && w.shapes[0].id() == 1;
        let chk = chk && w.shapes[1].id() == 2;

        if chk {
            Ok(())
        } else {
            loge!("test_chap_7_2", "world:{}", w);
            loge!(
                "test_chap_7_2",
                "w.shapes[0].material().color:{}",
                w.shapes[0].material().color
            );
            loge!("test_chap_7_2", "shape len : {}", w.shapes.len());
            loge!("test_chap_7_2", "shape id 1: {}", w.shapes[0].id());
            loge!("test_chap_7_2", "shape id 2: {}", w.shapes[1].id());
            loge!("test_chap_7_2", "xform[(0,0)]: {}", xform[(0, 0)]);
            loge!("test_chap_7_2", "xform[(1,1)]: {}", xform[(1, 1)]);
            loge!("test_chap_7_2", "xform[(2,2)]: {}", xform[(2, 2)]);
            Err("The default world".into())
        }
    }

    /// Chap 7 - Intersect a world with a ray
    #[test]
    fn test_chap_7_3() -> Result<(), String> {
        let w = World::default_world();

        logi!("test_chap_7_3", "world:{}", w);

        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs = w.intersect(&r);

        let chk = xs.count() == 4;
        let chk = chk && approx_eq(xs[0].t, 4.0);
        let chk = chk && approx_eq(xs[1].t, 4.5);
        let chk = chk && approx_eq(xs[2].t, 5.5);
        let chk = chk && approx_eq(xs[3].t, 6.0);
        if chk {
            Ok(())
        } else {
            Err("Intersect a world with a ray".into())
        }
    }

    /// Chap 7 - Precomputing the state of an intersection
    #[test]
    fn test_chap_7_4() -> Result<(), String> {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Sphere::new();
        let i = Intersection::new(4.0, shape.id());
        let comps = prepare_computations(i, &r, &shape as &dyn Shape);

        let chk = approx_eq(comps.t, i.t);
        let chk = chk && comps.object.id() == shape.id();
        let chk = chk && comps.point.approx_eq(Tuple::point(0.0, 0.0, -1.0));
        let chk = chk && comps.eyev.approx_eq(Tuple::vector(0.0, 0.0, -1.0));
        let chk = chk && comps.normalv.approx_eq(Tuple::vector(0.0, 0.0, -1.0));
        let chk = chk && !comps.inside;
        if chk {
            Ok(())
        } else {
            Err("Precomputing the state of an intersection".into())
        }
    }

    /// Chap x - The hit, when an intersection occurs on the outside
    #[test]
    fn test_chap_7_5() -> Result<(), String> {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Sphere::new();
        let i = Intersection::new(4.0, shape.id());
        let comps = prepare_computations(i, &r, &shape as &dyn Shape);

        let chk = !comps.inside;
        if chk {
            Ok(())
        } else {
            Err("The hit, when an intersection occurs on the outside".into())
        }
    }

    /// Chap x - The hit, when an intersection occurs on the inside
    #[test]
    fn test_chap_7_6() -> Result<(), String> {
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Sphere::new();
        let i = Intersection::new(1.0, shape.id());
        let comps = prepare_computations(i, &r, &shape as &dyn Shape);

        let chk = Tuple::point(0.0, 0.0, 1.0).approx_eq(comps.point);
        let chk = chk && Tuple::vector(0.0, 0.0, -1.0).approx_eq(comps.eyev);
        let chk = chk && comps.inside;
        // normal would have been (0, 0, 1), but is inverted!
        let chk = chk && Tuple::vector(0.0, 0.0, -1.0).approx_eq(comps.normalv);
        if chk {
            Ok(())
        } else {
            Err("The hit, when an intersection occurs on the inside".into())
        }
    }

    /// Chap 7 - Shading an intersection
    #[test]
    fn test_chap_7_7() -> Result<(), String> {
        let w = World::default_world();
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = w.shapes[0].as_ref();
        let i = Intersection::new(4.0, shape.id());
        let comps = prepare_computations(i, &r, shape as &dyn Shape);
        let c = w.shade_hit(&comps);

        let chk = c.approx_eq(Tuple::color(0.380661193081, 0.475826491351, 0.285495894811));
        if chk {
            Ok(())
        } else {
            loge!("test_chap_7_7", "c: {}", c);
            Err("Shading an intersection".into())
        }
    }

    /// Chap 7 - Shading an intersection from the inside
    #[test]
    fn test_chap_7_8() -> Result<(), String> {
        let mut w = World::default_world();

        w.light = Some(Light::point_light(
            Tuple::point(0.0, 0.25, 0.0),
            Tuple::color(1.0, 1.0, 1.0),
        ));

        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = w.shapes[1].as_ref();
        let i = Intersection::new(0.5, shape.id());
        let comps = prepare_computations(i, &r, shape as &dyn Shape);
        let c = w.shade_hit(&comps);

        let chk = c.approx_eq(Tuple::color(0.904984472083, 0.904984472083, 0.904984472083));
        if chk {
            Ok(())
        } else {
            loge!("test_chap_7_8", "c: {}", c);
            Err("Shading an intersection from the inside".into())
        }
    }

    /// Chap 7 - The color when a ray misses
    #[test]
    fn test_chap_7_9() -> Result<(), String> {
        let w = World::default_world();
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 1.0, 0.0));
        let c = w.color_at(&r);
        let chk = c.approx_eq(Tuple::color(0.0, 0.0, 0.0));
        if chk {
            Ok(())
        } else {
            Err("The color when a ray misses".into())
        }
    }
    /// Chap 7 - The color when a ray hits
    #[test]
    fn test_chap_7_10() -> Result<(), String> {
        let w = World::default_world();
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let c = w.color_at(&r);
        let chk = c.approx_eq(Tuple::color(0.380661193081, 0.475826491351, 0.285495894811));
        if chk {
            Ok(())
        } else {
            loge!("test_chap_7_10", "c: {}", c);
            Err("The color when a ray hits".into())
        }
    }

    /// Chap 7 - The color with an intersection behind the ray
    #[test]
    fn test_chap_7_11() -> Result<(), String> {
        let mut w = World::default_world();

        {
            let outer = w.shapes[0].as_mut();
            let mut mat = *outer.material();
            mat.ambient = 1.0;
            outer.set_material(mat);
        }

        {
            let inner = w.shapes[1].as_mut();
            let mut mat = *inner.material();
            mat.ambient = 1.0;
            inner.set_material(mat);
        }

        let r = Ray::new(Tuple::point(0.0, 0.0, 0.75), Tuple::vector(0.0, 0.0, -1.0));
        let c = w.color_at(&r);

        let inner = w.shapes[1].as_ref();
        let chk = c.approx_eq(inner.material().color);
        if chk {
            Ok(())
        } else {
            Err("The color with an intersection behind the ray".into())
        }
    }

    /// Chap 7 - The transformation matrix for the default orientation
    #[test]
    fn test_chap_7_12() -> Result<(), String> {
        let from = Tuple::point(0.0, 0.0, 0.0);
        let to = Tuple::point(0.0, 0.0, -1.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let t = view_transform(from, to, up);
        let chk = t.approx_eq(Matrix4::identity());
        if chk {
            Ok(())
        } else {
            Err("The transformation matrix for the default orientation".into())
        }
    }

    /// Chap 7 - A view transformation matrix looking in positive z direction
    #[test]
    fn test_chap_7_13() -> Result<(), String> {
        let from = Tuple::point(0.0, 0.0, 0.0);
        let to = Tuple::point(0.0, 0.0, 1.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let t = view_transform(from, to, up);
        let m = Matrix4::scaling(-1.0, 1.0, -1.0);
        let chk = t.approx_eq(m);
        if chk {
            Ok(())
        } else {
            loge!("test_chap_7_13", "\nt:{}\nm:{}", t, m);
            Err("A view transformation matrix looking in positive z direction".into())
        }
    }

    /// Chap 7 - The view transformation moves the world
    #[test]
    fn test_chap_7_14() -> Result<(), String> {
        let from = Tuple::point(0.0, 0.0, 8.0);
        let to = Tuple::point(0.0, 0.0, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let t = view_transform(from, to, up);
        let m = Matrix4::translation(0.0, 0.0, -8.0);
        let chk = t.approx_eq(m);
        if chk {
            Ok(())
        } else {
            Err("The view transformation moves the world".into())
        }
    }

    /// Chap 7 - An arbitrary view transformation
    #[test]
    fn test_chap_7_15() -> Result<(), String> {
        let from = Tuple::point(1.0, 3.0, 2.0);
        let to = Tuple::point(4.0, -2.0, 8.0);
        let up = Tuple::vector(1.0, 1.0, 0.0);
        let t = view_transform(from, to, up);
        let m = Matrix4::new([
            [
                -0.507092552837110,
                0.507092552837110,
                0.676123403782813,
                -2.366431913239846,
            ],
            [
                0.767715933859680,
                0.606091526731326,
                0.121218305346265,
                -2.828427124746189,
            ],
            [
                -0.358568582800318,
                0.597614304667197,
                -0.717137165600636,
                0.000000000000000,
            ],
            [
                0.000000000000000,
                0.000000000000000,
                0.000000000000000,
                1.000000000000000,
            ],
        ]);
        let chk = t.approx_eq(m);
        if chk {
            Ok(())
        } else {
            loge!("test_chap_7_15", "t:{}", t);
            Err("An arbitrary view transformation".into())
        }
    }

    /// Chap 7 - Rendering a world with a camera
    #[test]
    fn test_chap_7_22() -> Result<(), String> {
        let w = World::default_world();
        let from = Tuple::point(0.0, 0.0, -5.0);
        let to = Tuple::point(0.0, 0.0, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let transform = view_transform(from, to, up);

        // Option 1 - unmutable camera
        let c = Camera::new(11, 11, std::f64::consts::PI / 2.0).with_transform(transform);

        // Option 2 - mutable camera
        // let mut c = Camera::new(11, 11, std::f64::consts::PI / 2.0);
        // c.set_transform(&view_transform(from, to, up));

        let image = w.render(c);
        let rc = image.write_ppm("test_chap_7_22.ppm");
        let chk = rc.is_ok();

        if chk {
            Ok(())
        } else {
            Err("Rendering a world with a camera".into())
        }
    }

    /// Chap x - Chapter 7 Putting It  Together
    #[test]
    fn test_chap_7_23_putting_it_together() -> Result<(), String> {
        let mut world = World::new();
        let light = Light::point_light(
            Tuple::point(-10.0, 10.0, -10.0),
            Tuple::color(1.0, 1.0, 1.0),
        );
        world.light = Some(light);

        let mut material = Material::new();
        material.color = Tuple::color(1.0, 0.9, 0.9);

        let mut floor = Sphere::new();
        floor.set_transform(Matrix4::scaling(10.0, 0.01, 10.0));
        material.diffuse = 0.7;
        material.specular = 0.3;
        floor.set_material(material);

        let mut left_wall = Sphere::new();
        left_wall.set_transform(
            Matrix4::translation(0.0, 0.0, 5.0)
                * Matrix4::rotation_y(-std::f64::consts::PI / 4.0)
                * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(10.0, 0.01, 10.0),
        );
        left_wall.set_material(*floor.material());

        let mut right_wall = Sphere::new();
        right_wall.set_transform(
            Matrix4::translation(0.0, 0.0, 5.0)
                * Matrix4::rotation_y(std::f64::consts::PI / 4.0)
                * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(10.0, 0.01, 10.0),
        );
        right_wall.set_material(*floor.material());

        let mut middle = Sphere::new();
        middle.set_transform(Matrix4::translation(-0.5, 1.0, 0.5));
        material.color = Tuple::color(0.1, 1.0, 0.5);
        material.diffuse = 0.7;
        material.specular = 0.3;
        middle.set_material(material);

        let mut right = Sphere::new();
        right.set_transform(Matrix4::translation(1.5, 0.5, -0.5) * Matrix4::scaling(0.5, 0.5, 0.5));
        material.color = Tuple::color(0.5, 1.0, 0.1);
        material.diffuse = 0.7;
        material.specular = 0.3;
        right.set_material(material);

        let mut left = Sphere::new();
        left.set_transform(
            Matrix4::translation(-1.5, 0.33, -0.75) * Matrix4::scaling(0.33, 0.33, 0.33),
        );
        material.color = Tuple::color(1.0, 0.8, 0.1);
        material.diffuse = 0.7;
        material.specular = 0.3;
        left.set_material(material);

        let from = Tuple::point(0.0, 1.5, -5.0);
        let to = Tuple::point(0.0, 1.0, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let transform = view_transform(from, to, up);

        let camera = Camera::new(4096, 3192, std::f64::consts::PI / 3.0).with_transform(transform);

        world.add_shape(Box::new(floor));
        world.add_shape(Box::new(left_wall));
        world.add_shape(Box::new(right_wall));
        world.add_shape(Box::new(middle));
        world.add_shape(Box::new(left));
        world.add_shape(Box::new(right));

        let image = world.render(camera);
        let rc = image.write_ppm("test_chap_7_23_putting_it_together.ppm");
        let chk = rc.is_ok();

        if chk {
            Ok(())
        } else {
            Err("Chapter 7 Putting It  Together".into())
        }
    }

    /// Chap 8 - There is no shadow when nothing is collinear with point and light
    #[test]
    fn test_chap_8_2() -> Result<(), String> {
        let w = World::default_world();
        let p = Tuple::point(0.0, 10.0, 0.0);
        let is_shadowed = w.is_shadowed(p);

        let chk = !is_shadowed;
        if chk {
            Ok(())
        } else {
            Err("There is no shadow when nothing is collinear with point and light".into())
        }
    }

    /// Chap 8 - The shadow when an object is between the point and the light
    #[test]
    fn test_chap_8_3() -> Result<(), String> {
        let w = World::default_world();
        let p = Tuple::point(10.0, -10.0, 10.0);
        let is_shadowed = w.is_shadowed(p);

        let chk = is_shadowed;
        if chk {
            Ok(())
        } else {
            Err("The shadow when an object is between the point and the light".into())
        }
    }

    /// Chap 8 - There is no shadow when an object is behind the light
    #[test]
    fn test_chap_8_5() -> Result<(), String> {
        let w = World::default_world();
        let p = Tuple::point(0.0, 10.0, 0.0);
        let is_shadowed = w.is_shadowed(p);

        let chk = !is_shadowed;
        if chk {
            Ok(())
        } else {
            Err("There is no shadow when an object is behind the light".into())
        }
    }

    /// Chap 8 - There is no shadow when an object is behind the point
    #[test]
    fn test_chap_8_6() -> Result<(), String> {
        let w = World::default_world();
        let p = Tuple::point(0.0, 10.0, 0.0);
        let is_shadowed = w.is_shadowed(p);

        let chk = !is_shadowed;
        if chk {
            Ok(())
        } else {
            Err("There is no shadow when an object is behind the point".into())
        }
    }

    /// Chap 8 - shade_hit() is given an intersection in shadow
    #[test]
    fn test_chap_8_7() -> Result<(), String> {
        let mut w = World::new();
        let light = Light::point_light(Tuple::point(0.0, 0.0, -10.0), Tuple::color(1.0, 1.0, 1.0));
        w.set_light(light);
        let s1 = Sphere::new();
        w.add_shape(Box::new(s1));
        let mut s2 = Sphere::new();
        s2.set_transform(Matrix4::translation(0.0, 0.0, 10.0));
        w.add_shape(Box::new(s2));
        let r = Ray::new(Tuple::point(0.0, 0.0, 5.0), Tuple::vector(0.0, 0.0, 1.0));
        let i = Intersection::new(4.0, s2.id());
        let comps = prepare_computations(i, &r, &s2 as &dyn Shape);
        let c = w.shade_hit(&comps);

        let chk = c.approx_eq(Tuple::color(0.1, 0.1, 0.1));
        if chk {
            Ok(())
        } else {
            Err("shade_hit() is given an intersection in shadow".into())
        }
    }

    /// Chap 8 - The hit should offset the point
    #[test]
    fn test_chap_8_8() -> Result<(), String> {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let mut s1 = Sphere::new();
        s1.set_transform(Matrix4::translation(0.0, 0.0, 1.0));
        let i = Intersection::new(5.0, s1.id());
        let comps = prepare_computations(i, &r, &s1 as &dyn Shape);
        let chk = comps.over_point.z < -crate::math::EPSILON / 2.0;
        let chk = chk && comps.point.z > comps.over_point.z;

        if chk {
            Ok(())
        } else {
            Err("The hit should offset the point".into())
        }
    }
}
