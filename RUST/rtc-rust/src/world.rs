//! World
//!
//! This module defines what rays can hit.
//!

use rayon::prelude::*;

use crate::camera::Camera;
use crate::canvas::Canvas;
use crate::intersection::{Intersection, Intersections};
use crate::light::Light;
use crate::material::Material;
use crate::math::approx_eq;
use crate::matrix::Matrix4;
use crate::ray::Ray;
use crate::shape::Shape;
use crate::shapes::sphere::Sphere;
use crate::tuple::Tuple;
use crate::{log::*, tuple};
use std::fmt;
use std::sync::atomic::{AtomicUsize, Ordering};

pub struct Computations<'a> {
    pub t: f64,
    pub object: &'a dyn Shape,
    pub point: Tuple,
    pub over_point: Tuple,
    pub under_point: Tuple,
    pub eyev: Tuple,
    pub normalv: Tuple,
    pub reflectv: Tuple,
    pub inside: bool,
    pub n1: f64,
    pub n2: f64,
}

/// Encapsulating some precomputed information relating to the intersection.
///
/// Precomputes the following:
/// * the point (in world space) where the intersection occurred
/// * the eye vector (pointing back toward the eye or camera)
/// * the normal vector
///
pub fn prepare_computations_upto_chap10<'a>(
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
    let reflectv = ray.direction.reflect(normalv);
    let over_point = point + normalv * crate::math::EPSILON;
    let under_point = point - normalv * crate::math::EPSILON;

    Computations {
        t,
        object: shape,
        point,
        over_point,
        under_point,
        eyev,
        normalv,
        reflectv,
        inside,
        n1: 1.0,
        n2: 1.0,
    }
}

/// Schlick computation to approximate Fresnel's equation.
///
/// returns a number between 0 and 1, inclusive. This number is called the
/// reflectance and represents what fraction of the light is reflected, given the
/// surface information at the hit
///
pub fn schlick(comps: &Computations) -> f64 {
    // find the cosine of the angle between the eye and normal vectors
    let mut cos = comps.eyev.dot(comps.normalv);

    // total internal reflection can only occur if n1 > n2
    if comps.n1 > comps.n2 {
        let n = comps.n1 / comps.n2;
        let sin2_t = n * n * (1.0 - cos * cos);
        if sin2_t > 1.0 {
            return 1.0;
        }

        // compute cosine of theta_t using trig identity and
        // when n1 > n2, use cos(theta_t) instead
        cos = (1.0 - sin2_t).sqrt();
    }

    let r0 = ((comps.n1 - comps.n2) / (comps.n1 + comps.n2)).powf(2.0);

    r0 + (1.0 - r0) * (1.0 - cos).powf(5.0)
}

pub fn prepare_computations<'a>(
    intersection: Intersection,
    ray: &Ray,
    shapes: &'a [Box<dyn Shape>],
    xs: &Intersections,
) -> Computations<'a> {
    let shape_by_id = |id: usize| -> &dyn Shape {
        shapes
            .iter()
            .find(|s| s.id() == id)
            .expect("shape id must exist in shapes")
            .as_ref()
    };

    let shape = shape_by_id(intersection.object_id);
    let refractive_index_of =
        |id: usize| -> f64 { shape_by_id(id).data().material.refractive_index };

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

    let reflectv = ray.direction.reflect(normalv);
    let over_point = point + normalv * crate::math::EPSILON;
    let under_point = point - normalv * crate::math::EPSILON;

    let hit = &intersection;
    let mut n1 = 1.0;
    let mut n2 = 1.0;

    let mut containers = [0usize; 32]; // stack-allocated, no heap
    let last = |c: &[usize], len: usize| -> Option<usize> {
        if len == 0 { None } else { Some(c[len - 1]) }
    };
    let mut len = 0;

    for i in xs.iter() {
        let is_hit = approx_eq(i.t, hit.t) && i.object_id == hit.object_id;

        if is_hit {
            n1 = last(&containers, len).map_or(1.0, refractive_index_of);
        };

        // toggle membership: already inside => we're EXITING; otherwise ENTERING
        if let Some(pos) = containers[..len].iter().position(|&id| id == i.object_id) {
            containers.copy_within(pos + 1..len, pos); // shift left, preserve order
            len -= 1;
        } else {
            debug_assert!(len < containers.len(), "container overflow");
            containers[len] = i.object_id;
            len += 1;
        }

        // n2 = material the ray is ENTERING (last container, after the toggle)
        if is_hit {
            n2 = last(&containers, len).map_or(1.0, refractive_index_of);
            break; // (4) stop at the hit
        }
    }

    Computations {
        t,
        object: shape,
        point,
        over_point,
        under_point,
        eyev,
        normalv,
        reflectv,
        inside,
        n1,
        n2,
    }
}

pub struct World {
    pub shapes: Vec<Box<dyn Shape>>,
    pub light: Option<Light>,
    next_id: AtomicUsize,
}

impl World {
    /// Increment the shape id and add the shape to world.
    /// # Examples
    /// ```
    /// use rtc_rust::world::World;
    /// use rtc_rust::shapes::sphere::Sphere;
    ///
    /// let mut w = World::default_world();
    /// let ball = Sphere::new();
    /// let ball_id = w.add_shape(Box::new(ball));
    ///
    /// assert!(ball_id > 0);
    /// ```
    pub fn add_shape(&mut self, mut shape: Box<dyn Shape>) -> usize {
        let id = self.next_id.fetch_add(1, Ordering::Relaxed);
        shape.set_id(id);
        self.shapes.push(shape);
        id
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

    pub fn color_at(&self, ray: &Ray, remaining: i32) -> Tuple {
        let xs = self.intersect(ray);
        match xs.hit() {
            None => Tuple::color(0.0, 0.0, 0.0),
            Some(hit) => match self.shapes.iter().find(|s| s.id() == hit.object_id) {
                None => Tuple::color(0.0, 0.0, 0.0),
                Some(_shape) => {
                    let comps = prepare_computations(hit, ray, &self.shapes, &xs);
                    self.shade_hit(&comps, remaining)
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

    // pub fn shade_hit_old(&self, comps: &Computations) -> Tuple {
    //     let shadowed = self.is_shadowed(comps.over_point);
    //     match self.light {
    //         Some(light) => light.lighting_old(
    //             comps.object.material(),
    //             comps.over_point,
    //             comps.eyev,
    //             comps.normalv,
    //             shadowed,
    //         ),
    //         None => Tuple::color(0.0, 0.0, 0.0),
    //     }
    // }

    pub fn shade_hit(&self, comps: &Computations, remaining: i32) -> Tuple {
        let shadowed = self.is_shadowed(comps.over_point);
        match self.light {
            Some(light) => {
                let surface = light.lighting(
                    comps.object,
                    comps.over_point,
                    comps.eyev,
                    comps.normalv,
                    shadowed,
                );
                let reflected = self.reflected_color(comps, remaining);
                let refracted = self.refracted_color(comps, remaining);
                if comps.object.material().reflective > 0.0
                    && comps.object.material().transparency > 0.0
                {
                    let reflectance = schlick(comps);
                    return surface + reflected * reflectance + refracted * (1.0 - reflectance);
                }
                surface + reflected + refracted
            }
            None => Tuple::color(0.0, 0.0, 0.0),
        }
    }

    pub fn render_single(&self, camera: Camera) -> Canvas {
        let mut image = Canvas::new(camera.hsize, camera.vsize);

        for y in 0..camera.vsize {
            for x in 0..camera.hsize {
                let ray = camera.ray_for_pixel(x, y);
                let color = self.color_at(&ray, 10);
                image.write_pixel(x, y, color);
            }
        }

        image
    }

    pub fn reflected_color(&self, comps: &Computations, remaining: i32) -> Tuple {
        if remaining <= 0 {
            return tuple::colors::BLACK;
        }

        if approx_eq(comps.object.material().reflective, 0.0) {
            return tuple::colors::BLACK;
        }

        let reflect_ray = Ray::new(comps.over_point, comps.reflectv);
        let color = self.color_at(&reflect_ray, remaining - 1);
        color * comps.object.material().reflective
    }

    pub fn refracted_color(&self, comps: &Computations, remaining: i32) -> Tuple {
        if remaining <= 0 {
            return tuple::colors::BLACK;
        }

        if approx_eq(comps.object.material().transparency, 0.0) {
            return Tuple::color(0.0, 0.0, 0.0);
        }

        // Handle total internal reflection.
        // Find the ratio of first index of refraction to the second.
        // (Yup, this is inverted from the definition of Snell's Law.)
        let n_ratio = comps.n1 / comps.n2;

        // cos(theta_i) is the same as the dot product of the two vectors
        let cos_i = comps.eyev.dot(comps.normalv);

        // Find sin(theta_t)^2 via trigonometric identity
        let sin2_t = n_ratio * n_ratio * (1.0 - cos_i * cos_i);

        // Return black when there is total internal reflection.
        if sin2_t > 1.0 {
            return tuple::colors::BLACK;
        }

        // Find cos(theta_t) via trigonometric identity
        let cos_t = (1.0 - sin2_t).sqrt();

        // Compute the direction of the refracted ray
        let direction = comps.normalv * (n_ratio * cos_i - cos_t) - comps.eyev * n_ratio;

        // Create the refracted ray
        let refracted_ray = Ray::new(comps.under_point, direction);

        // Find the color of the refracted ray, making sure to multiply
        // by the transparency value to account for any opacity
        self.color_at(&refracted_ray, remaining - 1) * comps.object.material().transparency
    }

    pub fn render_parallel(&self, camera: Camera) -> Canvas {
        let width = camera.hsize;
        let height = camera.vsize;

        let mut pixels = vec![Tuple::color(0.0, 0.0, 0.0); width * height];
        pixels.par_iter_mut().enumerate().for_each(|(i, pixel)| {
            let x = i % width;
            let y = i / width;
            *pixel = self.color_at(&camera.ray_for_pixel(x, y), 10);
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
    use crate::math::{EPSILON, approx_eq};
    use crate::pattern::Pattern;
    use crate::patterns::blendedpattern::BlendedPattern;
    use crate::patterns::checkerspattern::CheckersPattern;
    use crate::patterns::gradientpattern::GradientPattern;
    use crate::patterns::nestedpattern::NestedPattern;
    use crate::patterns::ringpattern::RingPattern;
    use crate::patterns::stripepattern::*;
    use crate::patterns::testpattern::TestPattern;
    use crate::shape::Shape;
    use crate::shapes::cone::Cone;
    use crate::shapes::cube::Cube;
    use crate::shapes::cylinder::Cylinder;
    use crate::shapes::plane::Plane;
    use crate::tuple::colors::*;
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
        let comps = prepare_computations_upto_chap10(i, &r, &shape as &dyn Shape);

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
        let comps = prepare_computations_upto_chap10(i, &r, &shape as &dyn Shape);

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
        let comps = prepare_computations_upto_chap10(i, &r, &shape as &dyn Shape);

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
        let comps = prepare_computations_upto_chap10(i, &r, shape as &dyn Shape);
        let c = w.shade_hit(&comps, 10);

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
        let comps = prepare_computations_upto_chap10(i, &r, shape as &dyn Shape);
        let c = w.shade_hit(&comps, 10);

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
        let c = w.color_at(&r, 10);
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
        let c = w.color_at(&r, 10);
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
            let mut mat = outer.material().clone();
            mat.ambient = 1.0;
            outer.set_material(mat);
        }

        {
            let inner = w.shapes[1].as_mut();
            let mut mat = inner.material().clone();
            mat.ambient = 1.0;
            inner.set_material(mat);
        }

        let r = Ray::new(Tuple::point(0.0, 0.0, 0.75), Tuple::vector(0.0, 0.0, -1.0));
        let c = w.color_at(&r, 10);

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

    /// Chap 7 - Chapter 7 Putting It  Together
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
        floor.set_material(material.clone());

        let mut left_wall = Sphere::new();
        left_wall.set_transform(
            Matrix4::translation(0.0, 0.0, 5.0)
                * Matrix4::rotation_y(-std::f64::consts::PI / 4.0)
                * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(10.0, 0.01, 10.0),
        );
        left_wall.set_material(floor.material().clone());

        let mut right_wall = Sphere::new();
        right_wall.set_transform(
            Matrix4::translation(0.0, 0.0, 5.0)
                * Matrix4::rotation_y(std::f64::consts::PI / 4.0)
                * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(10.0, 0.01, 10.0),
        );
        right_wall.set_material(floor.material().clone());

        let mut middle = Sphere::new();
        middle.set_transform(Matrix4::translation(-0.5, 1.0, 0.5));
        material.color = Tuple::color(0.1, 1.0, 0.5);
        material.diffuse = 0.7;
        material.specular = 0.3;
        middle.set_material(material.clone());

        let mut right = Sphere::new();
        right.set_transform(Matrix4::translation(1.5, 0.5, -0.5) * Matrix4::scaling(0.5, 0.5, 0.5));
        material.color = Tuple::color(0.5, 1.0, 0.1);
        material.diffuse = 0.7;
        material.specular = 0.3;
        right.set_material(material.clone());

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

        let camera = Camera::new(100, 50, std::f64::consts::PI / 3.0).with_transform(transform);

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
        w.add_shape(Box::new(s2.clone()));
        let r = Ray::new(Tuple::point(0.0, 0.0, 5.0), Tuple::vector(0.0, 0.0, 1.0));
        let i = Intersection::new(4.0, s2.id());
        let comps = prepare_computations_upto_chap10(i, &r, &s2 as &dyn Shape);
        let c = w.shade_hit(&comps, 10);

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
        let comps = prepare_computations_upto_chap10(i, &r, &s1 as &dyn Shape);
        let chk = comps.over_point.z < -crate::math::EPSILON / 2.0;
        let chk = chk && comps.point.z > comps.over_point.z;

        if chk {
            Ok(())
        } else {
            Err("The hit should offset the point".into())
        }
    }

    /// Chap 9 - Chapter 9 Putting It  Together
    #[test]
    fn test_chap_9_6_putting_it_together() -> Result<(), String> {
        let mut world = World::new();
        let light = Light::point_light(
            Tuple::point(-10.0, 10.0, -10.0),
            Tuple::color(1.0, 1.0, 1.0),
        );
        world.light = Some(light);

        let mut material = Material::new();
        material.color = Tuple::color(1.0, 0.9, 0.9);

        let mut floor = Plane::new();
        floor.set_transform(Matrix4::scaling(10.0, 0.01, 10.0));
        material.diffuse = 0.7;
        material.specular = 0.3;
        floor.set_material(material.clone());

        let mut left_wall = Plane::new();
        left_wall.set_transform(
            Matrix4::translation(0.0, 0.0, 5.0)
                * Matrix4::rotation_y(-std::f64::consts::PI / 4.0)
                * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(10.0, 0.01, 10.0),
        );
        left_wall.set_material(floor.material().clone());

        let mut right_wall = Plane::new();
        right_wall.set_transform(
            Matrix4::translation(0.0, 0.0, 5.0)
                * Matrix4::rotation_y(std::f64::consts::PI / 4.0)
                * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(10.0, 0.01, 10.0),
        );
        right_wall.set_material(floor.material().clone());

        let mut middle = Sphere::new();
        middle.set_transform(Matrix4::translation(-0.5, 1.0, 0.5));
        material.color = Tuple::color(0.1, 1.0, 0.5);
        material.diffuse = 0.7;
        material.specular = 0.3;
        middle.set_material(material.clone());

        let mut right = Sphere::new();
        right.set_transform(Matrix4::translation(1.5, 0.5, -0.5) * Matrix4::scaling(0.5, 0.5, 0.5));
        material.color = Tuple::color(0.5, 1.0, 0.1);
        material.diffuse = 0.7;
        material.specular = 0.3;
        right.set_material(material.clone());

        let mut left = Sphere::new();
        left.set_transform(
            Matrix4::translation(-1.5, 0.33, -0.75) * Matrix4::scaling(0.33, 0.33, 0.33),
        );
        material.color = Tuple::color(1.0, 0.8, 0.1);
        material.diffuse = 0.7;
        material.specular = 0.3;
        left.set_material(material);

        let from = Tuple::point(0.0, 1.5, -50.0);
        let to = Tuple::point(0.0, 1.0, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let transform = view_transform(from, to, up);

        let camera = Camera::new(100, 50, std::f64::consts::PI / 1.1).with_transform(transform);

        world.add_shape(Box::new(floor));
        world.add_shape(Box::new(left_wall));
        world.add_shape(Box::new(right_wall));
        world.add_shape(Box::new(middle));
        world.add_shape(Box::new(left));
        world.add_shape(Box::new(right));

        let image = world.render(camera);
        let rc = image.write_ppm("test_chap_9_6_putting_it_together.ppm");
        let chk = rc.is_ok();

        if chk {
            Ok(())
        } else {
            Err("Chapter 7 Putting It  Together".into())
        }
    }

    /// Chap 10 - Stripes with an object transformation
    #[test]
    fn test_chap_10_6() -> Result<(), String> {
        let mut s1 = Sphere::new();
        s1.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));

        let pattern = StripePattern::new(WHITE, BLACK);
        let world_point = Tuple::point(1.5, 0.0, 0.0);

        let c = pattern.color_at_shape(&s1, world_point);

        let chk = c.approx_eq(WHITE);
        if chk {
            Ok(())
        } else {
            Err("Stripes with an object transformation".into())
        }
    }

    /// Chap 10 - Stripes with a pattern transformation
    #[test]
    fn test_chap_10_7() -> Result<(), String> {
        let s1 = Sphere::new();

        let mut pattern = StripePattern::new(WHITE, BLACK);
        pattern.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));

        let world_point = Tuple::point(1.5, 0.0, 0.0);

        let c = pattern.color_at_shape(&s1, world_point);

        let chk = c.approx_eq(WHITE);

        if chk {
            Ok(())
        } else {
            Err("Stripes with a pattern transformation".into())
        }
    }

    /// Chap 10 - Stripes with both an object and a pattern transformation
    #[test]
    fn test_chap_10_8() -> Result<(), String> {
        let mut s1 = Sphere::new();
        s1.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));

        let mut pattern = StripePattern::new(WHITE, BLACK);
        pattern.set_transform(Matrix4::translation(0.5, 0.0, 0.0));

        let world_point = Tuple::point(2.5, 0.0, 0.0);

        let c = pattern.color_at_shape(&s1, world_point);

        let chk = c.approx_eq(WHITE);
        if chk {
            Ok(())
        } else {
            Err("Stripes with both an object and a pattern transformation".into())
        }
    }

    /// Chap 10 - A pattern with an object transformation
    #[test]
    fn test_chap_10_11() -> Result<(), String> {
        let mut shape = Sphere::new();
        shape.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));
        let pattern = TestPattern::new();

        let world_point = Tuple::point(2.0, 3.0, 4.0);

        let c = pattern.color_at_shape(&shape, world_point);

        let chk = c.approx_eq(Tuple::color(1.0, 1.5, 2.0));
        if chk {
            Ok(())
        } else {
            Err("A pattern with an object transformation".into())
        }
    }

    /// Chap x - A pattern with a pattern transformation
    #[test]
    fn test_chap_10_12() -> Result<(), String> {
        let shape = Sphere::new();
        let mut pattern = TestPattern::new();
        pattern.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));

        let world_point = Tuple::point(2.0, 3.0, 4.0);

        let c = pattern.color_at_shape(&shape, world_point);

        let chk = c.approx_eq(Tuple::color(1.0, 1.5, 2.0));

        if chk {
            Ok(())
        } else {
            Err("A pattern with a pattern transformation".into())
        }
    }

    /// Chap 10 - A pattern with both an object and a pattern transformation
    #[test]
    fn test_chap_10_13() -> Result<(), String> {
        let mut shape = Sphere::new();
        shape.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));
        let mut pattern = TestPattern::new();
        pattern.set_transform(Matrix4::translation(0.5, 1.0, 1.5));

        let world_point = Tuple::point(2.5, 3.0, 3.5);

        let c = pattern.color_at_shape(&shape, world_point);

        let chk = c.approx_eq(Tuple::color(0.75, 0.5, 0.25));
        if chk {
            Ok(())
        } else {
            Err("A pattern with both an object and a pattern transformation".into())
        }
    }

    /// Chap 10 - A gradient linearly interpolates between colors
    #[test]
    fn test_chap_10_14() -> Result<(), String> {
        let pattern = GradientPattern::new(WHITE, BLACK);
        let chk_0 = pattern.color_at(Tuple::point(0.0, 0.0, 0.0));
        let chk_25 = pattern.color_at(Tuple::point(0.25, 0.0, 0.0));
        let chk_50 = pattern.color_at(Tuple::point(0.50, 0.0, 0.0));
        let chk_75 = pattern.color_at(Tuple::point(0.75, 0.0, 0.0));

        let chk = chk_0.approx_eq(WHITE);
        let chk = chk && chk_25.approx_eq(Tuple::color(0.75, 0.75, 0.75));
        let chk = chk && chk_50.approx_eq(Tuple::color(0.5, 0.5, 0.5));
        let chk = chk && chk_75.approx_eq(Tuple::color(0.25, 0.25, 0.25));
        if chk {
            Ok(())
        } else {
            loge!("test_chap_10_14", "chk_0: {}", chk_0);
            loge!("test_chap_10_14", "chk_25: {}", chk_25);
            loge!("test_chap_10_14", "chk_50: {}", chk_50);
            loge!("test_chap_10_14", "chk_75: {}", chk_75);
            Err("A gradient linearly interpolates between colors".into())
        }
    }

    /// Chap 10 - A ring should extend in both x and z
    #[test]
    fn test_chap_10_15() -> Result<(), String> {
        let pattern = RingPattern::new(WHITE, BLACK);
        let chk_0 = pattern.color_at(Tuple::point(0.0, 0.0, 0.0));
        let chk_1 = pattern.color_at(Tuple::point(1.0, 0.0, 0.0));
        let chk_2 = pattern.color_at(Tuple::point(0.0, 0.0, 1.0));
        let chk_3 = pattern.color_at(Tuple::point(0.708, 0.0, 0.708));

        let chk = chk_0.approx_eq(WHITE);
        let chk = chk && chk_1.approx_eq(BLACK);
        let chk = chk && chk_2.approx_eq(BLACK);
        let chk = chk && chk_3.approx_eq(BLACK);
        if chk {
            Ok(())
        } else {
            loge!("test_chap_10_15", "chk_0: {}", chk_0);
            Err("A ring should extend in both x and z".into())
        }
    }

    /// Chap 10 - Checkers should repeat in x
    #[test]
    fn test_chap_10_16() -> Result<(), String> {
        let pattern = CheckersPattern::new(WHITE, BLACK);
        let chk_0 = pattern.color_at(Tuple::point(0.0, 0.0, 0.0));
        let chk_1 = pattern.color_at(Tuple::point(0.99, 0.0, 0.0));
        let chk_2 = pattern.color_at(Tuple::point(1.01, 0.0, 0.0));

        let chk = chk_0.approx_eq(WHITE);
        let chk = chk && chk_1.approx_eq(WHITE);
        let chk = chk && chk_2.approx_eq(BLACK);
        if chk {
            Ok(())
        } else {
            Err("Checkers should repeat in x".into())
        }
    }

    /// Chap 10 - Checkers should repeat in y
    #[test]
    fn test_chap_10_17() -> Result<(), String> {
        let pattern = CheckersPattern::new(WHITE, BLACK);
        let chk_0 = pattern.color_at(Tuple::point(0.0, 0.0, 0.0));
        let chk_1 = pattern.color_at(Tuple::point(0.0, 0.99, 0.0));
        let chk_2 = pattern.color_at(Tuple::point(0.0, 1.01, 0.0));

        let chk = chk_0.approx_eq(WHITE);
        let chk = chk && chk_1.approx_eq(WHITE);
        let chk = chk && chk_2.approx_eq(BLACK);
        if chk {
            Ok(())
        } else {
            Err("Checkers should repeat in y".into())
        }
    }

    /// Chap 10 - Checkers should repeat in z
    #[test]
    fn test_chap_10_18() -> Result<(), String> {
        let pattern = CheckersPattern::new(WHITE, BLACK);
        let chk_0 = pattern.color_at(Tuple::point(0.0, 0.0, 0.0));
        let chk_1 = pattern.color_at(Tuple::point(0.0, 0.0, 0.99));
        let chk_2 = pattern.color_at(Tuple::point(0.0, 0.0, 1.01));

        let chk = chk_0.approx_eq(WHITE);
        let chk = chk && chk_1.approx_eq(WHITE);
        let chk = chk && chk_2.approx_eq(BLACK);
        if chk {
            Ok(())
        } else {
            Err("Checkers should repeat in z".into())
        }
    }

    /// Chap 10 - Chapter 10 Putting It  Together
    #[test]
    fn test_chap_10_19() -> Result<(), String> {
        let mut world = World::new();
        let light = Light::point_light(
            Tuple::point(-10.0, 10.0, -10.0),
            Tuple::color(1.0, 1.0, 1.0),
        );
        world.light = Some(light);

        let mut pattern = GradientPattern::new(WHITE, BLACK);
        pattern.set_transform(Matrix4::scaling(0.25, 0.25, 0.25));

        let mut material = Material::new();
        material.color = Tuple::color(1.0, 0.9, 0.9);
        material.pattern = Some(Box::new(pattern));
        material.reflective = 0.5;

        let mut floor = Plane::new();
        floor.set_transform(Matrix4::scaling(10.0, 0.01, 10.0));
        material.diffuse = 0.7;
        material.specular = 0.3;
        floor.set_material(material.clone());

        let mut left_wall = Plane::new();
        left_wall.set_transform(
            Matrix4::translation(0.0, 0.0, 5.0)
                * Matrix4::rotation_y(-std::f64::consts::PI / 4.0)
                * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(10.0, 0.01, 10.0),
        );
        left_wall.set_material(floor.material().clone());

        let mut right_wall = Plane::new();
        right_wall.set_transform(
            Matrix4::translation(0.0, 0.0, 5.0)
                * Matrix4::rotation_y(std::f64::consts::PI / 4.0)
                * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(10.0, 0.01, 10.0),
        );
        right_wall.set_material(floor.material().clone());

        let mut middle = Sphere::new();
        middle.set_transform(Matrix4::translation(-0.5, 1.0, 0.5));
        material.color = Tuple::color(0.1, 1.0, 0.5);
        material.diffuse = 0.7;
        material.specular = 0.3;
        middle.set_material(material.clone());

        let mut right = Sphere::new();
        right.set_transform(Matrix4::translation(1.5, 0.5, -0.5) * Matrix4::scaling(0.5, 0.5, 0.5));
        material.color = Tuple::color(0.5, 1.0, 0.1);
        material.diffuse = 0.7;
        material.specular = 0.3;
        right.set_material(material.clone());

        let mut left = Sphere::new();
        left.set_transform(
            Matrix4::translation(-1.5, 0.33, -0.75) * Matrix4::scaling(0.33, 0.33, 0.33),
        );
        material.color = Tuple::color(1.0, 0.8, 0.1);
        material.diffuse = 0.7;
        material.specular = 0.3;
        left.set_material(material);

        // VIEW TRANSFORM SETTINGS
        let from = Tuple::point(0.0, 1.5, -12.0);
        let to = Tuple::point(0.0, 1.0, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let transform = view_transform(from, to, up);

        let camera = Camera::new(100, 50, std::f64::consts::PI / 3.0).with_transform(transform);

        world.add_shape(Box::new(floor));
        world.add_shape(Box::new(left_wall));
        world.add_shape(Box::new(right_wall));
        world.add_shape(Box::new(middle));
        world.add_shape(Box::new(left));
        world.add_shape(Box::new(right));

        let image = world.render(camera);
        let rc = image.write_ppm("test_chap_10_19_putting_it_together.ppm");
        let chk = rc.is_ok();

        if chk {
            Ok(())
        } else {
            Err("Chapter 10 Putting It  Together".into())
        }
    }

    /// Chap 10 - Nested patterns showcase (4K, full rabbit hole)
    ///
    /// Renders a scene that exercises every pattern type AND nesting:
    /// - floor: checkers nested with a gradient (alternating tiles)
    /// - walls: stripes nested with rings
    /// - spheres: gradients, rings, stripes, and a doubly-nested pattern
    ///
    /// Rendered at 4K (3840x2160) via the parallel renderer.
    #[test]
    fn test_chap_10_20_nested_showcase() -> Result<(), String> {
        // --- colors -------------------------------------------------------
        let red = Tuple::color(0.9, 0.1, 0.1);
        let green = Tuple::color(0.1, 0.9, 0.2);
        let blue = Tuple::color(0.1, 0.2, 0.9);
        let cyan = Tuple::color(0.1, 0.9, 0.9);
        let magenta = Tuple::color(0.9, 0.1, 0.9);
        let yellow = Tuple::color(0.95, 0.85, 0.1);
        let orange = Tuple::color(1.0, 0.55, 0.0);

        let mut world = World::new();
        world.light = Some(Light::point_light(
            Tuple::point(-10.0, 10.0, -10.0),
            Tuple::color(1.0, 1.0, 1.0),
        ));

        // --- floor: checkers tiles, each tile filled by a gradient --------
        let mut floor_checkers = CheckersPattern::new(WHITE, BLACK);
        floor_checkers.set_transform(Matrix4::scaling(0.5, 0.5, 0.5));
        let mut floor_gradient = GradientPattern::new(blue, cyan);
        floor_gradient.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));
        let mut floor_pattern =
            NestedPattern::new(Box::new(floor_checkers), Box::new(floor_gradient));
        floor_pattern.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));

        let mut floor_mat = Material::new();
        floor_mat.pattern = Some(Box::new(floor_pattern));
        floor_mat.diffuse = 0.7;
        floor_mat.specular = 0.1;

        let mut floor = Plane::new();
        floor.set_material(floor_mat);

        // --- back wall: stripes alternating with rings --------------------
        let mut wall_stripes = StripePattern::new(magenta, WHITE);
        wall_stripes.set_transform(Matrix4::scaling(0.25, 0.25, 0.25));
        let mut wall_rings = RingPattern::new(yellow, orange);
        wall_rings.set_transform(Matrix4::scaling(0.5, 0.5, 0.5));
        let mut wall_pattern = NestedPattern::new(Box::new(wall_stripes), Box::new(wall_rings));
        wall_pattern.set_transform(Matrix4::rotation_y(std::f64::consts::PI / 6.0));

        let mut wall_mat = Material::new();
        wall_mat.pattern = Some(Box::new(wall_pattern));
        wall_mat.specular = 0.0;

        let mut back_wall = Plane::new();
        back_wall.set_transform(
            Matrix4::translation(0.0, 0.0, 8.0) * Matrix4::rotation_x(std::f64::consts::PI / 2.0),
        );
        back_wall.set_material(wall_mat);

        // --- middle sphere: doubly-nested (stripe-of-gradients vs ring) ---
        let mut inner_stripe = StripePattern::new(red, green);
        inner_stripe.set_transform(Matrix4::scaling(0.25, 0.25, 0.25));
        let inner_gradient = GradientPattern::new(yellow, magenta);
        // first nest: stripes alternating with a gradient
        let mut nest_a = NestedPattern::new(Box::new(inner_stripe), Box::new(inner_gradient));
        nest_a.set_transform(Matrix4::scaling(0.5, 0.5, 0.5));
        let mut ring_child = RingPattern::new(cyan, blue);
        ring_child.set_transform(Matrix4::scaling(0.3, 0.3, 0.3));
        // second nest: the whole thing above alternating with a ring
        let mut middle_pattern = NestedPattern::new(Box::new(nest_a), Box::new(ring_child));
        middle_pattern.set_transform(
            Matrix4::scaling(0.6, 0.6, 0.6) * Matrix4::rotation_z(std::f64::consts::PI / 4.0),
        );

        let mut middle_mat = Material::new();
        middle_mat.pattern = Some(Box::new(middle_pattern));
        middle_mat.diffuse = 0.7;
        middle_mat.specular = 0.3;

        let mut middle = Sphere::new();
        middle.set_transform(Matrix4::translation(-0.5, 1.0, 0.5));
        middle.set_material(middle_mat);

        // --- right sphere: rotated gradient -------------------------------
        let mut right_grad = GradientPattern::new(green, magenta);
        right_grad.set_transform(
            Matrix4::scaling(0.5, 0.5, 0.5) * Matrix4::rotation_y(std::f64::consts::PI / 4.0),
        );
        let mut right_mat = Material::new();
        right_mat.pattern = Some(Box::new(right_grad));
        right_mat.diffuse = 0.7;
        right_mat.specular = 0.3;
        let mut right = Sphere::new();
        right.set_transform(Matrix4::translation(1.5, 0.5, -0.5) * Matrix4::scaling(0.5, 0.5, 0.5));
        right.set_material(right_mat);

        // --- left sphere: fine rings --------------------------------------
        let mut left_rings = RingPattern::new(red, yellow);
        left_rings.set_transform(Matrix4::scaling(0.15, 0.15, 0.15));
        let mut left_mat = Material::new();
        left_mat.pattern = Some(Box::new(left_rings));
        left_mat.diffuse = 0.7;
        left_mat.specular = 0.3;
        let mut left = Sphere::new();
        left.set_transform(
            Matrix4::translation(-1.5, 0.33, -0.75) * Matrix4::scaling(0.33, 0.33, 0.33),
        );
        left.set_material(left_mat);

        // --- extra ball #1: nested stripe/ring, blue family ---------------
        let mut b1_stripe = StripePattern::new(blue, cyan);
        b1_stripe.set_transform(Matrix4::scaling(0.2, 0.2, 0.2));
        let mut b1_ring = RingPattern::new(WHITE, blue);
        b1_ring.set_transform(Matrix4::scaling(0.2, 0.2, 0.2));
        let mut b1_pattern = BlendedPattern::new(Box::new(b1_stripe), Box::new(b1_ring));
        b1_pattern.set_transform(Matrix4::rotation_z(std::f64::consts::PI / 3.0));
        let mut b1_mat = Material::new();
        b1_mat.pattern = Some(Box::new(b1_pattern));
        b1_mat.diffuse = 0.7;
        b1_mat.specular = 0.4;
        let mut ball1 = Sphere::new();
        ball1.set_transform(
            Matrix4::translation(2.6, 0.75, 1.2) * Matrix4::scaling(0.75, 0.75, 0.75),
        );
        ball1.set_material(b1_mat);

        // --- extra ball #2: warm gradient ---------------------------------
        let mut b2_grad = GradientPattern::new(orange, red);
        b2_grad.set_transform(Matrix4::scaling(0.5, 0.5, 0.5));
        let mut b2_mat = Material::new();
        b2_mat.pattern = Some(Box::new(b2_grad));
        b2_mat.diffuse = 0.8;
        b2_mat.specular = 0.5;
        b2_mat.shininess = 300.0;
        let mut ball2 = Sphere::new();
        ball2.set_transform(Matrix4::translation(-2.7, 0.5, 0.3) * Matrix4::scaling(0.5, 0.5, 0.5));
        ball2.set_material(b2_mat);

        // --- camera (4K) --------------------------------------------------
        let from = Tuple::point(0.0, 1.5, -12.0);
        let to = Tuple::point(0.0, 1.0, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let transform = view_transform(from, to, up);
        let camera = Camera::new(96, 54, std::f64::consts::PI / 3.0).with_transform(transform);

        world.add_shape(Box::new(floor));
        world.add_shape(Box::new(back_wall));
        world.add_shape(Box::new(middle));
        world.add_shape(Box::new(left));
        world.add_shape(Box::new(right));
        world.add_shape(Box::new(ball1));
        world.add_shape(Box::new(ball2));

        let image = world.render_parallel(camera);
        let rc = image.write_ppm("test_chap_10_20_nested_showcase.ppm");

        if rc.is_ok() {
            Ok(())
        } else {
            Err("Chapter 10 Nested patterns showcase".into())
        }
    }

    /// Chap 11 - Precomputing the reflection vector
    #[test]
    fn test_chap_11_2() -> Result<(), String> {
        let sqrt2_over_2 = (1.0 / std::f64::consts::FRAC_1_SQRT_2) / 2.0;

        let shape = Plane::new();
        let r = Ray::new(
            Tuple::point(0.0, 1.0, 1.0),
            Tuple::vector(0.0, -sqrt2_over_2, sqrt2_over_2),
        );
        let i = Intersection::new((2.0_f64).sqrt(), shape.id());
        let comps = prepare_computations_upto_chap10(i, &r, &shape);
        let chk = comps
            .reflectv
            .approx_eq(Tuple::vector(0.0, sqrt2_over_2, sqrt2_over_2));

        if chk {
            Ok(())
        } else {
            loge!("test_chap_11_2", "comps.reflectv:{}", comps.reflectv);
            Err("Precomputing the reflection vector ".into())
        }
    }

    /// Chap 11 - The reflected color for a nonreflective material
    #[test]
    fn test_chap_11_3() -> Result<(), String> {
        let mut w = World::default_world();

        let point = Tuple::point(0.0, 0.0, 0.0);
        let direction = Tuple::vector(0.0, 0.0, 1.0);
        let r = Ray::new(point, direction);
        w.shapes[1].data_mut().material.ambient = 1.0;
        let i = Intersection::new(1.0, w.shapes[1].as_ref().id());
        let comps = prepare_computations_upto_chap10(i, &r, w.shapes[1].as_ref());
        let color = w.reflected_color(&comps, 10);

        let chk = color.approx_eq(BLACK);
        if chk {
            Ok(())
        } else {
            Err("The reflected color for a nonreflective material".into())
        }
    }

    /// Chap 11 - The reflected color for a reflective material
    #[test]
    fn test_chap_11_4() -> Result<(), String> {
        let mut w = World::default_world();

        let mut material = Material::new();
        material.reflective = 0.5;

        let mut shape = Plane::new();
        shape.set_transform(Matrix4::translation(0.0, -1.0, 0.0));
        shape.set_material(material.clone());

        w.add_shape(Box::new(shape));

        let sqrt2_over_2 = (1.0 / std::f64::consts::FRAC_1_SQRT_2) / 2.0;
        let point = Tuple::point(0.0, 0.0, -3.0);
        let direction = Tuple::vector(0.0, -sqrt2_over_2, sqrt2_over_2);
        let r = Ray::new(point, direction);

        let shape = w.shapes.last().expect("world must have at least one shape");
        let i = Intersection::new(2.0_f64.sqrt(), shape.id());
        let comps = prepare_computations_upto_chap10(i, &r, shape.as_ref());
        let color = w.reflected_color(&comps, 10);

        let chk = color.approx_eq(Tuple::color(0.190330596701, 0.237913245876, 0.142747947526));
        if chk {
            Ok(())
        } else {
            loge!("test_chap_11_4", "color:{}", color);
            Err("The reflected color for a reflective material".into())
        }
    }

    /// Chap 11 - shade_hit() with a reflective material
    #[test]
    fn test_chap_11_5() -> Result<(), String> {
        let mut w = World::default_world();

        let mut material = Material::new();
        material.reflective = 0.5;

        let mut shape = Plane::new();
        shape.set_transform(Matrix4::translation(0.0, -1.0, 0.0));
        shape.set_material(material.clone());

        w.add_shape(Box::new(shape));

        let sqrt2_over_2 = (1.0 / std::f64::consts::FRAC_1_SQRT_2) / 2.0;
        let point = Tuple::point(0.0, 0.0, -3.0);
        let direction = Tuple::vector(0.0, -sqrt2_over_2, sqrt2_over_2);
        let r = Ray::new(point, direction);

        let shape = w.shapes.last().expect("world must have at least one shape");
        let i = Intersection::new(2.0_f64.sqrt(), shape.id());
        let comps = prepare_computations_upto_chap10(i, &r, shape.as_ref());
        let color = w.shade_hit(&comps, 10);

        let chk = color.approx_eq(Tuple::color(0.876755985652, 0.924338634827, 0.829173336477));

        if chk {
            Ok(())
        } else {
            loge!("test_chap_11_5", "color:{}", color);
            Err("shade_hit() with a reflective material".into())
        }
    }

    /// Chap 11 - color_at() with mutually reflective surfaces Test #6: Avoid Infinite Recursion
    /// Show that your code safely handles infinite recursion caused by two objects that mutually
    /// reflect rays between themselves. Create two parallel mirrors by positioning one plane above
    /// another and making them both reflective. Orient a ray so that it strikes one plane and
    /// bounces to the other. What will happen?
    #[test]
    fn test_chap_11_6() -> Result<(), String> {
        let mut w = World::new(); //World::default_world();

        let light = Light::point_light(Tuple::point(0.0, 0.0, 0.0), Tuple::color(1.0, 1.0, 1.0));
        w.set_light(light);

        let mut material = Material::new();
        material.reflective = 1.0;

        let mut lower = Plane::new();
        lower.set_transform(Matrix4::translation(0.0, -1.0, 0.0));
        lower.set_material(material.clone());

        w.add_shape(Box::new(lower));

        let mut upper = Plane::new();
        upper.set_transform(Matrix4::translation(0.0, 1.0, 0.0));
        upper.set_material(material.clone());

        w.add_shape(Box::new(upper));

        let point = Tuple::point(0.0, 0.0, 0.0);
        let direction = Tuple::vector(0.0, 1.0, 0.0);
        let r = Ray::new(point, direction);

        let color = w.color_at(&r, 10);
        let chk = color.x >= 1.0;

        if chk {
            Ok(())
        } else {
            loge!("test_chap_11_6", "color:{}", color);
            Err("shade_hit() with a reflective material".into())
        }
    }

    /// Chap 11 - Transparency and Refractive Index for the default material
    #[test]
    fn test_chap_11_7() -> Result<(), String> {
        let material = Material::new();
        let chk =
            approx_eq(material.transparency, 0.0) && approx_eq(material.refractive_index, 1.0);
        if chk {
            Ok(())
        } else {
            loge!("test_chap_11_7", "material:{}", material);
            Err("Transparency and Refractive Index for the default material".into())
        }
    }

    /// Chap 11 - A helper for producing a sphere with a glassy material
    #[test]
    fn test_chap_11_8() -> Result<(), String> {
        let s = Sphere::glass();
        let chk = s.transform().approx_eq(Matrix4::identity())
            && approx_eq(s.data.material.transparency, 1.0)
            && approx_eq(s.data.material.refractive_index, 1.5);
        if chk {
            Ok(())
        } else {
            Err("A helper for producing a sphere with a glassy material".into())
        }
    }

    /// Chap 11 - Finding n1 and n2 at various intersections
    #[test]
    fn test_chap_11_9() -> Result<(), String> {
        let mut a = Sphere::glass();
        a.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));
        a.data.material.refractive_index = 1.5;

        let mut b = Sphere::glass();
        b.set_transform(Matrix4::translation(0.0, 0.0, -0.25));
        b.data.material.refractive_index = 2.0;

        let mut c = Sphere::glass();
        c.set_transform(Matrix4::translation(0.0, 0.0, 0.25));
        c.data.material.refractive_index = 2.5;

        // capture ids BEFORE the moves below
        let (a_id, b_id, c_id) = (a.id(), b.id(), c.id());

        // now move each Sphere into a boxed trait object
        let shapes: [Box<dyn Shape>; 3] = [Box::new(a), Box::new(b), Box::new(c)];

        let r = Ray::new(Tuple::point(0.0, 0.0, -4.0), Tuple::vector(0.0, 0.0, 1.0));
        let mut xs = Intersections::new();
        xs.push(Intersection::new(2.0, a_id));
        xs.push(Intersection::new(2.75, b_id));
        xs.push(Intersection::new(3.25, c_id));
        xs.push(Intersection::new(4.75, b_id));
        xs.push(Intersection::new(5.25, c_id));
        xs.push(Intersection::new(6.0, a_id));

        let expected = [
            (1.0, 1.5),
            (1.5, 2.0),
            (2.0, 2.5),
            (2.5, 2.5),
            (2.5, 1.5),
            (1.5, 1.0),
        ];

        for (idx, (n1, n2)) in expected.iter().enumerate() {
            let comps = prepare_computations(xs[idx], &r, &shapes, &xs);
            if !approx_eq(comps.n1, *n1) || !approx_eq(comps.n2, *n2) {
                return Err(format!(
                    "idx {idx}: got ({}, {}), expected ({n1}, {n2})",
                    comps.n1, comps.n2
                ));
            }
        }
        Ok(())
    }

    /// Chap x - The under point is offset below the surface
    #[test]
    fn test_chap_11_10() -> Result<(), String> {
        let ray = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let mut shape = Sphere::glass();
        shape.set_transform(Matrix4::translation(0.0, 0.0, 1.0));

        let i = Intersection::new(5.0, shape.id());
        let mut xs = Intersections::new();
        xs.push(i);

        // now move each Sphere into a boxed trait object
        let shapes: [Box<dyn Shape>; 1] = [Box::new(shape)];

        let comps = prepare_computations(i, &ray, &shapes, &xs);

        let chk = comps.under_point.z > EPSILON / 2_f64 && comps.point.z < comps.under_point.z;

        if chk {
            Ok(())
        } else {
            Err(" The under point is offset below the surface".into())
        }
    }

    /// Chap 11 - The refracted color with an opaque surface
    #[test]
    fn test_chap_11_11() -> Result<(), String> {
        let w = World::default_world();
        let shape = &w.shapes[0];
        let ray = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));

        let mut xs = Intersections::new();
        xs.push(Intersection::new(4.0, shape.id()));
        xs.push(Intersection::new(6.0, shape.id()));

        let comps = prepare_computations(xs[0], &ray, &w.shapes, &xs);
        let c = w.refracted_color(&comps, 0);

        let chk = c.approx_eq(Tuple::color(0.0, 0.0, 0.0));
        if chk {
            Ok(())
        } else {
            Err("The refracted color with an opaque surface".into())
        }
    }

    /// Chap 11 - The refracted color at the maximum recursive depth
    #[test]
    fn test_chap_11_12() -> Result<(), String> {
        let mut w = World::default_world();
        let shape = w.shapes[0].as_mut();

        let mut m = Material::new();
        m.transparency = 1.0;
        m.refractive_index = 1.5;
        shape.set_material(m.clone());

        let ray = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));

        let mut xs = Intersections::new();
        xs.push(Intersection::new(4.0, shape.id()));
        xs.push(Intersection::new(6.0, shape.id()));

        let comps = prepare_computations(xs[0], &ray, &w.shapes, &xs);
        let c = w.refracted_color(&comps, 0);

        let chk = c.approx_eq(Tuple::color(0.0, 0.0, 0.0));
        if chk {
            Ok(())
        } else {
            Err("The refracted color at the maximum recursive depth".into())
        }
    }

    /// Chap x - The refracted color under total internal reflection
    #[test]
    fn test_chap_11_13() -> Result<(), String> {
        let mut w = World::default_world();
        let shape = w.shapes[0].as_mut();

        let mut m = Material::new();
        m.transparency = 1.0;
        m.refractive_index = 1.5;
        shape.set_material(m.clone());

        let sqrt2_over_2 = (1.0 / std::f64::consts::FRAC_1_SQRT_2) / 2.0;
        let ray = Ray::new(
            Tuple::point(0.0, 0.0, sqrt2_over_2),
            Tuple::vector(0.0, 1.0, 0.0),
        );

        let mut xs = Intersections::new();
        xs.push(Intersection::new(-sqrt2_over_2, shape.id()));
        xs.push(Intersection::new(sqrt2_over_2, shape.id()));

        // NOTE: this time you're inside the sphere, so you need
        // to look at the second intersection, xs[1], not xs[0]
        let comps = prepare_computations(xs[1], &ray, &w.shapes, &xs);
        let c = w.refracted_color(&comps, 5);

        let chk = c.approx_eq(Tuple::color(0.0, 0.0, 0.0));
        if chk {
            Ok(())
        } else {
            loge!("test_chap_11_13", "c:{}", c);
            Err("The refracted color under total internal reflection".into())
        }
    }

    /// Chap x - The refracted color with a refracted ray
    #[test]
    fn test_chap_11_14() -> Result<(), String> {
        let mut w = World::default_world();

        // Avoid mutating twice on the vector by
        // split once; `left` holds shapes[0], `right` holds shapes[1..]
        let (left, right) = w.shapes.split_at_mut(1);

        let a = left[0].as_mut();
        let b = right[0].as_mut();

        let mut m = Material::new();
        m.ambient = 1.0;
        m.pattern = Some(Box::new(TestPattern::new()));
        a.set_material(m.clone());

        let mut m = Material::new();
        m.transparency = 1.0;
        m.refractive_index = 1.5;
        b.set_material(m.clone());

        let ray = Ray::new(Tuple::point(0.0, 0.0, 0.1), Tuple::vector(0.0, 1.0, 0.0));

        let mut xs = Intersections::new();
        xs.push(Intersection::new(-0.9899, a.id()));
        xs.push(Intersection::new(-0.4899, b.id()));
        xs.push(Intersection::new(0.4899, b.id()));
        xs.push(Intersection::new(0.9899, a.id()));

        let comps = prepare_computations(xs[2], &ray, &w.shapes, &xs);
        let c = w.refracted_color(&comps, 5);

        let chk = c.approx_eq(Tuple::color(0.000000000000, 0.998884681786, 0.047216421860));
        if chk {
            Ok(())
        } else {
            loge!("test_chap_11_14", "refracted_color c: {}", c);
            Err("The refracted color with a refracted ray".into())
        }
    }

    /// Chap 11 - shade_hit() with a transparent material
    #[test]
    fn test_chap_11_15() -> Result<(), String> {
        let mut w = World::default_world();
        let mut floor = Plane::new();
        floor.set_transform(Matrix4::translation(0.0, -1.0, 0.0));

        let mut m = Material::new();
        m.transparency = 0.5;
        m.refractive_index = 1.5;
        floor.set_material(m.clone());

        let mut ball = Sphere::new();
        let mut m = Material::new();
        m.color = Tuple::color(1.0, 0.0, 0.0);
        m.ambient = 0.5;
        ball.set_transform(Matrix4::translation(0.0, -3.5, -0.5));
        ball.set_material(m.clone());

        let floor_id = w.add_shape(Box::new(floor));
        w.add_shape(Box::new(ball));

        let sqrt2_over_2 = (1.0 / std::f64::consts::FRAC_1_SQRT_2) / 2.0;
        let ray = Ray::new(
            Tuple::point(0.0, 0.0, -3.0),
            Tuple::vector(0.0, -sqrt2_over_2, sqrt2_over_2),
        );

        let mut xs = Intersections::new();
        xs.push(Intersection::new(2_f64.sqrt(), floor_id));

        let comps = prepare_computations(xs[0], &ray, &w.shapes, &xs);
        let color = w.shade_hit(&comps, 5);

        let chk = color.approx_eq(Tuple::color(0.936425388951, 0.686425388951, 0.686425388951));
        if chk {
            Ok(())
        } else {
            loge!("test_chap_11_15", "shade hit produced color: {}", color);
            Err("shade_hit() with a transparent material".into())
        }
    }

    /// Chap 11 - The Schlick approximation under total internal reflection
    #[test]
    fn test_chap_11_16() -> Result<(), String> {
        let mut w = World::default_world();
        let shape = Sphere::glass();
        let shape_id = w.add_shape(Box::new(shape));

        let sqrt2_over_2 = (1.0 / std::f64::consts::FRAC_1_SQRT_2) / 2.0;
        let ray = Ray::new(
            Tuple::point(0.0, 0.0, sqrt2_over_2),
            Tuple::vector(0.0, 1.0, 0.0),
        );

        let mut xs = Intersections::new();
        xs.push(Intersection::new(-sqrt2_over_2, shape_id));
        xs.push(Intersection::new(sqrt2_over_2, shape_id));

        let comps = prepare_computations(xs[1], &ray, &w.shapes, &xs);
        let reflectance = schlick(&comps);

        let chk = approx_eq(reflectance, 1.0);
        if chk {
            Ok(())
        } else {
            Err("The Schlick approximation under total internal reflection".into())
        }
    }

    /// Chap 11 - The Schlick approximation with a perpendicular viewing angle
    #[test]
    fn test_chap_11_17() -> Result<(), String> {
        let mut w = World::default_world();
        let shape = Sphere::glass();
        let shape_id = w.add_shape(Box::new(shape));

        let ray = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 1.0, 0.0));

        let mut xs = Intersections::new();
        xs.push(Intersection::new(-1.0, shape_id));
        xs.push(Intersection::new(1.0, shape_id));

        let comps = prepare_computations(xs[1], &ray, &w.shapes, &xs);
        let reflectance = schlick(&comps);

        let chk = approx_eq(reflectance, 0.04);
        if chk {
            Ok(())
        } else {
            loge!("test_chap_11_17", "reflectance: {}", reflectance);
            Err("The Schlick approximation with a perpendicular viewing angle".into())
        }
    }

    /// Chap 11 - The Schlick approximation with small angle and n2 > n1
    #[test]
    fn test_chap_11_18() -> Result<(), String> {
        let mut w = World::default_world();
        let shape = Sphere::glass();
        let shape_id = w.add_shape(Box::new(shape));

        let ray = Ray::new(Tuple::point(0.0, 0.99, -2.0), Tuple::vector(0.0, 0.0, 1.0));

        let mut xs = Intersections::new();
        xs.push(Intersection::new(1.8589, shape_id));

        let comps = prepare_computations(xs[0], &ray, &w.shapes, &xs);
        let reflectance = schlick(&comps);

        let chk = approx_eq(reflectance, 0.48873081012212183);
        if chk {
            Ok(())
        } else {
            loge!("test_chap_11_18", "reflectance: {}", reflectance);
            Err("The Schlick approximation with small angle and n2 > n1".into())
        }
    }

    /// Chap 11 - shade_hit() with a reflective, transparent material
    #[test]
    fn test_chap_11_19() -> Result<(), String> {
        let mut w = World::default_world();
        let mut floor = Plane::new();
        floor.set_transform(Matrix4::translation(0.0, -1.0, 0.0));

        let mut m = Material::new();
        m.reflective = 0.5;
        m.transparency = 0.5;
        m.refractive_index = 1.5;
        floor.set_material(m.clone());

        let mut ball = Sphere::new();
        let mut m = Material::new();
        m.color = Tuple::color(1.0, 0.0, 0.0);
        m.ambient = 0.5;
        ball.set_transform(Matrix4::translation(0.0, -3.5, -0.5));
        ball.set_material(m.clone());

        let floor_id = w.add_shape(Box::new(floor));
        w.add_shape(Box::new(ball));

        let sqrt2_over_2 = (1.0 / std::f64::consts::FRAC_1_SQRT_2) / 2.0;
        let ray = Ray::new(
            Tuple::point(0.0, 0.0, -3.0),
            Tuple::vector(0.0, -sqrt2_over_2, sqrt2_over_2),
        );

        let mut xs = Intersections::new();
        xs.push(Intersection::new(2_f64.sqrt(), floor_id));

        let comps = prepare_computations(xs[0], &ray, &w.shapes, &xs);
        let color = w.shade_hit(&comps, 5);

        let chk = color.approx_eq(Tuple::color(0.933915140526, 0.696434226271, 0.692430691343));

        if chk {
            Ok(())
        } else {
            loge!("test_chap_11_19", "color: {}", color);
            Err("shade_hit() with a reflective, transparent material".into())
        }
    }

    /// Chap 11 - Chapter 11 Putting It  Together
    #[test]
    fn test_chap_11_20() -> Result<(), String> {
        let mut world = World::new();
        let light = Light::point_light(
            Tuple::point(-10.0, 10.0, -10.0),
            Tuple::color(1.0, 1.0, 1.0),
        );
        world.light = Some(light);

        let mut pattern = GradientPattern::new(WHITE, BLACK);
        pattern.set_transform(Matrix4::scaling(0.25, 0.25, 0.25));

        let mut material = Material::new();
        material.color = Tuple::color(1.0, 0.9, 0.9);
        // material.pattern = Some(Box::new(pattern));
        material.reflective = 0.25;

        let mut floor = Plane::new();
        floor.set_transform(Matrix4::scaling(10.0, 0.01, 10.0));
        material.diffuse = 0.7;
        material.specular = 0.3;
        floor.set_material(material.clone());

        let mut left_wall = Plane::new();
        left_wall.set_transform(
            Matrix4::translation(0.0, 0.0, 5.0)
                * Matrix4::rotation_y(-std::f64::consts::PI / 4.0)
                * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(10.0, 0.01, 10.0),
        );
        left_wall.set_material(floor.material().clone());

        let mut right_wall = Plane::new();
        right_wall.set_transform(
            Matrix4::translation(0.0, 0.0, 5.0)
                * Matrix4::rotation_y(std::f64::consts::PI / 4.0)
                * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(10.0, 0.01, 10.0),
        );
        right_wall.set_material(floor.material().clone());

        let mut material = Material::new();
        material.color = Tuple::color(1.0, 0.9, 0.9);
        let mut pattern = CheckersPattern::new(WHITE, BLACK);
        pattern.set_transform(Matrix4::scaling(0.1, 0.1, 0.1));
        material.pattern = Some(Box::new(pattern));

        let mut middle = Sphere::new();
        middle.set_transform(Matrix4::translation(-0.5, 1.0, 0.5));
        material.color = Tuple::color(0.1, 1.0, 0.5);
        material.diffuse = 0.7;
        material.specular = 0.3;
        middle.set_material(material.clone());

        let mut right = Sphere::glass();
        right.set_transform(Matrix4::translation(1.5, 0.5, -0.5) * Matrix4::scaling(0.5, 0.5, 0.5));
        let mut material = Material::new();
        material.pattern = Some(Box::new(pattern));
        material.color = Tuple::color(0.82, 0.0, 0.0);
        material.transparency = 0.95;
        material.reflective = 1.0;
        material.shininess = 300.0;
        material.specular = 1.0;
        right.set_material(material.clone());

        let mut left = Sphere::new();
        left.set_transform(
            Matrix4::translation(-1.5, 0.33, -0.75) * Matrix4::scaling(0.33, 0.33, 0.33),
        );
        material.color = Tuple::color(1.0, 0.8, 0.1);
        material.diffuse = 0.7;
        material.specular = 0.3;
        left.set_material(material);

        // VIEW TRANSFORM SETTINGS
        let from = Tuple::point(0.0, 1.5, -12.0);
        let to = Tuple::point(0.0, 1.0, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let transform = view_transform(from, to, up);

        let camera = Camera::new(60, 40, std::f64::consts::PI / 3.0).with_transform(transform);

        world.add_shape(Box::new(floor));
        world.add_shape(Box::new(left_wall));
        world.add_shape(Box::new(right_wall));
        world.add_shape(Box::new(middle));
        world.add_shape(Box::new(left));
        world.add_shape(Box::new(right));

        let image = world.render(camera);
        let rc = image.write_ppm("test_chap_11_20_putting_it_together.ppm");
        let chk = rc.is_ok();

        if chk {
            Ok(())
        } else {
            Err("Chapter 11_20 Putting It  Together".into())
        }
    }

    /// Chap 12 - A ray intersects a cube
    #[test]
    fn test_chap_12_1() -> Result<(), String> {
        let c = Cube::new();

        // +x
        let r = Ray::new(Tuple::point(5.0, 0.5, 0.0), Tuple::vector(-1.0, 0.0, 0.0));
        let xs = c.local_intersect(&r);
        let chk = xs.count() == 2;
        let chk = chk && approx_eq(xs[0].t, 4.0) && approx_eq(xs[1].t, 6.0);

        // -x
        let r = Ray::new(Tuple::point(-5.0, 0.5, 0.0), Tuple::vector(1.0, 0.0, 0.0));
        let xs = c.local_intersect(&r);
        let chk = chk && approx_eq(xs[0].t, 4.0) && approx_eq(xs[1].t, 6.0);

        // +y
        let r = Ray::new(Tuple::point(0.5, 5.0, 0.0), Tuple::vector(0.0, -1.0, 0.0));
        let xs = c.local_intersect(&r);
        let chk = chk && approx_eq(xs[0].t, 4.0) && approx_eq(xs[1].t, 6.0);

        // -y
        let r = Ray::new(Tuple::point(0.5, -5.0, 0.0), Tuple::vector(0.0, 1.0, 0.0));
        let xs = c.local_intersect(&r);
        let chk = chk && approx_eq(xs[0].t, 4.0) && approx_eq(xs[1].t, 6.0);

        // +z
        let r = Ray::new(Tuple::point(0.5, 0.0, 5.0), Tuple::vector(0.0, 0.0, -1.0));
        let xs = c.local_intersect(&r);
        let chk = chk && approx_eq(xs[0].t, 4.0) && approx_eq(xs[1].t, 6.0);

        // -z
        let r = Ray::new(Tuple::point(0.5, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs = c.local_intersect(&r);
        let chk = chk && approx_eq(xs[0].t, 4.0) && approx_eq(xs[1].t, 6.0);

        // inside
        let r = Ray::new(Tuple::point(0.0, 0.5, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs = c.local_intersect(&r);
        let chk = chk && approx_eq(xs[0].t, -1.0) && approx_eq(xs[1].t, 1.0);

        if chk {
            Ok(())
        } else {
            Err("A ray intersects a cube".into())
        }
    }

    /// Chap 12 - A ray misses a cube
    #[test]
    fn test_chap_12_2() -> Result<(), String> {
        let c = Cube::new();

        let r = Ray::new(
            Tuple::point(-2.0, 0.0, 0.0),
            Tuple::vector(0.26730, 0.5345, 0.8018),
        );
        let xs = c.local_intersect(&r);
        let chk = xs.count() == 0;

        let r = Ray::new(
            Tuple::point(0.0, -2.0, 0.0),
            Tuple::vector(0.8018, 0.2673, 0.5345),
        );
        let xs = c.local_intersect(&r);
        let chk = chk && xs.count() == 0;

        let r = Ray::new(
            Tuple::point(0.0, 0.0, -2.0),
            Tuple::vector(0.5345, 0.8018, 0.2673),
        );
        let xs = c.local_intersect(&r);
        let chk = chk && xs.count() == 0;

        let r = Ray::new(Tuple::point(2.0, 0.0, 2.0), Tuple::vector(0.0, 0.0, -1.0));
        let xs = c.local_intersect(&r);
        let chk = chk && xs.count() == 0;

        let r = Ray::new(Tuple::point(0.0, 2.0, 2.0), Tuple::vector(0.0, -1.0, 0.0));
        let xs = c.local_intersect(&r);
        let chk = chk && xs.count() == 0;

        let r = Ray::new(Tuple::point(2.0, 2.0, 0.0), Tuple::vector(-1.0, 0.0, 0.0));
        let xs = c.local_intersect(&r);
        let chk = chk && xs.count() == 0;

        if chk {
            Ok(())
        } else {
            Err("A ray misses a cube".into())
        }
    }

    /// Chap 12 - The normal on the surface of a cube
    #[test]
    fn test_chap_12_3() -> Result<(), String> {
        let c = Cube::new();

        let p = Tuple::point(1.0, 0.5, -0.8);
        let normal = c.local_normal_at(p);
        let chk = Tuple::vector(1.0, 0.0, 0.0).approx_eq(normal);

        let p = Tuple::point(-1.0, -0.5, 0.9);
        let normal = c.local_normal_at(p);
        let chk = chk && Tuple::vector(-1.0, 0.0, 0.0).approx_eq(normal);

        let p = Tuple::point(-0.4, 1.0, -0.1);
        let normal = c.local_normal_at(p);
        let chk = chk && Tuple::vector(0.0, 1.0, 0.0).approx_eq(normal);

        let p = Tuple::point(0.3, -1.0, -0.7);
        let normal = c.local_normal_at(p);
        let chk = chk && Tuple::vector(0.0, -1.0, 0.0).approx_eq(normal);

        let p = Tuple::point(-0.6, 0.3, 1.0);
        let normal = c.local_normal_at(p);
        let chk = chk && Tuple::vector(0.0, 0.0, 1.0).approx_eq(normal);

        let p = Tuple::point(0.4, 0.4, -1.0);
        let normal = c.local_normal_at(p);
        let chk = chk && Tuple::vector(0.0, 0.0, -1.0).approx_eq(normal);

        let p = Tuple::point(1.0, 1.0, 1.0);
        let normal = c.local_normal_at(p);
        let chk = chk && Tuple::vector(1.0, 0.0, 0.0).approx_eq(normal);

        let p = Tuple::point(-1.0, -1.0, -1.0);
        let normal = c.local_normal_at(p);
        let chk = chk && Tuple::vector(-1.0, 0.0, 0.0).approx_eq(normal);

        if chk {
            Ok(())
        } else {
            Err("The normal on the surface of a cube".into())
        }
    }

    /// Chap 12 - Chapter 12 Putting It  Together
    #[test]
    fn test_chap_12_4() -> Result<(), String> {
        let mut world = World::new();
        let light = Light::point_light(
            Tuple::point(-10.0, 10.0, -10.0),
            Tuple::color(1.0, 1.0, 1.0),
        );
        world.light = Some(light);

        let mut pattern = GradientPattern::new(WHITE, BLACK);
        pattern.set_transform(Matrix4::scaling(0.25, 0.25, 0.25));

        let mut material = Material::new();
        material.color = Tuple::color(1.0, 0.9, 0.9);
        // material.pattern = Some(Box::new(pattern));
        material.reflective = 0.05;

        let mut floor = Plane::new();
        floor.set_transform(Matrix4::scaling(10.0, 0.01, 10.0));
        material.diffuse = 0.7;
        material.specular = 0.3;
        floor.set_material(material.clone());

        let mut left_wall = Plane::new();
        left_wall.set_transform(
            Matrix4::translation(0.0, 0.0, 5.0)
                * Matrix4::rotation_y(-std::f64::consts::PI / 4.0)
                * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(10.0, 0.01, 10.0),
        );
        left_wall.set_material(floor.material().clone());

        let mut right_wall = Plane::new();
        right_wall.set_transform(
            Matrix4::translation(0.0, 0.0, 5.0)
                * Matrix4::rotation_y(std::f64::consts::PI / 4.0)
                * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(10.0, 0.01, 10.0),
        );
        right_wall.set_material(floor.material().clone());

        let mut table_top = Cube::new();
        table_top.set_transform(
            Matrix4::translation(-1.0, 0.5, -6.0)
                // * Matrix4::rotation_y(std::f64::consts::PI / 4.0)
                * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(1.0, 1.0, 0.05),
        );
        let mut material = Material::new();
        material.color = Tuple::color(1.0, 0.0, 0.0);
        table_top.set_material(material.clone());

        let mut leg1 = Cube::new();
        leg1.set_transform(
            Matrix4::translation(-1.95, 0.0, -6.95)
                // * Matrix4::rotation_y(std::f64::consts::PI / 4.0)
                // * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(0.05, 0.5, 0.05),
        );
        let mut material = Material::new();
        material.color = Tuple::color(1.0, 1.0, 0.0);
        leg1.set_material(material.clone());

        let mut leg2 = leg1.clone();
        leg2.set_transform(
            Matrix4::translation(-0.05, 0.0, -6.95)
                // * Matrix4::rotation_y(std::f64::consts::PI / 4.0)
                // * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(0.05, 0.5, 0.05),
        );
        let mut material = Material::new();
        material.color = Tuple::color(0.0, 1.0, 1.0);
        leg2.set_material(material.clone());

        let mut x_pos: f64 = 1.5;
        let mut y_pos: f64 = 2.0;
        let mut z_pos: f64 = -1.95;
        let pos_incr: f64 = 0.2;
        let mut b = Cube::new();
        b.set_transform(
            Matrix4::translation(x_pos, 1.0, z_pos)
                // * Matrix4::rotation_y(std::f64::consts::PI / 4.0)
                // * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(0.05, 0.05, 0.05),
        );

        let mut vbb: Vec<Box<dyn Shape>> = vec![Box::new(b.clone())];
        loop {
            loop {
                x_pos += pos_incr;
                if x_pos > 3.0 {
                    x_pos = 1.5;
                    y_pos -= pos_incr / 2.0;
                    z_pos -= pos_incr;

                    let mut material = Material::new();
                    material.color = Tuple::color(x_pos / 10.0, y_pos / 10.0, y_pos / x_pos);
                    b.set_material(material.clone());

                    break;
                }

                b.set_transform(
                    Matrix4::translation(x_pos + pos_incr, y_pos, z_pos - pos_incr)
                // * Matrix4::rotation_y(std::f64::consts::PI / 4.0)
                // * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(0.05, 0.05, 0.05),
                );
                vbb.push(Box::new(b.clone()));
            }

            if z_pos < -4.0 {
                break;
            }
        }

        let from = Tuple::point(0.0, 1.5, -12.0);
        let to = Tuple::point(0.0, 1.0, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let transform = view_transform(from, to, up);

        let (display_x, display_y) = (60, 40);
        // let (display_x, display_y) = (3456, 2234);
        let camera =
            Camera::new(display_x, display_y, std::f64::consts::PI / 3.0).with_transform(transform);

        world.add_shape(Box::new(floor));
        world.add_shape(Box::new(left_wall));
        world.add_shape(Box::new(right_wall));
        world.add_shape(Box::new(table_top));
        world.add_shape(Box::new(leg1));
        world.add_shape(Box::new(leg2));
        for b_elem in vbb {
            world.add_shape(b_elem);
        }

        let image = world.render(camera);
        let rc = image.write_ppm("test_chap_12_4_putting_it_together.ppm");
        let chk = rc.is_ok();

        if chk {
            Ok(())
        } else {
            Err("Chapter 12_4 Putting It  Together".into())
        }
    }

    /// Chap 13 - A ray misses a cylinder
    #[test]
    fn test_chap_13_1() -> Result<(), String> {
        let cyl = Cylinder::new();

        let direction = Tuple::vector(0.0, 1.0, 0.0).normalize();
        let r = Ray::new(Tuple::point(1.0, 0.0, 0.0), direction);
        let xs = cyl.local_intersect(&r);
        let chk = xs.count() == 0;

        let direction = Tuple::vector(0.0, 1.0, 0.0).normalize();
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), direction);
        let xs = cyl.local_intersect(&r);
        let chk = chk && xs.count() == 0;

        let direction = Tuple::vector(1.0, 1.0, 1.0).normalize();
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), direction);
        let xs = cyl.local_intersect(&r);
        let chk = chk && xs.count() == 0;

        if chk {
            Ok(())
        } else {
            Err("A ray misses a cylinder".into())
        }
    }

    /// Chap 13 - A ray strikes a cylinder
    #[test]
    fn test_chap_13_2() -> Result<(), String> {
        let cyl = Cylinder::new();

        let direction = Tuple::vector(0.0, 0.0, 1.0).normalize();
        let r = Ray::new(Tuple::point(1.0, 0.0, -5.0), direction);
        let xs = cyl.local_intersect(&r);
        let chk = xs.count() == 2 && approx_eq(xs[0].t, 5.0) && approx_eq(xs[1].t, 5.0);

        let direction = Tuple::vector(0.0, 0.0, 1.0).normalize();
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), direction);
        let xs = cyl.local_intersect(&r);
        let chk = chk && xs.count() == 2 && approx_eq(xs[0].t, 4.0) && approx_eq(xs[1].t, 6.0);

        let direction = Tuple::vector(0.1, 1.0, 1.0).normalize();
        let r = Ray::new(Tuple::point(0.5, 0.0, -5.0), direction);
        let xs = cyl.local_intersect(&r);
        let chk = chk
            && xs.count() == 2
            && approx_eq(xs[0].t, 6.80798191702732)
            && approx_eq(xs[1].t, 7.088723439378861);

        if chk {
            Ok(())
        } else {
            loge!("test_chap_13_2", "hit t's:: t0:{} t1:{}", xs[0].t, xs[1].t);
            Err("A ray misses a cylinder".into())
        }
    }

    /// Chap 13 - Normal vector on a cylinder
    #[test]
    fn test_chap_13_3() -> Result<(), String> {
        let cyl = Cylinder::new();

        let n = cyl.local_normal_at(Tuple::point(1.0, 0.0, 0.0));
        let chk = n.approx_eq(Tuple::vector(1.0, 0.0, 0.0));

        let n = cyl.local_normal_at(Tuple::point(0.0, 5.0, -1.0));
        let chk = chk && n.approx_eq(Tuple::vector(0.0, 0.0, -1.0));

        let n = cyl.local_normal_at(Tuple::point(0.0, -2.0, 1.0));
        let chk = chk && n.approx_eq(Tuple::vector(0.0, 0.0, 1.0));

        let n = cyl.local_normal_at(Tuple::point(-1.0, 1.0, 0.0));
        let chk = chk && n.approx_eq(Tuple::vector(-1.0, 0.0, 0.0));

        if chk {
            Ok(())
        } else {
            Err("Normal vector on a cylinder".into())
        }
    }

    /// Chap 13 - The default minimum and maximum for a cylinder
    #[test]
    fn test_chap_13_4() -> Result<(), String> {
        let cyl = Cylinder::new();

        let chk = cyl.minimum == f64::NEG_INFINITY && cyl.maximum == f64::INFINITY;
        if chk {
            Ok(())
        } else {
            loge!("test_chap_13_4", "minimum: {}", cyl.minimum);
            Err("The default minimum and maximum for a cylinder".into())
        }
    }

    /// Chap 13 - Intersecting a constrained cylinder
    #[test]
    fn test_chap_13_5() -> Result<(), String> {
        let mut cyl = Cylinder::new();
        cyl.minimum = 1.0;
        cyl.maximum = 2.0;

        let origin = Tuple::point(0.0, 1.5, 0.0);
        let direction = Tuple::vector(0.1, 1.0, 0.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cyl.local_intersect(&r);
        let chk = xs.count() == 0;

        let origin = Tuple::point(0.0, 3.0, -5.0);
        let direction = Tuple::vector(0.0, 0.0, 1.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cyl.local_intersect(&r);
        let chk = chk && xs.count() == 0;

        let origin = Tuple::point(0.0, 0.0, -5.0);
        let direction = Tuple::vector(0.0, 0.0, 1.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cyl.local_intersect(&r);
        let chk = chk && xs.count() == 0;

        let origin = Tuple::point(0.0, 2.0, 0.0);
        let direction = Tuple::vector(0.0, 0.0, 1.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cyl.local_intersect(&r);
        let chk = chk && xs.count() == 0;

        let origin = Tuple::point(0.0, 1.0, 0.0);
        let direction = Tuple::vector(0.0, 0.0, 1.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cyl.local_intersect(&r);
        let chk = chk && xs.count() == 0;

        let origin = Tuple::point(0.0, 1.5, -2.0);
        let direction = Tuple::vector(0.0, 0.0, 1.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cyl.local_intersect(&r);
        let chk = chk && xs.count() == 2;

        if chk {
            Ok(())
        } else {
            Err("Intersecting a constrained cylinder".into())
        }
    }

    /// Chap 13 - The default closed value for a cylinder
    #[test]
    fn test_chap_13_6() -> Result<(), String> {
        let cyl = Cylinder::new();

        let chk = !cyl.closed;
        if chk {
            Ok(())
        } else {
            Err("The default closed value for a cylinder".into())
        }
    }

    /// Chap 13 - Intersecting the caps of a closed cylinder
    #[test]
    fn test_chap_13_7() -> Result<(), String> {
        let mut cyl = Cylinder::new();
        cyl.minimum = 1.0;
        cyl.maximum = 2.0;
        cyl.closed = true;

        // | point        | direction   | count |
        let origin = Tuple::point(0.0, 3.0, 0.0);
        let direction = Tuple::vector(0.0, -1.0, 0.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cyl.local_intersect(&r);
        let chk = xs.count() == 2;

        let origin = Tuple::point(0.0, 3.0, -2.0);
        let direction = Tuple::vector(0.0, -1.0, 2.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cyl.local_intersect(&r);
        let chk = chk && xs.count() == 2;

        let origin = Tuple::point(0.0, 4.0, -2.0); // corner case
        let direction = Tuple::vector(0.0, -1.0, 1.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cyl.local_intersect(&r);
        let chk = chk && xs.count() == 2;

        let origin = Tuple::point(0.0, 0.0, -2.0);
        let direction = Tuple::vector(0.0, 1.0, 2.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cyl.local_intersect(&r);
        let chk = chk && xs.count() == 2;

        let origin = Tuple::point(0.0, -1.0, -2.0); // corner case
        let direction = Tuple::vector(0.0, 1.0, 1.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cyl.local_intersect(&r);
        let chk = chk && xs.count() == 2;

        if chk {
            Ok(())
        } else {
            Err("Intersecting the caps of a closed cylinder".into())
        }
    }

    /// Chap 13 - The normal vector on a cylinder's end caps
    #[test]
    fn test_chap_13_8() -> Result<(), String> {
        let mut cyl = Cylinder::new();
        cyl.minimum = 1.0;
        cyl.maximum = 2.0;
        cyl.closed = true;

        let point = Tuple::point(0.0, 1.0, 0.0);
        let n = cyl.local_normal_at(point);
        let chk = n.approx_eq(Tuple::vector(0.0, -1.0, 0.0));

        let point = Tuple::point(0.5, 1.0, 0.0);
        let n = cyl.local_normal_at(point);
        let chk = chk && n.approx_eq(Tuple::vector(0.0, -1.0, 0.0));

        let point = Tuple::point(0.0, 1.0, 0.5);
        let n = cyl.local_normal_at(point);
        let chk = chk && n.approx_eq(Tuple::vector(0.0, -1.0, 0.0));

        let point = Tuple::point(0.0, 2.0, 0.0);
        let n = cyl.local_normal_at(point);
        let chk = chk && n.approx_eq(Tuple::vector(0.0, 1.0, 0.0));

        let point = Tuple::point(0.5, 2.0, 0.0);
        let n = cyl.local_normal_at(point);
        let chk = chk && n.approx_eq(Tuple::vector(0.0, 1.0, 0.0));

        let point = Tuple::point(0.0, 2.0, 0.5);
        let n = cyl.local_normal_at(point);
        let chk = chk && n.approx_eq(Tuple::vector(0.0, 1.0, 0.0));

        if chk {
            Ok(())
        } else {
            Err("The normal vector on a cylinder's end caps".into())
        }
    }

    /// Chap 13 - Intersecting a cone with a ray
    #[test]
    fn test_chap_13_9() -> Result<(), String> {
        let cone = Cone::new();

        let origin = Tuple::point(0.0, 0.0, -5.0);
        let direction = Tuple::vector(0.0, 0.0, 1.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cone.local_intersect(&r);
        let chk = xs.count() == 2 && approx_eq(xs[0].t, 5.0) && approx_eq(xs[1].t, 5.0);

        let origin = Tuple::point(0.0, 0.0, -5.0);
        let direction = Tuple::vector(1.0, 1.0, 1.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cone.local_intersect(&r);
        let chk = chk
            && xs.count() == 2
            && approx_eq(xs[0].t, 8.660254037844386)
            && approx_eq(xs[1].t, 8.660254037844386);

        let origin = Tuple::point(1.0, 1.0, -5.0);
        let direction = Tuple::vector(-0.5, -1.0, 1.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cone.local_intersect(&r);
        let chk = chk
            && xs.count() == 2
            && approx_eq(xs[0].t, 4.550055679356349)
            && approx_eq(xs[1].t, 49.449944320643645);

        if chk {
            Ok(())
        } else {
            if xs.count() > 1 {
                loge!("test_chap_13_9", "xs[0].t:{} xs[1].t:{}", xs[0].t, xs[1].t);
            }
            Err("A ray misses a cone".into())
        }
    }

    /// Chap x - Intersecting a cone with a ray parallel to one of its halves
    #[test]
    fn test_chap_13_10() -> Result<(), String> {
        let shape = Cone::new();

        let origin = Tuple::point(0.0, 0.0, -1.0);
        let direction = Tuple::vector(0.0, 1.0, 1.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = shape.local_intersect(&r);

        let chk = xs.count() == 1 && approx_eq(xs[0].t, 0.3535533905932738);
        if chk {
            Ok(())
        } else {
            if xs.count() > 0 {
                loge!("test_chap_13_10", "xs[0].t:{} ", xs[0].t);
            }
            Err("Intersecting a cone with a ray parallel to one of its halves".into())
        }
    }

    /// Chap x - Intersecting a cone's end caps
    #[test]
    fn test_chap_13_11() -> Result<(), String> {
        let mut shape = Cone::new();
        shape.minimum = -0.5;
        shape.maximum = 0.5;
        shape.closed = true;

        let origin = Tuple::point(0.0, 0.0, -5.0);
        let direction = Tuple::vector(0.0, 1.0, 0.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = shape.local_intersect(&r);
        let chk = xs.count() == 0;

        let origin = Tuple::point(0.0, 0.0, -0.25);
        let direction = Tuple::vector(0.0, 1.0, 1.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = shape.local_intersect(&r);
        let chk = chk && xs.count() == 2;

        let origin = Tuple::point(0.0, 0.0, -0.25);
        let direction = Tuple::vector(0.0, 1.0, 0.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = shape.local_intersect(&r);
        let chk = chk && xs.count() == 4;

        if chk {
            Ok(())
        } else {
            Err("Intersecting a cone's end caps".into())
        }
    }

    /// Chap x - Computing the normal vector on a cone
    #[test]
    fn test_chap_13_12() -> Result<(), String> {
        let shape = Cone::new();

        let point = Tuple::point(0.0, 0.0, 0.0);
        let n = shape.local_normal_at(point);
        let chk = n.approx_eq(Tuple::vector(0.0, 0.0, 0.0));

        let point = Tuple::point(1.0, 1.0, 1.0);
        let n = shape.local_normal_at(point);
        let chk = chk && n.approx_eq(Tuple::vector(1.0, -2.0_f64.sqrt(), 1.0));

        let point = Tuple::point(-1.0, -1.0, 0.0);
        let n = shape.local_normal_at(point);
        let chk = chk && n.approx_eq(Tuple::vector(-1.0, 1.0, 0.0));

        if chk {
            Ok(())
        } else {
            loge!("test_chap_13_12", "n:{}", n);
            Err("Computing the normal vector on a cone".into())
        }
    }
    /// Chap 13 - Putting It Together
    ///
    /// The book leaves this scene open-ended ("render cylinders and cones"), so
    /// this composes a showcase from the cylinder features built in this
    /// chapter: capped solids, a truncated open tube, and a transformed
    /// cylinder, following the same layout as the earlier "putting it together"
    /// tests. Renders to a PPM in the working directory.
    #[test]
    fn test_chap_13_putting_it_all_together() -> Result<(), String> {
        let mut world = World::new();

        let light = Light::point_light(
            Tuple::point(-10.0, 10.0, -10.0),
            Tuple::color(1.0, 1.0, 1.0),
        );
        world.light = Some(light);

        let mut floor_checkers = CheckersPattern::new(WHITE, BLACK);
        floor_checkers.set_transform(Matrix4::scaling(0.5, 0.5, 0.5));

        // Floor - a slightly reflective plane
        let mut floor = Plane::new();
        let mut floor_material = Material::new();
        floor_material.color = Tuple::color(0.8, 0.8, 0.85);
        floor_material.specular = 0.0;
        floor_material.reflective = 0.2;
        floor_material.pattern = Some(Box::new(floor_checkers));
        floor.set_material(floor_material);
        world.add_shape(Box::new(floor));

        // Tall capped cylinder (green)
        let mut tall = Cylinder::new();
        tall.minimum = 0.0;
        tall.maximum = 3.0;
        tall.closed = true;
        tall.set_transform(Matrix4::translation(-1.5, 0.0, 0.5) * Matrix4::scaling(0.5, 1.0, 0.5));
        let mut tall_material = Material::new();
        tall_material.color = Tuple::color(0.1, 0.8, 0.3);
        tall_material.diffuse = 0.7;
        tall_material.specular = 0.3;
        tall.set_material(tall_material);
        world.add_shape(Box::new(tall));

        // Capped drum (red, slightly reflective)
        let mut drum = Cylinder::new();
        drum.minimum = 0.0;
        drum.maximum = 2.5;
        drum.closed = true;
        drum.set_transform(
            Matrix4::translation(2.3, 0.0, -0.5) * Matrix4::scaling(0.29, 1.0, 0.29),
        );
        let mut drum_material = Material::new();
        drum_material.color = Tuple::color(0.9, 0.2, 0.2);
        drum_material.diffuse = 0.7;
        drum_material.specular = 0.3;
        drum_material.reflective = 0.2;
        drum.set_material(drum_material);
        world.add_shape(Box::new(drum));

        let mut cyl_radius = 0.39;
        let mut cyl_maximum = 2.0;
        let mut col_factor = 1.0;
        loop {
            // Thin open tube () - not closed, so you can see through it
            let mut tube = Cylinder::new();
            tube.minimum = 0.0;
            tube.maximum = cyl_maximum;
            tube.closed = false;
            tube.set_transform(
                Matrix4::translation(2.3, 0.0, -0.5)
                    * Matrix4::scaling(cyl_radius, 1.0, cyl_radius),
            );

            let mut mirror = Material::new();
            mirror.color = Tuple::color(0.0, 0.0, 0.0); // near-black base; reflection provides the look
            mirror.ambient = 0.0;
            mirror.diffuse = 0.0;
            mirror.specular = 1.0; // bright highlight where the light hits
            mirror.shininess = 300.0; // tight, sharp highlight (mirror-like, not matte)
            mirror.reflective = 0.94 * col_factor; // perfect mirror; 0.9 for "very polished but not perfect"

            // let mut tube_material = Material::new();
            // tube_material.color =
            //     Tuple::color(0.2 / col_factor, 0.4 / col_factor, 0.19 / col_factor);
            // // tube_material.diffuse = 0.7 * col_factor;
            // tube_material.transparency = 0.95;
            // tube_material.reflective = 1.0;
            // tube_material.shininess = 300.0 * col_factor;
            // tube_material.specular = 1.0 * col_factor;
            // tube.set_material(tube_material);
            tube.set_material(mirror.clone());
            world.add_shape(Box::new(tube));

            if cyl_radius > 1.5 {
                break;
            }

            cyl_radius += 0.4;
            cyl_maximum -= 0.4;
            col_factor *= 0.9;
        }

        // Thin open tube (blue) - not closed, so you can see through it
        let mut tube = Cylinder::new();
        tube.minimum = 0.0;
        tube.maximum = 2.0;
        tube.closed = false;
        tube.set_transform(Matrix4::translation(0.4, 0.0, 1.6) * Matrix4::scaling(0.3, 1.0, 0.3));
        let mut tube_material = Material::new();
        tube_material.color = Tuple::color(0.2, 0.4, 0.9);
        tube_material.diffuse = 0.7;
        tube_material.specular = 0.3;
        tube.set_material(tube_material);
        world.add_shape(Box::new(tube));

        // Tilted capped cylinder lying on its side (yellow)
        let mut tilted = Cylinder::new();
        tilted.minimum = 0.0;
        tilted.maximum = 2.0;
        tilted.closed = true;
        tilted.set_transform(
            Matrix4::translation(0.0, 0.25, -1.5)
                * Matrix4::rotation_y(std::f64::consts::PI / 4.0)
                * Matrix4::rotation_z(std::f64::consts::PI / 4.0)
                * Matrix4::scaling(0.25, 1.5, 0.25),
        );
        let mut tilted_material = Material::new();
        tilted_material.color = Tuple::color(0.9, 0.8, 0.1);
        tilted_material.diffuse = 0.7;
        tilted_material.specular = 0.3;
        tilted.set_material(tilted_material);
        world.add_shape(Box::new(tilted));

        // Thin cone ()
        let mut cone = Cone::new();
        cone.minimum = -1.4;
        cone.maximum = 1.4;
        cone.closed = false;
        cone.set_transform(
            Matrix4::translation(-2.0, 1.0, -1.6)
                * Matrix4::rotation_y(std::f64::consts::PI / 4.0)
                * Matrix4::rotation_z(std::f64::consts::PI / 4.0)
                * Matrix4::scaling(0.3, 1.0, 0.3),
        );
        let mut cone_material = Material::new();
        cone_material.color = Tuple::color(0.7, 0.9, 0.3);
        cone_material.diffuse = 0.7;
        cone_material.specular = 0.3;
        cone.set_material(cone_material);
        world.add_shape(Box::new(cone));

        // View transform / camera
        let from = Tuple::point(0.0, 2.5, -7.0);
        let to = Tuple::point(0.0, 1.0, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let transform = view_transform(from, to, up);
        // let (display_x, display_y) = (60, 40);
        let (display_x, display_y) = (3456, 2234);
        let camera =
            Camera::new(display_x, display_y, std::f64::consts::PI / 3.0).with_transform(transform);

        let image = world.render(camera);
        let rc = image.write_ppm("test_chap_13_putting_it_all_together.ppm");

        if rc.is_ok() {
            Ok(())
        } else {
            Err("Chapter 13 Putting It Together".into())
        }
    }
}
