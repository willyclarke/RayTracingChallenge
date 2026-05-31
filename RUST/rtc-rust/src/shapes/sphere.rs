//! Sphere defintion
//!

use crate::intersection::{Intersection, Intersections};
use crate::log::*;
use crate::ray::Ray;
use crate::shape::{Shape, ShapeData};
use crate::tuple::Tuple;
use std::fmt;

#[derive(Debug, Clone, Copy)]
pub struct Sphere {
    pub data: ShapeData,
}

impl Sphere {
    pub fn new() -> Self {
        Self {
            data: ShapeData::new(),
        }
    }
}

impl Default for Sphere {
    fn default() -> Self {
        Self::new()
    }
}

impl fmt::Display for Sphere {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "{}Sphere ID:{}{}{}{}",
            Color::Yellow,
            Color::Reset,
            Color::Green,
            self.data.id,
            Color::Reset
        )?;

        Ok(())
    }
}

impl Shape for Sphere {
    fn data(&self) -> &ShapeData {
        &self.data
    }
    fn data_mut(&mut self) -> &mut ShapeData {
        &mut self.data
    }

    fn local_intersect(&self, ray: &Ray) -> Intersections {
        let sphere_to_ray = ray.origin - Tuple::point(0.0, 0.0, 0.0);
        let a = ray.direction.dot(ray.direction);
        let b = 2.0 * ray.direction.dot(sphere_to_ray);
        let c = sphere_to_ray.dot(sphere_to_ray) - 1.0;

        let discriminant = b * b - 4.0 * a * c;
        let mut xs = Intersections::new();
        if discriminant < 0.0 {
            return xs;
        }

        let t1 = (-b - discriminant.sqrt()) / (2.0 * a);
        let t2 = (-b + discriminant.sqrt()) / (2.0 * a);
        xs.push(Intersection::new(t1, self.id()));
        xs.push(Intersection::new(t2, self.id()));
        xs
    }

    fn local_normal_at(&self, object_point: Tuple) -> Tuple {
        // sphere-specific normal logic goes here
        object_point - Tuple::point(0.0, 0.0, 0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::light::Light;
    // use crate::canvas::Canvas;
    use crate::material::Material;
    use crate::math::approx_eq;
    use crate::matrix::Matrix4;
    use crate::{logd, loge, logi, tuple::Tuple};

    /// Chap 5 - A ray intersects a sphere at two points
    #[test]
    fn test_chap_5_3() -> Result<(), String> {
        let ray = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = Sphere::new();

        let xs = s.intersect(&ray);
        let chk = xs.count() == 2_usize;
        let chk = chk && approx_eq(xs[0].t, 4.0);
        let chk = chk && approx_eq(xs[1].t, 6.0);

        if chk {
            Ok(())
        } else {
            logi!("test_chap_5_3", "chk:{:?}", chk);
            logd!("test_chap_5_3", "xs.count():{:?}", xs.count());
            loge!("test_chap_5_3", "chk:{:?}", chk);
            Err("A ray intersects a sphere at two points".into())
        }
    }

    /// Chap 5 -A ray intersects a sphere at a tangent
    #[test]
    fn test_chap_5_4() -> Result<(), String> {
        let ray = Ray::new(Tuple::point(0.0, 1.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = Sphere::new();

        let xs = s.intersect(&ray);

        let chk = xs.count() == 2_usize;
        let chk = chk && approx_eq(xs[0].t, 5.0);
        let chk = chk && approx_eq(xs[1].t, 5.0);

        if chk {
            Ok(())
        } else {
            loge!("test_chap_5_4", "num intersects:{}", xs.count());
            Err("A ray intersects a sphere at a tangent".into())
        }
    }

    /// Chap 5 - A ray misses a sphere
    #[test]
    fn test_chap_5_5() -> Result<(), String> {
        let ray = Ray::new(Tuple::point(0.0, 2.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = Sphere::new();

        let xs = s.intersect(&ray);

        let chk = xs.is_empty();
        let chk = chk && xs.count() == 0_usize;

        if chk {
            Ok(())
        } else {
            Err("A ray misses a sphere".into())
        }
    }

    /// Chap 5 - A ray originates inside a sphere
    #[test]
    fn test_chap_5_6() -> Result<(), String> {
        let ray = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = Sphere::new();

        let xs = s.intersect(&ray);

        let chk = xs.count() == 2_usize;
        let chk = chk && approx_eq(xs[0].t, -1.0);
        let chk = chk && approx_eq(xs[1].t, 1.0);
        if chk {
            Ok(())
        } else {
            Err("A ray originates inside a sphere".into())
        }
    }

    /// Chap 5 - A sphere is behind a ray
    #[test]
    fn test_chap_5_7() -> Result<(), String> {
        let ray = Ray::new(Tuple::point(0.0, 0.0, 5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = Sphere::new();

        let xs = s.intersect(&ray);

        let chk = xs.count() == 2_usize;
        let chk = chk && approx_eq(xs[0].t, -6.0);
        let chk = chk && approx_eq(xs[1].t, -4.0);
        if chk {
            Ok(())
        } else {
            Err("A sphere is behind a ray".into())
        }
    }

    /// Chap 5 - An intersection encapsulates t and object
    #[test]
    fn test_chap_5_8() -> Result<(), String> {
        let s = Sphere::new();
        let i = Intersection::new(3.5, s.id());
        let chk = approx_eq(i.t, 3.5);
        let chk = chk && i.object_id == s.id();
        if chk {
            Ok(())
        } else {
            Err("An intersection encapsulates t and object".into())
        }
    }

    /// Chap 5 - Aggregating intersections
    #[test]
    fn test_chap_5_9() -> Result<(), String> {
        let s = Sphere::new();
        let i1 = Intersection::new(1.0, s.id());
        let i2 = Intersection::new(2.0, s.id());
        let mut xs = Intersections::new();
        xs.push(i1);
        xs.push(i2);
        let chk = xs.count() == 2_usize;
        let chk = chk && approx_eq(1.0, xs[0].t);
        let chk = chk && approx_eq(2.0, xs[1].t);
        if chk {
            Ok(())
        } else {
            Err("Aggregating intersections".into())
        }
    }

    /// Chap 5 - Intersect sets the object on the intersection
    #[test]
    fn test_chap_5_10() -> Result<(), String> {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = Sphere::new();

        let xs = s.intersect(&r);

        let chk = xs.count() == 2;
        let chk = chk && xs[0].object_id == s.id();
        let chk = chk && xs[1].object_id == s.id();

        if chk {
            Ok(())
        } else {
            Err("Intersect sets the object on the intersection".into())
        }
    }

    /// Chap 5 - The hit, when all intersections have positive t
    #[test]
    fn test_chap_5_11() -> Result<(), String> {
        let s = Sphere::new();
        let i1 = Intersection::new(1.0, s.id());
        let i2 = Intersection::new(2.0, s.id());
        let mut xs = Intersections::new();
        xs.push(i1);
        xs.push(i2);
        let i = xs.hit().unwrap();
        let chk = i.object_id == i1.object_id;

        if chk {
            Ok(())
        } else {
            Err("The hit, when all intersections have positive t".into())
        }
    }

    /// Chap x - The hit, when some intersections have negative t
    #[test]
    fn test_chap_5_12() -> Result<(), String> {
        let s = Sphere::new();
        let i1 = Intersection::new(-1.0, s.id());
        let i2 = Intersection::new(1.0, s.id());
        let mut xs = Intersections::new();
        xs.push(i1);
        xs.push(i2);
        let i = xs.hit().unwrap();
        let chk = i.object_id == i2.object_id && approx_eq(i.t, 1.0);
        // logi!("test_chap_5_12", "i1:{:?}", i1);
        // logi!("test_chap_5_12", "i2:{:?}", i2);
        // logi!("test_chap_5_12", "i:{:?}", i);
        if chk {
            Ok(())
        } else {
            Err("The hit, when some intersections have negative t".into())
        }
    }

    /// Chap 5 - The hit, when all intersections have negative t
    #[test]
    fn test_chap_5_13() -> Result<(), String> {
        let s = Sphere::new();
        let i1 = Intersection::new(-2.0, s.id());
        let i2 = Intersection::new(-1.0, s.id());
        let mut xs = Intersections::new();
        xs.push(i1);
        xs.push(i2);
        logi!("test_chap_5_13", "i1:{:?}", i1);
        logi!("test_chap_5_13", "i2:{:?}", i2);
        let i = xs.hit();
        let chk = i.is_none();
        if chk {
            Ok(())
        } else {
            Err("The hit, when all intersections have negative t".into())
        }
    }

    /// Chap 5 - The hit is always the lowest nonnegative intersection
    #[test]
    fn test_chap_5_14() -> Result<(), String> {
        let s = Sphere::new();
        let i1 = Intersection::new(5.0, s.id());
        let i2 = Intersection::new(7.0, s.id());
        let i3 = Intersection::new(-3.0, s.id());
        let i4 = Intersection::new(2.0, s.id());
        let mut xs = Intersections::new();
        xs.push(i1);
        xs.push(i2);
        xs.push(i3);
        xs.push(i4);
        logi!("test_chap_5_14", "i1:{:?}", i1);
        logi!("test_chap_5_14", "i2:{:?}", i2);
        logi!("test_chap_5_14", "i3:{:?}", i3);
        logi!("test_chap_5_14", "i4:{:?}", i4);
        let i = xs.hit().unwrap();
        let chk = i.object_id == i4.object_id && approx_eq(i.t, 2.0);
        if chk {
            Ok(())
        } else {
            Err("The hit is always the lowest nonnegative intersection".into())
        }
    }

    /// Chap 5 - A sphere's default transformation
    #[test]
    fn test_chap_5_17() -> Result<(), String> {
        let s = Sphere::new();
        let chk = s.transform().approx_eq(Matrix4::identity());
        if chk {
            Ok(())
        } else {
            Err("A sphere's default transformation".into())
        }
    }

    /// Chap 5 - Changing a sphere's transformation
    #[test]
    fn test_chap_5_18() -> Result<(), String> {
        let mut s = Sphere::new();
        let t = Matrix4::translation(2.0, 3.0, 4.0);
        s.set_transform(t);
        let chk = s.transform().approx_eq(t);
        if chk {
            Ok(())
        } else {
            Err("Changing a sphere's transformation".into())
        }
    }

    /// Chap 5 - Intersecting a scaled sphere with a ray
    #[test]
    fn test_chap_5_19() -> Result<(), String> {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let mut s = Sphere::new();
        s.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));
        let xs = s.intersect(&r);
        let chk = xs.count() == 2;
        let chk = chk && approx_eq(xs[0].t, 3.0);
        // let chk = chk && approx_eq(xs[1].t, 7.0);
        if chk {
            Ok(())
        } else {
            Err("Intersecting a scaled sphere with a ray".into())
        }
    }

    /// Chap 6 - The normal on a sphere at a point on the x axis
    #[test]
    fn test_chap_6_1() -> Result<(), String> {
        let s = Sphere::new();
        let n = s.normal_at(Tuple::point(1.0, 0.0, 0.0));
        let chk = n.approx_eq(Tuple::vector(1.0, 0.0, 0.0));
        if chk {
            Ok(())
        } else {
            loge!("test_chap_6_1", "n:{}", n);
            Err("The normal on a sphere at a point on the x axis".into())
        }
    }

    /// Chap 6 - The normal on a sphere at a point on the y axis
    #[test]
    fn test_chap_6_2() -> Result<(), String> {
        let s = Sphere::new();
        let n = s.normal_at(Tuple::point(0.0, 1.0, 0.0));
        let chk = n.approx_eq(Tuple::vector(0.0, 1.0, 0.0));
        if chk {
            Ok(())
        } else {
            loge!("test_chap_6_2", "n:{}", n);
            Err("The normal on a sphere at a point on the y axis".into())
        }
    }

    /// Chap 6 - The normal on a sphere at a point on the z axis
    #[test]
    fn test_chap_6_3() -> Result<(), String> {
        let s = Sphere::new();
        let n = s.normal_at(Tuple::point(0.0, 0.0, 1.0));
        let chk = n.approx_eq(Tuple::vector(0.0, 0.0, 1.0));
        if chk {
            Ok(())
        } else {
            loge!("test_chap_6_3", "n:{}", n);
            Err("The normal on a sphere at a point on the z axis".into())
        }
    }

    /// Chap 6 - The normal on a sphere at a nonaxial point
    #[test]
    fn test_chap_6_4() -> Result<(), String> {
        let s = Sphere::new();
        let sqrt3o3 = 3_f64.sqrt() / 3_f64;
        let n = s.normal_at(Tuple::point(sqrt3o3, sqrt3o3, sqrt3o3));
        let chk = n.approx_eq(Tuple::vector(sqrt3o3, sqrt3o3, sqrt3o3));
        if chk {
            Ok(())
        } else {
            loge!("test_chap_6_4", "n:{}", n);
            Err("The normal on a sphere at a nonaxial point".into())
        }
    }

    /// Chap 6 - The normal is a normalized vector
    #[test]
    fn test_chap_6_5() -> Result<(), String> {
        let s = Sphere::new();
        let sqrt3o3 = 3_f64.sqrt() / 3_f64;
        let n = s.normal_at(Tuple::point(sqrt3o3, sqrt3o3, sqrt3o3));
        let chk = n.normalize().approx_eq(n);
        if chk {
            Ok(())
        } else {
            loge!("test_chap_6_5", "chk:{:?}", n);
            Err("The normal is a normalized vector".into())
        }
    }

    /// Chap 6 - Computing the normal on a translated sphere
    #[test]
    fn test_chap_6_6() -> Result<(), String> {
        let mut s = Sphere::new();
        s.set_transform(Matrix4::translation(0.0, 1.0, 0.0));
        let n = s.normal_at(Tuple::point(
            0.0,
            1.0 + std::f64::consts::FRAC_1_SQRT_2,
            -std::f64::consts::FRAC_1_SQRT_2,
        ));
        let chk = n.approx_eq(Tuple::vector(
            0.0,
            std::f64::consts::FRAC_1_SQRT_2,
            -std::f64::consts::FRAC_1_SQRT_2,
        ));
        if chk {
            Ok(())
        } else {
            loge!("test_chap_6_6", "chk:{:?}", n);
            Err("Computing the normal on a translated sphere".into())
        }
    }

    /// Chap 6 - Computing the normal on a transformed sphere
    #[test]
    fn test_chap_6_7() -> Result<(), String> {
        let mut s = Sphere::new();
        let m = Matrix4::scaling(1.0, 0.5, 1.0) * Matrix4::rotation_z(std::f64::consts::PI / 5.0);
        s.set_transform(m);
        let n = s.normal_at(Tuple::point(
            0.0,
            std::f64::consts::SQRT_2 / 2_f64,
            -std::f64::consts::SQRT_2 / 2_f64,
        ));
        let chk = n.approx_eq(Tuple::vector(0.0, 0.970142500145, -0.242535625036));
        if chk {
            Ok(())
        } else {
            loge!("test_chap_6_7", "chk:{:}", n);
            Err("Computing the normal on a transformed sphere".into())
        }
    }

    /// Chap 6 - Reflecting a vector approaching at 45°
    #[test]
    fn test_chap_6_8() -> Result<(), String> {
        let v = Tuple::vector(1.0, -1.0, 0.0);
        let n = Tuple::vector(0.0, 1.0, 0.0);
        let r = v.reflect(n);

        let chk = Tuple::vector(1.0, 1.0, 0.0).approx_eq(r);
        if chk {
            Ok(())
        } else {
            Err("Reflecting a vector approaching at 45°".into())
        }
    }

    /// Chap 6 - Reflecting a vector off a slanted surface
    #[test]
    fn test_chap_6_9() -> Result<(), String> {
        let sqrt2o2 = 2_f64.sqrt() / 2_f64;
        let v = Tuple::vector(0.0, -1.0, 0.0);
        let n = Tuple::vector(sqrt2o2, sqrt2o2, 0.0);
        let r = v.reflect(n);

        let chk = Tuple::vector(1.0, 0.0, 0.0).approx_eq(r);
        if chk {
            Ok(())
        } else {
            Err("Reflecting a vector off a slanted surface".into())
        }
    }

    /// Chap 6 - A point light has a position and intensity
    #[test]
    fn test_chap_6_10() -> Result<(), String> {
        let intensity = Tuple::color(1.0, 1.0, 1.0);
        let position = Tuple::point(0.0, 0.0, 0.0);
        let light = Light::point_light(position, intensity);
        let chk = position.approx_eq(light.position);
        let chk = chk && intensity.approx_eq(light.intensity);
        if chk {
            logi!("", "Light: {}.", light);
            Ok(())
        } else {
            Err("A point light has a position and intensity".into())
        }
    }

    /// Chap 6 - A sphere may be assigned a material
    #[test]
    fn test_chap_6_12() -> Result<(), String> {
        let mut s = Sphere::new();
        let mut m = Material::new();
        m.ambient = 1.0;
        s.set_material(m);

        let chk = m.approx_eq(*s.material());
        if chk {
            Ok(())
        } else {
            Err("A sphere may be assigned a material".into())
        }
    }
}
