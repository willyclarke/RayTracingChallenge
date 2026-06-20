//! Plane definition
//!
//! An infinite flat surface extending in the xz plane by default.
//! The normal of a plane is constant everywhere.
//!

use crate::intersection::{Intersection, Intersections};
use crate::ray::Ray;
use crate::shape::{Shape, ShapeData};
use crate::tuple::Tuple;

// #[derive(Debug, Clone, Copy)]
#[derive(Debug, Clone)]
pub struct Plane {
    pub data: ShapeData,
}

impl Plane {
    pub fn new() -> Self {
        Self {
            data: ShapeData::new(),
        }
    }
}

impl Default for Plane {
    fn default() -> Self {
        Self::new()
    }
}

impl Shape for Plane {
    fn data(&self) -> &ShapeData {
        &self.data
    }

    fn data_mut(&mut self) -> &mut ShapeData {
        &mut self.data
    }

    ///
    /// The normal of a plane is constant everywhere
    ///
    fn local_normal_at(&self, _point: Tuple) -> Tuple {
        Tuple::vector(0.0, 1.0, 0.0)
    }

    fn local_intersect(&self, ray: &Ray) -> Intersections {

        if ray.direction.y.abs() < crate::math::EPSILON {
            return Intersections::new();
        }

        let mut i = Intersections::new();
        let t = -ray.origin.y / ray.direction.y;
        let xs = Intersection::new(t, self.id());
        i.push(xs);
        i
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::approx_eq;
    use crate::ray::Ray;
    use crate::tuple::Tuple;

    /// Chap 9 - The normal of a plane is constant everywhere
    #[test]
    fn test_chap_9_1() -> Result<(), String> {
        let p = Plane::new();
        let n1 = p.local_normal_at(Tuple::point(0.0, 0.0, 0.0));
        let n2 = p.local_normal_at(Tuple::point(10.0, 0.0, -10.0));
        let n3 = p.local_normal_at(Tuple::point(-5.0, 0.0, 150.0));
        let expected = Tuple::vector(0.0, 1.0, 0.0);
        let chk = n1.approx_eq(expected) && n2.approx_eq(expected) && n3.approx_eq(expected);
        if chk {
            Ok(())
        } else {
            Err("The normal of a plane is constant everywhere".into())
        }
    }

    /// Chap 9 - Intersect with a ray parallel to the plane
    #[test]
    fn test_chap_9_2() -> Result<(), String> {
        let p = Plane::new();
        let r = Ray::new(Tuple::point(0.0, 10.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs = p.local_intersect(&r);
        let chk = xs.is_empty();
        if chk {
            Ok(())
        } else {
            Err("Intersect with a ray parallel to the plane".into())
        }
    }

    /// Chap 9 - Intersect with a coplanar ray
    #[test]
    fn test_chap_9_3() -> Result<(), String> {
        let p = Plane::new();
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs = p.local_intersect(&r);
        let chk = xs.is_empty();
        if chk {
            Ok(())
        } else {
            Err("Intersect with a coplanar ray".into())
        }
    }

    /// Chap 9 - A ray intersecting a plane from above
    #[test]
    fn test_chap_9_4() -> Result<(), String> {
        let p = Plane::new();
        let r = Ray::new(Tuple::point(0.0, 1.0, 0.0), Tuple::vector(0.0, -1.0, 0.0));
        let xs = p.local_intersect(&r);
        let chk = xs.count() == 1;
        let chk = chk && approx_eq(xs[0].t, 1.0);
        let chk = chk && xs[0].object_id == p.id();
        if chk {
            Ok(())
        } else {
            Err("A ray intersecting a plane from above".into())
        }
    }

    /// Chap 9 - A ray intersecting a plane from below
    #[test]
    fn test_chap_9_5() -> Result<(), String> {
        let p = Plane::new();
        let r = Ray::new(Tuple::point(0.0, -1.0, 0.0), Tuple::vector(0.0, 1.0, 0.0));
        let xs = p.local_intersect(&r);
        let chk = xs.count() == 1;
        let chk = chk && approx_eq(xs[0].t, 1.0);
        let chk = chk && xs[0].object_id == p.id();
        if chk {
            Ok(())
        } else {
            Err("A ray intersecting a plane from below".into())
        }
    }
}
