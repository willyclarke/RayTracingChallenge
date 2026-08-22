//! triangle.rs - add support for polygonal models with triangle's.

use crate::bounds::BoundingBox;
use crate::intersection::{Intersection, Intersections};
use crate::math::approx_eq;
use crate::shape::{Shape, ShapeData};
use crate::tuple::Tuple;

/// Möller–Trumbore ray/triangle intersection, shared by `Triangle` and
/// `TriangleUV`. Returns `Some((t, u, v))` on a hit, `None` on a miss.
pub(crate) fn moller_trumbore(
    p1: Tuple,
    e1: Tuple,
    e2: Tuple,
    ray: &crate::ray::Ray,
) -> Option<(f64, f64, f64)> {
    let dir_cross_e2 = ray.direction.cross(e2);
    let det = e1.dot(dir_cross_e2);

    if approx_eq(det.abs(), 0.0) {
        return None;
    }

    let f = 1.0 / det;
    let p1_to_origin = ray.origin - p1;
    let u = f * p1_to_origin.dot(dir_cross_e2);

    if !(0.0..=1.0).contains(&u) {
        return None;
    }

    let origin_cross_e1 = p1_to_origin.cross(e1);
    let v = f * ray.direction.dot(origin_cross_e1);

    if v < 0.0 || (u + v) > 1.0 {
        return None;
    }

    let t = f * e2.dot(origin_cross_e1);
    Some((t, u, v))
}

#[derive(Debug, Clone)]
pub struct Triangle {
    pub data: ShapeData,
    pub p1: Tuple,
    pub p2: Tuple,
    pub p3: Tuple,
    pub e1: Tuple,
    pub e2: Tuple,
    pub normal: Tuple,
}

impl Triangle {
    pub fn new(p1: Tuple, p2: Tuple, p3: Tuple) -> Self {
        let e1 = p2 - p1;
        let e2 = p3 - p1;
        let normal = (e2.cross(e1)).normalize();
        Self {
            data: ShapeData::new(),
            p1,
            p2,
            p3,
            e1,
            e2,
            normal,
        }
    }
}

impl Default for Triangle {
    fn default() -> Self {
        Self::new(
            Tuple::point(0.0, 1.0, 0.0),
            Tuple::point(-1.0, 0.0, 0.0),
            Tuple::point(1.0, 0.0, 0.0),
        )
    }
}

impl Shape for Triangle {
    /// Tight box over the three vertices: min/max per axis of p1, p2, p3.
    fn bounds(&self) -> BoundingBox {
        let mut bb = BoundingBox::empty();
        bb.add_point(self.p1);
        bb.add_point(self.p2);
        bb.add_point(self.p3);
        bb
    }

    fn local_normal_at(&self, _n: Tuple, _hit: Intersection) -> Tuple {
        self.normal
    }

    fn data(&self) -> &ShapeData {
        &self.data
    }

    fn data_mut(&mut self) -> &mut ShapeData {
        &mut self.data
    }

    fn local_intersect(&self, ray: &crate::ray::Ray) -> Intersections {
        let mut xs = Intersections::new();
        if let Some((t, u, v)) = moller_trumbore(self.p1, self.e1, self.e2, ray) {
            xs.push(Intersection::new_with_uv(t, self.id(), u, v));
        }
        xs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{loge, ray::Ray};

    /// Chap 15 - Constructing a triangle
    #[test]
    fn test_chap_15_1() -> Result<(), String> {
        let t = Triangle::new(
            Tuple::point(0.0, 1.0, 0.0),
            Tuple::point(-1.0, 0.0, 0.0),
            Tuple::point(1.0, 0.0, 0.0),
        );
        let chk = t.p1.approx_eq(Tuple::point(0.0, 1.0, 0.0));
        let chk = chk && t.p2.approx_eq(Tuple::point(-1.0, 0.0, 0.0));
        let chk = chk && t.p3.approx_eq(Tuple::point(1.0, 0.0, 0.0));
        let chk = chk && t.e1.approx_eq(Tuple::vector(-1.0, -1.0, 0.0));
        let chk = chk && t.e2.approx_eq(Tuple::vector(1.0, -1.0, 0.0));
        let chk = chk && t.normal.approx_eq(Tuple::vector(0.0, 0.0, -1.0));
        if chk {
            Ok(())
        } else {
            loge!("test_chap_15_1", "p1:{}", t.p1);
            loge!("test_chap_15_1", "p2:{}", t.p2);
            loge!("test_chap_15_1", "p3:{}", t.p3);
            loge!("test_chap_15_1", "e1:{}", t.e1);
            loge!("test_chap_15_1", "e2:{}", t.e2);
            loge!("test_chap_15_1", "normal:{}", t.normal);
            Err("Constructing a triangle".into())
        }
    }

    /// Chap 15 - Finding the normal on a triangle
    #[test]
    fn test_chap_15_2() -> Result<(), String> {
        let t = Triangle::new(
            Tuple::point(0.0, 1.0, 0.0),
            Tuple::point(-1.0, 0.0, 0.0),
            Tuple::point(1.0, 0.0, 0.0),
        );

        let n1 = t.local_normal_at_no_hit(Tuple::point(0.0, 0.5, 0.0));
        let n2 = t.local_normal_at_no_hit(Tuple::point(-0.5, 0.75, 0.0));
        let n3 = t.local_normal_at_no_hit(Tuple::point(0.5, 0.25, 0.0));
        let chk = n1.approx_eq(t.normal) && n2.approx_eq(t.normal) && n3.approx_eq(t.normal);
        if chk {
            Ok(())
        } else {
            Err("Finding the normal on a triangle".into())
        }
    }

    /// Chap 15 - Intersecting a ray parallel to the triangle
    #[test]
    fn test_chap_15_3() -> Result<(), String> {
        let t = Triangle::new(
            Tuple::point(0.0, 1.0, 0.0),
            Tuple::point(-1.0, 0.0, 0.0),
            Tuple::point(1.0, 0.0, 0.0),
        );

        let r = Ray::new(Tuple::point(0.0, -1.0, -2.0), Tuple::vector(0.0, 1.0, 0.0));

        let xs = t.local_intersect(&r);

        let chk = xs.is_empty();
        if chk {
            Ok(())
        } else {
            Err("Intersecting a ray parallel to the triangle".into())
        }
    }

    /// Chap 15 - A ray misses the p1-p3 edge
    #[test]
    fn test_chap_15_4() -> Result<(), String> {
        let t = Triangle::new(
            Tuple::point(0.0, 1.0, 0.0),
            Tuple::point(-1.0, 0.0, 0.0),
            Tuple::point(1.0, 0.0, 0.0),
        );

        let r = Ray::new(Tuple::point(0.0, -1.0, -2.0), Tuple::vector(0.0, 1.0, 0.0));

        let xs = t.local_intersect(&r);

        let chk = xs.is_empty();
        if chk {
            Ok(())
        } else {
            Err("A ray misses the p1-p3 edge".into())
        }
    }

    /// Chap 15 - A ray misses the p1-p2 edge
    #[test]
    fn test_chap_15_5() -> Result<(), String> {
        let t = Triangle::new(
            Tuple::point(0.0, 1.0, 0.0),
            Tuple::point(-1.0, 0.0, 0.0),
            Tuple::point(1.0, 0.0, 0.0),
        );

        let r = Ray::new(Tuple::point(-1.0, 1.0, -2.0), Tuple::vector(0.0, 0.0, 1.0));

        let xs = t.local_intersect(&r);

        let chk = xs.is_empty();
        if chk {
            Ok(())
        } else {
            Err("A ray misses the p1-p2 edge".into())
        }
    }

    /// Chap 15 - A ray misses the p2-p3 edge
    #[test]
    fn test_chap_15_6() -> Result<(), String> {
        let t = Triangle::new(
            Tuple::point(0.0, 1.0, 0.0),
            Tuple::point(-1.0, 0.0, 0.0),
            Tuple::point(1.0, 0.0, 0.0),
        );

        let r = Ray::new(Tuple::point(0.0, -1.0, -2.0), Tuple::vector(0.0, 0.0, 1.0));

        let xs = t.local_intersect(&r);

        let chk = xs.is_empty();
        if chk {
            Ok(())
        } else {
            Err("A ray misses the p2-p3 edge".into())
        }
    }

    /// Chap x - A ray strikes a triangle
    #[test]
    fn test_chap_15_7() -> Result<(), String> {
        let t = Triangle::new(
            Tuple::point(0.0, 1.0, 0.0),
            Tuple::point(-1.0, 0.0, 0.0),
            Tuple::point(1.0, 0.0, 0.0),
        );

        let r = Ray::new(Tuple::point(0.0, 0.5, -2.0), Tuple::vector(0.0, 0.0, 1.0));

        let xs = t.local_intersect(&r);

        let chk = xs.count() == 1;
        if chk {
            Ok(())
        } else {
            Err("A ray strikes a triangle".into())
        }
    }

    /// Chap 15 - A triangle has a bounding box
    #[test]
    fn test_chap_15_bounds() -> Result<(), String> {
        let t = Triangle::new(
            Tuple::point(-3.0, 7.0, 2.0),
            Tuple::point(6.0, 2.0, -4.0),
            Tuple::point(2.0, -1.0, -1.0),
        );
        let bb = t.bounds();
        let chk = bb.min.approx_eq(Tuple::point(-3.0, -1.0, -4.0))
            && bb.max.approx_eq(Tuple::point(6.0, 7.0, 2.0));
        if chk {
            Ok(())
        } else {
            Err("A triangle has a bounding box".into())
        }
    }
}
