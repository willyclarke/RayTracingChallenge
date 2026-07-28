//! Cube definition
//!
//! A generic cube with sides 1x1x1
//!

use crate::intersection::{Intersection, Intersections};
use crate::math::approx_eq;
use crate::ray::Ray;
use crate::shape::{Shape, ShapeData};
use crate::tuple::Tuple;

// #[derive(Debug, Clone, Copy)]
#[derive(Debug, Clone)]
pub struct Cube {
    pub data: ShapeData,
}

impl Cube {
    pub fn new() -> Self {
        Self {
            data: ShapeData::new(),
        }
    }
}

impl Default for Cube {
    fn default() -> Self {
        Self::new()
    }
}

pub fn check_axis(origin: f64, direction: f64) -> (f64, f64) {
    let tmin_numerator = -1.0 - origin; // cube spans -1..1 on each axis
    let tmax_numerator = 1.0 - origin;

    let (tmin, tmax) = if direction.abs() >= crate::math::EPSILON {
        (tmin_numerator / direction, tmax_numerator / direction)
    } else {
        // direction ~0: avoid 0/0; keep the sign via ±infinity
        (
            tmin_numerator * f64::INFINITY,
            tmax_numerator * f64::INFINITY,
        )
    };

    if tmin > tmax {
        (tmax, tmin)
    } else {
        (tmin, tmax)
    } // ensure min ≤ max
}

impl Shape for Cube {
    fn data(&self) -> &ShapeData {
        &self.data
    }

    fn data_mut(&mut self) -> &mut ShapeData {
        &mut self.data
    }

    fn local_normal_at(&self, _point: Tuple) -> Tuple {
        let maxc = _point.x.abs().max(_point.y.abs().max(_point.z.abs()));
        if approx_eq(maxc, _point.x.abs()) {
            return Tuple::vector(_point.x, 0.0, 0.0);
        }
        if approx_eq(maxc, _point.y.abs()) {
            return Tuple::vector(0.0, _point.y, 0.0);
        }
        Tuple::vector(0.0, 0.0, _point.z)
    }

    fn local_intersect(&self, ray: &Ray) -> Intersections {
        let (xtmin, xtmax) = check_axis(ray.origin.x, ray.direction.x);
        let (ytmin, ytmax) = check_axis(ray.origin.y, ray.direction.y);
        let (ztmin, ztmax) = check_axis(ray.origin.z, ray.direction.z);

        let tmin = xtmin.max(ytmin).max(ztmin); // latest entry
        let tmax = xtmax.min(ytmax).min(ztmax); // earliest exit

        let mut xs = Intersections::new();
        if tmin > tmax {
            return xs; // ray misses the cube
        }
        xs.push(Intersection::new(tmin, self.id()));
        xs.push(Intersection::new(tmax, self.id()));
        xs
    }
}
