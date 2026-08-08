//! bounds.rs
//! Implement axis aligned bounding box - AABB for the different shapes.

use crate::matrix::Matrix4;
use crate::ray::Ray;
use crate::tuple::Tuple;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundingBox {
    pub min: Tuple, // a point
    pub max: Tuple, // a point
}

fn check_axis(origin: f64, direction: f64, min: f64, max: f64) -> (f64, f64) {
    let tmin = (min - origin) / direction; // was hardcoded -1.0
    let tmax = (max - origin) / direction; // was hardcoded  1.0
    if tmin > tmax {
        (tmax, tmin)
    } else {
        (tmin, tmax)
    }
}

impl BoundingBox {
    pub fn add_point(&mut self, p: Tuple) {
        self.min.x = self.min.x.min(p.x);
        self.min.y = self.min.y.min(p.y);
        self.min.z = self.min.z.min(p.z);

        self.max.x = self.max.x.max(p.x);
        self.max.y = self.max.y.max(p.y);
        self.max.z = self.max.z.max(p.z);
    }

    pub fn add_box(&mut self, other: &BoundingBox) {
        self.add_point(other.min);
        self.add_point(other.max);
    }

    pub fn contains_point(&self, p: Tuple) -> bool {
        if self.min.x > p.x {
            return false;
        }
        if self.min.y > p.y {
            return false;
        }
        if self.min.z > p.z {
            return false;
        }

        if self.max.x < p.x {
            return false;
        }
        if self.max.y < p.y {
            return false;
        }
        if self.max.z < p.z {
            return false;
        }

        true
    }

    pub fn contains_box(&self, other: &BoundingBox) -> bool {
        self.contains_point(other.min) && self.contains_point(other.max)
    }

    pub fn empty() -> Self {
        Self {
            min: Tuple::point(f64::INFINITY, f64::INFINITY, f64::INFINITY),
            max: Tuple::point(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
        }
    }

    pub fn intersects(&self, ray: &Ray) -> bool {
        let (xtmin, xtmax) = check_axis(ray.origin.x, ray.direction.x, self.min.x, self.max.x);
        let (ytmin, ytmax) = check_axis(ray.origin.y, ray.direction.y, self.min.y, self.max.y);
        let (ztmin, ztmax) = check_axis(ray.origin.z, ray.direction.z, self.min.z, self.max.z);
        let tmin = xtmin.max(ytmin).max(ztmin);
        let tmax = xtmax.min(ytmax).min(ztmax);
        tmin <= tmax
    }

    pub fn new(min: Tuple, max: Tuple) -> Self {
        Self { min, max }
    }

    pub fn transform(&self, m: Matrix4) -> BoundingBox {
        // the 8 corners from all min/max combinations
        let corners = [
            (self.min.x, self.min.y, self.min.z),
            (self.min.x, self.min.y, self.max.z),
            (self.min.x, self.max.y, self.min.z),
            (self.min.x, self.max.y, self.max.z),
            (self.max.x, self.min.y, self.min.z),
            (self.max.x, self.min.y, self.max.z),
            (self.max.x, self.max.y, self.min.z),
            (self.max.x, self.max.y, self.max.z),
        ];
        let mut out = BoundingBox::empty();
        for c in corners {
            out.add_point(m * Tuple::point(c.0, c.1, c.2));
        }
        out
    }
}
