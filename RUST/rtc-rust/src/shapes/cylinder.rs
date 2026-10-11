//! Cylinder definition
//!
//! A cylinder of radius 1 centered on the y axis. By default it is infinite
//! (`minimum` = -∞, `maximum` = +∞) and open (`closed` = false). Set the
//! `minimum`/`maximum` bounds to truncate it along y, and `closed` to seal the
//! ends with flat caps.

use crate::bounds::BoundingBox;
use crate::intersection::{Intersection, Intersections};
use crate::math::{EPSILON, approx_eq};
use crate::ray::Ray;
use crate::shape::{Shape, ShapeData};
use crate::tuple::Tuple;

/// Test whether the ray at parameter `t` lands inside the cylinder's radius.
///
/// Used when checking end-cap hits: after intersecting the ray with a cap
/// plane, this confirms the point actually falls on the radius-1 cap disc.
/// Returns `true` when the point at `t` satisfies `x² + z² <= 1`.
pub fn check_cap(ray: &Ray, t: f64) -> bool {
    let x = ray.origin.x + t * ray.direction.x;
    let z = ray.origin.z + t * ray.direction.z;

    x * x + z * z <= 1.0
}

/// A unit-radius cylinder aligned with the y axis.
///
/// Truncated along y by `minimum`/`maximum` and optionally sealed with flat end
/// caps via `closed`. See the module docs for the defaults.
#[derive(Debug, Clone)]
pub struct Cylinder {
    pub data: ShapeData,
    /// Lower y bound; points at or below are not part of the wall.
    pub minimum: f64,
    /// Upper y bound; points at or above are not part of the wall.
    pub maximum: f64,
    /// Whether the cylinder has flat end caps at `minimum` and `maximum`.
    pub closed: bool,
}

impl Cylinder {
    /// Add any intersections of `ray` with the cylinder's end caps to `xs`.
    ///
    /// Does nothing for an open cylinder, or when the ray runs parallel to the
    /// caps (`direction.y ≈ 0`). Otherwise it intersects the ray with the cap
    /// planes at `y = minimum` and `y = maximum`, keeping only the hits that
    /// fall inside the radius-1 disc (see [`check_cap`]). Hits are appended to
    /// `xs`, which `Intersections::push` keeps sorted by `t`.
    pub fn intersect_caps(&self, ray: &Ray, xs: &mut Intersections) {
        // caps only matter if the cylinder is closed, and might possibly be intersected by the ray.
        if !self.closed || approx_eq(ray.direction.y, 0.0) {
            return;
        }

        // check for an intersection with the lower end cap by intersecting
        // the ray with the plane at y=cyl.minimum
        let t = (self.minimum - ray.origin.y) / ray.direction.y;
        if check_cap(ray, t) {
            xs.push(Intersection::new(t, self.id()));
        }

        // check for an intersection with the upper end cap by intersecting
        // the ray with the plane at y=cyl.maximum
        let t = (self.maximum - ray.origin.y) / ray.direction.y;
        if check_cap(ray, t) {
            xs.push(Intersection::new(t, self.id()));
        }
    }

    /// Create an infinite, open cylinder of radius 1 along the y axis.
    pub fn new() -> Self {
        Self {
            data: ShapeData::new(),
            minimum: f64::NEG_INFINITY,
            maximum: f64::INFINITY,
            closed: false,
        }
    }
}

impl Default for Cylinder {
    fn default() -> Self {
        Self::new()
    }
}

impl Shape for Cylinder {
    fn bounds(&self) -> BoundingBox {
        BoundingBox::new(
            Tuple::point(-1.0, self.minimum, -1.0),
            Tuple::point(1.0, self.maximum, 1.0),
        )
    }

    fn data(&self) -> &ShapeData {
        &self.data
    }

    fn data_mut(&mut self) -> &mut ShapeData {
        &mut self.data
    }

    /// Return the outward surface normal at `point`, in object space.
    ///
    /// A point within radius 1 and at (within `EPSILON` of) `maximum` is on the
    /// top cap, so the normal is `+y`; at `minimum` it is the bottom cap, `-y`.
    /// Everywhere else the point is on the wall, where the normal points
    /// radially outward with no y component.
    fn local_normal_at(&self, point: Tuple, _hit: Intersection) -> Tuple {
        // squared distance from the y axis
        let dist = point.x * point.x + point.z * point.z;

        if dist < 1.0 && point.y >= self.maximum - EPSILON {
            return Tuple::vector(0.0, 1.0, 0.0);
        }

        if dist < 1.0 && point.y <= self.minimum + EPSILON {
            return Tuple::vector(0.0, -1.0, 0.0);
        }

        Tuple::vector(point.x, 0.0, point.z)
    }

    /// Intersect `ray` with the cylinder, in object space.
    ///
    /// Solves the quadratic for the ray against the infinite radius-1 wall,
    /// keeps only the roots whose y lies within `minimum..maximum`, and adds any
    /// end-cap hits via [`Cylinder::intersect_caps`]. Returns 0–4 intersections
    /// (walls plus caps), in a collection kept sorted by `t`.
    fn local_intersect(&self, ray: &Ray) -> Intersections {
        let a = ray.direction.x * ray.direction.x + ray.direction.z * ray.direction.z;

        let mut xs = Intersections::new();

        self.intersect_caps(ray, &mut xs);

        // ray is parallel to the y axis
        if approx_eq(a, 0.0) {
            return xs;
        }

        let b = 2.0 * ray.origin.x * ray.direction.x + 2.0 * ray.origin.z * ray.direction.z;
        let c = ray.origin.x * ray.origin.x + ray.origin.z * ray.origin.z - 1.0;

        let disc = b * b - 4.0 * a * c;

        // ray does not intersect the cylinder
        if disc < 0.0 {
            return xs;
        }

        // Compute once use twice ;-)
        let disc_sqrt = disc.sqrt();
        let inv_a_mul_2 = 1.0 / (2.0 * a);

        let (t0, t1) = (
            (-b - disc_sqrt) * inv_a_mul_2,
            (-b + disc_sqrt) * inv_a_mul_2,
        );

        // No swap needed: a = dx² + dz² > 0 here (the a ≈ 0 case already
        // returned), so dividing by 2a preserves order and t0 <= t1 always.

        let y0 = ray.origin.y + t0 * ray.direction.y;
        if self.minimum < y0 && y0 < self.maximum {
            xs.push(Intersection::new(t0, self.id()));
        }

        let y1 = ray.origin.y + t1 * ray.direction.y;
        if self.minimum < y1 && y1 < self.maximum {
            xs.push(Intersection::new(t1, self.id()));
        }

        xs
    }
}
