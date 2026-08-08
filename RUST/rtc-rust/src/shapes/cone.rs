//! Cone definition
//!
//! A double-napped cone about the y axis (`x² + z² = y²`), so its radius at
//! height `y` equals `|y|`. Infinite and open by default; `minimum`/`maximum`
//! truncate it along y and `closed` adds flat end caps.

use crate::bounds::BoundingBox;
use crate::intersection::{Intersection, Intersections};
use crate::math::{EPSILON, approx_eq};
use crate::ray::Ray;
use crate::shape::{Shape, ShapeData};
use crate::tuple::Tuple;

/// Whether the ray at parameter `t` lands within the cap disc at height `y`.
///
/// A cone's radius at height `y` is `|y|`, so the point is inside the cap when
/// `x² + z² <= y²`. Used to validate end-cap hits.
pub fn check_cap(ray: &Ray, t: f64, y: f64) -> bool {
    let x = ray.origin.x + t * ray.direction.x;
    let z = ray.origin.z + t * ray.direction.z;

    x * x + z * z <= y * y
}

/// A double-napped cone about the y axis, truncated by `minimum`/`maximum` and
/// optionally sealed with flat end caps via `closed`.
#[derive(Debug, Clone)]
pub struct Cone {
    pub data: ShapeData,
    /// Lower y bound; points at or below are not part of the wall.
    pub minimum: f64,
    /// Upper y bound; points at or above are not part of the wall.
    pub maximum: f64,
    /// Whether the cone has flat end caps at `minimum` and `maximum`.
    pub closed: bool,
}

impl Cone {
    /// Append any intersections of `ray` with the cone's end caps to `xs`.
    ///
    /// No-op for an open cone or a ray parallel to the caps (`direction.y ≈ 0`).
    /// Otherwise intersects the cap planes at `y = minimum` and `y = maximum`,
    /// keeping hits that fall within the cap disc (see [`check_cap`]). `xs` is
    /// kept sorted by `t`.
    pub fn intersect_caps(&self, ray: &Ray, xs: &mut Intersections) {
        // caps only matter for a closed cone the ray could actually reach
        if !self.closed || approx_eq(ray.direction.y, 0.0) {
            return;
        }

        // lower cap: intersect the ray with the plane y = minimum
        let t = (self.minimum - ray.origin.y) / ray.direction.y;
        if check_cap(ray, t, self.minimum) {
            xs.push(Intersection::new(t, self.id()));
        }

        // upper cap: intersect the ray with the plane y = maximum
        let t = (self.maximum - ray.origin.y) / ray.direction.y;
        if check_cap(ray, t, self.maximum) {
            xs.push(Intersection::new(t, self.id()));
        }
    }

    /// Create an infinite, open cone about the y axis.
    pub fn new() -> Self {
        Self {
            data: ShapeData::new(),
            minimum: f64::NEG_INFINITY,
            maximum: f64::INFINITY,
            closed: false,
        }
    }
}

impl Default for Cone {
    fn default() -> Self {
        Self::new()
    }
}

impl Shape for Cone {
    /// A cone's radius at height y is |y| (that's the x²+z² = y² surface), so the box's x/z extent
    /// is the largest radius over the y-range = max(|min|, |max|)
    fn bounds(&self) -> BoundingBox {
        let limit = self.minimum.abs().max(self.maximum.abs());
        BoundingBox::new(
            Tuple::point(-limit, self.minimum, -limit),
            Tuple::point(limit, self.maximum, limit),
        )
    }

    fn data(&self) -> &ShapeData {
        &self.data
    }

    fn data_mut(&mut self) -> &mut ShapeData {
        &mut self.data
    }

    /// Outward surface normal at `point`, in object space.
    ///
    /// On an end cap the normal is `+y` (top) or `-y` (bottom); on the wall it
    /// points radially outward with a y component set by the cone's slope.
    fn local_normal_at(&self, point: Tuple) -> Tuple {
        // squared distance from the y axis
        let dist = point.x * point.x + point.z * point.z;

        if dist < self.maximum * self.maximum && point.y >= self.maximum - EPSILON {
            return Tuple::vector(0.0, 1.0, 0.0);
        }

        if dist < self.minimum * self.minimum && point.y <= self.minimum + EPSILON {
            return Tuple::vector(0.0, -1.0, 0.0);
        }

        let y = (point.x * point.x + point.z * point.z).sqrt();
        let y = if point.y > 0.0 { -y } else { y };

        Tuple::vector(point.x, y, point.z)
    }

    /// Intersect `ray` with the cone, in object space.
    ///
    /// Solves the cone quadratic against the wall, keeping roots whose `y` lies
    /// in `minimum..maximum`, plus end-cap hits via [`Cone::intersect_caps`]. A
    /// ray parallel to a cone half (`a ≈ 0`) uses the single linear root.
    /// Returns 0–4 intersections, kept sorted by `t`.
    fn local_intersect(&self, ray: &Ray) -> Intersections {
        let a = ray.direction.x * ray.direction.x - ray.direction.y * ray.direction.y
            + ray.direction.z * ray.direction.z;

        let mut xs = Intersections::new();

        self.intersect_caps(ray, &mut xs);

        let b = 2.0 * ray.origin.x * ray.direction.x - 2.0 * ray.origin.y * ray.direction.y
            + 2.0 * ray.origin.z * ray.direction.z;

        let c =
            ray.origin.x * ray.origin.x - ray.origin.y * ray.origin.y + ray.origin.z * ray.origin.z;

        // ray is parallel to one of the cone’s halves
        if approx_eq(a, 0.0) {
            if approx_eq(b, 0.0) {
                return xs;
            }

            // single point of intersection
            let t = -c / (2.0 * b);
            let y = ray.origin.y + t * ray.direction.y;
            if self.minimum < y && y < self.maximum {
                xs.push(Intersection::new(t, self.id()));
            }
            return xs; // ← stop here; the quadratic below doesn't apply
        }

        let disc = b * b - 4.0 * a * c;

        // ray does not intersect the cone
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

        // a = dx² - dy² + dz² can be negative for a cone (unlike the cylinder), so
        // t0 <= t1 is NOT guaranteed. That's fine: each root is range-checked and
        // pushed independently, and Intersections::push keeps xs sorted by t.

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
