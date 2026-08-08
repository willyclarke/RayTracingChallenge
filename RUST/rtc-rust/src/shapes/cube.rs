//! Cube definition
//!
//! An axis-aligned cube centered at the origin, spanning `-1..=1` on each axis
//! (the "unit cube", analogous to the unit sphere). Ray intersection uses the
//! slab method — see [`check_axis`].

use crate::bounds::BoundingBox;
use crate::intersection::{Intersection, Intersections};
use crate::math::approx_eq;
use crate::ray::Ray;
use crate::shape::{Shape, ShapeData};
use crate::tuple::Tuple;

/// An axis-aligned unit cube centered at the origin, spanning `-1..=1` on each
/// axis. Its transform (in `data`) scales, rotates, and positions it in the world.
#[derive(Debug, Clone)]
pub struct Cube {
    pub data: ShapeData,
}

impl Cube {
    /// Create a unit cube centered at the origin (spanning `-1..=1` per axis).
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

/// Compute the pair of `t` values where a ray crosses the two parallel
/// planes that bound the unit cube along a single axis.
///
/// This is one axis of the *slab method* for ray/AABB intersection: the cube
/// spans `-1..=1` on every axis, so for a given axis the ray enters the slab at
/// one plane and exits at the other. Call it once per axis (x, y, z); the cube
/// is hit only where all three `t` intervals overlap.
///
/// # Arguments
/// * `origin` — the ray origin's component for this axis
/// * `direction` — the ray direction's component for this axis
///
/// # Returns
/// `(tmin, tmax)`, ordered so that `tmin <= tmax`. When `direction` is `0` the
/// ray is parallel to this axis and the values are `±infinity`, which the
/// caller's `min`/`max` combination handles correctly.
///
/// # Examples
/// ```
/// use rtc_rust::shapes::cube::check_axis;
///
/// // Ray at x = -5 heading +x: enters the x-slab at t=4, exits at t=6.
/// let (tmin, tmax) = check_axis(-5.0, 1.0);
/// assert_eq!((tmin, tmax), (4.0, 6.0));
/// ```
pub fn check_axis(origin: f64, direction: f64) -> (f64, f64) {
    let tmin = (-1.0 - origin) / direction;
    let tmax = (1.0 - origin) / direction;
    if tmin > tmax {
        (tmax, tmin)
    } else {
        (tmin, tmax)
    }
}

impl Shape for Cube {
    fn bounds(&self) -> BoundingBox {
        BoundingBox::new(Tuple::point(-1.0, -1.0, -1.0), Tuple::point(1.0, 1.0, 1.0))
    }

    fn data(&self) -> &ShapeData {
        &self.data
    }

    fn data_mut(&mut self) -> &mut ShapeData {
        &mut self.data
    }

    /// Return the outward surface normal at `point`, in the cube's object space.
    ///
    /// A point on the cube's surface lies on the face whose axis has the
    /// largest-magnitude coordinate (e.g. a point with `x = 1` and `|x|` greatest
    /// is on the `+x` face). The normal is the unit vector along that axis, taking
    /// its sign from the coordinate — so `x = 1` gives `+x` and `x = -1` gives `-x`.
    ///
    /// `point` is assumed to already be in object space; the [`Shape`] machinery
    /// transforms world-space points before calling this.
    fn local_normal_at(&self, point: Tuple) -> Tuple {
        let maxc = point.x.abs().max(point.y.abs().max(point.z.abs()));
        if approx_eq(maxc, point.x.abs()) {
            return Tuple::vector(point.x, 0.0, 0.0);
        }
        if approx_eq(maxc, point.y.abs()) {
            return Tuple::vector(0.0, point.y, 0.0);
        }
        Tuple::vector(0.0, 0.0, point.z)
    }

    /// Intersect a ray with the unit cube, in the cube's own object space.
    ///
    /// Uses the *slab method with early termination*: the cube is the region
    /// `-1..=1` on all three axes, so a hit requires the ray's entry/exit `t`
    /// intervals for x, y, and z to all overlap. The axes are combined
    /// incrementally and the search bails as soon as the running interval is
    /// empty (`tmin > tmax`), skipping any remaining axes on a miss.
    ///
    /// Returns two intersections — the entry (`tmin`) and exit (`tmax`) `t`
    /// values, in ascending order — or an empty [`Intersections`] when the ray
    /// misses the cube. Called by the [`Shape`] machinery after the ray has
    /// already been transformed into object space.
    fn local_intersect(&self, ray: &Ray) -> Intersections {
        let mut tmin = f64::NEG_INFINITY;
        let mut tmax = f64::INFINITY;

        let axes = [
            (ray.origin.x, ray.direction.x),
            (ray.origin.y, ray.direction.y),
            (ray.origin.z, ray.direction.z),
        ];

        for (origin, direction) in axes {
            let (axis_min, axis_max) = check_axis(origin, direction);
            tmin = tmin.max(axis_min);
            tmax = tmax.min(axis_max);
            if tmin > tmax {
                return Intersections::new(); // ray misses — skip remaining axes
            }
        }

        let mut xs = Intersections::new();
        xs.push(Intersection::new(tmin, self.id()));
        xs.push(Intersection::new(tmax, self.id()));
        xs
    }
}
