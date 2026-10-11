//! Torus (book chapter 17, "Next Steps")
//!
//! A ring around the y axis: `major_radius` from the origin to the centre of
//! the tube, `minor_radius` for the tube itself. Implicit surface
//! `(|p|² − R² − r²)² = 4R²(r² − y²)`, so a ray meets it where a quartic in
//! `t` is zero — up to four hits.

use crate::bounds::BoundingBox;
use crate::intersection::{Intersection, Intersections};
use crate::math::EPSILON;
use crate::ray::Ray;
use crate::shape::{Shape, ShapeData};
use crate::tuple::Tuple;

#[derive(Debug, Clone)]
pub struct Torus {
    pub data: ShapeData,
    /// Distance from the axis to the centre of the tube.
    pub major_radius: f64,
    /// Radius of the tube.
    pub minor_radius: f64,
}

impl Torus {
    /// Default proportions: ring radius 1, tube radius 0.25.
    pub fn new() -> Self {
        Self::with_radii(1.0, 0.25)
    }

    pub fn with_radii(major_radius: f64, minor_radius: f64) -> Self {
        Self {
            data: ShapeData::new(),
            major_radius,
            minor_radius,
        }
    }

    /// Real roots of the ray/torus quartic, ascending.
    fn roots(&self, ray: &Ray) -> Vec<f64> {
        let o = ray.origin;
        let d = ray.direction;
        let (rr, r) = (self.major_radius, self.minor_radius);
        let dd = d.x * d.x + d.y * d.y + d.z * d.z;
        let od = o.x * d.x + o.y * d.y + o.z * d.z;
        let e = o.x * o.x + o.y * o.y + o.z * o.z - rr * rr - r * r;
        let four_rr2 = 4.0 * rr * rr;

        solve_quartic(
            dd * dd,
            4.0 * dd * od,
            4.0 * od * od + 2.0 * dd * e + four_rr2 * d.y * d.y,
            4.0 * od * e + 2.0 * four_rr2 * o.y * d.y,
            e * e - four_rr2 * (r * r - o.y * o.y),
        )
    }
}

impl Default for Torus {
    fn default() -> Self {
        Self::new()
    }
}

impl Shape for Torus {
    fn bounds(&self) -> BoundingBox {
        let outer = self.major_radius + self.minor_radius;
        BoundingBox::new(
            Tuple::point(-outer, -self.minor_radius, -outer),
            Tuple::point(outer, self.minor_radius, outer),
        )
    }

    fn data(&self) -> &ShapeData {
        &self.data
    }

    fn data_mut(&mut self) -> &mut ShapeData {
        &mut self.data
    }

    fn local_intersect(&self, ray: &Ray) -> Intersections {
        let mut xs = Intersections::new();
        for t in self.roots(ray) {
            xs.push(Intersection::new(t, self.id()));
        }
        xs
    }

    /// The normal points from the nearest point on the central ring to
    /// `point` (inward on the inner side of the tube, outward on the outer).
    fn local_normal_at(&self, point: Tuple, _hit: Intersection) -> Tuple {
        let axis_dist = (point.x * point.x + point.z * point.z).sqrt();
        let ring = if axis_dist < EPSILON {
            Tuple::point(self.major_radius, 0.0, 0.0)
        } else {
            Tuple::point(
                point.x / axis_dist * self.major_radius,
                0.0,
                point.z / axis_dist * self.major_radius,
            )
        };
        (point - ring).normalize()
    }
}

/// Real roots of `a t⁴ + b t³ + c t² + d t + e = 0`, ascending and
/// de-duplicated. Ferrari's method through the resolvent cubic, then a few
/// Newton steps on the original quartic to clean up the rounding.
pub fn solve_quartic(a: f64, b: f64, c: f64, d: f64, e: f64) -> Vec<f64> {
    if a.abs() < 1e-12 {
        return solve_cubic(b, c, d, e);
    }
    // Normalise to t⁴ + b t³ + c t² + d t + e.
    let (b, c, d, e) = (b / a, c / a, d / a, e / a);
    // Depress: t = y − b/4  →  y⁴ + p y² + q y + r.
    let b2 = b * b;
    let p = c - 3.0 * b2 / 8.0;
    let q = d - b * c / 2.0 + b2 * b / 8.0;
    let r = e - b * d / 4.0 + b2 * c / 16.0 - 3.0 * b2 * b2 / 256.0;
    let shift = -b / 4.0;

    let mut roots: Vec<f64> = Vec::with_capacity(4);
    if q.abs() < 1e-12 {
        // Biquadratic: y² = z.
        for z in solve_quadratic(1.0, p, r) {
            if z >= 0.0 {
                let s = z.sqrt();
                roots.push(shift + s);
                roots.push(shift - s);
            }
        }
    } else {
        // Resolvent cubic: m³ + p m² + (p²/4 − r) m − q²/8 = 0, take m > 0.
        let m = solve_cubic(1.0, p, p * p / 4.0 - r, -q * q / 8.0)
            .into_iter()
            .filter(|&m| m > 0.0)
            .fold(f64::NAN, f64::max);
        if m.is_nan() {
            return roots;
        }
        let sqrt_2m = (2.0 * m).sqrt();
        // y² ± sqrt(2m) y + (p/2 + m ∓ q/(2 sqrt(2m))) = 0
        for sign in [1.0, -1.0] {
            let k = p / 2.0 + m - sign * q / (2.0 * sqrt_2m);
            for y in solve_quadratic(1.0, sign * sqrt_2m, k) {
                roots.push(shift + y);
            }
        }
    }

    // Polish against the original (normalised) quartic.
    let f = |t: f64| (((t + b) * t + c) * t + d) * t + e;
    let df = |t: f64| ((4.0 * t + 3.0 * b) * t + 2.0 * c) * t + d;
    for t in roots.iter_mut() {
        for _ in 0..3 {
            let slope = df(*t);
            if slope.abs() < 1e-12 {
                break;
            }
            *t -= f(*t) / slope;
        }
    }

    roots.sort_by(|x, y| x.partial_cmp(y).unwrap());
    roots.dedup_by(|x, y| (*x - *y).abs() < 1e-9);
    roots
}

/// Real roots of `a t² + b t + c = 0`.
pub fn solve_quadratic(a: f64, b: f64, c: f64) -> Vec<f64> {
    if a.abs() < 1e-12 {
        return if b.abs() < 1e-12 {
            vec![]
        } else {
            vec![-c / b]
        };
    }
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return vec![];
    }
    let s = disc.sqrt();
    // Numerically stable form: avoid cancellation in −b ± s.
    let q = -0.5 * (b + b.signum() * s);
    if q == 0.0 {
        return vec![0.0];
    }
    vec![q / a, c / q]
}

/// Real roots of `a t³ + b t² + c t + d = 0` (Cardano / trigonometric).
pub fn solve_cubic(a: f64, b: f64, c: f64, d: f64) -> Vec<f64> {
    if a.abs() < 1e-12 {
        return solve_quadratic(b, c, d);
    }
    let (b, c, d) = (b / a, c / a, d / a);
    // Depress: t = y − b/3  →  y³ + p y + q.
    let p = c - b * b / 3.0;
    let q = 2.0 * b * b * b / 27.0 - b * c / 3.0 + d;
    let shift = -b / 3.0;
    let disc = q * q / 4.0 + p * p * p / 27.0;

    if disc > 0.0 {
        let s = disc.sqrt();
        let y = (-q / 2.0 + s).cbrt() + (-q / 2.0 - s).cbrt();
        vec![shift + y]
    } else if disc == 0.0 {
        let u = (-q / 2.0).cbrt();
        vec![shift + 2.0 * u, shift - u]
    } else {
        // Three real roots.
        let rho = (-p * p * p / 27.0).sqrt();
        let theta = (-q / (2.0 * rho)).clamp(-1.0, 1.0).acos();
        let m = 2.0 * rho.cbrt();
        (0..3)
            .map(|k| shift + m * ((theta + 2.0 * std::f64::consts::PI * k as f64) / 3.0).cos())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::loge;
    use crate::math::approx_eq;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-6
    }

    /// Chap 17 - The quartic solver finds all real roots
    #[test]
    fn test_chap_17_21() -> Result<(), String> {
        // (t−1)(t−2)(t−3)(t−4) = t⁴ − 10t³ + 35t² − 50t + 24
        let r = solve_quartic(1.0, -10.0, 35.0, -50.0, 24.0);
        let four = r.len() == 4
            && r.iter()
                .zip([1.0, 2.0, 3.0, 4.0])
                .all(|(a, b)| close(*a, b));
        // (t² + 1)(t − 1)(t + 2) = t⁴ + t³ − t² + t − 2
        let r2 = solve_quartic(1.0, 1.0, -1.0, 1.0, -2.0);
        let two = r2.len() == 2 && close(r2[0], -2.0) && close(r2[1], 1.0);
        // (t² + 1)(t² + 4): no real roots
        let none = solve_quartic(1.0, 0.0, 5.0, 0.0, 4.0).is_empty();
        // t⁴ − 5t² + 4 = (t²−1)(t²−4): biquadratic
        let r3 = solve_quartic(1.0, 0.0, -5.0, 0.0, 4.0);
        let bi = r3.len() == 4
            && r3
                .iter()
                .zip([-2.0, -1.0, 1.0, 2.0])
                .all(|(a, b)| close(*a, b));
        if four && two && none && bi {
            Ok(())
        } else {
            loge!("test_chap_17_21", "r:{:?} r2:{:?} r3:{:?}", r, r2, r3);
            Err("The quartic solver finds all real roots".into())
        }
    }

    /// Chap 17 - A ray through the tube of a torus hits it four times
    #[test]
    fn test_chap_17_22() -> Result<(), String> {
        let t = Torus::with_radii(2.0, 0.5);
        let r = Ray::new(Tuple::point(-5.0, 0.0, 0.0), Tuple::vector(1.0, 0.0, 0.0));
        let xs = t.local_intersect(&r);
        let expected = [2.5, 3.5, 6.5, 7.5];
        let chk = xs.count() == 4 && (0..4).all(|i| close(xs[i].t, expected[i]));
        if chk {
            Ok(())
        } else {
            let got: Vec<f64> = xs.iter().map(|i| i.t).collect();
            loge!("test_chap_17_22", "got:{:?}", got);
            Err("A ray through the tube of a torus hits it four times".into())
        }
    }

    /// Chap 17 - A ray down the axis passes through the hole; a ray beside
    /// the torus misses it
    #[test]
    fn test_chap_17_23() -> Result<(), String> {
        let t = Torus::with_radii(2.0, 0.5);
        let hole = Ray::new(Tuple::point(0.0, 5.0, 0.0), Tuple::vector(0.0, -1.0, 0.0));
        let beside = Ray::new(Tuple::point(-5.0, 0.0, 3.0), Tuple::vector(1.0, 0.0, 0.0));
        let above = Ray::new(Tuple::point(-5.0, 0.6, 0.0), Tuple::vector(1.0, 0.0, 0.0));
        let chk = t.local_intersect(&hole).is_empty()
            && t.local_intersect(&beside).is_empty()
            && t.local_intersect(&above).is_empty();
        if chk {
            Ok(())
        } else {
            Err("Rays through the hole or beside the torus must miss".into())
        }
    }

    /// Chap 17 - A ray starting inside the tube has two hits ahead of it
    #[test]
    fn test_chap_17_24() -> Result<(), String> {
        let t = Torus::with_radii(2.0, 0.5);
        let r = Ray::new(Tuple::point(2.0, 0.0, 0.0), Tuple::vector(1.0, 0.0, 0.0));
        let xs = t.local_intersect(&r);
        let hit = xs.hit().ok_or("no hit")?;
        if xs.count() == 4 && close(hit.t, 0.5) {
            Ok(())
        } else {
            let got: Vec<f64> = xs.iter().map(|i| i.t).collect();
            loge!("test_chap_17_24", "got:{:?}", got);
            Err("A ray starting inside the tube".into())
        }
    }

    /// Chap 17 - The normal on a torus
    #[test]
    fn test_chap_17_25() -> Result<(), String> {
        let t = Torus::with_radii(2.0, 0.5);
        let cases = [
            (Tuple::point(2.5, 0.0, 0.0), Tuple::vector(1.0, 0.0, 0.0)),
            (Tuple::point(2.0, 0.5, 0.0), Tuple::vector(0.0, 1.0, 0.0)),
            (Tuple::point(0.0, 0.0, -1.5), Tuple::vector(0.0, 0.0, 1.0)),
            (Tuple::point(0.0, -0.5, 2.0), Tuple::vector(0.0, -1.0, 0.0)),
        ];
        for (p, expected) in cases {
            let n = t.local_normal_at_no_hit(p);
            if !n.approx_eq(expected) {
                loge!("test_chap_17_25", "p:{} got:{}", p, n);
                return Err("The normal on a torus".into());
            }
        }
        Ok(())
    }

    /// Chap 17 - A torus has a bounding box
    #[test]
    fn test_chap_17_26() -> Result<(), String> {
        let t = Torus::with_radii(2.0, 0.5);
        let bb = t.bounds();
        let chk = bb.min.approx_eq(Tuple::point(-2.5, -0.5, -2.5))
            && bb.max.approx_eq(Tuple::point(2.5, 0.5, 2.5))
            && approx_eq(Torus::new().major_radius, 1.0);
        if chk {
            Ok(())
        } else {
            Err("A torus has a bounding box".into())
        }
    }
}
