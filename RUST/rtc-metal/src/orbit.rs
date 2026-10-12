//! Orbit camera maths: turns mouse drag/scroll into a view transform.

use rtc_rust::matrix::Matrix4;
use rtc_rust::tuple::Tuple;
use rtc_rust::world::view_transform;

/// Camera orbiting `to`. Assumes `up = (0, 1, 0)` (true of every scene in the
/// repo): yaw rotates around world +y, pitch is the elevation above the
/// xz-plane through `to`.
pub struct Orbit {
    to: Tuple,
    up: Tuple,
    radius: f64,
    yaw: f64,
    pitch: f64,
}

impl Orbit {
    pub fn new(from: Tuple, to: Tuple, up: Tuple) -> Self {
        let d = from - to;
        let radius = d.magnitude();
        Orbit {
            to,
            up,
            radius,
            yaw: d.x.atan2(-d.z),
            pitch: (d.y / radius).asin(),
        }
    }

    pub fn from(&self) -> Tuple {
        let (sy, cy) = self.yaw.sin_cos();
        let (sp, cp) = self.pitch.sin_cos();
        self.to + Tuple::vector(cp * sy, sp, -cp * cy) * self.radius
    }

    /// `dx`, `dy` in points; 0.01 rad per point. Dragging right moves the
    /// camera right; dragging down raises it.
    pub fn drag(&mut self, dx: f64, dy: f64) {
        let limit = 89f64.to_radians();
        self.yaw += dx * 0.01;
        self.pitch = (self.pitch + dy * 0.01).clamp(-limit, limit);
    }

    pub fn zoom(&mut self, dy: f64) {
        self.radius = (self.radius * (-0.1 * dy).exp()).max(1e-3);
    }

    pub fn transform(&self) -> Matrix4 {
        view_transform(self.from(), self.to, self.up)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn p(x: f64, y: f64, z: f64) -> Tuple {
        Tuple::point(x, y, z)
    }
    fn up() -> Tuple {
        Tuple::vector(0.0, 1.0, 0.0)
    }

    #[test]
    fn new_round_trips_from() {
        let (from, to) = (p(-3.0, 1.0, 2.5), p(0.0, 0.5, 0.0));
        let o = Orbit::new(from, to, up());
        assert!(o.from().approx_eq(from));
        assert!(o.transform().approx_eq(view_transform(from, to, up())));
    }

    #[test]
    fn drag_zero_is_identity() {
        let from = p(-3.0, 1.0, 2.5);
        let mut o = Orbit::new(from, p(0.0, 0.5, 0.0), up());
        o.drag(0.0, 0.0);
        assert!(o.from().approx_eq(from));
    }

    #[test]
    fn pitch_is_clamped() {
        let to = p(0.0, 0.5, 0.0);
        let mut o = Orbit::new(p(-3.0, 1.0, 2.5), to, up());
        for dy in [1e6, -1e6] {
            o.drag(0.0, dy);
            let d = o.from() - to;
            let elevation = (d.y / d.magnitude()).asin();
            assert!(elevation.abs() <= 89f64.to_radians() + 1e-9);
        }
    }

    #[test]
    fn zoom_keeps_positive_radius() {
        let to = p(0.0, 0.5, 0.0);
        let mut o = Orbit::new(p(-3.0, 1.0, 2.5), to, up());
        o.zoom(1e6);
        assert!((o.from() - to).magnitude() >= 1e-3 - 1e-12); // float slack on the unit vector
    }

    #[test]
    fn yaw_quarter_turn() {
        let mut o = Orbit::new(p(0.0, 0.0, -5.0), p(0.0, 0.0, 0.0), up());
        // Dragging right moves the camera right (+x as seen from -z looking +z).
        o.drag(PI / 2.0 / 0.01, 0.0);
        assert!(o.from().approx_eq(p(5.0, 0.0, 0.0)));
    }
}
