//! Orbit camera maths: turns mouse drag/scroll into a view transform.

use rtc_rust::matrix::Matrix4;
use rtc_rust::tuple::Tuple;
use rtc_rust::world::view_transform;

/// Camera orbiting `to` around world +y: yaw turns around the y axis through
/// `to`, pitch is the elevation above the xz-plane through `to`. The scene's
/// `up` is only passed on to `view_transform`, so a tilted `up` (e.g.
/// rtc-rust/scenes/cover.json's (-0.45, 1, 0)) is kept as a camera roll.
pub struct Orbit {
    to: Tuple,
    up: Tuple,
    radius: f64,
    /// Zooming out stops at 1e3 x the starting radius.
    max_radius: f64,
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
            // `max` keeps the clamp in `zoom` well-formed for tiny radii.
            max_radius: (1e3 * radius).max(1e-3),
            yaw: d.x.atan2(-d.z),
            pitch: (d.y / radius).asin(),
        }
    }

    pub fn from(&self) -> Tuple {
        let (sy, cy) = self.yaw.sin_cos();
        let (sp, cp) = self.pitch.sin_cos();
        self.to + Tuple::vector(cp * sy, sp, -cp * cy) * self.radius
    }

    /// Distance from the camera to the point it orbits.
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// `dx`, `dy` in points with +y = down (AppKit's mouse deltas); 0.01 rad
    /// per point. `Orbit` owns the sign convention, "grab the scene" (the
    /// surface follows the cursor, as in Blender or three.js OrbitControls):
    /// dragging right orbits the camera left, dragging down raises it.
    pub fn drag(&mut self, dx: f64, dy: f64) {
        let limit = 89f64.to_radians();
        self.yaw -= dx * 0.01;
        self.pitch = (self.pitch + dy * 0.01).clamp(-limit, limit);
    }

    /// Positive `dy` moves in: each unit scales the radius by exp(-0.1). The
    /// radius stays within [1e-3, 1e3 x the starting radius].
    pub fn zoom(&mut self, dy: f64) {
        self.radius = (self.radius * (-0.1 * dy).exp()).clamp(1e-3, self.max_radius);
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
        // Dragging right grabs the scene: the camera orbits left (-x as seen
        // from -z looking +z), so the surface follows the cursor.
        o.drag(PI / 2.0 / 0.01, 0.0);
        assert!(o.from().approx_eq(p(-5.0, 0.0, 0.0)));
    }

    #[test]
    fn drag_down_raises_camera() {
        let from = p(0.0, 0.0, -5.0);
        let mut o = Orbit::new(from, p(0.0, 0.0, 0.0), up());
        o.drag(0.0, 10.0);
        assert!(o.from().y > from.y);
    }

    #[test]
    fn zoom_step_is_exp_minus_tenth() {
        let to = p(0.0, 0.5, 0.0);
        let mut o = Orbit::new(p(-3.0, 1.0, 2.5), to, up());
        let radius = (o.from() - to).magnitude();
        o.zoom(1.0);
        let zoomed = (o.from() - to).magnitude();
        assert!((zoomed - radius * (-0.1f64).exp()).abs() < 1e-12);
    }

    #[test]
    fn radius_reports_the_clamped_zoom() {
        // The window scales the focal distance by the actual radius ratio,
        // which is 1 once zoom is clamped.
        let mut o = Orbit::new(p(0.0, 0.0, -5.0), p(0.0, 0.0, 0.0), up());
        assert!((o.radius() - 5.0).abs() < 1e-12);
        o.zoom(-1e6);
        let r = o.radius();
        o.zoom(-1.0);
        assert_eq!(o.radius(), r);
    }

    #[test]
    fn zoom_out_is_bounded() {
        let to = p(0.0, 0.5, 0.0);
        let mut o = Orbit::new(p(-3.0, 1.0, 2.5), to, up());
        let radius = (o.from() - to).magnitude();
        o.zoom(-1e6);
        let zoomed = (o.from() - to).magnitude();
        assert!(zoomed.is_finite());
        assert!(zoomed <= 1e3 * radius * (1.0 + 1e-12));
    }
}
