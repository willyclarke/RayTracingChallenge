//! Camera
//!
//! This module takes pictures of a scene.
//!
//! Map the three-dimensional scene onto a two-dimensional canvas.
//!
//! The camera is defined by the following four attributes:
//!
//! • hsize is the horizontal size (in pixels) of the canvas that the picture will be rendered to.
//!
//! • vsize is the canvas's vertical size (in pixels).
//!
//! • field_of_view is an angle that describes how much the camera can see. When the field of view
//! is small, the view will be "zoomed in," magnifying a smaller area of the scene.
//!
//! • transform is a matrix describing how the world should be oriented relative to the camera. This
//! is usually a view transformation like you implemented in the previous section.
//!

use crate::{matrix::Matrix4, ray::Ray, tuple::Tuple};

#[derive(Debug, Copy, Clone)]
pub struct Camera {
    pub hsize: usize,
    pub vsize: usize,
    pub field_of_view: f64,
    pub pixel_size: f64,
    pub half_width: f64,
    pub half_height: f64,
    transform: Matrix4,
    transform_inv: Matrix4,
    /// Anti-aliasing: pixels that differ from a neighbour by more than
    /// `aa_threshold` are re-rendered with an `aa_samples` x `aa_samples`
    /// sub-pixel grid. `1` = off (the default).
    pub aa_samples: usize,
    /// Max per-channel difference to a neighbour before a pixel counts as an
    /// edge for anti-aliasing.
    pub aa_threshold: f64,
    /// Focal blur: lens radius in world units. `0` = pinhole (the default).
    pub aperture: f64,
    /// Distance from the camera at which objects are in sharp focus.
    pub focal_distance: f64,
    /// Lens samples per pixel when `aperture > 0`.
    pub dof_samples: usize,
    /// Motion blur: rays per pixel spread over the shutter interval.
    /// `1` = shutter closed instantly at time 0 (the default).
    pub motion_samples: usize,
}

impl Camera {
    pub fn new(hsize: usize, vsize: usize, field_of_view: f64) -> Self {
        let half_view = (field_of_view / 2.0).tan();
        let aspect = hsize as f64 / vsize as f64;
        let half_width = if aspect >= 1.0 {
            half_view
        } else {
            half_view * aspect
        };
        let half_height = if aspect >= 1.0 {
            half_view / aspect
        } else {
            half_view
        };
        let pixel_size = (half_width * 2.0) / hsize as f64;
        Self {
            hsize,
            vsize,
            field_of_view,
            pixel_size,
            half_width,
            half_height,
            transform: Matrix4::identity(),
            transform_inv: Matrix4::identity(),
            aa_samples: 1,
            aa_threshold: 0.05,
            aperture: 0.0,
            focal_distance: 1.0,
            dof_samples: 1,
            motion_samples: 1,
        }
    }

    /// Enable motion blur: average `samples` rays per pixel at jittered
    /// times across the shutter interval `[0, 1]`.
    pub fn with_motion_blur(mut self, samples: usize) -> Self {
        self.motion_samples = samples.max(1);
        self
    }

    /// Does this camera need several rays per pixel (focal or motion blur)?
    pub fn is_distributed(&self) -> bool {
        self.aperture > 0.0 || self.motion_samples > 1
    }

    /// Enable focal blur: a lens of radius `aperture`, sharp at
    /// `focal_distance`, averaged over `samples` lens positions per pixel.
    pub fn with_focal_blur(mut self, aperture: f64, focal_distance: f64, samples: usize) -> Self {
        self.aperture = aperture;
        self.focal_distance = focal_distance;
        self.dof_samples = samples.max(1);
        self
    }

    /// Enable edge-detected supersampling with an `n` x `n` sub-pixel grid.
    pub fn with_antialias(mut self, n: usize) -> Self {
        self.aa_samples = n.max(1);
        self
    }

    /// Ray through the centre of pixel `(px, py)`.
    pub fn ray_for_pixel(&self, px: usize, py: usize) -> Ray {
        self.ray_for_subpixel(px, py, 0.5, 0.5)
    }

    /// Ray through pixel `(px, py)` at fractional offset `(dx, dy)` in
    /// `[0, 1)` from its top-left corner; `(0.5, 0.5)` is the centre.
    pub fn ray_for_subpixel(&self, px: usize, py: usize, dx: f64, dy: f64) -> Ray {
        //
        // the offset from the edge of the canvas to the sample point
        let xoffset = (px as f64 + dx) * self.pixel_size;
        let yoffset = (py as f64 + dy) * self.pixel_size;

        // the untransformed coordinates of the pixel in world space.
        // (remember that the camera looks toward -z, so +x is to the *left*.)
        let world_x = self.half_width - xoffset;
        let world_y = self.half_height - yoffset;

        // Using the camera matrix, transform the canvas point and the origin,
        // and then compute the ray's direction vector.
        // (remember that the canvas is at z=-1)
        let pixel = self.transform_inv * Tuple::point(world_x, world_y, -1.0);
        let origin = self.transform_inv * Tuple::point(0.0, 0.0, 0.0);
        let direction = (pixel - origin).normalize();
        Ray::new(origin, direction)
    }

    /// Ray from a point on the lens through the sub-pixel's point on the
    /// focal plane. `(lx, ly)` is the lens offset in units of the aperture
    /// radius (inside the unit disk); `(0, 0)` is the pinhole ray.
    pub fn ray_for_lens(&self, px: usize, py: usize, dx: f64, dy: f64, lx: f64, ly: f64) -> Ray {
        let pinhole = self.ray_for_subpixel(px, py, dx, dy);
        if self.aperture == 0.0 {
            return pinhole;
        }
        // Where this pixel's pinhole ray is in focus.
        let focal_point = pinhole.position(self.focal_distance);
        // The lens lies in the camera's z = 0 plane.
        let origin = self.transform_inv * Tuple::point(lx * self.aperture, ly * self.aperture, 0.0);
        Ray::new(origin, (focal_point - origin).normalize())
    }

    /// Deterministic per-pixel sample `i` of `n`: a stratified point in the
    /// unit disk (lens), a jittered sub-pixel offset, and a stratified
    /// shutter time. Hashing the pixel coordinates decorrelates neighbours
    /// so the blur is noise, not ghosts.
    pub fn lens_sample(px: usize, py: usize, i: usize, n: usize) -> (f64, f64, f64, f64, f64) {
        let hash = |i: usize, k: u64| -> f64 {
            // small integer hash -> [0, 1)
            let mut x = (px as u64).wrapping_mul(0x9E37_79B9)
                ^ (py as u64).wrapping_mul(0x85EB_CA6B)
                ^ (i as u64).wrapping_mul(0xC2B2_AE35)
                ^ k.wrapping_mul(0x27D4_EB2F);
            x ^= x >> 15;
            x = x.wrapping_mul(0x2C1B_3C6D);
            x ^= x >> 12;
            x = x.wrapping_mul(0x297A_2D39);
            x ^= x >> 15;
            (x & 0xFF_FFFF) as f64 / 16_777_216.0
        };
        // per-sample jitter, and per-pixel (sample-independent) permutation
        let h = |k: u64| hash(i, k);
        let hp = |k: u64| hash(usize::MAX, k);
        // stratify the radius over the samples; random angle
        let r = ((i as f64 + h(1)) / n as f64).sqrt();
        let theta = 2.0 * std::f64::consts::PI * h(2);
        // stratify time too, visiting the strata in a per-pixel permuted
        // order (stride coprime to n) so time isn't correlated with radius
        let stride = {
            let gcd = |mut a: usize, mut b: usize| {
                while b != 0 {
                    (a, b) = (b, a % b);
                }
                a
            };
            let want = 1 + (hp(5) * n as f64) as usize;
            (want..want + n).find(|&k| gcd(k, n) == 1).unwrap_or(1)
        };
        let offset = (hp(7) * n as f64) as usize;
        let slot = (i * stride + offset) % n;
        let time = (slot as f64 + h(6)) / n as f64;
        (r * theta.cos(), r * theta.sin(), h(3), h(4), time)
    }

    pub fn transform(&self) -> Matrix4 {
        self.transform
    }

    pub fn set_transform(&mut self, transform: &Matrix4) {
        self.transform = *transform;
        self.transform_inv = self.transform.inverse().unwrap()
    }

    pub fn with_transform(mut self, transform: Matrix4) -> Self {
        self.transform = transform;
        self.transform_inv = transform.inverse().unwrap();
        self
    }
}

impl Default for Camera {
    fn default() -> Self {
        Self::new(160, 120, std::f64::consts::PI / 2.0)
    }
}

#[cfg(test)]
mod tests {
    use crate::camera::Camera;
    use crate::loge;
    use crate::math::approx_eq;
    use crate::matrix::Matrix4;
    use crate::ray::Ray;
    use crate::tuple::Tuple;

    /// Chap 7 - Constructing a camera
    #[test]
    fn test_chap_7_16() -> Result<(), String> {
        let c = Camera::default();
        let chk = c.hsize == 160;
        let chk = chk && c.vsize == 120;
        let chk = chk && approx_eq(c.field_of_view, std::f64::consts::PI / 2.0);
        let chk = chk && c.transform().approx_eq(Matrix4::identity());
        if chk {
            Ok(())
        } else {
            Err("Constructing a camera".into())
        }
    }

    /// Chap 7 - The pixel size for a horizontal canvas
    #[test]
    fn test_chap_7_17() -> Result<(), String> {
        let c = Camera::new(200, 125, std::f64::consts::PI / 2.0);
        let chk = approx_eq(c.pixel_size, 0.01);
        if chk {
            Ok(())
        } else {
            Err("The pixel size for a horizontal canvas".into())
        }
    }

    /// Chap 7 - The pixel size for a vertical canvas
    #[test]
    fn test_chap_7_18() -> Result<(), String> {
        let c = Camera::new(125, 200, std::f64::consts::PI / 2.0);
        let chk = approx_eq(c.pixel_size, 0.01);
        if chk {
            Ok(())
        } else {
            Err("The pixel size for a vertical canvas".into())
        }
    }

    /// Chap 7 - Constructing a ray through the center of the canvas
    #[test]
    fn test_chap_7_19() -> Result<(), String> {
        let c = Camera::new(201, 101, std::f64::consts::PI / 2.0);
        let r = c.ray_for_pixel(100, 50);

        let r_expect = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, -1.0));
        let chk = r.origin.approx_eq(r_expect.origin);
        let chk = chk && r.direction.approx_eq(r_expect.direction);
        if chk {
            Ok(())
        } else {
            loge!("test_chap_7_19", "Ray r_expect:{}", r_expect);
            loge!("test_chap_7_19", "Ray r       :{}", r);
            Err("Constructing a ray through the center of the canvas".into())
        }
    }

    /// Chap 7 - Constructing a ray through a corner of the canvas
    #[test]
    fn test_chap_7_20() -> Result<(), String> {
        let c = Camera::new(201, 101, std::f64::consts::PI / 2.0);
        let r = c.ray_for_pixel(0, 0);

        let r_expect = Ray::new(
            Tuple::point(0.0, 0.0, 0.0),
            Tuple::vector(0.6651864261, 0.3325932131, -0.6685123583),
        );
        let chk = r.origin.approx_eq(r_expect.origin);
        let chk = chk && r.direction.approx_eq(r_expect.direction);
        if chk {
            Ok(())
        } else {
            loge!("test_chap_7_20", "Ray r_expect:{}", r_expect);
            loge!("test_chap_7_20", "Ray r       :{}", r);
            Err("Constructing a ray through a corner of the canvas".into())
        }
    }

    /// Chap 7 - Constructing a ray when the camera is transformed
    #[test]
    fn test_chap_7_21() -> Result<(), String> {
        let mut c = Camera::new(201, 101, std::f64::consts::PI / 2.0);
        let m =
            Matrix4::rotation_y(std::f64::consts::PI / 4.0) * Matrix4::translation(0.0, -2.0, 5.0);
        c.set_transform(&m);

        let sqrt2_over_2 = (1.0 / std::f64::consts::FRAC_1_SQRT_2) / 2.0;
        let r_expect = Ray::new(
            Tuple::point(0.0, 2.0, -5.0),
            Tuple::vector(sqrt2_over_2, 0.0, -sqrt2_over_2),
        );

        let r = c.ray_for_pixel(100, 50);

        let chk = r.origin.approx_eq(r_expect.origin);
        let chk = chk && r.direction.approx_eq(r_expect.direction);
        if chk {
            Ok(())
        } else {
            loge!("test_chap_7_20", "Ray r_expect:{}", r_expect);
            loge!("test_chap_7_20", "Ray r       :{}", r);
            Err("Constructing a ray when the camera is transformed".into())
        }
    }

    /// Chap 17 - A sub-pixel ray at (0.5, 0.5) is the pixel-centre ray
    #[test]
    fn test_chap_17_1() -> Result<(), String> {
        let c = Camera::new(201, 101, std::f64::consts::PI / 2.0);
        let centre = c.ray_for_pixel(100, 50);
        let sub = c.ray_for_subpixel(100, 50, 0.5, 0.5);
        let chk = centre.origin.approx_eq(sub.origin) && centre.direction.approx_eq(sub.direction);
        if chk {
            Ok(())
        } else {
            Err("A sub-pixel ray at (0.5, 0.5) is the pixel-centre ray".into())
        }
    }

    /// Chap 17 - Sub-pixel rays tile the pixel: the far corner of one pixel
    /// is the near corner of the next
    #[test]
    fn test_chap_17_2() -> Result<(), String> {
        let c = Camera::new(201, 101, std::f64::consts::PI / 2.0);
        let a = c.ray_for_subpixel(10, 20, 1.0, 1.0);
        let b = c.ray_for_subpixel(11, 21, 0.0, 0.0);
        let chk = a.direction.approx_eq(b.direction);
        if chk {
            Ok(())
        } else {
            loge!("test_chap_17_2", "a:{} b:{}", a, b);
            Err("Sub-pixel rays tile the pixel".into())
        }
    }

    /// Chap 17 - Anti-aliasing is off by default and enabled with a grid size
    #[test]
    fn test_chap_17_3() -> Result<(), String> {
        let c = Camera::new(10, 10, std::f64::consts::PI / 2.0);
        let aa = Camera::new(10, 10, std::f64::consts::PI / 2.0).with_antialias(4);
        let chk = c.aa_samples == 1 && aa.aa_samples == 4;
        if chk {
            Ok(())
        } else {
            Err("Anti-aliasing is off by default and enabled with a grid size".into())
        }
    }

    /// Chap 17 - With a pinhole camera the lens ray is the sub-pixel ray,
    /// whatever the lens offset
    #[test]
    fn test_chap_17_11() -> Result<(), String> {
        let c = Camera::new(201, 101, std::f64::consts::PI / 2.0);
        let a = c.ray_for_subpixel(50, 25, 0.3, 0.6);
        let b = c.ray_for_lens(50, 25, 0.3, 0.6, 0.7, -0.2);
        if a.origin.approx_eq(b.origin) && a.direction.approx_eq(b.direction) {
            Ok(())
        } else {
            Err("Pinhole lens ray must equal the sub-pixel ray".into())
        }
    }

    /// Chap 17 - Every lens ray for a pixel passes through the same point on
    /// the focal plane, and starts on the lens
    #[test]
    fn test_chap_17_12() -> Result<(), String> {
        // A level camera: the book's view_transform only yields an
        // orthonormal basis when `up` is perpendicular to the view direction.
        let from = Tuple::point(1.0, 2.0, -5.0);
        let to = Tuple::point(0.0, 2.0, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let c = Camera::new(201, 101, std::f64::consts::PI / 3.0)
            .with_transform(crate::world::view_transform(from, to, up))
            .with_focal_blur(0.5, 4.0, 8);
        let pinhole = c.ray_for_pixel(60, 40);
        let focal_point = pinhole.position(4.0);
        for (lx, ly) in [(0.0, 0.0), (1.0, 0.0), (0.0, -1.0), (0.6, 0.6), (-0.3, 0.9)] {
            let r = c.ray_for_lens(60, 40, 0.5, 0.5, lx, ly);
            // origin is on the lens disk, centred on the camera
            let lens_offset = (r.origin - from).magnitude();
            let expected_offset = 0.5 * (lx * lx + ly * ly).sqrt();
            if !approx_eq(lens_offset, expected_offset) {
                loge!(
                    "test_chap_17_12",
                    "lens offset {lens_offset} != {expected_offset}"
                );
                return Err("Lens ray must start on the lens".into());
            }
            // and the ray passes through the focal point
            let t = (focal_point - r.origin).magnitude();
            if !r.position(t).approx_eq(focal_point) {
                loge!(
                    "test_chap_17_12",
                    "ray {} misses focal point {}",
                    r,
                    focal_point
                );
                return Err("Lens rays must converge on the focal point".into());
            }
        }
        Ok(())
    }

    /// Chap 17 - Lens samples lie in the unit disk and differ between pixels
    #[test]
    fn test_chap_17_13() -> Result<(), String> {
        let mut distinct = false;
        for i in 0..16 {
            let (lx, ly, dx, dy, time) = Camera::lens_sample(3, 7, i, 16);
            let (mx, my, _, _, _) = Camera::lens_sample(4, 7, i, 16);
            if lx * lx + ly * ly > 1.0 + 1e-12
                || !(0.0..1.0).contains(&dx)
                || !(0.0..1.0).contains(&dy)
                || !(0.0..1.0).contains(&time)
            {
                return Err("Lens sample outside the unit disk / pixel".into());
            }
            if (lx, ly) != (mx, my) {
                distinct = true;
            }
        }
        if distinct {
            Ok(())
        } else {
            Err("Neighbouring pixels must get different lens samples".into())
        }
    }

    /// Chap 17 - Shutter times are stratified: with n samples, every 1/n
    /// slice of the shutter interval gets exactly one sample
    #[test]
    fn test_chap_17_17() -> Result<(), String> {
        let n = 8;
        let mut slots = [0usize; 8];
        for i in 0..n {
            let (_, _, _, _, time) = Camera::lens_sample(11, 5, i, n);
            slots[(time * n as f64) as usize] += 1;
        }
        if slots.iter().all(|&c| c == 1) {
            Ok(())
        } else {
            loge!("test_chap_17_17", "slots:{:?}", slots);
            Err("Shutter times must be stratified".into())
        }
    }
}
