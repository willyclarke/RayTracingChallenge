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
        }
    }

    pub fn ray_for_pixel(&self, px: usize, py: usize) -> Ray {
        //
        // the offset from the edge of the canvas to the pixel's center
        let xoffset = (px as f64 + 0.5) * self.pixel_size;
        let yoffset = (py as f64 + 0.5) * self.pixel_size;

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
}
