//! Ray definition and implementation
//!

// use crate::canvas::Canvas;
use crate::log::*;
// use crate::math::approx_eq;
use crate::matrix::Matrix4;
use crate::tuple::Tuple;
use std::fmt;

#[derive(Debug, Copy, Clone)]
pub struct Ray {
    pub origin: Tuple,
    pub direction: Tuple,
    /// Moment within the shutter interval `[0, 1]` this ray samples (motion
    /// blur, book chapter 17). `0` unless set with `with_time`.
    pub time: f64,
}

impl Ray {
    pub fn new(origin: Tuple, direction: Tuple) -> Self {
        Self {
            origin,
            direction,
            time: 0.0,
        }
    }

    pub fn with_time(mut self, time: f64) -> Self {
        self.time = time;
        self
    }

    pub fn transform(&self, m: Matrix4) -> Self {
        Self {
            origin: m * self.origin,
            direction: m * self.direction,
            time: self.time,
        }
    }

    /// The same ray shifted by `-offset`: how a ray looks to a shape that has
    /// moved by `offset`.
    pub fn shifted(&self, offset: Tuple) -> Self {
        Self {
            origin: self.origin - offset,
            direction: self.direction,
            time: self.time,
        }
    }

    pub fn position(&self, t: f64) -> Tuple {
        self.origin + self.direction * t
    }
}
impl fmt::Display for Ray {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "{}Ori{}:[{:.10}, {:.10}, {:.10}]. {}Dir{}:[{:.10}, {:.10}, {:.10}]",
            Color::Green,
            Color::Reset,
            self.origin.x,
            self.origin.y,
            self.origin.z,
            Color::Yellow,
            Color::Reset,
            self.direction.x,
            self.direction.y,
            self.direction.z
        )?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    // use core::f64;

    use super::*;
    // use crate::canvas::Canvas;
    // use crate::log::*;
    // use crate::math::approx_eq;
    // use crate::tuple::Tuple;
    // use std::fmt;
    use crate::{logd, loge, logi, tuple::Tuple};

    /// Chap 5 - Creating and querying a ray
    #[test]
    fn test_chap_5_1() -> Result<(), String> {
        let origin = Tuple::point(1.0, 2.0, 3.0);
        let direction = Tuple::vector(4.0, 5.0, 6.0);
        let ray = Ray::new(origin, direction);
        logi!("test_chap_5_1", "ray:{}", ray);
        logd!("test_chap_5_1", "ray:{}", ray);
        loge!("test_chap_5_1", "ray:{}", ray);

        let chk = ray.origin.approx_eq(origin);
        let chk = chk && ray.direction.approx_eq(direction);
        if chk {
            Ok(())
        } else {
            logi!("test_chap_5_1", "ray:{}", ray);
            logd!("test_chap_5_1", "ray:{}", ray);
            loge!("test_chap_5_1", "ray:{}", ray);
            Err("Creating and querying a ray".into())
        }
    }

    /// Chap 5 - Computing a point from a distance
    #[test]
    fn test_chap_5_2() -> Result<(), String> {
        let ray = Ray::new(Tuple::point(2.0, 3.0, 4.0), Tuple::vector(1.0, 0.0, 0.0));
        let chk = ray.position(0.0).approx_eq(Tuple::point(2.0, 3.0, 4.0));
        let chk = chk && ray.position(1.0).approx_eq(Tuple::point(3.0, 3.0, 4.0));
        let chk = chk && ray.position(-1.0).approx_eq(Tuple::point(1.0, 3.0, 4.0));
        let chk = chk && ray.position(2.5).approx_eq(Tuple::point(4.5, 3.0, 4.0));
        if chk {
            Ok(())
        } else {
            Err("Computing a point from a distance".into())
        }
    }

    /// Chap 5 - Translating a ray
    #[test]
    fn test_chap_5_15() -> Result<(), String> {
        let r = Ray::new(Tuple::point(1.0, 2.0, 3.0), Tuple::vector(0.0, 1.0, 0.0));
        let m = Matrix4::translation(3.0, 4.0, 5.0);
        let r2 = r.transform(m);
        let chk = r2.origin.approx_eq(Tuple::point(4.0, 6.0, 8.0));
        let chk = chk && r2.direction.approx_eq(Tuple::vector(0.0, 1.0, 0.0));
        if chk {
            Ok(())
        } else {
            Err("Translating a ray".into())
        }
    }

    /// Chap 5 - Scaling a ray
    #[test]
    fn test_chap_5_16() -> Result<(), String> {
        let r = Ray::new(Tuple::point(1.0, 2.0, 3.0), Tuple::vector(0.0, 1.0, 0.0));
        let m = Matrix4::scaling(2.0, 3.0, 4.0);
        let r2 = r.transform(m);
        let chk = r2.origin.approx_eq(Tuple::point(2.0, 6.0, 12.0));
        let chk = chk && r2.direction.approx_eq(Tuple::vector(0.0, 3.0, 0.0));
        if chk {
            Ok(())
        } else {
            Err("Scaling a ray".into())
        }
    }

    /// Chap 17 - A ray's time defaults to 0 and is kept by transforms
    #[test]
    fn test_chap_17_18() -> Result<(), String> {
        let r = Ray::new(Tuple::point(1.0, 2.0, 3.0), Tuple::vector(0.0, 1.0, 0.0));
        let moved = r
            .with_time(0.75)
            .transform(Matrix4::translation(3.0, 4.0, 5.0));
        let shifted = moved.shifted(Tuple::vector(1.0, 0.0, 0.0));
        let chk = r.time == 0.0
            && moved.time == 0.75
            && shifted.time == 0.75
            && shifted.origin.approx_eq(Tuple::point(3.0, 6.0, 8.0));
        if chk {
            Ok(())
        } else {
            Err("A ray's time defaults to 0 and is kept by transforms".into())
        }
    }
}
