//! Tuple math utilities for the ray tracer.
//!
//! This module defines the core geometric types used
//! throughout the rendering pipeline.
//!

use crate::math::approx_eq;

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Tuple {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub w: f64,
}

/// A 4D tuple used to represent points and vectors.
///
/// A tuple with `w = 1.0` represents a point,
/// while `w = 0.0` represents a vector.
impl Tuple {
    /// Compare each element of the tuple with the other.
    ///
    /// # Returns true when delta is less than epsilon for each element.
    ///
    /// # Examples
    /// ```
    /// use rtc_rust::tuple::Tuple;
    ///
    /// let t = Tuple::vector(1.0, 2.0, 3.0);
    /// let chk = t.approx_eq(Tuple::new(1.0, 2.0, 3.0, 0.0));
    /// assert!(chk);
    /// ```
    pub fn approx_eq(self, other: Tuple) -> bool {
        approx_eq(self.x, other.x)
            && approx_eq(self.y, other.y)
            && approx_eq(self.z, other.z)
            && approx_eq(self.w, other.w)
    }

    /// Create a color Tuple
    ///
    /// Assign red, green and blue to tuple vars x, y, z. w=alpha (default 0.0)
    ///
    /// # Examples
    /// ```
    /// # use rtc_rust::tuple::Tuple;
    ///
    /// let col = Tuple::color(0.1, 0.1, 0.1);
    ///
    /// assert!(col.approx_eq(Tuple::color(0.1, 0.1, 0.1)));
    /// assert!(col.normalize().magnitude() > 0.0);
    /// ```
    pub fn color(red: f64, blue: f64, green: f64) -> Self {
        Self {
            x: red,
            y: blue,
            z: green,
            w: 0.0,
        }
    }

    pub fn red(self) -> f64 {
        self.x
    }

    pub fn green(self) -> f64 {
        self.y
    }

    pub fn blue(self) -> f64 {
        self.z
    }

    pub fn alpha(self) -> f64 {
        self.w
    }

    /// Cross product between vectors
    ///
    ///
    /// # Examples
    /// ```
    /// use rtc_rust::tuple::Tuple;
    ///
    /// let a = Tuple::vector(1.0, 0.0, 0.0);
    /// let b = Tuple::vector(0.0, 1.0, 0.0);
    /// let a_x_b = a.cross(b);
    ///
    /// assert!(a_x_b.approx_eq(Tuple::vector(0.0, 0.0, 1.0)));
    /// ```
    pub fn cross(self, other: Tuple) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
            w: 0.0,
        }
    }

    /// Dot product
    ///
    /// return s.x * o.x + s.y * o.y + s.z * o.z + s.w * o.w
    ///
    /// # Examples
    /// ```
    /// # use rtc_rust::tuple::Tuple;
    ///
    /// let result = Tuple::new(1.0, 0.0, 0.0, 0.0).dot(Tuple::new(1.0, 0.0, 0.0, 0.0));
    ///
    /// assert!(result == 1.0);
    /// ```
    pub fn dot(self, other: Tuple) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z + self.w * other.w
    }

    /// Normalize the Tuple
    ///
    /// Divide each element with the magnitude
    ///
    /// # Examples
    /// ```
    /// use rtc_rust::tuple::Tuple;
    ///
    /// let v_n = Tuple::vector(4.0, 0.0, 0.0).normalize();
    ///
    /// assert!(v_n.approx_eq(Tuple::vector(1.0, 0.0, 0.0)));
    /// ```
    pub fn normalize(self) -> Self {
        let m = self.magnitude();
        Self {
            x: self.x / m,
            y: self.y / m,
            z: self.z / m,
            w: self.w / m,
        }
    }

    /// Computes the magnitude (length) of the tuple.
    ///
    /// This uses the Euclidean norm:
    /// sqrt(x^2 + y^2 + z^2 + w^2)
    ///
    /// # Examples
    /// ```
    /// use rtc_rust::tuple::Tuple;
    ///
    /// let t = Tuple::vector(1.0, 2.0, 3.0);
    /// let m = t.magnitude();
    /// assert!(m > 0.0);
    /// ```
    pub fn magnitude(self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z + self.w * self.w).sqrt()
    }

    /// CTOR - all four elements.
    ///
    /// Create a new tuple from x, y, z, w f64 inputs.
    ///
    /// # Examples
    /// ```
    /// use rtc_rust::tuple::Tuple;
    ///
    /// let t = Tuple::new(1.0, 2.0, 3.0, 4.0);
    /// assert!(t.w > 0.0);
    /// ```
    pub fn new(x: f64, y: f64, z: f64, w: f64) -> Self {
        Self { x, y, z, w }
    }

    /// point  CTOR - set w to 1.0
    ///
    /// When ```w``` is exactly 1.0_f64 the tuple is a point.
    ///
    /// # Examples
    /// ```
    /// use rtc_rust::tuple::Tuple;
    ///
    /// let v = Tuple::point(1.0, 2.0, 3.0);
    ///
    /// assert!(v.w == 1.0_f64);
    /// ```
    pub fn point(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z, w: 1.0 }
    }

    /// Check if tuple is a point.
    ///
    /// When ```w``` is exactly 1.0_f64 the tuple is a point.
    ///
    /// # Examples
    /// ```
    /// use rtc_rust::tuple::Tuple;
    /// let p = Tuple::point(1.0, 2.0, 3.0);
    /// assert!(p.is_point());
    /// ```
    pub fn is_point(self) -> bool {
        approx_eq(self.w, 1.0)
    }

    /// vector CTOR - set w to 0.0
    ///
    /// When ```w``` is exactly 0.0_f64 the tuple is a vector.
    ///
    /// # Examples
    /// ```
    /// use rtc_rust::tuple::Tuple;
    ///
    /// let v = Tuple::vector(1.0, 2.0, 3.0);
    ///
    /// assert!(v.w == 0.0_f64);
    /// ```
    pub fn vector(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z, w: 0.0 }
    }

    pub fn is_vector(self) -> bool {
        approx_eq(self.w, 0.0)
    }
}

use std::fmt;
use std::ops::Add;
use std::ops::Div;
use std::ops::Mul;
use std::ops::Neg;
use std::ops::Sub;

impl Add for Tuple {
    type Output = Tuple;

    fn add(self, rhs: Tuple) -> Tuple {
        Tuple {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
            w: self.w + rhs.w,
        }
    }
}

impl fmt::Display for Tuple {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if approx_eq(self.w, 1.0) {
            write!(f, "Point({:.3}, {:.3}, {:.3})", self.x, self.y, self.z)
        } else if approx_eq(self.w, 0.0) {
            write!(f, "Vector({:.3}, {:.3}, {:.3})", self.x, self.y, self.z)
        } else {
            write!(
                f,
                "Tuple({:.3}, {:.3}, {:.3}, {:.3})",
                self.x, self.y, self.z, self.w
            )
        }
    }
}

impl Div<f64> for Tuple {
    type Output = Tuple;

    fn div(self, rhs: f64) -> Tuple {
        Tuple {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
            w: self.w / rhs,
        }
    }
}

impl Mul<f64> for Tuple {
    type Output = Tuple;

    fn mul(self, rhs: f64) -> Tuple {
        Tuple {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
            w: self.w * rhs,
        }
    }
}

impl Mul<Tuple> for f64 {
    type Output = Tuple;

    fn mul(self, rhs: Tuple) -> Tuple {
        Tuple {
            x: rhs.x * self,
            y: rhs.y * self,
            z: rhs.z * self,
            w: rhs.w * self,
        }
    }
}

impl Mul<Tuple> for Tuple {
    type Output = Tuple;

    /// Blend two tuples - create a new tuple with Hadamard product.
    ///
    /// Each component is multplied by rhs's component.
    ///
    fn mul(self, rhs: Tuple) -> Tuple {
        Tuple {
            x: rhs.x * self.x,
            y: rhs.y * self.y,
            z: rhs.z * self.z,
            w: rhs.w * self.w,
        }
    }
}

impl Neg for Tuple {
    type Output = Tuple;

    fn neg(self) -> Tuple {
        Tuple {
            x: -self.x,
            y: -self.y,
            z: -self.z,
            w: -self.w,
        }
    }
}

impl Sub for Tuple {
    type Output = Tuple;

    fn sub(self, rhs: Tuple) -> Tuple {
        Tuple {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
            w: self.w - rhs.w,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::log::{logd, loge, logi};

    #[test]
    fn test_chap_1_01() {
        let p = Tuple::point(4.0, -4.0, 3.0);
        assert_eq!(p.w, 1.0);
    }

    /// Chap1 - A tuple with w=1.0 is a point
    #[test]
    fn test_chap_1_02() -> Result<(), String> {
        let a = Tuple::point(4.3, -4.2, 3.1);
        let approx_ok = a.approx_eq(Tuple::point(4.3, -4.2, 3.1));

        if a.is_point() && !a.is_vector() && approx_ok {
            Ok(())
        } else {
            logd!("test_chap_1_02", "w:{}", a.w);
            loge!("test_chap_1_02", "w:{}", a.w);
            logi!("test_chap_1_02", "w:{}", a.w);
            Err("Chap1 - A tuple with w=1.0 is a point".into())
        }
    }

    /// Chap1 - A tuple with w=0 is a vector
    #[test]
    fn test_chap_1_03() -> Result<(), String> {
        let a = Tuple::vector(4.3, -4.2, 3.1);
        let approx_ok = a.approx_eq(Tuple::vector(4.3, -4.2, 3.1));

        if !a.is_point() && a.is_vector() && approx_ok {
            Ok(())
        } else {
            Err("Chap1 - A tuple with w=0 is a vector".into())
        }
    }

    /// Chap1 - point() creates tuples with w=1
    #[test]
    fn test_chap_1_04() -> Result<(), String> {
        let p = Tuple::point(4.0, -4.0, 3.0);
        let approx_ok = p.approx_eq(Tuple::point(4.0, -4.0, 3.0));
        if p.is_point() && !p.is_vector() && approx_ok {
            Ok(())
        } else {
            Err("Chap1 - point() creates tuples with w=1".into())
        }
    }

    /// Chap 1 - Adding two tuples
    #[test]
    fn test_chap_1_05() -> Result<(), String> {
        let a1 = Tuple::point(3.0, -2.0, 5.0);
        let a2 = Tuple::vector(-2.0, 3.0, 1.0);
        let sum = a1 + a2;
        let chk = sum.approx_eq(Tuple::point(1.0, 1.0, 6.0));
        if chk {
            Ok(())
        } else {
            Err("Chap 1 - Adding two tuples".into())
        }
    }

    /// Chap 1 - Subtracting two points
    #[test]
    fn test_chap_1_06() -> Result<(), String> {
        let a1 = Tuple::point(3.0, 2.0, 1.0);
        let a2 = Tuple::point(5.0, 6.0, 7.0);
        let sub = a1 - a2;
        let chk = sub.approx_eq(Tuple::vector(-2.0, -4.0, -6.0));
        if chk {
            Ok(())
        } else {
            Err("Subtracting two points".into())
        }
    }

    /// Chap 1 - Subtracting a vector from a point
    #[test]
    fn test_chap_1_07() -> Result<(), String> {
        let p = Tuple::point(3.0, 2.0, 1.0);
        let v = Tuple::vector(5.0, 6.0, 7.0);
        logi!("test_chap_1_07", "p:{}. v:{}", p, v);
        let sub = p - v;
        let chk = sub.approx_eq(Tuple::point(-2.0, -4.0, -6.0));
        if chk {
            Ok(())
        } else {
            Err("Subtracting a vector from a point".into())
        }
    }

    /// Chap 1 -Subtracting two vectors
    #[test]
    fn test_chap_1_08() -> Result<(), String> {
        let v1 = Tuple::vector(3.0, 2.0, 1.0);
        let v2 = Tuple::vector(5.0, 6.0, 7.0);

        logi!("test_chap_1_08", "v1:{}", v1);

        let sub = v1 - v2;
        let chk = sub.approx_eq(Tuple::vector(-2.0, -4.0, -6.0));
        if chk {
            Ok(())
        } else {
            Err("Subtracting two vectors".into())
        }
    }

    /// Chap x - Subtracting a vector from the zero vector
    #[test]
    fn test_chap_1_09() -> Result<(), String> {
        let zero = Tuple::vector(0.0, 0.0, 0.0);
        let v = Tuple::vector(1.0, -2.0, 3.0);
        let chk = (zero - v).approx_eq(Tuple::vector(-1.0, 2.0, -3.0));
        if chk {
            Ok(())
        } else {
            Err("Subtracting a vector from the zero vector".into())
        }
    }

    /// Chap x - Negating a tuple
    #[test]
    fn test_chap_1_10() -> Result<(), String> {
        let a = Tuple::new(1.0, -2.0, 3.0, -4.0);
        let neg_a = -a;
        let chk = neg_a.approx_eq(Tuple::new(-1.0, 2.0, -3.0, 4.0));
        if chk {
            Ok(())
        } else {
            Err("Negating a tuple".into())
        }
    }

    /// Chap x - Multiplying a tuple by a scalar
    #[test]
    fn test_chap_1_11() -> Result<(), String> {
        let a = Tuple::new(1.0, -2.0, 3.0, -4.0);
        let b = a * 3.5;
        let c = 3.5 * a;
        // logd!(
        //     "Mul<f64>",
        //     "x = {:?}, y = {:?}, z = {:?}, w = {:?}",
        //     b.x,
        //     b.y,
        //     b.z,
        //     b.w
        // );
        // loge!(
        //     "Mul<f64>",
        //     "x = {:?}, y = {:?}, z = {:?}, w = {:?}",
        //     b.x,
        //     b.y,
        //     b.z,
        //     b.w
        // );
        // logi!(
        //     "Mul<Tuple>",
        //     "x = {:?}, y = {:?}, z = {:?}, w = {:?}",
        //     c.x,
        //     c.y,
        //     c.z,
        //     c.w
        // );
        let chk = b.approx_eq(Tuple::new(3.5, -7.0, 10.5, -14.0));
        let chk = chk && c.approx_eq(Tuple::new(3.5, -7.0, 10.5, -14.0));
        if chk {
            Ok(())
        } else {
            Err("Multiplying a tuple by a scalar".into())
        }
    }

    /// Chap 1 - Multiplying a tuple by a fraction
    #[test]
    fn test_chap_1_12() -> Result<(), String> {
        let a = Tuple::new(1.0, -2.0, 3.0, -4.0);
        let b = 0.5 * a;
        let chk = b.approx_eq(Tuple::new(0.5, -1.0, 1.5, -2.0));
        if chk {
            Ok(())
        } else {
            Err("Multiplying a tuple by a fraction".into())
        }
    }

    /// Chap 1 - Dividing a tuple by a scalar
    #[test]
    fn test_chap_1_13() -> Result<(), String> {
        let a = Tuple::new(1.0, -2.0, 3.0, -4.0);
        let b = a / 2.0;

        logi!(
            "Mul<Tuple> b",
            "x = {:?}, y = {:?}, z = {:?}, w = {:?}. b:{}",
            b.x,
            b.y,
            b.z,
            b.w,
            b
        );

        let chk = b.approx_eq(Tuple::new(0.5, -1.0, 1.5, -2.0));
        if chk {
            Ok(())
        } else {
            Err("Dividing a tuple by a scalar".into())
        }
    }

    /// Chap 1 - Computing the magnitude of vector(1, 0, 0)
    #[test]
    fn test_chap_1_14() -> Result<(), String> {
        let v = Tuple::vector(1.0, 0.0, 0.0);
        let chk = approx_eq(v.magnitude(), 1.0);
        if chk {
            Ok(())
        } else {
            Err("Computing the magnitude of vector(1, 0, 0)".into())
        }
    }

    /// Chap 1 - Computing the magnitude of vector(0, 1, 0)
    #[test]
    fn test_chap_1_15() -> Result<(), String> {
        let v = Tuple::vector(0.0, 1.0, 0.0);
        let chk = approx_eq(v.magnitude(), 1.0);
        if chk {
            Ok(())
        } else {
            Err("Computing the magnitude of vector(0, 1, 0)".into())
        }
    }

    /// Chap 1 - Computing the magnitude of vector(0, 0, 1)
    #[test]
    fn test_chap_1_16() -> Result<(), String> {
        let v = Tuple::vector(0.0, 0.0, 1.0);
        let chk = approx_eq(v.magnitude(), 1.0);
        if chk {
            Ok(())
        } else {
            Err("Computing the magnitude of vector(0, 0, 1)".into())
        }
    }

    /// Chap 1 - Computing the magnitude of vector(1, 2, 3)
    #[test]
    fn test_chap_1_17() -> Result<(), String> {
        let v = Tuple::vector(1.0, 2.0, 3.0);
        let chk = approx_eq(v.magnitude(), (14.0_f64).sqrt());
        if chk {
            Ok(())
        } else {
            Err("Computing the magnitude of vector(1, 2, 3)".into())
        }
    }

    /// Chap 1 - Computing the magnitude of vector(-1, -2, -3)
    #[test]
    fn test_chap_1_18() -> Result<(), String> {
        let v = Tuple::vector(-1.0, -2.0, -3.0);
        let chk = approx_eq(v.magnitude(), (14.0_f64).sqrt());
        if chk {
            Ok(())
        } else {
            Err("Computing the magnitude of vector(-1, -2, -3)".into())
        }
    }

    /// Chap 1 - Normalizing vector(4, 0, 0) gives (1, 0, 0)
    #[test]
    fn test_chap_1_19() -> Result<(), String> {
        let v = Tuple::vector(4.0, 0.0, 0.0);
        let chk = v.normalize().approx_eq(Tuple::vector(1.0, 0.0, 0.0));
        if chk {
            Ok(())
        } else {
            Err("Normalizing vector(4, 0, 0) gives (1, 0, 0)".into())
        }
    }

    /// Chap 1 -Normalizing vector(1, 2, 3)
    #[test]
    fn test_chap_1_20() -> Result<(), String> {
        let v = Tuple::vector(1.0, 2.0, 3.0);
        let v_n = v.normalize();
        let sqrt_14 = (14.0_f64).sqrt();
        let chk = v_n.approx_eq(Tuple::vector(1.0 / sqrt_14, 2.0 / sqrt_14, 3.0 / sqrt_14));
        if chk {
            Ok(())
        } else {
            Err("Normalizing vector(1, 2, 3)".into())
        }
    }

    /// Chap 1 -The magnitude of a normalized vector
    #[test]
    fn test_chap_1_21() -> Result<(), String> {
        let v = Tuple::vector(1.0, 2.0, 3.0);
        let chk = approx_eq(1.0, v.normalize().magnitude());
        if chk {
            Ok(())
        } else {
            Err("The magnitude of a normalized vector".into())
        }
    }

    /// Chap x - The dot product of two tuples
    #[test]
    fn test_chap_1_22() -> Result<(), String> {
        let a = Tuple::vector(1.0, 2.0, 3.0);
        let b = Tuple::vector(2.0, 3.0, 4.0);
        let result = a.dot(b);
        let chk = approx_eq(result, 20.0);
        if chk {
            Ok(())
        } else {
            Err("The dot product of two tuples".into())
        }
    }

    /// Chap 1 - The cross product of two vectors
    #[test]
    fn test_chap_1_23() -> Result<(), String> {
        let a = Tuple::vector(1.0, 2.0, 3.0);
        let b = Tuple::vector(2.0, 3.0, 4.0);
        let a_x_b = a.cross(b);
        let b_x_a = b.cross(a);
        let chk = Tuple::vector(-1.0, 2.0, -1.0).approx_eq(a_x_b);
        let chk = chk && Tuple::vector(1.0, -2.0, 1.0).approx_eq(b_x_a);
        if chk {
            Ok(())
        } else {
            Err("The cross product of two vectors".into())
        }
    }

    #[derive(Debug, Copy, Clone, PartialEq)]
    pub struct Projectile {
        pub position: Tuple,
        pub velocity: Tuple,
    }

    impl Projectile {
        pub fn new(p: Tuple, v: Tuple) -> Self {
            Self {
                position: p,
                velocity: v,
            }
        }

        pub fn tick(self, env: Environment) -> Projectile {
            Self {
                position: self.position + self.velocity,
                velocity: self.velocity + env.gravity + env.wind,
            }
        }
    }

    #[derive(Debug, Copy, Clone, PartialEq)]
    pub struct Environment {
        pub gravity: Tuple,
        pub wind: Tuple,
    }

    impl Environment {
        pub fn new(g: Tuple, w: Tuple) -> Self {
            Self {
                gravity: g,
                wind: w,
            }
        }
    }

    /// Chap 1 - Putting It Together
    #[test]
    #[ignore]
    fn test_chap_1_24_putting_it_together() -> Result<(), String> {
        let start_pos = Tuple::point(0.0, 1.0, 0.0);
        let start_vel = Tuple::vector(1.0, 1.0, 0.0);
        let projectile = Projectile::new(start_pos, start_vel);

        let gravity = Tuple::vector(0.0, -0.1, 0.0);
        let wind = Tuple::vector(-0.01, 0.0, 0.0);
        let environment = Environment::new(gravity, wind);

        let mut projectile = projectile.tick(environment);

        while projectile.position.y > 0.0 {
            logi!("projectile", "y:{}", projectile.position.y);
            projectile = projectile.tick(environment);
        }
        logi!("projectile", "final y:{}", projectile.position.y);

        let chk = projectile.position.y < 0.0;
        if chk {
            Ok(())
        } else {
            Err("Putting It Together".into())
        }
    }

    /// Chap 2 - Adding colors
    #[test]
    fn test_chap_2_1() -> Result<(), String> {
        let c1 = Tuple::color(0.9, 0.6, 0.75);
        let c2 = Tuple::color(0.7, 0.1, 0.25);
        let c3 = c1 + c2;
        let chk = c3.approx_eq(Tuple::color(1.6, 0.7, 1.0));
        let chk = chk && approx_eq(c3.red(), 1.6);
        let chk = chk && approx_eq(c3.green(), 0.7);
        let chk = chk && approx_eq(c3.blue(), 1.0);
        if chk {
            Ok(())
        } else {
            loge!(
                "Adding colors",
                "Result c3: r:{} g:{} b:{}",
                c3.x,
                c3.y,
                c3.z
            );
            Err("Adding colors".into())
        }
    }

    /// Chap 2 - Subtracting colors
    #[test]
    fn test_chap_2_2() -> Result<(), String> {
        let c1 = Tuple::color(0.9, 0.6, 0.75);
        let c2 = Tuple::color(0.7, 0.1, 0.25);
        let c3 = c1 - c2;
        let chk = c3.approx_eq(Tuple::color(0.2, 0.5, 0.5));
        let chk = chk && approx_eq(c3.red(), 0.2);
        let chk = chk && approx_eq(c3.green(), 0.5);
        let chk = chk && approx_eq(c3.blue(), 0.5);
        if chk {
            Ok(())
        } else {
            loge!(
                "Subtracting colors",
                "Result c3: r:{} g:{} b:{}",
                c3.x,
                c3.y,
                c3.z
            );
            Err("Subtracting colors".into())
        }
    }

    /// Chap 2 - Multiplying a color by a scalar
    #[test]
    fn test_chap_2_3() -> Result<(), String> {
        let c = Tuple::color(0.2, 0.3, 0.4) * 2.0;
        let chk = c.approx_eq(Tuple::color(0.4, 0.6, 0.8));
        let chk = chk && approx_eq(c.red(), 0.4);
        let chk = chk && approx_eq(c.green(), 0.6);
        let chk = chk && approx_eq(c.blue(), 0.80);
        if chk {
            Ok(())
        } else {
            loge!(
                "Multiplying a color by a scalar",
                "Result c: r:{} g:{} b:{}",
                c.x,
                c.y,
                c.z
            );
            Err("Multiplying a color by a scalar".into())
        }
    }

    /// Chap 2 - Multiplying colors
    #[test]
    fn test_chap_2_4() -> Result<(), String> {
        let c1 = Tuple::color(1.0, 0.2, 0.4);
        let c2 = Tuple::color(0.9, 1.0, 0.1);
        let c3 = c1.mul(c2);
        let chk = c3.approx_eq(Tuple::color(0.9, 0.2, 0.04));
        if chk {
            Ok(())
        } else {
            loge!(
                "Multiplying colors",
                "Result c3: r:{} g:{} b:{}",
                c3.red(),
                c3.green(),
                c3.blue()
            );
            Err("Multiplying colors".into())
        }
    }
}
