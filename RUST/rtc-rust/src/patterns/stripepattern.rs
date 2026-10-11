//! Stripe pattern definition
//!
//! Colors modulo 2 in x.

use crate::{
    matrix::Matrix4,
    pattern::{Pattern, PatternData},
    tuple::{
        Tuple,
        colors::{BLACK, WHITE},
    },
};

#[derive(Debug, Clone, Copy)]
pub struct StripePattern {
    pub data: PatternData,
    pub a: Tuple,
    pub b: Tuple,
}

impl StripePattern {
    pub fn new(a: Tuple, b: Tuple) -> Self {
        Self {
            data: PatternData::new(),
            a,
            b,
        }
    }
}

impl Default for StripePattern {
    fn default() -> Self {
        Self::new(WHITE, BLACK)
    }
}

impl Pattern for StripePattern {
    fn clone_box(&self) -> Box<dyn Pattern> {
        Box::new(*self)
    }

    fn color_at(&self, point: crate::tuple::Tuple) -> crate::tuple::Tuple {
        if (point.x.floor() as i32).rem_euclid(2) == 0 {
            return self.a;
        }
        self.b
    }

    fn data(&self) -> &PatternData {
        &self.data
    }

    fn data_mut(&mut self) -> &mut PatternData {
        &mut self.data
    }

    fn set_transform(&mut self, m: Matrix4) {
        self.data_mut().set_transform(m);
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
    use crate::shape::Shape;
    use crate::shapes::sphere::Sphere;
    use crate::{
        logd, loge, logi,
        tuple::{
            // Tuple,
            colors::{BLACK, WHITE},
        },
    };

    /// Chap 10 - Creating a stripe pattern
    #[test]
    fn test_chap_10_1() -> Result<(), String> {
        let pattern = StripePattern::new(WHITE, BLACK);

        let chk = pattern.a.approx_eq(WHITE);
        let chk = chk && pattern.b.approx_eq(BLACK);

        if chk {
            Ok(())
        } else {
            logi!("test_chap_10_1", "pattern: {}", pattern.a);
            logd!("test_chap_10_1", "pattern: {}", pattern.b);
            loge!("test_chap_10_1", "pattern: {}", pattern.data.transform());
            loge!("test_chap_10_1", "WHITE: {}", WHITE);
            loge!("test_chap_10_1", "BLACK: {}", BLACK);
            Err("Creating a stripe pattern".into())
        }
    }

    /// Chap 10 - A stripe pattern is constant in y
    #[test]
    fn test_chap_10_2() -> Result<(), String> {
        let pattern = StripePattern::new(WHITE, BLACK);

        let chk = pattern
            .color_at(Tuple::point(0.0, 0.0, 0.0))
            .approx_eq(WHITE);
        let chk = chk
            && pattern
                .color_at(Tuple::point(0.0, 1.0, 0.0))
                .approx_eq(WHITE);
        let chk = chk
            && pattern
                .color_at(Tuple::point(0.0, 2.0, 0.0))
                .approx_eq(WHITE);

        if chk {
            Ok(())
        } else {
            logi!("test_chap_10_2", "pattern: {}", pattern.a);
            logd!("test_chap_10_2", "pattern: {}", pattern.b);
            loge!("test_chap_10_2", "pattern: {}", pattern.data.transform());
            loge!("test_chap_10_2", "WHITE: {}", WHITE);
            loge!("test_chap_10_2", "BLACK: {}", BLACK);
            Err("Creating a stripe pattern".into())
        }
    }

    /// Chap 10 - A stripe pattern is constant in z
    #[test]
    fn test_chap_10_3() -> Result<(), String> {
        let pattern = StripePattern::new(WHITE, BLACK);

        let chk = pattern
            .color_at(Tuple::point(0.0, 0.0, 0.0))
            .approx_eq(WHITE);
        let chk = chk
            && pattern
                .color_at(Tuple::point(0.0, 0.0, 1.0))
                .approx_eq(WHITE);
        let chk = chk
            && pattern
                .color_at(Tuple::point(0.0, 0.0, 2.0))
                .approx_eq(WHITE);

        if chk {
            Ok(())
        } else {
            logi!("test_chap_10_3", "pattern: {}", pattern.a);
            logd!("test_chap_10_3", "pattern: {}", pattern.b);
            loge!("test_chap_10_3", "pattern: {}", pattern.data.transform());
            loge!("test_chap_10_3", "WHITE: {}", WHITE);
            loge!("test_chap_10_3", "BLACK: {}", BLACK);
            Err("Creating a stripe pattern".into())
        }
    }

    /// Chap 10 - A stripe pattern alternates in x
    #[test]
    fn test_chap_10_4() -> Result<(), String> {
        let pattern = StripePattern::new(WHITE, BLACK);

        let chk = pattern
            .color_at(Tuple::point(0.0, 0.0, 0.0))
            .approx_eq(WHITE);
        let chk = chk
            && pattern
                .color_at(Tuple::point(0.9, 0.0, 0.0))
                .approx_eq(WHITE);
        let chk = chk
            && pattern
                .color_at(Tuple::point(1.0, 0.0, 0.0))
                .approx_eq(BLACK);
        let chk = chk
            && pattern
                .color_at(Tuple::point(-0.1, 0.0, 0.0))
                .approx_eq(BLACK);
        let chk = chk
            && pattern
                .color_at(Tuple::point(-1.0, 0.0, 0.0))
                .approx_eq(BLACK);
        let chk = chk
            && pattern
                .color_at(Tuple::point(-1.1, 0.0, 0.0))
                .approx_eq(WHITE);

        if chk {
            Ok(())
        } else {
            logi!("test_chap_10_4", "pattern: {}", pattern.a);
            logd!("test_chap_10_4", "pattern: {}", pattern.b);
            loge!("test_chap_10_4", "pattern: {}", pattern.data.transform());
            loge!("test_chap_10_4", "WHITE: {}", WHITE);
            loge!("test_chap_10_4", "BLACK: {}", BLACK);
            loge!(
                "test_chap_10_4",
                "color_at: {}: color: {}",
                Tuple::point(-1.1, 0.0, 0.0),
                pattern.color_at(Tuple::point(-1.1, 0.0, 0.0))
            );
            Err("Creating a stripe pattern".into())
        }
    }

    /// Chap 10 - Stripes with an object transformation
    #[test]
    fn test_chap_10_6() -> Result<(), String> {
        let mut shape = Sphere::new();
        shape.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));
        let pattern = StripePattern::new(WHITE, BLACK);
        let c = pattern.color_at_shape(&shape, Tuple::point(1.5, 0.0, 0.0));
        if c.approx_eq(WHITE) {
            Ok(())
        } else {
            loge!("test_chap_10_6", "expected WHITE, got: {}", c);
            Err("Stripes with an object transformation".into())
        }
    }

    /// Chap 10 - Stripes with a pattern transformation
    #[test]
    fn test_chap_10_7() -> Result<(), String> {
        let shape = Sphere::new();
        let mut pattern = StripePattern::new(WHITE, BLACK);
        pattern.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));
        let c = pattern.color_at_shape(&shape, Tuple::point(1.5, 0.0, 0.0));
        if c.approx_eq(WHITE) {
            Ok(())
        } else {
            loge!("test_chap_10_7", "expected WHITE, got: {}", c);
            Err("Stripes with a pattern transformation".into())
        }
    }

    /// Chap 10 - Stripes with both an object and a pattern transformation
    #[test]
    fn test_chap_10_8() -> Result<(), String> {
        let mut shape = Sphere::new();
        shape.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));
        let mut pattern = StripePattern::new(WHITE, BLACK);
        pattern.set_transform(Matrix4::translation(0.5, 0.0, 0.0));
        let c = pattern.color_at_shape(&shape, Tuple::point(2.5, 0.0, 0.0));
        if c.approx_eq(WHITE) {
            Ok(())
        } else {
            loge!("test_chap_10_8", "expected WHITE, got: {}", c);
            Err("Stripes with both an object and a pattern transformation".into())
        }
    }

    /// Chap 10 - The default pattern transformation
    #[test]
    fn test_chap_10_9() -> Result<(), String> {
        let pattern = PatternData::new();

        let chk = pattern.transform().approx_eq(Matrix4::identity());
        if chk {
            Ok(())
        } else {
            Err("The default pattern transformation".into())
        }
    }

    /// Chap 10 - Assigning a transformation
    #[test]
    fn test_chap_10_10() -> Result<(), String> {
        let mut pattern = PatternData::new();
        pattern.set_transform(Matrix4::translation(1.0, 2.0, 3.0));

        let chk = pattern
            .transform()
            .approx_eq(Matrix4::translation(1.0, 2.0, 3.0));
        if chk {
            Ok(())
        } else {
            Err("Assigning a transformation".into())
        }
    }
}
