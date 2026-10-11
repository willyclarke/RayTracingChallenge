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
pub struct GradientPattern {
    pub data: PatternData,
    pub a: Tuple,
    pub b: Tuple,
}

impl GradientPattern {
    pub fn new(a: Tuple, b: Tuple) -> Self {
        Self {
            data: PatternData::new(),
            a,
            b,
        }
    }
}

impl Default for GradientPattern {
    fn default() -> Self {
        Self::new(WHITE, BLACK)
    }
}

impl Pattern for GradientPattern {
    fn clone_box(&self) -> Box<dyn Pattern> {
        Box::new(*self)
    }

    /// A blending function to create the gradient pattern.
    /// This is a function that takes two values and interpolates the values
    /// between them. A basic linear interpolation looks like this:
    ///
    /// color(point, ca, cb) = ca + (cb - ca) * (px - floor(px))
    ///
    /// The two colors a and b are set during construction.
    ///
    fn color_at(&self, point: crate::tuple::Tuple) -> crate::tuple::Tuple {
        self.a + (self.b - self.a) * (point.x - point.x.floor())
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
