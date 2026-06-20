//! Ring pattern definition
//!
//! Colors modulo 2 in x and z to make circle/ring.

use crate::{
    matrix::Matrix4,
    pattern::{Pattern, PatternData},
    shape::Shape,
    tuple::{
        Tuple,
        colors::{BLACK, WHITE},
    },
};

#[derive(Debug, Clone, Copy)]
pub struct RingPattern {
    pub data: PatternData,
    pub a: Tuple,
    pub b: Tuple,
}

impl RingPattern {
    pub fn new(a: Tuple, b: Tuple) -> Self {
        Self {
            data: PatternData::new(),
            a,
            b,
        }
    }
}

impl Default for RingPattern {
    fn default() -> Self {
        Self::new(WHITE, BLACK)
    }
}

impl Pattern for RingPattern {
    fn clone_box(&self) -> Box<dyn Pattern> {
        Box::new(*self)
    }

    ///
    /// Create a circle pattern
    ///
    /// color(point, ca, cb) =  color a when sqrt(ca^2 + cb^2) MOD 2 == 0
    ///                         else select color b
    ///
    /// The two colors a and b are set during construction.
    ///
    fn color_at(&self, point: crate::tuple::Tuple) -> crate::tuple::Tuple {
        if (((point.x * point.x + point.z * point.z).sqrt()).floor() as i32 % 2) == 0 {
            return self.a;
        }
        self.b
    }

    /// This method resembles the stripe_at_object from the RTC book.
    ///
    /// Multiplies world_point by the inverse of the transform to go to object space.
    /// And then multiplies the point in object space by the patterns inverse
    /// transform to go to pattern space.
    /// Then return result followed by getting the color at the pattern_point.
    ///
    fn color_at_shape(&self, shape: &dyn Shape, world_point: Tuple) -> crate::tuple::Tuple {
        let object_point = *shape.transform_inv() * world_point;
        let pattern_point = *self.data.transform_inv() * object_point;
        self.color_at(pattern_point)
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
