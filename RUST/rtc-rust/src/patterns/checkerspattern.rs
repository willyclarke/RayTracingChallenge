//! Checkers pattern definition
//!
//! Colors modulo 2 in x and z to make checkers/ring.

use crate::{
    matrix::Matrix4,
    pattern::{Pattern, PatternData},
    tuple::{
        Tuple,
        colors::{BLACK, WHITE},
    },
};

#[derive(Debug, Clone, Copy)]
pub struct CheckersPattern {
    pub data: PatternData,
    pub a: Tuple,
    pub b: Tuple,
}

impl CheckersPattern {
    pub fn new(a: Tuple, b: Tuple) -> Self {
        Self {
            data: PatternData::new(),
            a,
            b,
        }
    }
}

impl Default for CheckersPattern {
    fn default() -> Self {
        Self::new(WHITE, BLACK)
    }
}

impl Pattern for CheckersPattern {
    fn clone_box(&self) -> Box<dyn Pattern> {
        Box::new(*self)
    }

    ///
    /// Create a checkers pattern
    ///
    /// color(point, ca, cb) =  color a when abs(p_x) + abs(p_y) + abs(p_z) MOD 2 == 0
    ///                         else select color b
    ///
    /// The two colors a and b are set during construction.
    ///
    fn color_at(&self, point: crate::tuple::Tuple) -> crate::tuple::Tuple {
        if ((point.x.floor() + point.y.floor() + point.z.floor()) as i32 % 2) == 0 {
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
