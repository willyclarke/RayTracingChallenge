//! Test pattern definition
//!
//! These tests assume the test pattern’s concrete function is defined like this:
//! pattern at(pattern,point) = color(pointx,pointy,pointz)
//! In other words, it takes the given point and returns a new color where the
//! color’s red/green/blue components are set to the point’s x/y/z components.
//! You can then use the color to see that the point was transformed!

use crate::{
    matrix::Matrix4,
    pattern::{Pattern, PatternData},
    shape::Shape,
    tuple::Tuple,
};

#[derive(Debug, Clone, Copy)]
pub struct TestPattern {
    pub data: PatternData,
}

impl TestPattern {
    pub fn new() -> Self {
        Self {
            data: PatternData::new(),
        }
    }
}

impl Default for TestPattern {
    fn default() -> Self {
        Self::new()
    }
}

impl Pattern for TestPattern {
    fn clone_box(&self) -> Box<dyn Pattern> {
        Box::new(*self)
    }

    fn color_at(&self, point: crate::tuple::Tuple) -> crate::tuple::Tuple {
        Tuple::color(point.x, point.y, point.z)
    }

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
