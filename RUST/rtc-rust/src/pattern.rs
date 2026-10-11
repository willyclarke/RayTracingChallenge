//! Pattern definition
//! Abstract data structure for stripes, gradients, rings and checkers.

use crate::matrix::Matrix4;
use crate::shape::Shape;
use crate::tuple::Tuple;

#[derive(Debug, Clone, Copy)]
pub struct PatternData {
    transform: Matrix4, // private for protection. Use set_transform to update these two.
    transform_inv: Matrix4,
}

impl PatternData {
    pub fn new() -> Self {
        Self {
            transform: Matrix4::identity(),
            transform_inv: Matrix4::identity(),
        }
    }

    pub fn transform(&self) -> &Matrix4 {
        &self.transform
    }
    pub fn transform_inv(&self) -> &Matrix4 {
        &self.transform_inv
    }

    pub fn set_transform(&mut self, m: Matrix4) {
        self.transform = m;
        self.transform_inv = m.inverse().unwrap_or(Matrix4::identity());
    }
}

impl Default for PatternData {
    /// Construct a `PatternData` with no transformation applied
    /// (identity transform and its inverse).
    ///
    /// Implementing `Default` gives the type the conventional,
    /// argument-free constructor that the rest of Rust expects:
    /// - lets generic code that requires `T: Default` construct one with no arguments
    /// - documents the "zero value" — a pattern living in object space,
    ///   untransformed, so its colours map directly onto the unit pattern.
    ///
    /// It simply delegates to [`PatternData::new`].
    fn default() -> Self {
        Self::new()
    }
}

pub trait Pattern: std::fmt::Debug + Send + Sync {
    fn data(&self) -> &PatternData;
    fn data_mut(&mut self) -> &mut PatternData;

    fn color_at(&self, point: Tuple) -> Tuple;
    fn clone_box(&self) -> Box<dyn Pattern>;
    fn set_transform(&mut self, m: Matrix4);

    // provided: default body — patterns inherit unless they override
    fn color_at_local(&self, point: Tuple) -> Tuple {
        let local = *self.data().transform_inv() * point;
        self.color_at(local)
    }

    ///
    /// Multiplies world_point by the inverse of the transform to go to object space.
    /// And then multiplies the point in object space by the patterns inverse
    /// transform to go to pattern space.
    /// Then return result followed by getting the color at the pattern_point.
    ///
    fn color_at_shape(&self, shape: &dyn Shape, world_point: Tuple) -> Tuple {
        let object_point = *shape.transform_inv() * world_point;
        self.color_at_local(object_point)
    }
}

impl Clone for Box<dyn Pattern> {
    fn clone(&self) -> Self {
        self.clone_box()
    }
}
