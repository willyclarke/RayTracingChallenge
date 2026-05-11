//! Shape defintion
//! Abstract data structure for sphere, cubes etc.
//!

use crate::intersection::Intersections;
use crate::matrix::Matrix4;
use std::sync::atomic::{AtomicUsize, Ordering};

static NEXT_ID: AtomicUsize = AtomicUsize::new(1);

#[derive(Debug, Clone, Copy)]
pub struct ShapeData {
    pub id: usize,
    pub transform: Matrix4,
    pub transform_inv: Matrix4,
}

impl ShapeData {
    pub fn new() -> Self {
        Self {
            id: NEXT_ID.fetch_add(1, Ordering::Relaxed),
            transform: Matrix4::identity(),
            transform_inv: Matrix4::identity(),
        }
    }
}

impl Default for ShapeData {
    fn default() -> Self {
        Self::new()
    }
}

pub trait Shape {
    fn data(&self) -> &ShapeData;
    fn data_mut(&mut self) -> &mut ShapeData;

    fn id(&self) -> usize {
        self.data().id
    }

    fn transform(&self) -> &Matrix4 {
        &self.data().transform
    }

    fn transform_inv(&self) -> &Matrix4 {
        &self.data().transform_inv
    }

    fn set_transform(&mut self, m: Matrix4) {
        self.data_mut().transform = m;
        self.data_mut().transform_inv = m.inverse().unwrap_or(Matrix4::identity());
    }

    fn intersect(&self, ray: &crate::ray::Ray) -> Intersections;
    fn local_intersect(&self, ray: &crate::ray::Ray) -> Intersections;
    fn local_normal_at(&self, point: crate::tuple::Tuple) -> crate::tuple::Tuple;
}
