//! Shape defintion
//! Abstract data structure for sphere, cubes etc.
//!

use crate::matrix::Matrix4;
use crate::{intersection::Intersections, material::Material};
use std::sync::atomic::{AtomicUsize, Ordering};

static NEXT_ID: AtomicUsize = AtomicUsize::new(1);

#[derive(Debug, Clone, Copy)]
pub struct ShapeData {
    pub id: usize,
    pub transform: Matrix4,
    pub transform_inv: Matrix4,
    pub material: Material,
}

impl ShapeData {
    pub fn new() -> Self {
        Self {
            id: NEXT_ID.fetch_add(1, Ordering::Relaxed),
            transform: Matrix4::identity(),
            transform_inv: Matrix4::identity(),
            material: Material::new(),
        }
    }

    ///
    /// Reset the NEXT_ID to the value of 1.
    ///
    pub fn reset() {
        NEXT_ID.store(1, Ordering::Relaxed);
    }
}

impl Default for ShapeData {
    fn default() -> Self {
        Self::new()
    }
}

pub trait Shape: Send + Sync {
    fn data(&self) -> &ShapeData;
    fn data_mut(&mut self) -> &mut ShapeData;

    fn id(&self) -> usize {
        self.data().id
    }

    fn set_id(&mut self, id: usize) {
        self.data_mut().id = id;
    }

    fn material(&self) -> &Material {
        &self.data().material
    }

    fn set_material(&mut self, m: Material) {
        self.data_mut().material = m;
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
    fn normal_at(&self, point: crate::tuple::Tuple) -> crate::tuple::Tuple;
    fn local_normal_at(&self, point: crate::tuple::Tuple) -> crate::tuple::Tuple;
}
