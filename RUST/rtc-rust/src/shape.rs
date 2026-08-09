//! Shape defintion
//! Abstract data structure for sphere, cubes etc.
//!

use crate::bounds::BoundingBox;
use crate::matrix::Matrix4;
use crate::ray::Ray;
use crate::tuple::Tuple;
use crate::{intersection::Intersections, material::Material};
use std::sync::atomic::{AtomicUsize, Ordering};

static NEXT_ID: AtomicUsize = AtomicUsize::new(1);

#[derive(Debug, Clone)]
pub struct ShapeData {
    pub id: usize,
    pub transform: Matrix4,
    pub transform_inv: Matrix4,
    pub material: Material,
    pub parent: Option<usize>, // None = root; Some(id) = enclosing group
}

impl ShapeData {
    /// Accessor to children
    /// default: leaf
    /// override as needed in group's
    pub fn children(&self) -> Option<&[usize]> {
        None
    }

    pub fn new() -> Self {
        Self {
            id: NEXT_ID.fetch_add(1, Ordering::Relaxed),
            transform: Matrix4::identity(),
            transform_inv: Matrix4::identity(),
            material: Material::new(),
            parent: None,
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
    fn add_child_id(&mut self, _id: usize) {}

    fn bounds(&self) -> BoundingBox {
        BoundingBox::new(
            Tuple::point(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
            Tuple::point(f64::INFINITY, f64::INFINITY, f64::INFINITY),
        )
    }

    fn children(&self) -> Option<&[usize]> {
        None
    }

    fn data(&self) -> &ShapeData;
    fn data_mut(&mut self) -> &mut ShapeData;

    fn id(&self) -> usize {
        self.data().id
    }

    fn intersect(&self, ray: &Ray) -> Intersections {
        let local_ray = Ray::new(
            self.data().transform_inv * ray.origin,
            self.data().transform_inv * ray.direction,
        );
        self.local_intersect(&local_ray)
    }

    fn material(&self) -> &Material {
        &self.data().material
    }

    fn normal_at(&self, world_point: Tuple) -> Tuple {
        let local_point = self.data().transform_inv * world_point;
        // move to local coordinates by use of inverse matrix
        let local_normal = self.local_normal_at(local_point);
        let mut world_normal = self.data().transform_inv.transpose() * local_normal;
        world_normal.w = 0_f64;
        world_normal.normalize()
    }

    fn set_bounds(&mut self, _bb: BoundingBox) {}

    fn set_children(&mut self, _ids: Vec<usize>) {}

    fn set_id(&mut self, id: usize) {
        self.data_mut().id = id;
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

    fn local_intersect(&self, ray: &crate::ray::Ray) -> Intersections;
    fn local_normal_at(&self, point: crate::tuple::Tuple) -> crate::tuple::Tuple;
}
