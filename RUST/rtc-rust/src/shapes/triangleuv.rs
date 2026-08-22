//! triangleuv.rs - add support for polygonal models with triangle's
//! but make them smooth..

use crate::bounds::BoundingBox;
use crate::intersection::{Intersection, Intersections};
use crate::shape::{Shape, ShapeData};
use crate::shapes::triangle::moller_trumbore;
use crate::tuple::Tuple;

#[derive(Debug, Clone)]
pub struct TriangleUV {
    pub data: ShapeData,
    pub p1: Tuple,
    pub p2: Tuple,
    pub p3: Tuple,
    pub e1: Tuple,
    pub e2: Tuple,
    pub normal: Tuple,
    pub n1: Tuple,
    pub n2: Tuple,
    pub n3: Tuple,
}

impl TriangleUV {
    pub fn new(p1: Tuple, p2: Tuple, p3: Tuple, n1: Tuple, n2: Tuple, n3: Tuple) -> Self {
        let e1 = p2 - p1;
        let e2 = p3 - p1;
        let normal = (e2.cross(e1)).normalize();
        Self {
            data: ShapeData::new(),
            p1,
            p2,
            p3,
            e1,
            e2,
            normal,
            n1,
            n2,
            n3,
        }
    }
}

impl Default for TriangleUV {
    fn default() -> Self {
        Self::new(
            Tuple::point(0.0, 1.0, 0.0),
            Tuple::point(-1.0, 0.0, 0.0),
            Tuple::point(1.0, 0.0, 0.0),
            Tuple::vector(0.0, 1.0, 0.0),
            Tuple::vector(-1.0, 0.0, 0.0),
            Tuple::vector(1.0, 0.0, 0.0),
        )
    }
}

impl Shape for TriangleUV {
    /// Tight box over the three vertices: min/max per axis of p1, p2, p3.
    fn bounds(&self) -> BoundingBox {
        let mut bb = BoundingBox::empty();
        bb.add_point(self.p1);
        bb.add_point(self.p2);
        bb.add_point(self.p3);
        bb
    }

    fn local_normal_at_no_hit(&self, _n: Tuple) -> Tuple {
        self.normal
    }

    fn local_normal_at(&self, _n: Tuple, hit: Intersection) -> Tuple {
        self.n2 * hit.u + self.n3 * hit.v + self.n1 * (1.0 - hit.u - hit.v)
    }

    fn data(&self) -> &ShapeData {
        &self.data
    }

    fn data_mut(&mut self) -> &mut ShapeData {
        &mut self.data
    }

    fn local_intersect(&self, ray: &crate::ray::Ray) -> Intersections {
        let mut xs = Intersections::new();
        if let Some((t, u, v)) = moller_trumbore(self.p1, self.e1, self.e2, ray) {
            xs.push(Intersection::new_with_uv(t, self.id(), u, v));
        }
        xs
    }
}
