//! csg.rs - Constructive Solid Geometry: combine two shapes (or subtrees)
//! with union, intersection, or difference.
//!
//! A `Csg` node is a two-child group: `children[0]` is the left shape,
//! `children[1]` the right. Ray hits are always on the leaves; the CSG
//! operation only decides which of the children's intersections survive
//! (see `intersection_allowed` and `World::filter_intersections`).

use crate::bounds::BoundingBox;
use crate::intersection::{Intersection, Intersections};
use crate::ray::Ray;
use crate::shape::{Shape, ShapeData};
use crate::tuple::Tuple;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CsgOperation {
    Union,
    Intersection,
    Difference,
}

/// The CSG rule table: does an intersection survive, given which child was
/// hit (`lhit`) and whether the ray is currently inside the left (`inl`)
/// and right (`inr`) shapes?
///
/// ```
/// use rtc_rust::shapes::csg::{CsgOperation, intersection_allowed};
///
/// // union keeps hits on the outside of both shapes
/// assert!(intersection_allowed(CsgOperation::Union, true, false, false));
/// assert!(!intersection_allowed(CsgOperation::Union, true, false, true));
/// ```
pub fn intersection_allowed(op: CsgOperation, lhit: bool, inl: bool, inr: bool) -> bool {
    match op {
        CsgOperation::Union => (lhit && !inr) || (!lhit && !inl),
        CsgOperation::Intersection => (lhit && inr) || (!lhit && inl),
        CsgOperation::Difference => (lhit && !inr) || (!lhit && inl),
    }
}

#[derive(Debug, Clone)]
pub struct Csg {
    pub data: ShapeData,
    pub operation: CsgOperation,
    /// [left, right] once both children are added via `World::add_child`.
    pub children: Vec<usize>,
    pub bounds: BoundingBox,
}

impl Csg {
    pub fn new(operation: CsgOperation) -> Self {
        Self {
            data: ShapeData::new(),
            operation,
            children: Vec::new(),
            bounds: BoundingBox::empty(),
        }
    }

    pub fn left(&self) -> usize {
        self.children[0]
    }

    pub fn right(&self) -> usize {
        self.children[1]
    }
}

impl Shape for Csg {
    fn add_child_id(&mut self, id: usize) {
        debug_assert!(
            self.children.len() < 2,
            "a CSG node has exactly two children"
        );
        self.children.push(id);
    }

    fn bounds(&self) -> BoundingBox {
        self.bounds
    }

    fn children(&self) -> Option<&[usize]> {
        Some(&self.children)
    }

    fn csg_operation(&self) -> Option<CsgOperation> {
        Some(self.operation)
    }

    fn data(&self) -> &ShapeData {
        &self.data
    }

    fn data_mut(&mut self) -> &mut ShapeData {
        &mut self.data
    }

    fn local_intersect(&self, _ray: &Ray) -> Intersections {
        Intersections::new()
    }

    fn local_normal_at(&self, _object_point: Tuple, _hit: Intersection) -> Tuple {
        unreachable!("a CSG node has no local normal; only its leaves do")
    }

    fn set_bounds(&mut self, bb: BoundingBox) {
        self.bounds = bb;
    }

    fn set_children(&mut self, ids: Vec<usize>) {
        self.children = ids;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::loge;
    use crate::matrix::Matrix4;

    /// Chap 17 - CSG is created with an operation and two shapes
    #[test]
    fn test_chap_17_1() -> Result<(), String> {
        let c = Csg::new(CsgOperation::Union);
        let chk = c.operation == CsgOperation::Union
            && c.children.is_empty()
            && c.transform().approx_eq(Matrix4::identity());
        if chk {
            Ok(())
        } else {
            loge!("test_chap_17_1", "operation:{:?}", c.operation);
            Err("CSG is created with an operation and two shapes".into())
        }
    }

    /// Chap 17 - Evaluating the rule for a CSG operation
    #[test]
    fn test_chap_17_2() -> Result<(), String> {
        use CsgOperation::*;
        // (op, lhit, inl, inr, expected) — the book's scenario outline
        let table = [
            (Union, true, true, true, false),
            (Union, true, true, false, true),
            (Union, true, false, true, false),
            (Union, true, false, false, true),
            (Union, false, true, true, false),
            (Union, false, true, false, false),
            (Union, false, false, true, true),
            (Union, false, false, false, true),
            (Intersection, true, true, true, true),
            (Intersection, true, true, false, false),
            (Intersection, true, false, true, true),
            (Intersection, true, false, false, false),
            (Intersection, false, true, true, true),
            (Intersection, false, true, false, true),
            (Intersection, false, false, true, false),
            (Intersection, false, false, false, false),
            (Difference, true, true, true, false),
            (Difference, true, true, false, true),
            (Difference, true, false, true, false),
            (Difference, true, false, false, true),
            (Difference, false, true, true, true),
            (Difference, false, true, false, true),
            (Difference, false, false, true, false),
            (Difference, false, false, false, false),
        ];
        for (op, lhit, inl, inr, expected) in table {
            let got = intersection_allowed(op, lhit, inl, inr);
            if got != expected {
                loge!(
                    "test_chap_17_2",
                    "op:{:?} lhit:{} inl:{} inr:{} -> {} (expected {})",
                    op,
                    lhit,
                    inl,
                    inr,
                    got,
                    expected
                );
                return Err("Evaluating the rule for a CSG operation".into());
            }
        }
        Ok(())
    }
}
