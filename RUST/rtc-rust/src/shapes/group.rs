//! Group of shapes
//!

use crate::intersection::Intersections;
use crate::log::*;
use crate::ray::Ray;
use crate::shape::{Shape, ShapeData};
use crate::tuple::Tuple;
use std::fmt;

// #[derive(Debug, Clone, Copy)]
#[derive(Debug, Clone)]
pub struct Group {
    pub data: ShapeData,
    pub children: Vec<usize>,
}

impl Group {
    pub fn new() -> Self {
        Self {
            data: ShapeData::new(),
            children: [].to_vec(),
        }
    }
}

impl Default for Group {
    fn default() -> Self {
        Self::new()
    }
}

impl fmt::Display for Group {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "{}Group ID:{}{}{}{}",
            Color::Yellow,
            Color::Reset,
            Color::Green,
            self.data.id,
            Color::Reset
        )?;

        Ok(())
    }
}

impl Shape for Group {
    fn add_child_id(&mut self, id: usize) {
        self.children.push(id);
    }

    fn children(&self) -> Option<&[usize]> {
        Some(&self.children)
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

    fn local_normal_at(&self, _object_point: Tuple) -> Tuple {
        unreachable!("a group has no local normal; only its leaves do")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::loge;
    use crate::matrix::Matrix4;

    /// Chap 14 - Creating a new group
    #[test]
    fn test_chap_14_1() -> Result<(), String> {
        let group = Group::new();
        let chk = group.transform().approx_eq(Matrix4::identity());
        let chk = chk && group.children.is_empty();

        if chk {
            Ok(())
        } else {
            loge!("test_chap_14_1", "Group transform:{}", group.transform());
            Err("Creating a new group".into())
        }
    }
}
