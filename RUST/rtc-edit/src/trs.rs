//! A shape's transform list read as scale, rotation and translation, the
//! three rows of boxes the editor shows, and written back.

use rtc_rust::scene::TransformOp;

/// Scale, then rotate about x, y and z (radians), then translate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Trs {
    pub scale: [f64; 3],
    pub rotate: [f64; 3],
    pub translate: [f64; 3],
}

impl Default for Trs {
    fn default() -> Self {
        Trs {
            scale: [1.0; 3],
            rotate: [0.0; 3],
            translate: [0.0; 3],
        }
    }
}

impl Trs {
    /// The list as a `Trs`, or `None` unless its ops come in the order
    /// scale, rotate x, rotate y, rotate z, translate, each at most once.
    pub fn from_ops(ops: &[TransformOp]) -> Option<Trs> {
        let mut t = Trs::default();
        // The stage the next op must come after; it only moves forward.
        let mut next = 0;
        for op in ops {
            let stage = match *op {
                TransformOp::Scale(s) => {
                    t.scale = s;
                    0
                }
                TransformOp::RotateX(r) => {
                    t.rotate[0] = r;
                    1
                }
                TransformOp::RotateY(r) => {
                    t.rotate[1] = r;
                    2
                }
                TransformOp::RotateZ(r) => {
                    t.rotate[2] = r;
                    3
                }
                TransformOp::Translate(v) => {
                    t.translate = v;
                    4
                }
                TransformOp::Shear(_) | TransformOp::Matrix(_) => return None,
            };
            if stage < next {
                return None;
            }
            next = stage + 1;
        }
        Some(t)
    }

    /// The list in that order, leaving out the steps that do nothing.
    pub fn to_ops(&self) -> Vec<TransformOp> {
        let [rx, ry, rz] = self.rotate;
        let mut ops = vec![];
        if self.scale != [1.0; 3] {
            ops.push(TransformOp::Scale(self.scale));
        }
        for (r, op) in [
            (rx, TransformOp::RotateX as fn(f64) -> TransformOp),
            (ry, TransformOp::RotateY),
            (rz, TransformOp::RotateZ),
        ] {
            if r != 0.0 {
                ops.push(op(r));
            }
        }
        if self.translate != [0.0; 3] {
            ops.push(TransformOp::Translate(self.translate));
        }
        ops
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use TransformOp::*;
    use rtc_rust::scene::{SceneDescription, ShapeDescription, ShapeKind};
    use std::fs;

    #[test]
    fn an_empty_list_is_the_default() {
        assert_eq!(Trs::from_ops(&[]), Some(Trs::default()));
        assert!(Trs::default().to_ops().is_empty());
    }

    #[test]
    fn a_full_list_reads_and_writes_back() {
        let ops = vec![
            Scale([2.0, 2.0, 2.0]),
            RotateX(0.1),
            RotateY(0.2),
            RotateZ(0.3),
            Translate([1.0, 2.0, 3.0]),
        ];
        let t = Trs::from_ops(&ops).unwrap();
        assert_eq!(
            t,
            Trs {
                scale: [2.0, 2.0, 2.0],
                rotate: [0.1, 0.2, 0.3],
                translate: [1.0, 2.0, 3.0],
            }
        );
        assert_eq!(t.to_ops(), ops);
    }

    #[test]
    fn a_partial_list_fills_in_the_rest() {
        let ops = vec![RotateY(0.5), Translate([0.0, 1.0, 0.0])];
        let t = Trs::from_ops(&ops).unwrap();
        assert_eq!(
            t,
            Trs {
                scale: [1.0; 3],
                rotate: [0.0, 0.5, 0.0],
                translate: [0.0, 1.0, 0.0],
            }
        );
        assert_eq!(t.to_ops(), ops);
    }

    #[test]
    fn no_op_steps_are_left_out_when_written() {
        let t = Trs {
            scale: [1.0; 3],
            rotate: [0.0; 3],
            translate: [4.0, 0.0, 0.0],
        };
        assert_eq!(t.to_ops(), [Translate([4.0, 0.0, 0.0])]);
    }

    #[test]
    fn lists_that_dont_fit_give_none() {
        let lists = [
            vec![Shear([1.0, 0.0, 0.0, 0.0, 0.0, 0.0])],
            vec![Matrix([[1.0, 0.0, 0.0, 0.0]; 4])],
            vec![Translate([1.0, 0.0, 0.0]), Scale([2.0, 2.0, 2.0])],
            vec![RotateX(0.1), RotateX(0.2)],
            vec![RotateZ(0.1), RotateX(0.2)],
        ];
        for ops in lists {
            assert_eq!(Trs::from_ops(&ops), None, "{ops:?}");
        }
    }

    fn all_shapes<'a>(shape: &'a ShapeDescription, out: &mut Vec<&'a ShapeDescription>) {
        out.push(shape);
        match &shape.shape {
            ShapeKind::Group { children } => children.iter().for_each(|c| all_shapes(c, out)),
            ShapeKind::Csg { left, right, .. } => {
                all_shapes(left, out);
                all_shapes(right, out);
            }
            _ => {}
        }
    }

    fn is_no_op(op: &TransformOp) -> bool {
        matches!(op, Scale(s) if *s == [1.0; 3])
            || matches!(op, RotateX(r) | RotateY(r) | RotateZ(r) if *r == 0.0)
            || matches!(op, Translate(t) if *t == [0.0; 3])
    }

    #[test]
    fn every_scene_shape_reads_without_panicking() {
        let mut fits = 0;
        for dir in ["../rtc-metal/scenes", "../rtc-rust/scenes"] {
            for entry in fs::read_dir(dir).unwrap() {
                let path = entry.unwrap().path();
                if path.extension().is_none_or(|e| e != "json") {
                    continue;
                }
                let desc = SceneDescription::from_path(&path).unwrap();
                let mut shapes = vec![];
                desc.shapes.iter().for_each(|s| all_shapes(s, &mut shapes));
                for shape in shapes {
                    let Some(t) = Trs::from_ops(&shape.transform) else {
                        continue;
                    };
                    fits += 1;
                    if !shape.transform.iter().any(is_no_op) {
                        assert_eq!(t.to_ops(), shape.transform, "{}", path.display());
                    }
                }
            }
        }
        assert!(fits > 0);
    }
}
