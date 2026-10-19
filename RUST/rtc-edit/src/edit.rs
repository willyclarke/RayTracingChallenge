//! Making a new scene and adding shapes and lights to one, as pure edits of
//! a `SceneDescription`.

use crate::tree::{Node, shape_at};
use rtc_rust::scene::{
    CameraDescription, LightDescription, MaterialDescription, SceneDescription, ShapeDescription,
    ShapeKind, TransformOp,
};
use std::collections::HashSet;

/// Something the editor can add.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Item {
    Sphere,
    Cube,
    Plane,
    Cylinder,
    Cone,
    Light,
}

impl Item {
    /// The item an "add" key names.
    pub fn from_key(c: char) -> Option<Item> {
        match c {
            's' => Some(Item::Sphere),
            'c' => Some(Item::Cube),
            'p' => Some(Item::Plane),
            'y' => Some(Item::Cylinder),
            'o' => Some(Item::Cone),
            'l' => Some(Item::Light),
            _ => None,
        }
    }
}

fn point_light() -> LightDescription {
    LightDescription::Point {
        position: [-10., 10., -10.],
        intensity: [1., 1., 1.],
    }
}

fn shape(kind: ShapeKind, transform: Vec<TransformOp>) -> ShapeDescription {
    ShapeDescription {
        shape: kind,
        transform,
        material: Some(MaterialDescription {
            color: Some([0.8; 3]),
            ..Default::default()
        }),
        motion: None,
        comment: None,
    }
}

/// A scene with a camera, a light and a floor.
pub fn new_scene() -> SceneDescription {
    SceneDescription {
        comment: None,
        camera: CameraDescription {
            width: 800,
            height: 600,
            field_of_view: std::f64::consts::FRAC_PI_3,
            from: [0., 1.5, -5.],
            to: [0., 1., 0.],
            up: [0., 1., 0.],
            antialias: None,
            focal_blur: None,
            motion_blur: None,
        },
        light: Some(point_light()),
        lights: None,
        settings: Default::default(),
        shapes: vec![ShapeDescription {
            shape: ShapeKind::Plane,
            transform: vec![],
            material: None,
            motion: None,
            comment: None,
        }],
    }
}

/// Adds `item` to `desc` next to `selected` and returns the node it is at.
/// A shape goes at the end of a selected group, else right after the
/// selected shape; a CSG's children are fixed, so a shape directly in one
/// goes after the (outermost) CSG. Without a selected shape it goes at the
/// end of the scene. A light goes at the end of the lights.
pub fn add(desc: &mut SceneDescription, selected: Option<&Node>, item: Item) -> Node {
    let new = match item {
        Item::Light => {
            let mut lights = match (desc.light.take(), desc.lights.take()) {
                (Some(l), _) => vec![l],
                (None, Some(lights)) => lights,
                (None, None) => vec![],
            };
            lights.push(point_light());
            let at = lights.len() - 1;
            desc.lights = Some(lights);
            return Node::Light(at);
        }
        Item::Sphere => shape(
            ShapeKind::Sphere,
            vec![TransformOp::Translate([0., 1., 0.])],
        ),
        Item::Cube => shape(ShapeKind::Cube, vec![TransformOp::Translate([0., 1., 0.])]),
        Item::Plane => shape(ShapeKind::Plane, vec![]),
        Item::Cylinder => shape(
            ShapeKind::Cylinder {
                min: Some(0.),
                max: Some(1.),
                closed: true,
            },
            vec![],
        ),
        Item::Cone => shape(
            ShapeKind::Cone {
                min: Some(0.),
                max: Some(1.),
                closed: true,
            },
            vec![],
        ),
    };
    let mut path = match selected {
        Some(Node::Shape(path)) if shape_at(desc, path).is_some() => path.clone(),
        _ => {
            desc.shapes.push(new);
            return Node::Shape(vec![desc.shapes.len() - 1]);
        }
    };
    if let Some(ShapeKind::Group { .. }) = shape_at(desc, &path).map(|s| &s.shape) {
        let ShapeKind::Group { children } = &mut shape_at_mut(desc, &path).shape else {
            unreachable!()
        };
        children.push(new);
        path.push(children.len() - 1);
        return Node::Shape(path);
    }
    while path.len() > 1 && is_csg(shape_at(desc, &path[..path.len() - 1])) {
        path.pop();
    }
    let at = path.pop().expect("a shape path is not empty") + 1;
    let siblings = if path.is_empty() {
        &mut desc.shapes
    } else {
        match &mut shape_at_mut(desc, &path).shape {
            ShapeKind::Group { children } => children,
            _ => unreachable!("only a group or the scene holds a list of shapes"),
        }
    };
    siblings.insert(at, new);
    path.push(at);
    Node::Shape(path)
}

/// Why `node` can't be deleted, if it can't.
pub fn can_delete(desc: &SceneDescription, node: &Node) -> Result<(), String> {
    match node {
        Node::Shape(path) if shape_at(desc, path).is_some() => {
            if path.len() > 1 && is_csg(shape_at(desc, &path[..path.len() - 1])) {
                Err("a CSG needs both its shapes".to_string())
            } else {
                Ok(())
            }
        }
        Node::Light(i) if *i < light_count(desc) => {
            if light_count(desc) == 1 {
                Err("a scene needs a light".to_string())
            } else {
                Ok(())
            }
        }
        _ => Err("nothing to delete".to_string()),
    }
}

fn light_count(desc: &SceneDescription) -> usize {
    match (&desc.light, &desc.lights) {
        (Some(_), _) => 1,
        (None, Some(lights)) => lights.len(),
        (None, None) => 0,
    }
}

/// Deletes `node` and returns the node to select afterwards: the next
/// sibling, else the previous one, else the parent.
pub fn delete(desc: &mut SceneDescription, node: &Node) -> Result<Node, String> {
    can_delete(desc, node)?;
    match node {
        Node::Light(i) => {
            let lights = desc.lights.as_mut().expect("more than one light is a list");
            lights.remove(*i);
            Ok(Node::Light(if *i < lights.len() { *i } else { *i - 1 }))
        }
        Node::Shape(path) => {
            let (&at, parent) = path.split_last().expect("a shape path is not empty");
            let siblings = if parent.is_empty() {
                &mut desc.shapes
            } else {
                match &mut shape_at_mut(desc, parent).shape {
                    ShapeKind::Group { children } => children,
                    _ => unreachable!("only a group or the scene holds a list of shapes"),
                }
            };
            siblings.remove(at);
            let mut next = parent.to_vec();
            if at < siblings.len() {
                next.push(at);
            } else if at > 0 {
                next.push(at - 1);
            } else if parent.is_empty() {
                return Ok(Node::Shapes);
            }
            Ok(Node::Shape(next))
        }
        _ => unreachable!("can_delete accepts only shapes and lights"),
    }
}

/// The nodes that must be unfolded for `node` to show: its enclosing groups
/// and the `Shapes` or `Lights` heading.
pub fn ancestors(node: &Node) -> Vec<Node> {
    match node {
        Node::Light(_) => vec![Node::Lights],
        Node::Shape(p) => {
            let mut out = vec![Node::Shapes];
            out.extend((1..p.len()).map(|n| Node::Shape(p[..n].to_vec())));
            out
        }
        _ => vec![],
    }
}

/// Moves the folded shapes after an insert at `at` (`by = 1`) or a delete at
/// `at` (`by = -1`): the later siblings, and what is inside them, change
/// index. The folds of a deleted shape and of what was inside it go.
pub fn shift(folded: &mut HashSet<Node>, at: &[usize], by: isize) {
    let Some((&last, parent)) = at.split_last() else {
        return;
    };
    let depth = parent.len();
    *folded = folded
        .drain()
        .filter_map(|node| match node {
            Node::Shape(mut p) if p.len() > depth && p.starts_with(parent) => {
                if by < 0 && p[depth] == last {
                    return None;
                }
                if p[depth] >= last {
                    p[depth] = p[depth]
                        .checked_add_signed(by)
                        .expect("an index stays valid");
                }
                Some(Node::Shape(p))
            }
            other => Some(other),
        })
        .collect();
}

fn is_csg(shape: Option<&ShapeDescription>) -> bool {
    matches!(shape.map(|s| &s.shape), Some(ShapeKind::Csg { .. }))
}

/// The shape a `Node::Shape` path names, which must exist.
fn shape_at_mut<'a>(desc: &'a mut SceneDescription, path: &[usize]) -> &'a mut ShapeDescription {
    let (first, rest) = path.split_first().expect("a shape path is not empty");
    rest.iter().fold(&mut desc.shapes[*first], |shape, &i| {
        match &mut shape.shape {
            ShapeKind::Group { children } => &mut children[i],
            ShapeKind::Csg { left, right, .. } => {
                if i == 0 {
                    left
                } else {
                    right
                }
            }
            _ => panic!("no child {i}"),
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;
    use std::path::Path;

    #[test]
    fn ancestors_lists_what_hides_a_node() {
        assert_eq!(ancestors(&Node::Light(2)), vec![Node::Lights]);
        assert_eq!(
            ancestors(&Node::Shape(vec![1, 0, 3])),
            vec![Node::Shapes, Node::Shape(vec![1]), Node::Shape(vec![1, 0])]
        );
        assert!(ancestors(&Node::Camera).is_empty());
    }

    fn scene() -> SceneDescription {
        let v = json!({
            "camera": {"width": 10, "height": 5, "field_of_view": 1.0,
                       "from": [0, 0, -5], "to": [0, 0, 0], "up": [0, 1, 0]},
            "lights": [
                {"type": "point", "position": [0, 0, -5], "intensity": [1, 1, 1]},
                {"type": "point", "position": [5, 0, -5], "intensity": [1, 1, 1]}
            ],
            "shapes": [
                {"type": "sphere"},
                {"type": "group", "children": [{"type": "cube"}, {"type": "plane"}]},
                {"type": "csg", "operation": "difference",
                 "left": {"type": "group", "children": [{"type": "sphere"}]},
                 "right": {"type": "sphere"}}
            ]
        });
        SceneDescription::from_json(&v.to_string()).unwrap()
    }

    fn single_light_scene() -> SceneDescription {
        let mut desc = scene();
        desc.light = desc.lights.take().map(|mut l| l.remove(0));
        desc
    }

    fn kind(desc: &SceneDescription, path: &[usize]) -> ShapeKind {
        shape_at(desc, path).unwrap().shape.clone()
    }

    fn children_len(desc: &SceneDescription, path: &[usize]) -> usize {
        match kind(desc, path) {
            ShapeKind::Group { children } => children.len(),
            k => panic!("{k:?}"),
        }
    }

    #[test]
    fn the_new_scene_builds_and_reads_back_equal() {
        let d = new_scene();
        assert!(d.build(Path::new(".")).is_ok());
        // serde_json's float parsing can land one ulp off pi/3.
        let mut back = SceneDescription::from_json(&d.to_json()).unwrap();
        assert!((back.camera.field_of_view - d.camera.field_of_view).abs() < 1e-12);
        back.camera.field_of_view = d.camera.field_of_view;
        assert_eq!(back, d);
        assert_eq!(d.camera.width, 800);
        assert_eq!(d.camera.height, 600);
        assert_eq!(d.camera.field_of_view, std::f64::consts::FRAC_PI_3);
        assert_eq!(d.camera.from, [0., 1.5, -5.]);
        assert_eq!(d.camera.to, [0., 1., 0.]);
        assert_eq!(d.camera.up, [0., 1., 0.]);
        assert_eq!(d.light, Some(point_light()));
        assert_eq!(d.lights, None);
        assert_eq!(d.shapes.len(), 1);
        let plane = &d.shapes[0];
        assert_eq!(plane.shape, ShapeKind::Plane);
        assert!(plane.transform.is_empty() && plane.material.is_none());
    }

    #[test]
    fn a_shape_goes_after_the_selected_top_level_shape() {
        let mut d = scene();
        let at = add(&mut d, Some(&Node::Shape(vec![0])), Item::Cube);
        assert_eq!(at, Node::Shape(vec![1]));
        assert_eq!(kind(&d, &[1]), ShapeKind::Cube);
        assert_eq!(d.shapes.len(), 4);
    }

    #[test]
    fn a_shape_goes_after_the_selected_shape_inside_a_group() {
        let mut d = scene();
        let at = add(&mut d, Some(&Node::Shape(vec![1, 0])), Item::Sphere);
        assert_eq!(at, Node::Shape(vec![1, 1]));
        assert_eq!(kind(&d, &[1, 1]), ShapeKind::Sphere);
        assert_eq!(children_len(&d, &[1]), 3);
    }

    #[test]
    fn a_shape_goes_at_the_end_of_a_selected_group() {
        let mut d = scene();
        let at = add(&mut d, Some(&Node::Shape(vec![1])), Item::Cube);
        assert_eq!(at, Node::Shape(vec![1, 2]));
        assert_eq!(kind(&d, &[1, 2]), ShapeKind::Cube);
    }

    #[test]
    fn a_shape_inside_a_csg_goes_after_the_csg() {
        let mut d = scene();
        let before = d.shapes[2].clone();
        let at = add(&mut d, Some(&Node::Shape(vec![2, 1])), Item::Cube);
        assert_eq!(at, Node::Shape(vec![3]));
        assert_eq!(kind(&d, &[3]), ShapeKind::Cube);
        assert_eq!(d.shapes[2], before);
    }

    #[test]
    fn a_group_inside_a_csg_still_takes_a_shape() {
        let mut d = scene();
        let at = add(&mut d, Some(&Node::Shape(vec![2, 0])), Item::Cube);
        assert_eq!(at, Node::Shape(vec![2, 0, 1]));
        assert_eq!(kind(&d, &[2, 0, 1]), ShapeKind::Cube);
    }

    #[test]
    fn with_no_shape_selected_a_shape_goes_at_the_end() {
        for selected in [
            None,
            Some(Node::Camera),
            Some(Node::Lights),
            Some(Node::Light(0)),
            Some(Node::Shapes),
        ] {
            let mut d = scene();
            let at = add(&mut d, selected.as_ref(), Item::Plane);
            assert_eq!(at, Node::Shape(vec![3]), "{selected:?}");
            assert_eq!(kind(&d, &[3]), ShapeKind::Plane);
        }
    }

    #[test]
    fn added_shapes_have_their_defaults_and_build() {
        let translate = vec![TransformOp::Translate([0., 1., 0.])];
        let cap = |closed_kind: fn(Option<f64>, Option<f64>, bool) -> ShapeKind| {
            closed_kind(Some(0.), Some(1.), true)
        };
        let cases = [
            (Item::Sphere, ShapeKind::Sphere, translate.clone()),
            (Item::Cube, ShapeKind::Cube, translate),
            (Item::Plane, ShapeKind::Plane, vec![]),
            (
                Item::Cylinder,
                cap(|min, max, closed| ShapeKind::Cylinder { min, max, closed }),
                vec![],
            ),
            (
                Item::Cone,
                cap(|min, max, closed| ShapeKind::Cone { min, max, closed }),
                vec![],
            ),
        ];
        for (item, expected, transform) in cases {
            let mut d = scene();
            let Node::Shape(path) = add(&mut d, None, item) else {
                panic!("{item:?}")
            };
            let added = shape_at(&d, &path).unwrap();
            assert_eq!(added.shape, expected, "{item:?}");
            assert_eq!(added.transform, transform, "{item:?}");
            assert_eq!(
                added.material,
                Some(MaterialDescription {
                    color: Some([0.8; 3]),
                    ..Default::default()
                }),
                "{item:?}"
            );
            assert!(d.build(Path::new(".")).is_ok(), "{item:?}");
        }
    }

    #[test]
    fn adding_a_light_to_a_single_light_makes_a_list() {
        let mut d = single_light_scene();
        let first = d.light.clone().unwrap();
        let at = add(&mut d, None, Item::Light);
        assert_eq!(at, Node::Light(1));
        assert_eq!(d.light, None);
        assert_eq!(d.lights, Some(vec![first, point_light()]));
    }

    #[test]
    fn a_light_goes_at_the_end_of_the_lights() {
        let mut d = scene();
        let at = add(&mut d, Some(&Node::Light(0)), Item::Light);
        assert_eq!(at, Node::Light(2));
        assert_eq!(d.lights.as_ref().unwrap().len(), 3);
        assert_eq!(d.lights.unwrap()[2], point_light());
    }

    #[test]
    fn item_keys() {
        let keys: Vec<_> = "scpyol".chars().map(Item::from_key).collect();
        assert_eq!(
            keys,
            [
                Some(Item::Sphere),
                Some(Item::Cube),
                Some(Item::Plane),
                Some(Item::Cylinder),
                Some(Item::Cone),
                Some(Item::Light),
            ]
        );
        assert_eq!(Item::from_key('x'), None);
    }

    fn folds(nodes: &[Node]) -> HashSet<Node> {
        nodes.iter().cloned().collect()
    }

    fn sh(path: &[usize]) -> Node {
        Node::Shape(path.to_vec())
    }

    #[test]
    fn deleting_a_shape_in_a_group_selects_the_next_one() {
        let mut d = scene();
        assert_eq!(delete(&mut d, &sh(&[1, 0])), Ok(sh(&[1, 0])));
        assert_eq!(children_len(&d, &[1]), 1);
        assert_eq!(kind(&d, &[1, 0]), ShapeKind::Plane);
    }

    #[test]
    fn deleting_the_last_shape_of_a_level_selects_the_previous_one() {
        let mut d = scene();
        assert_eq!(delete(&mut d, &sh(&[2])), Ok(sh(&[1])));
        assert_eq!(d.shapes.len(), 2);
    }

    #[test]
    fn deleting_an_only_child_selects_its_group() {
        let mut d = scene();
        d.shapes[1] = ShapeDescription {
            shape: ShapeKind::Group {
                children: vec![d.shapes[0].clone()],
            },
            ..d.shapes[1].clone()
        };
        assert_eq!(delete(&mut d, &sh(&[1, 0])), Ok(sh(&[1])));
        assert_eq!(children_len(&d, &[1]), 0);
        d.shapes.truncate(1);
        assert_eq!(delete(&mut d, &sh(&[0])), Ok(Node::Shapes));
        assert!(d.shapes.is_empty());
    }

    #[test]
    fn deleting_a_group_removes_its_children() {
        let mut d = scene();
        assert_eq!(delete(&mut d, &sh(&[1])), Ok(sh(&[1])));
        assert_eq!(d.shapes.len(), 2);
    }

    #[test]
    fn a_csg_child_cant_be_deleted() {
        let mut d = scene();
        let before = d.clone();
        let msg = "a CSG needs both its shapes".to_string();
        assert_eq!(can_delete(&d, &sh(&[2, 0])), Err(msg.clone()));
        assert_eq!(delete(&mut d, &sh(&[2, 0])), Err(msg));
        assert_eq!(d, before);
    }

    #[test]
    fn the_last_light_cant_be_deleted() {
        let mut d = single_light_scene();
        let before = d.clone();
        let msg = "a scene needs a light".to_string();
        assert_eq!(can_delete(&d, &Node::Light(0)), Err(msg.clone()));
        assert_eq!(delete(&mut d, &Node::Light(0)), Err(msg));
        assert_eq!(d, before);
    }

    #[test]
    fn deleting_a_light_selects_the_next_or_previous() {
        let mut d = scene();
        assert_eq!(delete(&mut d, &Node::Light(1)), Ok(Node::Light(0)));
        assert_eq!(d.lights.as_ref().unwrap().len(), 1);
        let mut d = scene();
        assert_eq!(delete(&mut d, &Node::Light(0)), Ok(Node::Light(0)));
        assert_eq!(d.lights.as_ref().unwrap().len(), 1);
    }

    #[test]
    fn heading_rows_cant_be_deleted() {
        for node in [Node::Camera, Node::Settings, Node::Lights, Node::Shapes] {
            let mut d = scene();
            let before = d.clone();
            assert!(can_delete(&d, &node).is_err(), "{node:?}");
            assert!(delete(&mut d, &node).is_err(), "{node:?}");
            assert_eq!(d, before);
        }
        let mut d = scene();
        assert!(delete(&mut d, &sh(&[9])).is_err());
    }

    #[test]
    fn a_folded_group_stays_folded_after_an_insert_before_it() {
        let mut f = folds(&[sh(&[1]), sh(&[1, 0])]);
        shift(&mut f, &[0], 1);
        assert_eq!(f, folds(&[sh(&[2]), sh(&[2, 0])]));
    }

    #[test]
    fn a_folded_group_stays_folded_after_a_delete_before_it() {
        let mut f = folds(&[sh(&[2])]);
        shift(&mut f, &[0], -1);
        assert_eq!(f, folds(&[sh(&[1])]));
    }

    #[test]
    fn deleting_a_folded_group_drops_its_folds() {
        let mut f = folds(&[sh(&[1]), sh(&[1, 0]), sh(&[2])]);
        shift(&mut f, &[1], -1);
        assert_eq!(f, folds(&[sh(&[1])]));
    }

    #[test]
    fn folds_elsewhere_dont_move() {
        let mut f = folds(&[Node::Lights, sh(&[0]), sh(&[1, 0])]);
        let before = f.clone();
        shift(&mut f, &[1, 1], 1);
        assert_eq!(f, before);
    }
}
