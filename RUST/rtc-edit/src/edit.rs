//! Making a new scene and adding shapes and lights to one, as pure edits of
//! a `SceneDescription`.

use crate::tree::{Node, shape_at};
use rtc_rust::scene::{
    CameraDescription, LightDescription, MaterialDescription, SceneDescription, ShapeDescription,
    ShapeKind, TransformOp,
};

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
}
