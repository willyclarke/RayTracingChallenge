//! The boxes the inspector shows for one item, and writing what is typed in
//! them back into the scene.

use crate::tree::{Node, shape_at};
use crate::trs::Trs;
use rtc_rust::scene::{
    LightDescription, MaterialDescription, SceneDescription, ShapeDescription, ShapeKind,
};

/// One box. The `usize` is the axis, 0 to 2.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Field {
    Position(usize),
    Rotation(usize),
    Scale(usize),
    Min,
    Max,
    Color(usize),
    Ambient,
    Diffuse,
    Specular,
    Shininess,
    Reflective,
    Transparency,
    RefractiveIndex,
    LightPosition(usize),
    LightIntensity(usize),
}

const MATERIAL: [Field; 10] = [
    Field::Color(0),
    Field::Color(1),
    Field::Color(2),
    Field::Ambient,
    Field::Diffuse,
    Field::Specular,
    Field::Shininess,
    Field::Reflective,
    Field::Transparency,
    Field::RefractiveIndex,
];

/// The boxes for `node` with their starting text, in display order; empty
/// for items that can't be edited.
pub fn fields(desc: &SceneDescription, node: &Node) -> Vec<(Field, String)> {
    match node {
        Node::Shape(path) => match shape_at(desc, path) {
            Some(shape) if editable(shape) => shape_fields(shape),
            _ => vec![],
        },
        Node::Light(i) => match light_at(desc, *i) {
            Some(LightDescription::Point {
                position,
                intensity,
            }) => (0..3)
                .map(|a| (Field::LightPosition(a), number(position[a])))
                .chain((0..3).map(|a| (Field::LightIntensity(a), number(intensity[a]))))
                .collect(),
            _ => vec![],
        },
        _ => vec![],
    }
}

fn shape_fields(shape: &ShapeDescription) -> Vec<(Field, String)> {
    let mut out = vec![];
    if let Some(t) = Trs::from_ops(&shape.transform) {
        for a in 0..3 {
            out.push((Field::Position(a), number(t.translate[a])));
        }
        for a in 0..3 {
            let degrees = (t.rotate[a].to_degrees() * 1e6).round() / 1e6;
            out.push((Field::Rotation(a), number(degrees)));
        }
        for a in 0..3 {
            out.push((Field::Scale(a), number(t.scale[a])));
        }
    }
    if let ShapeKind::Cylinder { min, max, .. } | ShapeKind::Cone { min, max, .. } = &shape.shape {
        out.push((Field::Min, optional(*min)));
        out.push((Field::Max, optional(*max)));
    }
    let m = shape.material.clone().unwrap_or_default();
    for field in MATERIAL {
        let value = match field {
            Field::Color(a) => m.color.map(|c| c[a]),
            _ => *material_value(&mut m.clone(), field),
        };
        out.push((field, optional(value)));
    }
    out
}

/// Parse `text` and write it into `field` of `node`. On `Err` the scene is
/// unchanged.
pub fn set(
    desc: &mut SceneDescription,
    node: &Node,
    field: Field,
    text: &str,
) -> Result<(), String> {
    let value = parse(text)?;
    match node {
        Node::Shape(path) => {
            let shape = shape_at_mut(desc, path)
                .filter(|s| editable(s))
                .ok_or("this item can't be edited")?;
            set_shape(shape, field, value)
        }
        Node::Light(i) => {
            let Some(LightDescription::Point {
                position,
                intensity,
            }) = light_at_mut(desc, *i)
            else {
                return Err("this light can't be edited".into());
            };
            let v = needed(value)?;
            match field {
                Field::LightPosition(a) => position[a] = v,
                Field::LightIntensity(a) => intensity[a] = v,
                _ => return Err(no_such(field)),
            }
            Ok(())
        }
        _ => Err("this item can't be edited".into()),
    }
}

fn set_shape(shape: &mut ShapeDescription, field: Field, value: Option<f64>) -> Result<(), String> {
    match field {
        Field::Position(a) | Field::Rotation(a) | Field::Scale(a) => {
            let mut t = Trs::from_ops(&shape.transform).ok_or("the transform has no fields")?;
            let v = needed(value)?;
            match field {
                Field::Position(_) => t.translate[a] = v,
                Field::Rotation(_) => t.rotate[a] = v.to_radians(),
                _ if v == 0.0 => return Err("a scale of 0 can't be inverted".into()),
                _ => t.scale[a] = v,
            }
            shape.transform = t.to_ops();
        }
        Field::Min | Field::Max => {
            let (ShapeKind::Cylinder { min, max, .. } | ShapeKind::Cone { min, max, .. }) =
                &mut shape.shape
            else {
                return Err(no_such(field));
            };
            *(if field == Field::Min { min } else { max }) = value;
        }
        Field::LightPosition(_) | Field::LightIntensity(_) => return Err(no_such(field)),
        _ => {
            let mut m = shape.material.clone().unwrap_or_default();
            if let Field::Color(a) = field {
                // An unset color is the default material's white.
                let mut c = m.color.unwrap_or([1.0; 3]);
                c[a] = needed(value)?;
                m.color = Some(c);
            } else {
                *material_value(&mut m, field) = value;
            }
            shape.material = (m != MaterialDescription::default()).then_some(m);
        }
    }
    Ok(())
}

/// A cylinder's or cone's `closed`; `None` for anything else.
pub fn closed(desc: &SceneDescription, node: &Node) -> Option<bool> {
    let Node::Shape(path) = node else { return None };
    match shape_at(desc, path)?.shape {
        ShapeKind::Cylinder { closed, .. } | ShapeKind::Cone { closed, .. } => Some(closed),
        _ => None,
    }
}

/// Set a cylinder's or cone's `closed`; anything else is left alone.
pub fn set_closed(desc: &mut SceneDescription, node: &Node, value: bool) {
    let Node::Shape(path) = node else { return };
    if let Some(ShapeDescription {
        shape: ShapeKind::Cylinder { closed, .. } | ShapeKind::Cone { closed, .. },
        ..
    }) = shape_at_mut(desc, path)
    {
        *closed = value;
    }
}

/// The material's optional number behind a scalar field.
fn material_value(m: &mut MaterialDescription, field: Field) -> &mut Option<f64> {
    match field {
        Field::Ambient => &mut m.ambient,
        Field::Diffuse => &mut m.diffuse,
        Field::Specular => &mut m.specular,
        Field::Shininess => &mut m.shininess,
        Field::Reflective => &mut m.reflective,
        Field::Transparency => &mut m.transparency,
        Field::RefractiveIndex => &mut m.refractive_index,
        _ => unreachable!("{field:?} is not a material number"),
    }
}

fn editable(shape: &ShapeDescription) -> bool {
    matches!(
        shape.shape,
        ShapeKind::Sphere
            | ShapeKind::Cube
            | ShapeKind::Plane
            | ShapeKind::Cylinder { .. }
            | ShapeKind::Cone { .. }
    )
}

/// `None` for an empty box; otherwise a finite number.
fn parse(text: &str) -> Result<Option<f64>, String> {
    let text = text.trim();
    if text.is_empty() {
        return Ok(None);
    }
    match text.parse::<f64>() {
        Ok(v) if v.is_finite() => Ok(Some(v)),
        _ => Err(format!("not a number: {text}")),
    }
}

fn needed(value: Option<f64>) -> Result<f64, String> {
    value.ok_or_else(|| "this box needs a number".into())
}

fn no_such(field: Field) -> String {
    format!("this item has no {field:?}")
}

fn number(v: f64) -> String {
    format!("{v}")
}

fn optional(v: Option<f64>) -> String {
    v.map(number).unwrap_or_default()
}

/// The lights, from `light` or `lights`, as `tree` lists them.
fn light_at(desc: &SceneDescription, i: usize) -> Option<&LightDescription> {
    match (&desc.light, &desc.lights) {
        (Some(light), _) => (i == 0).then_some(light),
        (None, Some(lights)) => lights.get(i),
        (None, None) => None,
    }
}

fn light_at_mut(desc: &mut SceneDescription, i: usize) -> Option<&mut LightDescription> {
    match (&mut desc.light, &mut desc.lights) {
        (Some(light), _) => (i == 0).then_some(light),
        (None, Some(lights)) => lights.get_mut(i),
        (None, None) => None,
    }
}

fn children_mut(shape: &mut ShapeDescription) -> Vec<&mut ShapeDescription> {
    match &mut shape.shape {
        ShapeKind::Group { children } => children.iter_mut().collect(),
        ShapeKind::Csg { left, right, .. } => vec![left.as_mut(), right.as_mut()],
        _ => vec![],
    }
}

fn shape_at_mut<'a>(
    desc: &'a mut SceneDescription,
    path: &[usize],
) -> Option<&'a mut ShapeDescription> {
    let (first, rest) = path.split_first()?;
    let mut shape = desc.shapes.get_mut(*first)?;
    for &i in rest {
        shape = children_mut(shape).into_iter().nth(i)?;
    }
    Some(shape)
}

#[cfg(test)]
mod tests {
    use super::*;
    use Field::*;
    use rtc_rust::scene::TransformOp;
    use serde_json::json;
    use std::f64::consts::FRAC_PI_2;
    use std::path::Path;

    fn scene() -> SceneDescription {
        let v = json!({
            "camera": {"width": 10, "height": 5, "field_of_view": 1.0,
                       "from": [0, 0, -5], "to": [0, 0, 0], "up": [0, 1, 0]},
            "lights": [
                {"type": "point", "position": [0, 0, -5], "intensity": [1, 1, 1]},
                {"type": "area", "corner": [0, 0, 0], "uvec": [1, 0, 0], "usteps": 2,
                 "vvec": [0, 1, 0], "vsteps": 2, "intensity": [1, 1, 1]}
            ],
            "shapes": [
                {"type": "sphere",
                 "transform": [{"scale": [2, 2, 2]}, {"translate": [0, 1, 0]}],
                 "material": {"color": [1, 0.2, 0.2]}},
                {"type": "cylinder", "min": -1, "max": 1, "closed": true},
                {"type": "group", "children": [{"type": "cube"}]},
                {"type": "torus"},
                {"type": "obj", "file": "missing.obj"},
                {"type": "sphere", "transform": [{"shear": [1, 0, 0, 0, 0, 0]}]},
                {"type": "sphere",
                 "transform": [{"rotate_x": 0.7854321}, {"rotate_z": 0.123456789}]}
            ]
        });
        SceneDescription::from_json(&v.to_string()).unwrap()
    }

    fn shape(i: usize) -> Node {
        Node::Shape(vec![i])
    }

    fn text(desc: &SceneDescription, node: &Node, field: Field) -> String {
        fields(desc, node)
            .into_iter()
            .find(|(f, _)| *f == field)
            .map(|(_, t)| t)
            .unwrap_or_else(|| panic!("no {field:?}"))
    }

    fn material(desc: &SceneDescription, i: usize) -> Option<MaterialDescription> {
        desc.shapes[i].material.clone()
    }

    #[test]
    fn a_sphere_lists_transform_and_material() {
        let desc = scene();
        let got = fields(&desc, &shape(0));
        let first: Vec<(Field, &str)> = got[..9].iter().map(|(f, t)| (*f, t.as_str())).collect();
        assert_eq!(
            first,
            [
                (Position(0), "0"),
                (Position(1), "1"),
                (Position(2), "0"),
                (Rotation(0), "0"),
                (Rotation(1), "0"),
                (Rotation(2), "0"),
                (Scale(0), "2"),
                (Scale(1), "2"),
                (Scale(2), "2"),
            ]
        );
        assert_eq!(text(&desc, &shape(0), Color(0)), "1");
        assert_eq!(text(&desc, &shape(0), Color(1)), "0.2");
        assert_eq!(text(&desc, &shape(0), Color(2)), "0.2");
        assert_eq!(text(&desc, &shape(0), Ambient), "");
        let rest: Vec<Field> = got[9..].iter().map(|(f, _)| *f).collect();
        assert_eq!(
            rest,
            [
                Color(0),
                Color(1),
                Color(2),
                Ambient,
                Diffuse,
                Specular,
                Shininess,
                Reflective,
                Transparency,
                RefractiveIndex
            ]
        );
    }

    #[test]
    fn a_cylinder_lists_its_bounds() {
        let desc = scene();
        assert_eq!(text(&desc, &shape(1), Min), "-1");
        assert_eq!(text(&desc, &shape(1), Max), "1");
        assert_eq!(closed(&desc, &shape(1)), Some(true));
        assert_eq!(closed(&desc, &shape(0)), None);
    }

    #[test]
    fn a_point_light_lists_position_and_intensity() {
        let desc = scene();
        let got = fields(&desc, &Node::Light(0));
        let want = [
            (LightPosition(0), "0"),
            (LightPosition(1), "0"),
            (LightPosition(2), "-5"),
            (LightIntensity(0), "1"),
            (LightIntensity(1), "1"),
            (LightIntensity(2), "1"),
        ];
        let got: Vec<(Field, &str)> = got.iter().map(|(f, t)| (*f, t.as_str())).collect();
        assert_eq!(got, want);
    }

    #[test]
    fn a_shape_whose_transform_doesnt_fit_lists_material_only() {
        let desc = scene();
        assert_eq!(fields(&desc, &shape(5))[0].0, Color(0));
    }

    #[test]
    fn items_you_cant_edit_list_nothing() {
        let desc = scene();
        for node in [
            shape(2),
            shape(3),
            shape(4),
            Node::Light(1),
            Node::Camera,
            Node::Settings,
            Node::Lights,
            Node::Shapes,
        ] {
            assert!(fields(&desc, &node).is_empty(), "{node:?}");
        }
    }

    #[test]
    fn setting_position_and_color_writes_the_scene() {
        let mut desc = scene();
        set(&mut desc, &shape(0), Position(0), "3").unwrap();
        set(&mut desc, &shape(0), Color(1), "0.5").unwrap();
        assert_eq!(
            desc.shapes[0].transform,
            [
                TransformOp::Scale([2.0, 2.0, 2.0]),
                TransformOp::Translate([3.0, 1.0, 0.0])
            ]
        );
        assert_eq!(material(&desc, 0).unwrap().color, Some([1.0, 0.5, 0.2]));
        desc.shapes.truncate(1);
        assert!(desc.build(Path::new(".")).is_ok());
    }

    #[test]
    fn rotation_is_typed_in_degrees() {
        let mut desc = scene();
        set(&mut desc, &shape(0), Rotation(1), "90").unwrap();
        assert!(
            desc.shapes[0]
                .transform
                .contains(&TransformOp::RotateY(FRAC_PI_2))
        );
        assert_eq!(text(&desc, &shape(0), Rotation(1)), "90");
    }

    #[test]
    fn editing_one_rotation_keeps_the_others_exact() {
        let mut desc = scene();
        set(&mut desc, &shape(6), Rotation(1), "10").unwrap();
        let t = &desc.shapes[6].transform;
        assert_eq!(t[0], TransformOp::RotateX(0.7854321));
        assert_eq!(t[2], TransformOp::RotateZ(0.123456789));
    }

    #[test]
    fn an_empty_material_box_unsets_the_value() {
        let mut desc = scene();
        set(&mut desc, &shape(0), Ambient, "0.3").unwrap();
        assert_eq!(material(&desc, 0).unwrap().ambient, Some(0.3));
        set(&mut desc, &shape(0), Ambient, "").unwrap();
        assert_eq!(material(&desc, 0).unwrap().ambient, None);
    }

    #[test]
    fn setting_a_value_creates_a_material() {
        let mut desc = scene();
        set(&mut desc, &shape(1), Reflective, "0.5").unwrap();
        let want = MaterialDescription {
            reflective: Some(0.5),
            ..Default::default()
        };
        assert_eq!(material(&desc, 1), Some(want));
    }

    #[test]
    fn clearing_the_last_material_value_removes_the_material() {
        let mut desc = scene();
        set(&mut desc, &shape(1), Reflective, "0.5").unwrap();
        set(&mut desc, &shape(1), Reflective, "").unwrap();
        assert_eq!(material(&desc, 1), None);
    }

    #[test]
    fn set_rejects_text_that_isnt_a_finite_number() {
        let mut desc = scene();
        let before = desc.clone();
        for bad in ["abc", "inf", "NaN", "1e400", ""] {
            assert!(
                set(&mut desc, &shape(0), Position(0), bad).is_err(),
                "{bad}"
            );
            assert_eq!(desc, before, "{bad}");
        }
    }

    #[test]
    fn a_zero_scale_is_rejected() {
        let mut desc = scene();
        let before = desc.clone();
        assert!(set(&mut desc, &shape(0), Scale(2), "0").is_err());
        assert_eq!(desc, before);
    }

    #[test]
    fn a_field_the_item_doesnt_have_is_rejected() {
        let mut desc = scene();
        let before = desc.clone();
        assert!(set(&mut desc, &shape(0), Min, "1").is_err());
        assert!(set(&mut desc, &shape(5), Position(0), "1").is_err());
        assert!(set(&mut desc, &Node::Light(1), LightPosition(0), "1").is_err());
        assert_eq!(desc, before);
    }

    #[test]
    fn a_shape_in_a_group_is_editable() {
        let mut desc = scene();
        set(&mut desc, &Node::Shape(vec![2, 0]), Position(0), "5").unwrap();
        let ShapeKind::Group { children } = &desc.shapes[2].shape else {
            panic!("not a group");
        };
        assert_eq!(
            children[0].transform,
            [TransformOp::Translate([5.0, 0.0, 0.0])]
        );
    }

    #[test]
    fn a_single_light_is_editable() {
        let mut desc = scene();
        desc.light = desc.lights.take().map(|mut l| l.remove(0));
        set(&mut desc, &Node::Light(0), LightPosition(1), "7").unwrap();
        assert_eq!(
            desc.light,
            Some(LightDescription::Point {
                position: [0.0, 7.0, -5.0],
                intensity: [1.0, 1.0, 1.0]
            })
        );
    }

    #[test]
    fn set_closed_toggles_the_caps() {
        let mut desc = scene();
        set_closed(&mut desc, &shape(1), false);
        assert_eq!(closed(&desc, &shape(1)), Some(false));
    }

    #[test]
    fn whitespace_around_a_number_is_accepted() {
        let mut desc = scene();
        set(&mut desc, &shape(1), Max, " 2 ").unwrap();
        assert_eq!(text(&desc, &shape(1), Max), "2");
    }
}
