//! The boxes the inspector shows for one item, and writing what is typed in
//! them back into the scene.

use crate::tree::{Node, shape_at};
use crate::trs::Trs;
use rtc_rust::scene::{
    Antialias, CameraDescription, FocalBlur, LightDescription, MaterialDescription,
    SceneDescription, ShapeDescription, ShapeKind,
};
use std::mem::discriminant;

/// One box. The `usize` is the axis, 0 to 2.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Field {
    Comment,
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
    CameraFrom(usize),
    CameraTo(usize),
    CameraUp(usize),
    FieldOfView,
    Width,
    Height,
    Antialias,
    FNumber,
    FocalDistance,
    FocalSamples,
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
        Node::Shape(path) => {
            let Some(shape) = shape_at(desc, path) else {
                return vec![];
            };
            let comment = match &shape.comment {
                None => (Field::Comment, String::new()),
                Some(serde_json::Value::String(s)) => (Field::Comment, s.clone()),
                // Not text: no box, so it is never overwritten.
                Some(_) => {
                    return if editable(shape) {
                        shape_fields(shape)
                    } else {
                        vec![]
                    };
                }
            };
            let mut out = vec![comment];
            if editable(shape) {
                out.extend(shape_fields(shape));
            }
            out
        }
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
        Node::Camera => camera_fields(&desc.camera),
        _ => vec![],
    }
}

fn camera_fields(c: &CameraDescription) -> Vec<(Field, String)> {
    let mut out = vec![];
    for a in 0..3 {
        out.push((Field::CameraFrom(a), number(c.from[a])));
    }
    for a in 0..3 {
        out.push((Field::CameraTo(a), number(c.to[a])));
    }
    for a in 0..3 {
        out.push((Field::CameraUp(a), number(c.up[a])));
    }
    out.push((Field::FieldOfView, number(degrees(c.field_of_view))));
    out.push((Field::Width, c.width.to_string()));
    out.push((Field::Height, c.height.to_string()));
    out.push((
        Field::Antialias,
        c.antialias
            .as_ref()
            .map(|a| a.n.to_string())
            .unwrap_or_default(),
    ));
    let blur = c.focal_blur.as_ref();
    out.push((Field::FNumber, blur.map(f_number).unwrap_or_default()));
    out.push((
        Field::FocalDistance,
        optional(blur.map(|b| b.focal_distance)),
    ));
    out.push((
        Field::FocalSamples,
        blur.map(|b| b.samples.to_string()).unwrap_or_default(),
    ));
    out
}

/// The f-number a focal blur's aperture (the lens radius) stands for, rounded
/// to 6 decimals; `∞` for a lens of radius 0.
fn f_number(blur: &FocalBlur) -> String {
    if blur.aperture == 0.0 {
        return "∞".into();
    }
    number((blur.focal_distance / (2.0 * blur.aperture) * 1e6).round() / 1e6)
}

/// The camera's lens radius, when it has focal blur.
pub fn lens_radius(desc: &SceneDescription) -> Option<f64> {
    desc.camera.focal_blur.as_ref().map(|b| b.aperture)
}

/// Radians as the degrees a box shows, rounded to 6 decimals.
fn degrees(radians: f64) -> f64 {
    (radians.to_degrees() * 1e6).round() / 1e6
}

fn shape_fields(shape: &ShapeDescription) -> Vec<(Field, String)> {
    let mut out = vec![];
    if let Some(t) = Trs::from_ops(&shape.transform) {
        for a in 0..3 {
            out.push((Field::Position(a), number(t.translate[a])));
        }
        for a in 0..3 {
            out.push((Field::Rotation(a), number(degrees(t.rotate[a]))));
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
    if let Node::Camera = node {
        return set_camera(&mut desc.camera, field, text);
    }
    if let (Node::Shape(path), Field::Comment) = (node, field) {
        let shape = shape_at_mut(desc, path).ok_or("this item can't be edited")?;
        shape.comment = (!text.trim().is_empty()).then(|| text.into());
        return Ok(());
    }
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

/// Focus `blur` at `distance`, keeping its f-number, so the aperture follows
/// the focus.
pub fn set_focal_distance(blur: &mut FocalBlur, distance: f64) {
    let f = blur.focal_distance / (2.0 * blur.aperture);
    if f.is_finite() && f > 0.0 {
        blur.aperture = distance / (2.0 * f);
    }
    blur.focal_distance = distance;
}

/// `set`, except that the focal blur an emptied f-number removes goes to
/// `removed`, and the next f-number typed while there is none puts its focal
/// distance and samples back, so retyping the f-number doesn't lose them.
pub fn set_remembering(
    desc: &mut SceneDescription,
    node: &Node,
    field: Field,
    text: &str,
    removed: &mut Option<FocalBlur>,
) -> Result<(), String> {
    if !(matches!(node, Node::Camera) && field == Field::FNumber) {
        return set(desc, node, field, text);
    }
    let c = &mut desc.camera;
    let before = c.focal_blur.clone();
    if before.is_none() {
        c.focal_blur = removed.clone();
    }
    match set_camera(c, field, text) {
        Ok(()) => {
            if c.focal_blur.is_some() {
                *removed = None;
            } else if before.is_some() {
                *removed = before;
            }
            Ok(())
        }
        Err(e) => {
            c.focal_blur = before;
            Err(e)
        }
    }
}

/// The samples of a focal blur made by typing an aperture; its focal distance
/// is from the eye to the point looked at.
const FOCAL_SAMPLES: usize = 32;

fn set_camera(c: &mut CameraDescription, field: Field, text: &str) -> Result<(), String> {
    match field {
        Field::CameraFrom(a) => c.from[a] = needed(parse(text)?)?,
        Field::CameraTo(a) => c.to[a] = needed(parse(text)?)?,
        Field::CameraUp(a) => c.up[a] = needed(parse(text)?)?,
        Field::FieldOfView => {
            let v = needed(parse(text)?)?;
            if v <= 0.0 || v >= 180.0 {
                return Err("the field of view must be above 0 and below 180".into());
            }
            c.field_of_view = v.to_radians();
        }
        Field::Width => c.width = whole(text)?,
        Field::Height => c.height = whole(text)?,
        Field::Antialias => {
            c.antialias = match text.trim() {
                "" => None,
                _ => Some(Antialias { n: whole(text)? }),
            }
        }
        Field::FNumber => match parse(text)? {
            None => c.focal_blur = None,
            Some(f) if f <= 0.0 => return Err("the f-number must be above 0".into()),
            Some(f) => match &mut c.focal_blur {
                Some(blur) => blur.aperture = blur.focal_distance / (2.0 * f),
                None => {
                    let distance = c
                        .from
                        .iter()
                        .zip(&c.to)
                        .map(|(f, t)| (t - f).powi(2))
                        .sum::<f64>()
                        .sqrt();
                    c.focal_blur = Some(FocalBlur {
                        aperture: distance / (2.0 * f),
                        focal_distance: distance,
                        samples: FOCAL_SAMPLES,
                    });
                }
            },
        },
        Field::FocalDistance | Field::FocalSamples => {
            let Some(blur) = &mut c.focal_blur else {
                return Err("type an f-number first".into());
            };
            if field == Field::FocalDistance {
                let distance = needed(parse(text)?)?;
                if distance <= 0.0 {
                    return Err("the focal distance must be above 0".into());
                }
                set_focal_distance(blur, distance);
            } else {
                blur.samples = whole(text)?;
            }
        }
        _ => return Err(no_such(field)),
    }
    Ok(())
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

/// Whether `new` is `old` with one run of text removed: a cut, the only
/// change a box makes while Cmd is held that isn't a shortcut's letter
/// typed into it (paste arrives separately).
pub fn is_cut(old: &str, new: &str) -> bool {
    if new.len() >= old.len() {
        return false;
    }
    let same = old
        .char_indices()
        .zip(new.chars())
        .find(|((_, a), b)| a != b)
        .map_or(new.len(), |((i, _), _)| i);
    old[same..].ends_with(&new[same..])
}

/// The box an arrow moves the shading to, `right` along a row and `down`
/// across rows of the grid the inspector draws; it stops at the edges.
/// A row is one variant's fields; going up or down keeps the axis, or
/// takes the row's last box when the row is shorter.
pub fn next_box(order: &[Field], from: Field, right: isize, down: isize) -> Field {
    // The rows, each a run of fields of one variant, in display order.
    let mut grid: Vec<Vec<Field>> = vec![];
    for &f in order {
        match grid.last_mut() {
            Some(row) if discriminant(&row[0]) == discriminant(&f) => row.push(f),
            _ => grid.push(vec![f]),
        }
    }
    let Some((r, c)) = grid
        .iter()
        .enumerate()
        .find_map(|(r, row)| row.iter().position(|f| *f == from).map(|c| (r, c)))
    else {
        return from;
    };
    let clamp = |i: isize, len: usize| i.clamp(0, len as isize - 1) as usize;
    let row = &grid[clamp(r as isize + down, grid.len())];
    row[clamp(c as isize + right, row.len())]
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

/// A whole number of at least 1.
fn whole(text: &str) -> Result<usize, String> {
    match text.trim().parse::<usize>() {
        Ok(n) if n >= 1 => Ok(n),
        _ => Err(format!("not a whole number of at least 1: {}", text.trim())),
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
    use rtc_rust::scene::{FocalBlur, TransformOp};
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
        let first: Vec<(Field, &str)> = got[1..10].iter().map(|(f, t)| (*f, t.as_str())).collect();
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
        let rest: Vec<Field> = got[10..].iter().map(|(f, _)| *f).collect();
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
        assert_eq!(fields(&desc, &shape(5))[1].0, Color(0));
    }

    #[test]
    fn items_you_cant_edit_list_nothing() {
        let desc = scene();
        for node in [Node::Light(1), Node::Settings, Node::Lights, Node::Shapes] {
            assert!(fields(&desc, &node).is_empty(), "{node:?}");
        }
    }

    fn camera(desc: &SceneDescription) -> &rtc_rust::scene::CameraDescription {
        &desc.camera
    }

    fn set_camera(desc: &mut SceneDescription, field: Field, text: &str) -> Result<(), String> {
        set(desc, &Node::Camera, field, text)
    }

    #[test]
    fn the_camera_lists_its_boxes() {
        let desc = scene();
        let c = |f| text(&desc, &Node::Camera, f);
        assert_eq!(c(CameraFrom(2)), "-5");
        assert_eq!(c(Width), "10");
        assert_eq!(c(FieldOfView), "57.29578");
        for f in [Antialias, FNumber, FocalDistance, FocalSamples] {
            assert_eq!(c(f), "", "{f:?}");
        }
        let order: Vec<Field> = fields(&desc, &Node::Camera)
            .into_iter()
            .map(|(f, _)| f)
            .collect();
        assert_eq!(
            order,
            [
                CameraFrom(0),
                CameraFrom(1),
                CameraFrom(2),
                CameraTo(0),
                CameraTo(1),
                CameraTo(2),
                CameraUp(0),
                CameraUp(1),
                CameraUp(2),
                FieldOfView,
                Width,
                Height,
                Antialias,
                FNumber,
                FocalDistance,
                FocalSamples
            ]
        );
    }

    #[test]
    fn field_of_view_is_typed_in_degrees() {
        let mut desc = scene();
        set_camera(&mut desc, FieldOfView, "90").unwrap();
        assert_eq!(camera(&desc).field_of_view, FRAC_PI_2);
    }

    #[test]
    fn editing_another_camera_box_keeps_the_field_of_view_exact() {
        let mut desc = scene();
        set_camera(&mut desc, Width, "20").unwrap();
        assert_eq!(camera(&desc).field_of_view, 1.0);
        assert_eq!(camera(&desc).width, 20);
    }

    #[test]
    fn camera_vectors_are_edited_per_axis() {
        let mut desc = scene();
        set_camera(&mut desc, CameraTo(1), "2").unwrap();
        set_camera(&mut desc, CameraUp(0), "0").unwrap();
        assert_eq!(camera(&desc).to, [0.0, 2.0, 0.0]);
        assert_eq!(camera(&desc).up, [0.0, 1.0, 0.0]);
    }

    #[test]
    fn camera_values_out_of_range_are_rejected() {
        let mut desc = scene();
        let before = desc.clone();
        for (field, bad) in [
            (FieldOfView, "0"),
            (FieldOfView, "180"),
            (Width, "0"),
            (Width, "2.5"),
            (Height, "-1"),
            (Antialias, "0"),
            (CameraFrom(0), ""),
        ] {
            assert!(
                set_camera(&mut desc, field, bad).is_err(),
                "{field:?} {bad}"
            );
            assert_eq!(desc, before, "{field:?} {bad}");
        }
    }

    #[test]
    fn typing_n_adds_antialias_and_emptying_it_removes_it() {
        let mut desc = scene();
        set_camera(&mut desc, Antialias, "3").unwrap();
        assert_eq!(
            camera(&desc).antialias,
            Some(rtc_rust::scene::Antialias { n: 3 })
        );
        assert_eq!(text(&desc, &Node::Camera, Antialias), "3");
        set_camera(&mut desc, Antialias, "").unwrap();
        assert_eq!(camera(&desc).antialias, None);
    }

    #[test]
    fn an_f_number_creates_focal_blur() {
        let mut desc = scene();
        set_camera(&mut desc, FNumber, "8").unwrap();
        assert_eq!(
            camera(&desc).focal_blur,
            Some(FocalBlur {
                aperture: 0.3125,
                focal_distance: 5.0,
                samples: 32
            })
        );
        assert_eq!(text(&desc, &Node::Camera, FNumber), "8");
        assert_eq!(text(&desc, &Node::Camera, FocalDistance), "5");
        assert_eq!(text(&desc, &Node::Camera, FocalSamples), "32");
        assert_eq!(lens_radius(&desc), Some(0.3125));
    }

    #[test]
    fn an_f_number_with_from_equal_to_to_gives_distance_zero() {
        let mut desc = scene();
        desc.camera.to = desc.camera.from;
        set_camera(&mut desc, FNumber, "8").unwrap();
        let blur = camera(&desc).focal_blur.clone().unwrap();
        assert_eq!((blur.focal_distance, blur.aperture), (0.0, 0.0));
    }

    #[test]
    fn emptying_the_f_number_removes_focal_blur() {
        let mut desc = scene();
        set_camera(&mut desc, FNumber, "8").unwrap();
        set_camera(&mut desc, FNumber, "").unwrap();
        assert_eq!(camera(&desc).focal_blur, None);
        for f in [FNumber, FocalDistance, FocalSamples] {
            assert_eq!(text(&desc, &Node::Camera, f), "", "{f:?}");
        }
        assert_eq!(lens_radius(&desc), None);
    }

    #[test]
    fn emptying_then_retyping_the_f_number_keeps_focal_distance_and_samples() {
        let mut desc = scene();
        desc.camera.focal_blur = Some(FocalBlur {
            aperture: 3.5 / 16.0,
            focal_distance: 3.5,
            samples: 8,
        });
        let mut removed = None;
        let mut put = |desc: &mut SceneDescription, text| {
            set_remembering(desc, &Node::Camera, FNumber, text, &mut removed)
        };
        put(&mut desc, "").unwrap();
        assert_eq!(camera(&desc).focal_blur, None);
        put(&mut desc, "abc").unwrap_err();
        assert_eq!(camera(&desc).focal_blur, None);
        put(&mut desc, "").unwrap();
        put(&mut desc, "16").unwrap();
        assert_eq!(
            camera(&desc).focal_blur,
            Some(FocalBlur {
                aperture: 3.5 / 32.0,
                focal_distance: 3.5,
                samples: 8
            })
        );
        assert_eq!(text(&desc, &Node::Camera, FocalDistance), "3.5");
        assert_eq!(text(&desc, &Node::Camera, FocalSamples), "8");
    }

    #[test]
    fn retyping_without_a_remembered_blur_derives_it() {
        let mut desc = scene();
        set_remembering(&mut desc, &Node::Camera, FNumber, "8", &mut None).unwrap();
        assert_eq!(camera(&desc).focal_blur.clone().unwrap().samples, 32);
    }

    #[test]
    fn changing_the_focal_distance_keeps_the_f_number() {
        let mut desc = scene();
        set_camera(&mut desc, FNumber, "8").unwrap();
        set_camera(&mut desc, FocalDistance, "10").unwrap();
        assert_eq!(camera(&desc).focal_blur.as_ref().unwrap().aperture, 0.625);
        assert_eq!(text(&desc, &Node::Camera, FNumber), "8");
    }

    #[test]
    fn the_f_number_must_be_above_zero() {
        let mut desc = scene();
        let before = desc.clone();
        for bad in ["0", "-2"] {
            assert!(set_camera(&mut desc, FNumber, bad).is_err(), "{bad}");
            assert_eq!(desc, before, "{bad}");
        }
    }

    #[test]
    fn an_f_number_of_zero_shows_an_infinite_f_number() {
        let mut desc = scene();
        desc.camera.focal_blur = Some(FocalBlur {
            aperture: 0.0,
            focal_distance: 5.0,
            samples: 8,
        });
        assert_eq!(text(&desc, &Node::Camera, FNumber), "∞");
    }

    #[test]
    fn focal_distance_and_samples_edit_the_existing_focal_blur() {
        let mut desc = scene();
        set_camera(&mut desc, FNumber, "8").unwrap();
        set_camera(&mut desc, FocalDistance, "7.5").unwrap();
        set_camera(&mut desc, FocalSamples, "8").unwrap();
        let blur = camera(&desc).focal_blur.clone().unwrap();
        assert_eq!((blur.focal_distance, blur.samples), (7.5, 8));
        assert!(set_camera(&mut desc, FocalSamples, "0").is_err());
    }

    #[test]
    fn typing_a_focal_distance_through_zero_keeps_the_f_number() {
        let mut desc = scene();
        set_camera(&mut desc, FNumber, "8").unwrap();
        let before = desc.clone();
        assert!(set_camera(&mut desc, FocalDistance, "0").is_err());
        assert!(set_camera(&mut desc, FocalDistance, "-1").is_err());
        assert_eq!(desc, before);
        set_camera(&mut desc, FocalDistance, "0.5").unwrap();
        assert_eq!(camera(&desc).focal_blur.as_ref().unwrap().aperture, 0.03125);
        assert_eq!(text(&desc, &Node::Camera, FNumber), "8");
    }

    #[test]
    fn focal_distance_needs_a_number_while_focal_blur_exists() {
        let mut desc = scene();
        set_camera(&mut desc, FNumber, "8").unwrap();
        let before = desc.clone();
        assert!(set_camera(&mut desc, FocalDistance, "").is_err());
        assert!(set_camera(&mut desc, FocalSamples, "").is_err());
        assert_eq!(desc, before);
    }

    #[test]
    fn focal_distance_and_samples_do_nothing_without_focal_blur() {
        let mut desc = scene();
        let before = desc.clone();
        assert!(set_camera(&mut desc, FocalDistance, "3").is_err());
        assert!(set_camera(&mut desc, FocalSamples, "3").is_err());
        assert_eq!(desc, before);
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

    fn cylinder_order() -> Vec<Field> {
        fields(&scene(), &shape(1))
            .into_iter()
            .map(|(f, _)| f)
            .collect()
    }

    #[test]
    fn next_box_moves_along_a_vector_and_stops_at_its_end() {
        let order = cylinder_order();
        assert_eq!(next_box(&order, Position(0), 1, 0), Position(1));
        assert_eq!(next_box(&order, Position(2), 1, 0), Position(2));
        assert_eq!(next_box(&order, Position(0), -1, 0), Position(0));
    }

    #[test]
    fn next_box_keeps_the_axis_going_down_and_clamps_it_on_a_single_box() {
        let order = cylinder_order();
        assert_eq!(next_box(&order, Position(2), 0, 1), Rotation(2));
        assert_eq!(next_box(&order, Scale(2), 0, 1), Min);
        assert_eq!(next_box(&order, Max, 0, 1), Color(0));
        assert_eq!(next_box(&order, Ambient, 0, -1), Color(0));
    }

    #[test]
    fn only_a_cut_counts_as_typing_while_cmd_is_held() {
        assert!(is_cut("1.5", "1"));
        assert!(is_cut("1.5", "15"));
        assert!(is_cut("1.5", ""));
        assert!(!is_cut("1.5", "1.5l"));
        assert!(!is_cut("1.5", "s"));
        assert!(!is_cut("1.5", "1.s"));
        assert!(!is_cut("1.5", "1.5"));
    }

    #[test]
    fn next_box_stops_at_the_first_and_last_row() {
        let order = cylinder_order();
        assert_eq!(next_box(&order, Comment, 0, -1), Comment);
        assert_eq!(next_box(&order, Position(1), 0, -1), Comment);
        assert_eq!(next_box(&order, Comment, 0, 1), Position(0));
        assert_eq!(next_box(&order, RefractiveIndex, 0, 1), RefractiveIndex);
    }

    fn commented_scene() -> SceneDescription {
        let v = json!({
            "camera": {"width": 10, "height": 5, "field_of_view": 1.0,
                       "from": [0, 0, -5], "to": [0, 0, 0], "up": [0, 1, 0]},
            "light": {"type": "point", "position": [0, 0, -5], "intensity": [1, 1, 1]},
            "shapes": [
                {"type": "sphere", "comment": "ball"},
                {"type": "sphere"},
                {"type": "group", "comment": "stack", "children": [{"type": "cube"}]},
                {"type": "csg", "operation": "union",
                 "left": {"type": "sphere"}, "right": {"type": "cube"}},
                {"type": "sphere", "comment": ["a", "b"]}
            ]
        });
        SceneDescription::from_json(&v.to_string()).unwrap()
    }

    #[test]
    fn a_shape_lists_its_comment_first() {
        let desc = commented_scene();
        assert_eq!(fields(&desc, &shape(0))[0], (Comment, "ball".to_string()));
        assert_eq!(fields(&desc, &shape(1))[0], (Comment, String::new()));
    }

    #[test]
    fn a_group_and_a_csg_list_only_a_comment_box() {
        let desc = commented_scene();
        assert_eq!(
            fields(&desc, &shape(2)),
            vec![(Comment, "stack".to_string())]
        );
        assert_eq!(fields(&desc, &shape(3)), vec![(Comment, String::new())]);
    }

    #[test]
    fn a_non_string_comment_gets_no_box() {
        let desc = commented_scene();
        assert!(fields(&desc, &shape(4)).iter().all(|(f, _)| *f != Comment));
    }

    #[test]
    fn typing_a_comment_sets_it() {
        let mut desc = commented_scene();
        set(&mut desc, &shape(0), Comment, "red ball").unwrap();
        assert_eq!(desc.shapes[0].comment, Some(json!("red ball")));
        set(&mut desc, &shape(2), Comment, " a  stack ").unwrap();
        assert_eq!(desc.shapes[2].comment, Some(json!(" a  stack ")));
    }

    #[test]
    fn an_empty_comment_removes_it() {
        let mut desc = commented_scene();
        for text in ["", "   "] {
            set(&mut desc, &shape(0), Comment, "x").unwrap();
            set(&mut desc, &shape(0), Comment, text).unwrap();
            assert_eq!(desc.shapes[0].comment, None);
        }
        let saved = serde_json::to_string(&desc.shapes[0]).unwrap();
        assert!(!saved.contains("comment"), "{saved}");
    }

    #[test]
    fn lights_and_the_camera_have_no_comment_box() {
        let desc = commented_scene();
        for node in [Node::Light(0), Node::Camera] {
            assert!(fields(&desc, &node).iter().all(|(f, _)| *f != Comment));
        }
    }
}
