//! The scene as rows for the editor's tree, and the JSON its inspector
//! shows for one item.

use rtc_rust::scene::{LightDescription, SceneDescription, Settings, ShapeDescription, ShapeKind};
use serde_json::Value;
use std::collections::HashSet;

/// One item in the tree. `Lights` and `Shapes` are the section rows. A
/// `Shape` path goes top-level index first; in a group the child index
/// follows, and in a CSG 0 is left and 1 is right.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum Node {
    Camera,
    Settings,
    Lights,
    Light(usize),
    Shapes,
    Shape(Vec<usize>),
}

#[derive(Debug, Clone, PartialEq)]
pub struct Row {
    pub node: Node,
    pub depth: usize,
    pub label: String,
    pub foldable: bool,
}

/// The tree's rows, top to bottom; a folded node keeps its row and hides
/// everything under it.
pub fn rows(desc: &SceneDescription, folded: &HashSet<Node>) -> Vec<Row> {
    let mut rows = vec![row(Node::Camera, 0, "camera".into(), false)];
    if desc.settings != Settings::default() {
        rows.push(row(Node::Settings, 0, "settings".into(), false));
    }
    let lights = lights(desc);
    rows.push(row(
        Node::Lights,
        0,
        format!("lights ({})", lights.len()),
        true,
    ));
    if !folded.contains(&Node::Lights) {
        for (i, light) in lights.iter().enumerate() {
            rows.push(row(
                Node::Light(i),
                1,
                type_tag(serde_json::to_value(light)),
                false,
            ));
        }
    }
    rows.push(row(
        Node::Shapes,
        0,
        format!("shapes ({})", desc.shapes.len()),
        true,
    ));
    if !folded.contains(&Node::Shapes) {
        for (i, shape) in desc.shapes.iter().enumerate() {
            shape_rows(shape, vec![i], folded, &mut rows);
        }
    }
    rows
}

fn row(node: Node, depth: usize, label: String, foldable: bool) -> Row {
    Row {
        node,
        depth,
        label,
        foldable,
    }
}

/// The scene's lights, from `light` or `lights`.
fn lights(desc: &SceneDescription) -> Vec<&LightDescription> {
    match (&desc.light, &desc.lights) {
        (Some(light), _) => vec![light],
        (None, Some(lights)) => lights.iter().collect(),
        (None, None) => vec![],
    }
}

/// The `"type"` tag of an item serialized to JSON.
fn type_tag(item: serde_json::Result<Value>) -> String {
    match item {
        Ok(Value::Object(m)) => m
            .get("type")
            .and_then(Value::as_str)
            .unwrap_or("?")
            .to_string(),
        _ => "?".to_string(),
    }
}

/// A shape's children: a group's in order, a CSG's left then right.
fn children(shape: &ShapeDescription) -> Vec<&ShapeDescription> {
    match &shape.shape {
        ShapeKind::Group { children } => children.iter().collect(),
        ShapeKind::Csg { left, right, .. } => vec![left, right],
        _ => vec![],
    }
}

fn shape_rows(
    shape: &ShapeDescription,
    path: Vec<usize>,
    folded: &HashSet<Node>,
    rows: &mut Vec<Row>,
) {
    let base = match &shape.shape {
        ShapeKind::Group { children } => format!("group ({} children)", children.len()),
        ShapeKind::Csg { operation, .. } => {
            format!(
                "csg {}",
                serde_json::to_value(operation)
                    .ok()
                    .and_then(|v| v.as_str().map(String::from))
                    .unwrap_or_default()
            )
        }
        _ => type_tag(serde_json::to_value(shape)),
    };
    let label = match &shape.comment {
        Some(Value::String(c)) => format!("{base} — {c}"),
        _ => base,
    };
    let kids = children(shape);
    let node = Node::Shape(path.clone());
    let open = !folded.contains(&node);
    rows.push(row(node, path.len(), label, !kids.is_empty()));
    if open {
        for (i, kid) in kids.into_iter().enumerate() {
            let mut p = path.clone();
            p.push(i);
            shape_rows(kid, p, folded, rows);
        }
    }
}

/// The shape a `Node::Shape` path names.
fn shape_at<'a>(desc: &'a SceneDescription, path: &[usize]) -> Option<&'a ShapeDescription> {
    let (first, rest) = path.split_first()?;
    rest.iter().try_fold(desc.shapes.get(*first)?, |shape, &i| {
        children(shape).get(i).copied()
    })
}

/// The inspector's text for `node`: the item as pretty JSON, without a
/// group's or CSG's children; empty for the section rows.
pub fn json(desc: &SceneDescription, node: &Node) -> String {
    let value = match node {
        Node::Camera => serde_json::to_value(&desc.camera),
        Node::Settings => serde_json::to_value(&desc.settings),
        Node::Light(i) => match lights(desc).get(*i) {
            Some(light) => serde_json::to_value(light),
            None => return String::new(),
        },
        Node::Shape(path) => match shape_at(desc, path) {
            Some(shape) => serde_json::to_value(shape),
            None => return String::new(),
        },
        Node::Lights | Node::Shapes => return String::new(),
    };
    let Ok(mut value) = value else {
        return String::new();
    };
    if let (Node::Shape(_), Value::Object(m)) = (node, &mut value) {
        for key in ["children", "left", "right"] {
            m.remove(key);
        }
    }
    serde_json::to_string_pretty(&value).unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::{Value, json};

    fn scene() -> SceneDescription {
        let v = json!({
            "camera": {"width": 10, "height": 5, "field_of_view": 1.0,
                       "from": [0, 0, -5], "to": [0, 0, 0], "up": [0, 1, 0]},
            "lights": [
                {"type": "point", "position": [0, 0, -5], "intensity": [1, 1, 1]},
                {"type": "point", "position": [5, 0, -5], "intensity": [1, 1, 1]}
            ],
            "settings": {"path_tracing": {"samples": 4, "depth": 2}},
            "shapes": [
                {"type": "sphere", "comment": "ball"},
                {"type": "group", "children": [{"type": "cube"}, {"type": "plane"}]},
                {"type": "csg", "operation": "difference",
                 "left": {"type": "sphere"}, "right": {"type": "sphere"}}
            ]
        });
        SceneDescription::from_json(&v.to_string()).unwrap()
    }

    fn labels(rows: &[Row]) -> Vec<(&str, usize)> {
        rows.iter().map(|r| (r.label.as_str(), r.depth)).collect()
    }

    #[test]
    fn rows_list_the_scene_in_order() {
        let rows = rows(&scene(), &HashSet::new());
        assert_eq!(
            labels(&rows),
            [
                ("camera", 0),
                ("settings", 0),
                ("lights (2)", 0),
                ("point", 1),
                ("point", 1),
                ("shapes (3)", 0),
                ("sphere — ball", 1),
                ("group (2 children)", 1),
                ("cube", 2),
                ("plane", 2),
                ("csg difference", 1),
                ("sphere", 2),
                ("sphere", 2),
            ]
        );
        assert_eq!(rows[11].node, Node::Shape(vec![2, 0]));
        assert_eq!(rows[12].node, Node::Shape(vec![2, 1]));
        let foldable: Vec<&str> = rows
            .iter()
            .filter(|r| r.foldable)
            .map(|r| r.label.as_str())
            .collect();
        assert_eq!(
            foldable,
            [
                "lights (2)",
                "shapes (3)",
                "group (2 children)",
                "csg difference"
            ]
        );
    }

    #[test]
    fn folding_shapes_hides_every_shape() {
        let rows = rows(&scene(), &HashSet::from([Node::Shapes]));
        assert_eq!(rows.last().unwrap().label, "shapes (3)");
    }

    #[test]
    fn folding_a_group_hides_its_children() {
        let rows = rows(&scene(), &HashSet::from([Node::Shape(vec![1])]));
        let labels: Vec<&str> = rows.iter().map(|r| r.label.as_str()).collect();
        assert!(
            !labels.contains(&"cube") && !labels.contains(&"plane"),
            "{labels:?}"
        );
        assert!(labels.contains(&"csg difference"));
    }

    #[test]
    fn default_settings_have_no_row() {
        let mut desc = scene();
        desc.settings = Default::default();
        let rows = rows(&desc, &HashSet::new());
        assert!(rows.iter().all(|r| r.node != Node::Settings));
    }

    #[test]
    fn a_single_light_is_listed() {
        let mut desc = scene();
        desc.light = desc.lights.take().map(|mut l| l.remove(0));
        let rows = rows(&desc, &HashSet::new());
        assert_eq!(labels(&rows)[2..4], [("lights (1)", 0), ("point", 1)]);
    }

    #[test]
    fn json_leaves_out_a_groups_children() {
        let v: Value = serde_json::from_str(&json(&scene(), &Node::Shape(vec![1]))).unwrap();
        assert_eq!(v["type"], "group");
        assert!(v.get("children").is_none(), "{v}");
    }

    #[test]
    fn json_of_a_nested_shape() {
        let v: Value = serde_json::from_str(&json(&scene(), &Node::Shape(vec![2, 1]))).unwrap();
        assert_eq!(v["type"], "sphere");
    }
}
