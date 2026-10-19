//! Asset paths (`.obj` meshes, `.ppm` textures) in a scene description:
//! rewritten for another folder on Save As, or made absolute for the
//! preview file.

use rtc_rust::scene::{
    PatternDescription, PatternKind, SceneDescription, ShapeDescription, ShapeKind, UvDescription,
};
use std::path::{Component, Path, PathBuf};

/// Call `f` on the `file` of every `obj` shape and every `uv_image`, however
/// deeply nested in groups, CSG operands and patterns.
pub fn for_each_asset(desc: &mut SceneDescription, mut f: impl FnMut(&mut String)) {
    for shape in &mut desc.shapes {
        shape_assets(shape, &mut f);
    }
}

fn shape_assets(shape: &mut ShapeDescription, f: &mut impl FnMut(&mut String)) {
    match &mut shape.shape {
        ShapeKind::Obj { file } => f(file),
        ShapeKind::Group { children } => {
            for child in children {
                shape_assets(child, f);
            }
        }
        ShapeKind::Csg { left, right, .. } => {
            shape_assets(left, f);
            shape_assets(right, f);
        }
        _ => {}
    }
    if let Some(pattern) = shape.material.as_mut().and_then(|m| m.pattern.as_mut()) {
        pattern_assets(pattern, f);
    }
}

fn pattern_assets(pattern: &mut PatternDescription, f: &mut impl FnMut(&mut String)) {
    match &mut pattern.pattern {
        PatternKind::Nested { a, b } | PatternKind::Blended { a, b } => {
            pattern_assets(a, f);
            pattern_assets(b, f);
        }
        PatternKind::Perturbed { pattern, .. } => pattern_assets(pattern, f),
        PatternKind::TextureMap { uv, .. } => uv_assets(uv, f),
        PatternKind::CubeMap {
            left,
            front,
            right,
            back,
            up,
            down,
        } => {
            for uv in [left, front, right, back, up, down] {
                uv_assets(uv, f);
            }
        }
        _ => {}
    }
}

fn uv_assets(uv: &mut UvDescription, f: &mut impl FnMut(&mut String)) {
    if let UvDescription::UvImage { file } = uv {
        f(file);
    }
}

/// Rewrite every relative asset path with `new`, on a copy; `desc` changes
/// only if every path succeeds.
fn rewrite(
    desc: &mut SceneDescription,
    mut new: impl FnMut(&Path) -> Result<PathBuf, String>,
) -> Result<(), String> {
    let mut copy = desc.clone();
    let mut error = None;
    for_each_asset(&mut copy, |file| {
        if error.is_some() || Path::new(file.as_str()).is_absolute() {
            return;
        }
        match new(Path::new(file.as_str())) {
            Ok(path) => *file = path.to_string_lossy().into_owned(),
            Err(e) => error = Some(e),
        }
    });
    match error {
        Some(e) => Err(e),
        None => {
            *desc = copy;
            Ok(())
        }
    }
}

/// `from_dir.join(file)`, canonicalized; the error names `file` as written.
fn resolve(from_dir: &Path, file: &Path) -> Result<PathBuf, String> {
    from_dir
        .join(file)
        .canonicalize()
        .map_err(|e| format!("{}: {e}", file.display()))
}

/// Rewrite relative asset paths so they reach the same files from `to_dir`
/// instead of `from_dir`. Absolute paths are left alone. On error nothing
/// changes, and the message names the path as written.
pub fn rebase(desc: &mut SceneDescription, from_dir: &Path, to_dir: &Path) -> Result<(), String> {
    let to_dir = to_dir
        .canonicalize()
        .map_err(|e| format!("{}: {e}", to_dir.display()))?;
    rewrite(desc, |file| {
        Ok(relative(&resolve(from_dir, file)?, &to_dir))
    })
}

/// The path from folder `dir` to `target`, both canonical: one `..` per
/// component of `dir` past their common start, then the rest of `target`.
fn relative(target: &Path, dir: &Path) -> PathBuf {
    let target: Vec<Component> = target.components().collect();
    let dir: Vec<Component> = dir.components().collect();
    let common = target.iter().zip(&dir).take_while(|(a, b)| a == b).count();
    let mut out = PathBuf::new();
    for _ in common..dir.len() {
        out.push("..");
    }
    for c in &target[common..] {
        out.push(c);
    }
    out
}

/// Make every relative asset path absolute, resolved against `from_dir`.
pub fn absolutize(desc: &mut SceneDescription, from_dir: &Path) -> Result<(), String> {
    let from_dir = if from_dir.as_os_str().is_empty() {
        Path::new(".")
    } else {
        from_dir
    };
    rewrite(desc, |file| resolve(from_dir, file))
}

/// The folder a scene file is in; a bare file name is in `.`.
pub fn dir_of(path: &Path) -> &Path {
    match path.parent() {
        Some(p) if !p.as_os_str().is_empty() => p,
        _ => Path::new("."),
    }
}

/// A scene's file name, or `untitled.json` for a scene with no file yet.
pub fn scene_name(path: Option<&Path>) -> String {
    path.and_then(Path::file_name)
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| "untitled.json".into())
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::{Value, json};
    use std::fs;
    use std::path::PathBuf;

    /// An empty folder of the test's own under the temp dir.
    fn fresh(name: &str) -> PathBuf {
        let root =
            std::env::temp_dir().join(format!("rtc_edit_paths_{name}_{}", std::process::id()));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(&root).unwrap();
        root.canonicalize().unwrap()
    }

    /// Create an empty file (and its folders) at `root/rel`.
    fn touch(root: &Path, rel: &str) {
        let path = root.join(rel);
        fs::create_dir_all(path.parent().unwrap()).unwrap();
        fs::write(path, "").unwrap();
    }

    fn scene(shapes: Value) -> SceneDescription {
        let v = json!({
            "camera": {"width": 10, "height": 5, "field_of_view": 1.0,
                       "from": [0, 0, -5], "to": [0, 0, 0], "up": [0, 1, 0]},
            "light": {"type": "point", "position": [0, 0, -5], "intensity": [1, 1, 1]},
            "shapes": shapes
        });
        SceneDescription::from_json(&v.to_string()).unwrap()
    }

    fn obj(file: &str) -> Value {
        json!({"type": "obj", "file": file})
    }

    fn assets(desc: &mut SceneDescription) -> Vec<String> {
        let mut out = vec![];
        for_each_asset(desc, |f| out.push(f.clone()));
        out
    }

    /// The number of keys named "comment" anywhere in `v`.
    fn comment_keys(v: &Value) -> usize {
        match v {
            Value::Object(m) => m
                .iter()
                .map(|(k, x)| usize::from(k == "comment") + comment_keys(x))
                .sum(),
            Value::Array(a) => a.iter().map(comment_keys).sum(),
            _ => 0,
        }
    }

    #[test]
    fn rebase_into_a_child_folder() {
        let root = fresh("child");
        touch(&root, "models/m.obj");
        fs::create_dir_all(root.join("out")).unwrap();
        let mut desc = scene(json!([obj("models/m.obj")]));
        rebase(&mut desc, &root, &root.join("out")).unwrap();
        assert_eq!(assets(&mut desc), ["../models/m.obj"]);
    }

    #[test]
    fn rebase_into_the_parent_folder() {
        let root = fresh("parent");
        touch(&root, "scenes/models/m.obj");
        let mut desc = scene(json!([obj("models/m.obj")]));
        rebase(&mut desc, &root.join("scenes"), &root).unwrap();
        assert_eq!(assets(&mut desc), ["scenes/models/m.obj"]);
    }

    #[test]
    fn rebase_into_a_sibling_folder() {
        let root = fresh("sibling");
        touch(&root, "models/m.obj");
        fs::create_dir_all(root.join("a")).unwrap();
        fs::create_dir_all(root.join("b")).unwrap();
        let mut desc = scene(json!([obj("../models/m.obj")]));
        rebase(&mut desc, &root.join("a"), &root.join("b")).unwrap();
        assert_eq!(assets(&mut desc), ["../models/m.obj"]);
    }

    #[test]
    fn rebase_leaves_absolute_paths_alone() {
        let root = fresh("absolute");
        touch(&root, "m.obj");
        fs::create_dir_all(root.join("out")).unwrap();
        let abs = root.join("m.obj").to_string_lossy().into_owned();
        let mut desc = scene(json!([obj(&abs)]));
        rebase(&mut desc, &root, &root.join("out")).unwrap();
        assert_eq!(assets(&mut desc), [abs]);
    }

    #[test]
    fn rebase_reaches_nested_assets() {
        let root = fresh("nested");
        for f in ["m.obj", "t1.ppm", "t2.ppm", "t3.ppm"] {
            touch(&root, f);
        }
        fs::create_dir_all(root.join("out")).unwrap();
        let image = |f: &str| json!({"type": "uv_image", "file": f});
        let checkers = json!({"type": "uv_checkers", "width": 2, "height": 2,
                              "a": [0, 0, 0], "b": [1, 1, 1]});
        let texture = json!({"type": "texture_map", "mapping": "spherical", "uv": image("t1.ppm")});
        let cube_map = json!({"type": "cube_map", "left": checkers, "front": checkers,
                              "right": checkers, "back": checkers, "up": checkers,
                              "down": image("t2.ppm")});
        let perturbed = json!({"type": "perturbed", "amplitude": 0.1, "frequency": 1.0,
            "pattern": {"type": "texture_map", "mapping": "planar", "uv": image("t3.ppm")}});
        let with = |p: Value| json!({"type": "sphere", "material": {"pattern": p}});
        let mut desc = scene(json!([
            {"type": "csg", "operation": "union",
             "left": {"type": "sphere"},
             "right": {"type": "group", "children": [obj("m.obj")]}},
            with(texture), with(cube_map), with(perturbed)
        ]));
        rebase(&mut desc, &root, &root.join("out")).unwrap();
        let mut got = assets(&mut desc);
        got.sort();
        assert_eq!(got, ["../m.obj", "../t1.ppm", "../t2.ppm", "../t3.ppm"]);
    }

    #[test]
    fn rebase_names_a_missing_asset_and_changes_nothing() {
        let root = fresh("missing");
        touch(&root, "here.obj");
        fs::create_dir_all(root.join("out")).unwrap();
        let mut desc = scene(json!([obj("here.obj"), obj("missing.obj")]));
        let before = desc.clone();
        let err = rebase(&mut desc, &root, &root.join("out")).unwrap_err();
        assert!(err.contains("missing.obj"), "{err}");
        assert_eq!(desc, before);
    }

    #[test]
    fn rebasing_every_scene_into_its_own_folder_changes_nothing() {
        let mut checked = 0;
        for dir in ["../rtc-metal/scenes", "../rtc-rust/scenes"] {
            for entry in fs::read_dir(dir).unwrap() {
                let path = entry.unwrap().path();
                if path.extension().is_none_or(|x| x != "json") {
                    continue;
                }
                let text = fs::read_to_string(&path).unwrap();
                let mut desc = SceneDescription::from_json(&text).unwrap();
                let before = desc.clone();
                rebase(&mut desc, Path::new(dir), Path::new(dir)).unwrap();
                assert_eq!(desc, before, "{}", path.display());
                let saved: Value = serde_json::from_str(&desc.to_json()).unwrap();
                let original: Value = serde_json::from_str(&text).unwrap();
                assert_eq!(
                    comment_keys(&saved),
                    comment_keys(&original),
                    "comments lost in {}",
                    path.display()
                );
                checked += 1;
            }
        }
        assert!(checked >= 20, "only {checked} scene files checked");
    }

    #[test]
    fn absolutize_gives_a_scene_that_loads_from_anywhere() {
        let root = fresh("absolutize");
        let dir = Path::new("../rtc-metal/scenes");
        let mut desc = SceneDescription::from_path(dir.join("chapter15.json")).unwrap();
        absolutize(&mut desc, dir).unwrap();
        assert!(assets(&mut desc).iter().all(|f| Path::new(f).is_absolute()));
        let saved = root.join("chapter15.json");
        desc.save(&saved).unwrap();
        let (world, _) = rtc_rust::scene::load(&saved).unwrap();
        assert!(world.shapes.len() > 2, "teapots missing");
    }

    #[test]
    fn an_empty_folder_means_the_current_one() {
        assert_eq!(dir_of(Path::new("x.json")), Path::new("."));
        assert_eq!(dir_of(Path::new("a/x.json")), Path::new("a"));
    }

    #[test]
    fn scene_name_is_the_file_name_or_untitled() {
        assert_eq!(scene_name(Some(Path::new("a/b/c.json"))), "c.json");
        assert_eq!(scene_name(None), "untitled.json");
    }
}
