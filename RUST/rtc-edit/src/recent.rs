//! The recent-scenes list behind File ▸ Open Recent, kept in the settings
//! file under `recent_scenes`, newest first.

use rtc_config::{Map, RECENT_SCENES, Value};
use std::path::{Path, PathBuf};

/// The most scenes the list keeps (the macOS default).
pub const MAX: usize = 10;

/// `map[recent_scenes]` as paths: items that aren't non-empty strings are
/// skipped, and at most [`MAX`] are kept; a missing or non-array key is an
/// empty list.
pub fn from_config(map: &Map<String, Value>) -> Vec<PathBuf> {
    map.get(RECENT_SCENES)
        .and_then(Value::as_array)
        .map(|items| {
            items
                .iter()
                .filter_map(Value::as_str)
                .filter(|s| !s.is_empty())
                .take(MAX)
                .map(PathBuf::from)
                .collect()
        })
        .unwrap_or_default()
}

/// The list as the settings file stores it.
pub fn to_value(list: &[PathBuf]) -> Value {
    list.iter()
        .map(|p| Value::from(p.to_string_lossy().into_owned()))
        .collect()
}

/// The list with `path` first, made absolute (the working directory in front
/// of a relative path; symlinks are left alone), without any other copy of
/// it, and at most [`MAX`] long.
pub fn add(list: &[PathBuf], path: &Path) -> Vec<PathBuf> {
    let path = std::path::absolute(path).unwrap_or_else(|_| path.to_path_buf());
    std::iter::once(path.clone())
        .chain(list.iter().filter(|&p| *p != path).cloned())
        .take(MAX)
        .collect()
}

/// The list without `path`.
pub fn remove(list: &[PathBuf], path: &Path) -> Vec<PathBuf> {
    list.iter().filter(|p| *p != path).cloned().collect()
}

/// One menu label per entry: its file name, with ` — <folder>` added when
/// another entry has the same file name. A path with no file name shows
/// whole.
pub fn labels(list: &[PathBuf]) -> Vec<String> {
    let name = |p: &PathBuf| {
        p.file_name().map_or_else(
            || p.display().to_string(),
            |n| n.to_string_lossy().into_owned(),
        )
    };
    let names: Vec<String> = list.iter().map(name).collect();
    list.iter()
        .zip(&names)
        .map(|(p, n)| {
            let shared = names.iter().filter(|m| *m == n).count() > 1;
            match p.parent().and_then(Path::file_name) {
                Some(folder) if shared => format!("{n} — {}", folder.to_string_lossy()),
                _ => n.clone(),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rtc_config::json;

    fn paths(list: &[&str]) -> Vec<PathBuf> {
        list.iter().map(PathBuf::from).collect()
    }

    fn object(v: Value) -> Map<String, Value> {
        v.as_object().unwrap().clone()
    }

    #[test]
    fn add_puts_a_new_path_first() {
        let list = add(&paths(&["/a/x.json"]), Path::new("/b/y.json"));
        assert_eq!(list, paths(&["/b/y.json", "/a/x.json"]));
    }

    #[test]
    fn add_moves_a_duplicate_to_the_front() {
        let list = add(&paths(&["/a", "/b", "/c"]), Path::new("/c"));
        assert_eq!(list, paths(&["/c", "/a", "/b"]));
    }

    #[test]
    fn add_keeps_at_most_max() {
        let old: Vec<PathBuf> = (0..10)
            .map(|i| PathBuf::from(format!("/s/{i}.json")))
            .collect();
        let list = add(&old, Path::new("/s/new.json"));
        assert_eq!(list.len(), MAX);
        assert_eq!(list[0], Path::new("/s/new.json"));
        assert!(!list.contains(&PathBuf::from("/s/9.json")));
    }

    #[test]
    fn add_makes_relative_paths_absolute_and_dedupes() {
        let a = add(&[], Path::new("examples/x.json"));
        assert!(a[0].is_absolute(), "{}", a[0].display());
        assert!(a[0].ends_with("examples/x.json"));
        assert_eq!(add(&a, Path::new("examples/x.json")).len(), 1);
        assert_eq!(add(&a, &a[0]).len(), 1);
    }

    #[test]
    fn remove_drops_the_entry_and_keeps_order() {
        let list = paths(&["/a", "/b", "/c"]);
        assert_eq!(remove(&list, Path::new("/b")), paths(&["/a", "/c"]));
        assert_eq!(remove(&paths(&["/a"]), Path::new("/z")), paths(&["/a"]));
    }

    #[test]
    fn from_config_reads_the_array() {
        let map = object(json!({"recent_scenes": ["/a.json", "/b.json"]}));
        assert_eq!(from_config(&map), paths(&["/a.json", "/b.json"]));
    }

    #[test]
    fn from_config_skips_bad_items_and_caps() {
        for v in [
            json!({}),
            json!({"recent_scenes": 3}),
            json!({"recent_scenes": {"a": 1}}),
        ] {
            assert!(from_config(&object(v.clone())).is_empty(), "{v}");
        }
        let map = object(json!({"recent_scenes": ["/a", 7, "", null, "/b"]}));
        assert_eq!(from_config(&map), paths(&["/a", "/b"]));
        let many: Vec<String> = (0..15).map(|i| format!("/s/{i}.json")).collect();
        let list = from_config(&object(json!({ "recent_scenes": many })));
        assert_eq!(list.len(), MAX);
        assert_eq!(list[0], Path::new("/s/0.json"));
        assert_eq!(list[9], Path::new("/s/9.json"));
    }

    #[test]
    fn to_value_round_trips() {
        let list = paths(&["/a/x.json", "/b/y.json"]);
        let map = object(json!({ "recent_scenes": to_value(&list) }));
        assert_eq!(from_config(&map), list);
    }

    #[test]
    fn labels_are_file_names_when_unique() {
        assert_eq!(
            labels(&paths(&["/a/x.json", "/b/y.json"])),
            ["x.json", "y.json"]
        );
    }

    #[test]
    fn labels_add_the_folder_to_shared_names() {
        let list = paths(&[
            "/p/examples/untitled.json",
            "/p/rtc-edit/untitled.json",
            "/p/pencil.json",
        ]);
        assert_eq!(
            labels(&list),
            [
                "untitled.json — examples",
                "untitled.json — rtc-edit",
                "pencil.json"
            ]
        );
    }

    #[test]
    fn labels_show_a_nameless_path_whole() {
        assert_eq!(labels(&paths(&["/"])), ["/"]);
    }
}
