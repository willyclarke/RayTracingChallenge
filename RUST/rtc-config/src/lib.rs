//! The settings file rtc-edit and rtc-metal share:
//! `~/Library/Application Support/rtc-edit/config.json`, one JSON object
//! whose keys are all optional. A missing, unreadable or non-object file
//! reads as empty, and the next write replaces it. Keys a program doesn't
//! know are kept on every write.
//!
//! Both programs write it. Each [`set`] re-reads the file, changes one key
//! and renames a temporary file into place, so two writes at the same
//! instant lose at most one key's new value, never the file.

use std::fs;
use std::path::{Path, PathBuf};

pub use serde_json::{Map, Value, json};

/// The folder the last scene was opened or saved in.
pub const SCENE_FOLDER: &str = "scene_folder";
/// The folder the last image was saved in.
pub const IMAGE_FOLDER: &str = "image_folder";
/// The editor window's `{x, y, width, height}` when it last closed.
pub const EDITOR_FRAME: &str = "editor_frame";

/// The file under the home folder `home`.
pub fn path_in(home: &Path) -> PathBuf {
    home.join("Library/Application Support/rtc-edit/config.json")
}

/// The file under `$HOME`; `None` without one.
pub fn path() -> Option<PathBuf> {
    home().map(|h| path_in(&h))
}

fn home() -> Option<PathBuf> {
    std::env::var_os("HOME")
        .filter(|h| !h.is_empty())
        .map(PathBuf::from)
}

/// The file's object; empty when it is missing, unreadable or not an object.
pub fn load(path: &Path) -> Map<String, Value> {
    fs::read_to_string(path)
        .ok()
        .and_then(|text| serde_json::from_str::<Value>(&text).ok())
        .and_then(|v| match v {
            Value::Object(map) => Some(map),
            _ => None,
        })
        .unwrap_or_default()
}

/// Set `key` to `value` in the file as it is now, making its folder; the
/// new file is written beside it and renamed into place.
pub fn set(path: &Path, key: &str, value: Value) -> Result<(), String> {
    let mut map = load(path);
    map.insert(key.to_string(), value);
    let text = serde_json::to_string_pretty(&map).map_err(|e| e.to_string())?;
    if let Some(dir) = path.parent() {
        fs::create_dir_all(dir).map_err(|e| e.to_string())?;
    }
    let mut tmp = path.as_os_str().to_owned();
    tmp.push(format!(".{}.tmp", std::process::id()));
    let tmp = PathBuf::from(tmp);
    fs::write(&tmp, text + "\n").map_err(|e| e.to_string())?;
    fs::rename(&tmp, path).map_err(|e| {
        let _ = fs::remove_file(&tmp);
        e.to_string()
    })
}

/// `map[key]` as a folder when it is a non-empty string, else
/// `<home>/rtc/<default>`, else the current folder (`.`).
pub fn folder_in(
    map: &Map<String, Value>,
    key: &str,
    default: &str,
    home: Option<&Path>,
) -> PathBuf {
    match map.get(key).and_then(Value::as_str) {
        Some(dir) if !dir.is_empty() => PathBuf::from(dir),
        _ => home.map_or_else(|| PathBuf::from("."), |h| h.join("rtc").join(default)),
    }
}

/// [`folder_in`] under `$HOME`.
pub fn folder(map: &Map<String, Value>, key: &str, default: &str) -> PathBuf {
    folder_in(map, key, default, home().as_deref())
}

/// Make `dir`, ignoring errors: a panel given a missing folder opens
/// somewhere else.
pub fn ensure_dir(dir: &Path) {
    let _ = fs::create_dir_all(dir);
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    /// A fresh, empty folder for one test.
    fn scratch(name: &str) -> PathBuf {
        let dir =
            std::env::temp_dir().join(format!("rtc-config-test-{}-{name}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn object(v: Value) -> Map<String, Value> {
        v.as_object().unwrap().clone()
    }

    #[test]
    fn load_reads_missing_bad_and_non_object_files_as_empty() {
        let dir = scratch("load");
        let file = dir.join("config.json");
        assert!(load(&file).is_empty());
        for text in ["not json", "[]", "{}"] {
            fs::write(&file, text).unwrap();
            assert!(load(&file).is_empty(), "{text}");
        }
    }

    #[test]
    fn set_keeps_the_other_keys_and_unknown_keys() {
        let dir = scratch("keep");
        let file = dir.join("config.json");
        fs::write(&file, r#"{"mine": 1, "scene_folder": "/a"}"#).unwrap();
        set(&file, IMAGE_FOLDER, json!("/b")).unwrap();
        assert_eq!(
            load(&file),
            object(json!({"mine": 1, "scene_folder": "/a", "image_folder": "/b"}))
        );
    }

    #[test]
    fn set_makes_the_folder_and_leaves_no_temporary_file() {
        let dir = scratch("make").join("deep");
        let file = dir.join("config.json");
        set(&file, SCENE_FOLDER, json!("/s")).unwrap();
        assert!(file.is_file());
        let names: Vec<_> = fs::read_dir(&dir)
            .unwrap()
            .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
            .collect();
        assert_eq!(names, ["config.json"]);
    }

    #[test]
    fn set_replaces_bad_json() {
        let dir = scratch("bad");
        let file = dir.join("config.json");
        fs::write(&file, "oops").unwrap();
        set(&file, SCENE_FOLDER, json!("/s")).unwrap();
        assert_eq!(load(&file), object(json!({"scene_folder": "/s"})));
    }

    #[test]
    fn the_file_lives_in_application_support() {
        assert_eq!(
            path_in(Path::new("/Users/me")),
            Path::new("/Users/me/Library/Application Support/rtc-edit/config.json")
        );
    }

    #[test]
    fn folder_falls_back_to_rtc_under_home() {
        let home = Some(Path::new("/Users/me"));
        let map = object(json!({"scene_folder": "/s"}));
        assert_eq!(
            folder_in(&map, SCENE_FOLDER, "scenes", home),
            Path::new("/s")
        );
        for v in [
            json!({}),
            json!({"scene_folder": ""}),
            json!({"scene_folder": 3}),
        ] {
            assert_eq!(
                folder_in(&object(v.clone()), SCENE_FOLDER, "scenes", home),
                Path::new("/Users/me/rtc/scenes"),
                "{v}"
            );
        }
        assert_eq!(
            folder_in(&map, IMAGE_FOLDER, "renders", None),
            Path::new(".")
        );
    }
}
