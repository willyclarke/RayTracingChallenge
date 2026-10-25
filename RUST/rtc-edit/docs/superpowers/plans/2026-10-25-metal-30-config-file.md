# Config File Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** One JSON file in Application Support, read and written by rtc-edit and rtc-metal, remembering the scene folder, the image folder and the editor's frame.

**Architecture:** A new one-file crate `rtc-config` owns the file (path, load, set with re-read + atomic rename, folder defaults). rtc-edit uses it for Open…/Save As… and the window frame (saved from the `atexit` hook); rtc-metal uses it for its save panel. Pure helpers (`panel_dir`, `frame_from`, `fits`) carry the logic that tests can reach.

**Tech Stack:** Rust 2024, serde_json, iced 0.14, rfd 0.15, objc2-app-kit 0.3.2.

**Spec:** `rtc-edit/docs/superpowers/specs/2026-10-25-metal-30-config-file-design.md`

## Global Constraints

- File: `$HOME/Library/Application Support/rtc-edit/config.json`; keys `scene_folder`, `image_folder`, `editor_frame` (`{x, y, width, height}`).
- Defaults: `$HOME/rtc/scenes`, `$HOME/rtc/renders`, the top-left at 900×640.
- A missing, unreadable or non-object file reads as empty; unknown keys survive every write.
- Temporary file `config.json.<pid>.tmp`, renamed into place.
- Error text: `couldn't save settings: <error>` (status line); `rtc-edit: couldn't save settings: <error>` / `rtc-metal: couldn't save settings: <error>` on stderr.
- `dir()` in rtc-edit (asset resolution, preview base) does not change.
- Tests never touch the real `$HOME` file; no test sets environment variables.
- Commit prefixes ENHANCEMENT/FIX/DOC; every message ends with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. Do not push.

## Review Focus

1. A scene with a path: Save As… must still rebase assets from `dir()` (the scene's folder), not from `panel_dir` — an untitled scene rebases from the current folder even though the panel starts in `~/rtc/scenes`. Task 2 pins `panel_dir`; the reviewer checks `save_as` keeps `dir` for `rebase`.
2. ⌘Q and the close button both write the frame: the frame must reach `clean_up` through a static, since `Editor` is gone at `atexit`. By eye in Task 3.
3. The viewer's unsaved-image alert on quit writes through `write_image` too, so `image_folder` is stored on that path as well (Task 4 stores inside `write_image`).
4. A config file with `editor_frame` of strings or negative size must not panic or open a zero-size window (Task 3's `frame_from` tests).
5. The viewer and the editor writing within moments of each other keep both keys (Task 1's `set` keeps other keys; by eye in Task 5).

---

### Task 1: rtc-config

**Files:**
- Create: `RUST/rtc-config/Cargo.toml`, `RUST/rtc-config/src/lib.rs`

**Interfaces:**
- Produces (spec's API, plus two testable cores — Ruling: `path`/`folder` read `$HOME`, and tests must not set env vars, so each wraps an `_in` version taking the home folder):
  - `pub const SCENE_FOLDER, IMAGE_FOLDER, EDITOR_FRAME: &str`
  - `pub fn path_in(home: &Path) -> PathBuf`; `pub fn path() -> Option<PathBuf>`
  - `pub fn load(path: &Path) -> Map<String, Value>`
  - `pub fn set(path: &Path, key: &str, value: Value) -> Result<(), String>`
  - `pub fn folder_in(map: &Map<String, Value>, key: &str, default: &str, home: Option<&Path>) -> PathBuf`; `pub fn folder(map, key, default) -> PathBuf`
  - `pub fn ensure_dir(dir: &Path)`
  - Re-exports `pub use serde_json::{Map, Value, json};`

- [ ] **Step 1: Write the crate skeleton and failing tests**

`Cargo.toml`: package `rtc-config`, version 0.1.0, edition 2024, dependency `serde_json = "1"`. Tests use a per-test folder `std::env::temp_dir().join(format!("rtc-config-test-{}-{name}", std::process::id()))`, removed at the start of the test.

```rust
#[test] fn load_reads_missing_bad_and_non_object_files_as_empty()
// missing file; "not json"; "[]"; "{}" → load(...).is_empty()
#[test] fn set_keeps_the_other_keys_and_unknown_keys()
// write {"mine": 1, "scene_folder": "/a"}; set(IMAGE_FOLDER, json!("/b"))
// → load == {"mine":1, "scene_folder":"/a", "image_folder":"/b"}
#[test] fn set_makes_the_folder_and_leaves_no_temporary_file()
// path = dir/"deep"/"config.json" (dir missing); set ok; file exists;
// no entry in dir/"deep" ends in ".tmp"
#[test] fn set_replaces_bad_json()
// write "oops"; set(SCENE_FOLDER, json!("/s")) → load == {"scene_folder":"/s"}
#[test] fn the_file_lives_in_application_support()
// path_in(Path::new("/Users/me")) == "/Users/me/Library/Application Support/rtc-edit/config.json"
#[test] fn folder_falls_back_to_rtc_under_home()
// {"scene_folder":"/s"} → "/s"; missing, "" and 3 → "/Users/me/rtc/scenes"
// home None → Path::new(".")
```

- [ ] **Step 2: Run them to verify they fail**

Run: `cd RUST/rtc-config && cargo test`
Expected: compile errors for the missing functions.

- [ ] **Step 3: Implement the functions in `src/lib.rs`**

`path()` = `std::env::var_os("HOME")` non-empty → `path_in`. `set`: `load`, insert, `serde_json::to_string_pretty`, `create_dir_all(parent)`, write `config.json.<pid>.tmp` beside it, `rename`; errors as `e.to_string()`. `ensure_dir` ignores errors. Crate doc comment states the file, its keys, and the two-writer rule from the spec.

- [ ] **Step 4: Run the tests**

Run: `cd RUST/rtc-config && cargo test && cargo clippy --all-targets -- -D warnings`
Expected: 6 passed; no warnings.

- [ ] **Step 5: Commit**

```bash
git add RUST/rtc-config/Cargo.toml RUST/rtc-config/src/lib.rs RUST/rtc-config/Cargo.lock
git commit -m "ENHANCEMENT: rtc-config, one settings file for rtc-edit and rtc-metal (milestone 30)"
```

### Task 2: the editor's scene folder

**Files:**
- Modify: `RUST/rtc-edit/Cargo.toml` (`rtc-config = { path = "../rtc-config" }`), `RUST/rtc-edit/build.rs` (rerun on `../rtc-config/src`), `RUST/rtc-edit/src/paths.rs`, `RUST/rtc-edit/src/main.rs` (`open`, `save_as`)

**Interfaces:**
- Consumes: Task 1's `load`, `set`, `folder_in`, `folder`, `ensure_dir`, `SCENE_FOLDER`, `path`.
- Produces: `paths::panel_dir_in(scene: Option<&Path>, map: &Map<String, Value>, home: Option<&Path>) -> PathBuf` and `paths::panel_dir(scene, map) -> PathBuf` (wraps with `$HOME`); `fn config() -> Map<String, Value>` and `Editor::remember(&mut self, key: &str, value: Value)` in main.rs (Task 3 uses `config()`; it writes from `clean_up`, not through `remember`, since the status line is gone by then).

- [ ] **Step 1: Write the failing test in `paths.rs`**

```rust
#[test]
fn panels_start_in_the_scenes_folder_then_the_remembered_one() {
    let home = Some(Path::new("/Users/me"));
    let mut map = Map::new();
    assert_eq!(panel_dir_in(None, &map, home), Path::new("/Users/me/rtc/scenes"));
    map.insert("scene_folder".into(), json!("/s"));
    assert_eq!(panel_dir_in(None, &map, home), Path::new("/s"));
    assert_eq!(panel_dir_in(Some(Path::new("/x/a.json")), &map, home), Path::new("/x"));
}
```

- [ ] **Step 2: Run it to verify it fails**

Run: `cd RUST/rtc-edit && cargo test --lib panels_start`
Expected: compile error, `panel_dir_in` not found.

- [ ] **Step 3: Implement**

- `panel_dir_in`: `Some(p)` → `dir_of(p)`, else `folder_in(map, SCENE_FOLDER, "scenes", home)`.
- `open`: start the panel at `panel_dir(self.path.as_deref(), &config())` after `ensure_dir`; on a successful load, `self.remember(SCENE_FOLDER, json!(dir_of(&path)))`.
- `save_as`: the panel starts at `panel_dir(...)`; `rebase` keeps `dir = self.dir()`; on `Ok`, `remember(SCENE_FOLDER, ...)` after the status is set, so a failure's text replaces "saved".
- `config()` = `rtc_config::path().map_or_else(Map::new, |p| load(&p))`; `remember` calls `set` when `path()` is `Some`, and on `Err(e)` sets `status = format!("couldn't save settings: {e}")`.

- [ ] **Step 4: Run the suite**

Run: `cd RUST/rtc-edit && cargo test --lib && cargo clippy --all-targets -- -D warnings`
Expected: all pass (156 + 1); no warnings.

- [ ] **Step 5: Commit**

```bash
git add RUST/rtc-edit/Cargo.toml RUST/rtc-edit/Cargo.lock RUST/rtc-edit/build.rs RUST/rtc-edit/src/paths.rs RUST/rtc-edit/src/main.rs
git commit -m "ENHANCEMENT: rtc-edit's Open… and Save As… start in the remembered scene folder (milestone 30)"
```

### Task 3: the editor's frame

**Files:**
- Modify: `RUST/rtc-edit/src/windows.rs`, `RUST/rtc-edit/src/main.rs` (`boot`, `clean_up`, the `Message::Window` arm, a new `FRAME` static)

**Interfaces:**
- Consumes: Task 1's `load`, `set`, `EDITOR_FRAME`, `path`, `json`.
- Produces: `windows::frame_from(map: &Map<String, Value>) -> Option<Rectangle>`, `windows::fits(frame: Rectangle, screen: Option<Size>) -> bool`, `windows::frame_json(frame: Rectangle) -> Value`; `Message::Fit(Option<Size>)`.

- [ ] **Step 1: Write the failing tests in `windows.rs`**

```rust
#[test]
fn a_saved_frame_needs_four_numbers_and_a_size() {
    let map = |v: Value| Map::from_iter([("editor_frame".to_string(), v)]);
    let full = json!({"x": 10, "y": 25, "width": 900, "height": 640});
    assert_eq!(frame_from(&map(full.clone())),
        Some(Rectangle::new(Point::new(10.0, 25.0), Size::new(900.0, 640.0))));
    assert_eq!(frame_from(&map(frame_json(frame_from(&map(full)).unwrap()))).map(|r| r.x), Some(10.0));
    assert_eq!(frame_from(&Map::new()), None);
    for bad in [json!({"x": 0, "y": 0, "width": 900}), json!({"x": "0", "y": 0, "width": 900, "height": 640}),
                json!({"x": 0, "y": 0, "width": 0, "height": 640}), json!({"x": 0, "y": 0, "width": -5, "height": 640}),
                json!([1, 2, 3, 4])] {
        assert_eq!(frame_from(&map(bad.clone())), None, "{bad}");
    }
}

#[test]
fn a_frame_fits_only_inside_its_screen() {
    let screen = Some(Size::new(1440.0, 900.0));
    let at = |x, y, w, h| Rectangle::new(Point::new(x, y), Size::new(w, h));
    assert!(fits(at(0.0, 25.0, 900.0, 640.0), screen));
    assert!(fits(at(540.0, 260.0, 900.0, 640.0), screen)); // touching both edges
    assert!(!fits(at(600.0, 25.0, 900.0, 640.0), screen));
    assert!(!fits(at(0.0, 300.0, 900.0, 640.0), screen));
    assert!(!fits(at(-1.0, 25.0, 900.0, 640.0), screen));
    assert!(!fits(at(0.0, 25.0, 900.0, 640.0), None));
}
```

- [ ] **Step 2: Run them to verify they fail**

Run: `cd RUST/rtc-edit && cargo test --lib frame`
Expected: compile errors for `frame_from`, `fits`, `frame_json`.

- [ ] **Step 3: Implement**

- `windows.rs`: the three functions; numbers via `Value::as_f64`, cast to `f32`.
- `main.rs`:
  - `static FRAME: Mutex<Option<Rectangle>> = Mutex::new(None);`, set whenever `e.frame` changes in the Opened / Moved / Resized arms (one helper `track(frame)`).
  - `boot`: `let saved = windows::frame_from(&config());` — open at its position and size, else as today; `e.frame` starts from it.
  - Opened, first time (the `menus.is_none()` branch): also return `window::monitor_size(e.main).map(Message::Fit)`.
  - `Message::Fit(screen)`: when `!fits(e.frame, screen)`, `Task::batch([window::move_to(e.main, Point::ORIGIN), window::resize(e.main, Size::new(WIDTH, HEIGHT))])`.
  - `clean_up`: if `FRAME` holds a frame and `path()` is `Some`, `set(EDITOR_FRAME, frame_json(..))`; `Err(e)` → `eprintln!("rtc-edit: couldn't save settings: {e}")`. Before the viewer is dropped, so a slow viewer can't delay it out.

- [ ] **Step 4: Run the suite**

Run: `cd RUST/rtc-edit && cargo test --lib && cargo clippy --all-targets -- -D warnings`
Expected: all pass; no warnings.

- [ ] **Step 5: Commit**

```bash
git add RUST/rtc-edit/src/windows.rs RUST/rtc-edit/src/main.rs
git commit -m "ENHANCEMENT: rtc-edit reopens where it was, back at the top-left when that no longer fits (milestone 30)"
```

### Task 4: the viewer's image folder

**Files:**
- Modify: `RUST/rtc-metal/Cargo.toml` (`rtc-config = { path = "../rtc-config" }`), `RUST/rtc-metal/src/window.rs` (`save_path`, its two callers, `write_image`)

**Interfaces:**
- Consumes: Task 1's `path`, `load`, `set`, `folder`, `ensure_dir`, `IMAGE_FOLDER`.
- Produces: `fn image_dir() -> PathBuf` (`path()` → `folder(load, IMAGE_FOLDER, "renders")` after `ensure_dir`; no `path()` → `current_dir`, else `.`); `save_path(mtm, name, dir: &Path)`.

- [ ] **Step 1: Change `save_path` to take `dir: &Path`**

Both callers pass `&image_dir()`. Update the `ask` doc comment ("in the current directory" → "in the image folder").

- [ ] **Step 2: Store the folder after a write**

In `write_image`, on `Ok`: when `path()` is `Some`, `set(IMAGE_FOLDER, json!(parent of path))`; `Err(e)` → `eprintln!("rtc-metal: couldn't save settings: {e}")`.

Ruling (planned): no new unit test in rtc-metal — `save_path` and `write_image` need AppKit panels and a GPU, and the logic they use (`folder`, `set`) is tested in Task 1. The by-eye checks in Task 5 cover it.

- [ ] **Step 3: Run the suite**

Run: `cd RUST/rtc-metal && cargo test --lib && cargo clippy --all-targets -- -D warnings`
Expected: all pass as before; no warnings.

- [ ] **Step 4: Commit**

```bash
git add RUST/rtc-metal/Cargo.toml RUST/rtc-metal/Cargo.lock RUST/rtc-metal/src/window.rs
git commit -m "ENHANCEMENT: rtc-metal's Save Image… starts in the remembered image folder (milestone 30)"
```

### Task 5: docs and the eye check

**Files:**
- Modify: `RUST/rtc-edit/README.md`, `RUST/rtc-metal/README.md`

- [ ] **Step 1: README text**

rtc-edit: a "Settings" section — the file's path, the three keys with the example from the spec, the `~/rtc/scenes` / `~/rtc/renders` defaults, "delete the file to forget"; milestone list: `**30:** a config file: the scene and image folders and the editor's window.` (not yet "done") and `**31 (planned):** an app bundle.`; the layout table gets `rtc-config`. rtc-metal: Save Image… starts in the remembered image folder, `~/rtc/renders` at first.

- [ ] **Step 2: Build release binaries**

Run: `cd RUST/rtc-metal && cargo build --release && cd ../rtc-edit && cargo build --release`
Expected: both build.

- [ ] **Step 3: Commit**

```bash
git add RUST/rtc-edit/README.md RUST/rtc-metal/README.md
git commit -m "DOC: the settings file in rtc-edit's and rtc-metal's READMEs (milestone 30)"
```

The by-eye list from the spec's Tests section goes to the user after the final review.
