# Milestone 30: a config file that remembers folders and the window

Branch `metal/30-config-file`, merged into develop with `--no-ff` when
done. Milestone 29 added the About panel; 31 is planned as an app bundle.

## Goal

Every panel opens in the folder rtc-edit or rtc-metal was started from,
and the editor always opens at the top-left at 900×640. Started from
Finder (milestone 31's bundle) that folder is `/`. The two programs should
remember where scenes and images were last saved and where the editor's
window was, in one file both of them read and write.

Success:
- On the first run, Open… and Save As… for an untitled scene start in
  `~/rtc/scenes`, and the viewer's Save Image… in `~/rtc/renders`; each
  folder is made when its panel is about to open.
- After a scene or an image is saved elsewhere, the next panel (this run
  or a later one) starts there.
- The editor reopens where it was and at its size, unless that no longer
  fits its screen.
- The viewer started on its own from the terminal behaves the same.

## Out of scope

- File ▸ Open Recent, the viewer's path, the viewer's default samples and
  traversal (planned after 31).
- The help window's frame and the viewer's frame.
- Saving the frame while the window moves; it is saved when the editor
  closes.
- Frames on a second screen (see Decisions: the frame).
- A settings window; the file is edited by hand if at all.
- `NSUserDefaults` (needs a bundle id; revisit in 31 if useful).

## Decisions

### The file

- `~/Library/Application Support/rtc-edit/config.json` (`$HOME`); no file
  until the first write.
- A JSON object; every key optional:

  ```json
  {
    "scene_folder": "/Users/me/scenes",
    "image_folder": "/Users/me/renders",
    "editor_frame": { "x": 0, "y": 25, "width": 900, "height": 640 }
  }
  ```

- A missing key means the default: `~/rtc/scenes`, `~/rtc/renders`, the
  top-left at 900×640.
- A missing, unreadable or non-object file reads as empty, never an error;
  the next write replaces it.
- Keys a program doesn't know are kept on every write.

### rtc-config (new crate, `RUST/rtc-config`)

One file, `src/lib.rs`, depending only on `serde_json`:

```rust
pub const SCENE_FOLDER: &str = "scene_folder";
pub const IMAGE_FOLDER: &str = "image_folder";
pub const EDITOR_FRAME: &str = "editor_frame";

/// $HOME/Library/Application Support/rtc-edit/config.json; None without $HOME.
pub fn path() -> Option<PathBuf>;
/// The file's object; empty when missing, unreadable or not an object.
pub fn load(path: &Path) -> Map<String, Value>;
/// Re-read, set `key`, write a temporary file beside it and rename it
/// into place; makes the folder.
pub fn set(path: &Path, key: &str, value: Value) -> Result<(), String>;
/// `map[key]` as a path when it is a non-empty string, else `$HOME/rtc/<default>`.
pub fn folder(map: &Map<String, Value>, key: &str, default: &str) -> PathBuf;
/// `create_dir_all`, errors ignored: a panel given a missing folder
/// opens somewhere else.
pub fn ensure_dir(dir: &Path);
```

- `default` is `scenes` or `renders`; without `$HOME`, `folder` gives the
  current folder.
- Two programs writing at once: each `set` re-reads first, so a clash
  needs two writes in the same instant, and then loses one key's new
  value, never the file (the rename is atomic).
- The temporary file is `config.json.<pid>.tmp`, so the two programs never
  share one.

### The editor (rtc-edit)

- The scene folder:
  - `panel_dir(scene: Option<&Path>, map) -> PathBuf`: the scene's folder
    when it has a path, else `folder(map, SCENE_FOLDER, "scenes")`.
  - Open… and Save As… start there; `ensure_dir` before the panel.
  - A successful Open… or Save As… stores the chosen file's folder as
    `scene_folder`. `rtc-edit scene.json` doesn't.
  - `dir()` (where an untitled scene's relative asset paths resolve, and
    the preview's base) is unchanged: still the current folder.
- The frame:
  - `frame_from(map) -> Option<Rectangle>`: `editor_frame` with four
    finite numbers, width and height within 200×150 to 8192 points;
    anything else is `None` (review fix: a hand-edited `1e39` or `1` was
    taken).
  - Boot opens the editor at it, or at the top-left at 900×640.
  - On main's first Opened, `window::monitor_size(main)`; when
    `fits(frame, screen)` is false, `window::move_to(main, ORIGIN)` and
    `window::resize(main, 900×640)`.
  - `fits(frame: Rectangle, screen: Option<Size>) -> bool`: `x ≥ 0`,
    `y ≥ 0`, `x + width ≤ screen.width`, `y + height ≤ screen.height`;
    `None` doesn't fit. `fits` and `frame_from` live in `windows.rs`.
  - Global coordinates against one screen's size: a frame on a second
    screen to the right doesn't fit and returns to the top-left, the same
    limit `help_position` has.
  - Saved when the editor closes: the close button and ⌘Q (which skips
    `Drop`) both reach `clean_up`, the `atexit` hook, which writes
    `editor_frame` from a static the frame tracking keeps up to date, as
    `CHILDREN` is.
- A failed `set` puts `couldn't save settings: <error>` in the status line;
  at quit it goes to stderr as `rtc-edit: couldn't save settings: <error>`.

### The viewer (rtc-metal)

- `save_path(mtm, name, dir: &Path)` takes the starting folder; both
  callers (File ▸ Save Image… and the Save button of the unsaved-image
  alert) pass `folder(load(path), IMAGE_FOLDER, "renders")` after
  `ensure_dir`.
- After `write_image` succeeds, the image's folder, made absolute against
  the current folder, is stored as `image_folder` (review fix: `-o out.ppm`
  stored `""`, and `-o renders/x.ppm` a relative folder other programs
  would read as somewhere else). A cancelled panel or failed write stores nothing.
- A failed `set` prints `rtc-metal: couldn't save settings: <error>` to
  stderr; no alert.
- No `path()` (no `$HOME`): panels start in the current folder, as today,
  and nothing is stored.

### Docs

- rtc-edit's README: the file's location and keys, the `~/rtc` defaults;
  milestone 30 in the list and 31 (the bundle) planned. rtc-metal's README:
  where Save Image… starts.

## Tests

- rtc-config, each on a temporary folder:
  - `load` of a missing file, bad JSON, `[]` and `{}` is empty.
  - `set` keeps other keys and unknown keys; makes the folder; replaces
    bad JSON; leaves no `.tmp` file.
  - `path()` ends in `Library/Application Support/rtc-edit/config.json`.
  - `folder`: a string value; missing, empty or non-string gives
    `$HOME/rtc/<default>`.
- rtc-edit:
  - `panel_dir`: the scene's folder with a path; else `scene_folder`;
    else `~/rtc/scenes`.
  - `frame_from`: a full frame; missing, partial, non-number, zero width
    give `None`.
  - `fits`: inside yes; off the right, off the bottom, negative x, no
    screen no.
- rtc-metal: `save_path` is a panel; covered by eye.
- By eye:
  - With no file, Open… opens `~/rtc/scenes` and Save Image… `~/rtc/renders`.
  - Save a scene and an image elsewhere, relaunch: the panels start there.
  - Move and resize the editor, ⌘Q, relaunch: it comes back. The same
    with the close button.
  - Set an off-screen `editor_frame` by hand: the editor comes back at the
    top-left at 900×640.
  - Save a scene, then an image: both keys are in the file.

## Errors

- No `$HOME`: no file; defaults are the current folder; nothing stored.
- Bad file: read as empty, replaced on the next write.
- A write fails: the status line (editor), or stderr (at quit, and the
  viewer).
- A default folder can't be made: the panel opens wherever macOS puts it.
