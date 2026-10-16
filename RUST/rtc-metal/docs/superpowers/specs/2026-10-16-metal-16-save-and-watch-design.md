# rtc-metal milestone 16: saving scenes, and `--window --watch`

Branch `metal/16-save-and-watch`, merged into develop with `--no-ff` when
done. This is the first step towards a scene editor, which comes as two
programs: an iced editor (tree and inspector, milestones 17–19), and the
existing Metal window as its viewer. The editor writes the scene to a
file on every edit, and the viewer reloads that file. This milestone
builds the two foundations: rtc-rust can save scenes, and `--window` can
watch its scene file.

## Goal

- rtc-rust writes a `SceneDescription` back to JSON that loads into an
  equal description, and reads like a hand-written file.
- `rtc-metal scene.json --window --watch` reloads the scene whenever the
  file changes and starts refining again, keeping the window's size,
  traversal and sample target.

Success:
- Every scene in `rtc-rust/scenes/` and `rtc-metal/scenes/` round-trips:
  load → `to_json` → load gives an equal `SceneDescription`.
- Editing a watched scene in a text editor and saving updates the window
  within about half a second, also while it is paused at its target.
- All existing tests pass, and every existing render is unchanged.

## Out of scope

- The editor itself (milestones 17–19).
- Rewriting relative `.obj`/`.ppm` paths for another folder. Save As and
  the editor's preview file need it, so it comes with them in milestone 17.
- Keeping the exact text of a file: `1` becomes `1.0`, key order follows
  the structs, and unknown keys are already rejected by the loader.
- Blending an edit into the image instead of restarting. Every reload
  restarts from one sample, by decision.
- File-system event APIs (FSEvents, kqueue). Polling is enough at this
  rate.

## Decisions

- **`Serialize` and `PartialEq` on every description type** in
  `rtc-rust/src/scene.rs`, next to the existing `Deserialize`. Keys stay as
  they are, since the same serde attributes (`tag`, `rename_all`,
  `flatten`) drive both directions.
- **Absent and default values are left out**, with
  `skip_serializing_if`:
  - every `Option` field when `None`;
  - `transform` when empty, `closed` when `false`;
  - `settings` when it equals `Settings::default()`.
  `shapes` is always written, even when empty.
- **What the file used is kept**: `comment` (any JSON value), and whether
  it had `"light"` or `"lights"`.
- **API:**
  - `SceneDescription::to_json(&self) -> String`, pretty-printed by serde_json;
  - `SceneDescription::save(&self, path) -> Result<(), SceneError>`.
- **`save` replaces the file in one step**: it writes `<file>.tmp` in the
  same folder, then renames it over the target. A watcher never reads a
  half-written scene from us. Text editors that write in place can still
  produce one; that case is in Error handling.
- **`--watch` is a flag**, valid only with `--window`, and it watches the
  scene file given on the command line.
- **Polling:** an `NSTimer` on the main thread checks the file's
  modification time every 0.25 s. It is a timer, not part of `draw`,
  because `draw` stops running once the image reaches its target.
- **What a reload keeps and replaces:**
  - kept: the window's current size, traversal and sample target;
  - replaced: the scene (shapes, lights, settings).
- **Camera on reload:**
  - the description's camera unchanged → the current orbit and zoom are
    kept;
  - the camera changed → the window takes it as the new home view (as
    space would) at the current size.
  This compares `CameraDescription`s, so moving a sphere never throws
  away your orbit.
- **Accumulation restarts** on every successful reload, and the
  stopwatch starts from zero.

## rtc-rust

`src/scene.rs`:

```rust
#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct SceneDescription { ... }   // likewise every description type

impl SceneDescription {
    /// The scene as pretty-printed JSON; `from_json` reads it back into an
    /// equal description.
    pub fn to_json(&self) -> String;
    /// Write `to_json` to `path`, replacing the file in one step (a
    /// temporary file in the same folder, then a rename).
    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<(), SceneError>;
}
```

The module doc, `README.md` (two places) and CLAUDE.md's `scene`
paragraph stop calling the format load-only, and mention `to_json` and
`save`.

## rtc-metal

- **`src/watch.rs` (new):** `Watch { path, modified }`.
  - `Watch::new(path)` records the file's current modification time.
  - `changed(&mut self) -> bool` is true once per new modification time.
    A missing or unreadable file counts as unchanged, so a file being
    replaced is picked up on the next tick.
  - Testable without a window.
- **`src/window.rs`:**
  - `Renderer::reload(&mut self, desc: &SceneDescription, base_dir: &Path) -> Result<(), String>`
    builds the world, makes the `GpuScene` at the current size and
    traversal, applies the camera rule above, and rebuilds the
    accumulator. On error, nothing in the renderer changes.
  - `Renderer` keeps the loaded `CameraDescription` to compare against.
  - `run` takes an `Option<Watch>`. When given, it schedules the timer.
    The timer calls `checkWatchedFile:` on the view: if the file changed,
    load it, `reload`, then resume drawing.
- **`src/main.rs`:** parses `--watch`, and errors with the usage text
  without `--window`. `--help`, the module doc and README document it.

## Error handling

The window never exits on a bad file.

- **The file doesn't parse, or doesn't build** (half-written, a typo, a
  missing `.obj`): the old scene and image stay. The overlay and title
  show `reload failed: <message>` until the next successful reload. The
  message also goes to stderr.
- **The scene isn't supported on the GPU** (a spot light, more than 8
  lights): the same path, with `GpuScene::from_world`'s message.
- **The file is deleted:** nothing happens until it reappears.

## Testing

rtc-rust:
- `test_scene_10`: every `scenes/*.json` round-trips into an equal
  description.
- `test_scene_11`: a minimal scene (camera, one light, a sphere with no
  transform or material) serializes to exactly the expected JSON, with no
  `transform`, `settings`, `closed` or `null` in it.
- `test_scene_12`: a scene using `"lights"` and a `comment` keeps both
  through a round trip.
- `test_scene_13`: `save` writes a file that `load` builds into the same
  world (shape count, light count, camera size), and leaves no `.tmp`
  file behind.

rtc-metal:
- Every `rtc-metal/scenes/*.json` round-trips.
- `Watch::changed`: false at first, true once after the file is
  rewritten, false again until the next change, false while the file is
  missing.
- `Renderer::reload`:
  - with the same camera, the orbit and size are kept, and the frames
    restart at 0;
  - with a changed camera, it becomes the home view at the current size;
  - with a broken scene, it returns `Err` and the old scene and frames
    are untouched;
  - the traversal is kept.
- By eye: run `--window --watch` on a copy of a scene. Move a sphere in a
  text editor and save, while the window is refining and while it is
  paused; orbit first and check the view is kept; break the JSON and see
  the overlay message.

## Risks

- **`#[serde(flatten)]` with internally tagged enums**
  (`ShapeDescription.shape`, `PatternDescription`) is the part of serde
  most likely to serialize differently from what it reads. The round
  trip over every scene file catches that.
- **A text editor saving in two steps** (truncate, then write) can be
  read in between. That gives a reload error that heals on the next tick.
  It is acceptable, and the editor of milestone 17 uses `save`, which
  never does this.
- **A reload is synchronous on the main thread.** Rebuilding a large
  scene (the dice, the teapots) takes tens to hundreds of milliseconds
  and freezes the window that long. That is fine for editing; a
  background build is a later concern if ever needed.
