# rtc-edit, the read-only scene editor — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `rtc-edit scene.json` shows the scene as a foldable tree with a
read-only JSON inspector, starts `rtc-metal --window --watch` on a
preview copy, and saves the scene (Save, and Save As with asset paths
rebased).

**Architecture:**
- **New crate `RUST/rtc-edit`:** a library holding everything that can be
  tested without a window (`paths`, `tree`, `preview`, `viewer`), plus a
  binary `src/main.rs` with the iced application.
- **The viewer is a child process.** The editor talks to it only through
  the preview file.
- **rtc-rust** gains one field: `ShapeDescription.comment`.

**Tech Stack:** Rust 2024; serde_json; iced 0.14 (`iced::application`,
`keyboard::listen`, `window::close_requests`); rfd 0.15 (`FileDialog`);
std::process.

**Spec:** `RUST/rtc-edit/docs/superpowers/specs/2026-10-16-metal-17-edit-tree-design.md`

## Global Constraints

- **Branch** `metal/17-edit-tree`. Never push. Commits end with
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Commit subjects:**
  - rtc-rust: `ENHANCEMENT: rtc-rust …`;
  - rtc-edit: end with `(milestone 17)`;
  - prefixes `ENHANCEMENT:` / `TEST:` / `DOCS:` / `FIX:`.
- **Test names:** rtc-rust continues `test_scene_N` from 14. rtc-edit
  uses descriptive snake_case names.
- **No Cargo workspace.** rtc-edit has its own `target/`, ignored in
  the top-level `.gitignore` like the other crates'.
- **Local assets:** `rtc-rust/models/teapot*.obj` are git-ignored
  downloads. Tests that load chapter15 rely on them, as rtc-metal's
  already do.
- **Dependencies of rtc-edit:**
  - `rtc-rust = { path = "../rtc-rust" }`;
  - `serde_json = "1"`;
  - `iced = "0.14"`;
  - `rfd = "0.15"`;
  - nothing else; the relative-path computation is written by hand.
- **Preview file:** `<temp dir>/rtc-edit-<pid>/<scene file name>`, with
  every asset path absolute.
- **Viewer:** `$RTC_METAL` if set, otherwise
  `<CARGO_MANIFEST_DIR>/../rtc-metal/target/release/rtc-metal`, run as
  `rtc-metal <preview> --window --watch`.
- **Keys:** Cmd-S Save, Cmd-Shift-S Save As, Cmd-R show viewer.
- **Title:** `rtc-edit — <file name>`.
- **Commands:**
  - rtc-rust (from `RUST/rtc-rust`): `cargo test --lib`,
    `cargo clippy --lib`;
  - rtc-edit (from `RUST/rtc-edit`): `cargo test`, `cargo clippy --all-targets`.
- **Formatting:** run `cargo fmt` in each crate touched before committing.
- **YAGNI:** no editing, no Open/New, no unsaved-changes prompt, no menu
  bar.

## Review Focus

1. **A scene given as a bare file name** (`rtc-edit chapter15.json`):
   `Path::parent()` is `""`. Expected: it is treated as `.`. Task 2's
   `an_empty_folder_means_the_current_one` pins it for `absolutize`;
   `main.rs` must use the same rule.
2. **Save As cancelled.** Expected: nothing changes, no error in the
   status line. Checked by eye in Task 5.
3. **The viewer exits at once**, for example because rtc-metal rejects a
   scene with a spot light. Expected: Cmd-R starts it again rather than
   believing it is still running. Task 4's
   `running_is_false_after_the_child_exits` pins it.
4. **Save As onto the scene's own folder or its own file.** Expected: the
   asset paths are unchanged, and the file is rewritten. Task 2's
   `rebasing_every_scene_into_its_own_folder_changes_nothing` pins the
   paths.
5. **Two editors open at once.** Expected: each has its own preview
   folder (the pid is in the name), and quitting one doesn't remove the
   other's. This follows from the folder name; no test.

---

### Task 1: rtc-rust keeps shape comments

**Files:**
- Modify: `RUST/rtc-rust/src/scene.rs` (`ShapeDescription`, tests
  `test_scene_10`, new `test_scene_14`)
- Modify: `RUST/rtc-rust/README.md` (the shape keys list)

**Interfaces:**
- Produces: `ShapeDescription { …, pub comment: Option<serde_json::Value> }`.
  It is the last field, with `#[serde(default)]` and
  `#[serde(skip_serializing_if = "Option::is_none")]`. Every
  `ShapeDescription { … }` literal in the crate gets `comment: None`.

- [ ] **Step 1: Write the failing test**

```rust
/// Scene saving - a shape's comment survives a round trip and is ignored by build
#[test]
fn test_scene_14() -> Result<(), String> {
    let json = r#"{"camera":{"width":10,"height":5,"field_of_view":1.0,"from":[0,0,-5],"to":[0,0,0],"up":[0,1,0]},
        "light":{"type":"point","position":[0,0,-5],"intensity":[1,1,1]},
        "shapes":[{"type":"sphere","comment":"the ball"}]}"#;
    let desc = SceneDescription::from_json(json).map_err(|e| e.to_string())?;
    let again = SceneDescription::from_json(&desc.to_json()).map_err(|e| e.to_string())?;
    // assert: again.shapes[0].comment == Some(json!("the ball"))
    // assert: again.build(Path::new(".")) succeeds with 1 shape
}
```

Also extend `test_scene_10`: for each file, the number of keys named
`"comment"` anywhere in the original JSON equals the number in
`to_json()`'s output. Count them with a small recursive helper over
`serde_json::Value`, `fn comment_keys(v: &serde_json::Value) -> usize`.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cargo test --lib test_scene_1` (from `RUST/rtc-rust`)
Expected: `test_scene_14` fails to compile, because there is no
`comment` field on `ShapeDescription`. The extended `test_scene_10`
passes even before the fix, because rtc-rust's own scenes have only
top-level comments. It guards that case; the shape comments live in
rtc-metal's scenes, which Task 2's all-scenes test covers.

- [ ] **Step 3: Add the field**, as in Interfaces; `build` ignores it.
  Add a `comment` line to the README's shape keys ("free-form notes,
  kept when saving").

- [ ] **Step 4: Run the tests**

Run: `cargo test --lib && cargo clippy --lib`
Expected: 353 passed, 7 ignored; clippy clean.

- [ ] **Step 5: Commit**
`ENHANCEMENT: rtc-rust keeps shape comments when saving`

---

### Task 2: the rtc-edit crate and asset paths

**Files:**
- Create: `RUST/rtc-edit/Cargo.toml`. In this task the only dependencies
  are `rtc-rust` and `serde_json`; `iced` and `rfd` arrive in Task 5.
- Create: `RUST/rtc-edit/src/lib.rs` (`pub mod paths;` and a module doc)
- Create: `RUST/rtc-edit/src/paths.rs`
- Modify: `.gitignore` (add `/RUST/rtc-edit/target` after
  `/RUST/rtc-metal/*.ppm`)

**Interfaces:**
- Consumes: Task 1's `ShapeDescription.comment`.
- Produces, in `rtc_edit::paths`:
  - `pub fn for_each_asset(desc: &mut SceneDescription, f: impl FnMut(&mut String))`;
  - `pub fn rebase(desc: &mut SceneDescription, from_dir: &Path, to_dir: &Path) -> Result<(), String>`;
  - `pub fn absolutize(desc: &mut SceneDescription, from_dir: &Path) -> Result<(), String>`;
  - `pub fn dir_of(path: &Path) -> &Path`: the parent folder, with `""`
    read as `"."`.

**What they do:**
- **`for_each_asset`** visits the `file` of every `ShapeKind::Obj` and
  every `UvDescription::UvImage`. It recurses through:
  - `Group.children` and `Csg.left`/`right`;
  - `material.pattern`;
  - `Nested.a`/`b`, `Blended.a`/`b` and `Perturbed.pattern`;
  - `TextureMap.uv`;
  - all six `CubeMap` faces.
- **`rebase`** leaves a path alone when it is absolute. Otherwise:
  - canonicalize `from_dir.join(path)`; on error, return `Err` naming the
    path as written;
  - canonicalize `to_dir`;
  - drop their common leading components;
  - write one `..` for each component left in `to_dir`, then the rest of
    the asset's path, joined with `/`.

  `rebase` either rewrites every path or returns `Err` and leaves `desc`
  unchanged: work on a clone and assign it back on success.
- **`absolutize`** rewrites each relative path as
  `canonicalize(from_dir.join(path))`, with an empty `from_dir` read as
  `.`. Failures are reported in the same way.

- [ ] **Step 1: Write the failing tests** in `paths.rs`.

Each test builds a fresh folder under `std::env::temp_dir()`, named with
the test's name, and creates empty asset files there
(`models/m.obj`, `tex/t.ppm`).

- **`rebase_into_a_child_folder`:** `models/m.obj` from `root` to
  `root/out` becomes `../models/m.obj`.
- **`rebase_into_the_parent_folder`:** `models/m.obj` from `root/scenes`
  (asset at `root/scenes/models/m.obj`) to `root` becomes
  `scenes/models/m.obj`.
- **`rebase_into_a_sibling_folder`:** `../models/m.obj` from `root/a` to
  `root/b` stays `../models/m.obj`.
- **`rebase_leaves_absolute_paths_alone`.**
- **`rebase_reaches_nested_assets`:** one scene with:
  - an `obj` inside a group inside a CSG's `right`;
  - `uv_image`s in a texture map, in a cube map's `down` face, and in a
    texture map wrapped in a `perturbed` pattern.

  Every one is rewritten to `../…`.
- **`rebase_names_a_missing_asset_and_changes_nothing`:** `Err`
  contains `missing.obj`, and `desc` is unchanged.
- **`rebasing_every_scene_into_its_own_folder_changes_nothing`:** for
  every `*.json` in `../rtc-metal/scenes` and `../rtc-rust/scenes`:
  - `rebase(from = to = its folder)` leaves an equal description;
  - its `to_json` keeps as many `"comment"` keys as the file has (a
    local `comment_keys` helper, like Task 1's);
  - at least 20 files are checked.
- **`absolutize_gives_a_scene_that_loads_from_anywhere`:** load
  `../rtc-metal/scenes/chapter15.json`, absolutize it, save it into the
  temp folder, then `rtc_rust::scene::load` that file. It builds with
  the teapots.
- **`an_empty_folder_means_the_current_one`:** `dir_of("x.json")` is
  `"."`.

- [ ] **Step 2: Run them to verify they fail**

Run: `cargo test` (from `RUST/rtc-edit`)
Expected: compile errors for the missing functions. Add stubs that
return `Ok(())` or do nothing, then each test fails on its assertion,
except `rebase_leaves_absolute_paths_alone` and
`rebasing_every_scene_into_its_own_folder_changes_nothing`, which pass
against the stubs. Note that, and keep them: their job is guarding
against changes later.

- [ ] **Step 3: Implement `paths.rs`** as described above.

- [ ] **Step 4: Run the tests**

Run: `cargo test && cargo clippy --all-targets`
Expected: all pass; clippy clean.

- [ ] **Step 5: Commit**
`ENHANCEMENT: rtc-edit crate rebases a scene's asset paths (milestone 17)`

---

### Task 3: the tree

**Files:**
- Create: `RUST/rtc-edit/src/tree.rs`; add `pub mod tree;` to `lib.rs`.

**Interfaces:**
- Produces, in `rtc_edit::tree`:
  - `#[derive(Debug, Clone, PartialEq, Eq, Hash)] pub enum Node { Camera, Settings, Lights, Light(usize), Shapes, Shape(Vec<usize>) }`.
    `Lights` and `Shapes` are the section rows. A `Shape` path goes
    top-level index first; in a group the child index follows, and in a
    CSG 0 is left and 1 is right.
  - `#[derive(Debug, Clone, PartialEq)] pub struct Row { pub node: Node, pub depth: usize, pub label: String, pub foldable: bool }`.
  - `pub fn rows(desc: &SceneDescription, folded: &HashSet<Node>) -> Vec<Row>`.
  - `pub fn json(desc: &SceneDescription, node: &Node) -> String`.

**Rows and labels:**
- The order is:
  - `camera`, at depth 0;
  - `settings`, at depth 0, only when it is not `Settings::default()`;
  - `lights (n)`, at depth 0 and foldable, then one row per light at
    depth 1 labelled `point`, `area` or `spot`. The lights come from
    `light` or from `lights`;
  - `shapes (n)`, at depth 0 and foldable, then the shapes from depth 1,
    each group's or CSG's children one level deeper.
- **A shape's label** is a base, then ` — <comment>` when the comment is
  a JSON string. The base is:
  - for a group: `group (<n> children)`;
  - for a CSG: `csg <union|intersection|difference>`;
  - for anything else: the `"type"` tag (`sphere`, `obj`, …).
- **`foldable`** is true for the two sections, groups and CSGs.
- **A folded node** keeps its own row but hides everything under it.

**`json`:**
- `Camera`, `Settings`, `Light(i)` and `Shape(path)` give
  `serde_json::to_string_pretty` of that item.
- A group's `"children"` and a CSG's `"left"`/`"right"` are removed from
  the `serde_json::Value` first.
- `Lights` and `Shapes` give an empty string.

- [ ] **Step 1: Write the failing tests**

All of them use one scene from a `json!` literal: two point lights, path
tracing settings, and these shapes:
- a sphere with comment `"ball"`;
- a group of a cube and a plane;
- a CSG difference of two spheres.

The tests:
- **`rows_list_the_scene_in_order`:** the labels and depths are exactly:
  - `camera` 0, `settings` 0, `lights (2)` 0, `point` 1, `point` 1;
  - `shapes (3)` 0, `sphere — ball` 1, `group (2 children)` 1, `cube` 2,
    `plane` 2;
  - `csg difference` 1, `sphere` 2, `sphere` 2.

  The nodes of the last two rows are `Shape(vec![2, 0])` and
  `Shape(vec![2, 1])`.
- **`folding_shapes_hides_every_shape`:** with `Node::Shapes` folded,
  the last row is `shapes (3)`.
- **`folding_a_group_hides_its_children`:** with `Shape(vec![1])`
  folded, `cube` and `plane` are gone, and `csg difference` remains.
- **`default_settings_have_no_row`.**
- **`a_single_light_is_listed`:** a scene with `"light"` gives
  `lights (1)` and `point`.
- **`json_leaves_out_a_groups_children`:** `json(Shape(vec![1]))` parses
  to an object with `"type": "group"` and no `"children"`.
- **`json_of_a_nested_shape`:** `json(Shape(vec![2, 1]))` has
  `"type": "sphere"`.

- [ ] **Step 2: Run them to verify they fail**

Run: `cargo test tree`
Expected: compile errors, then assertion failures against stubs.

- [ ] **Step 3: Implement `tree.rs`.**

- [ ] **Step 4: Run the tests**

Run: `cargo test && cargo clippy --all-targets`
Expected: all pass; clippy clean.

- [ ] **Step 5: Commit**
`ENHANCEMENT: rtc-edit lists a scene as tree rows (milestone 17)`

---

### Task 4: the preview file and the viewer process

**Files:**
- Create: `RUST/rtc-edit/src/preview.rs`, `RUST/rtc-edit/src/viewer.rs`;
  add both modules to `lib.rs`.

**Interfaces:**
- Consumes: `paths::absolutize` and `paths::dir_of` (Task 2).
- Produces:
  - **`preview::Preview`:**
    - `pub fn new(root: &Path, scene_name: &OsStr) -> Result<Preview, String>`
      creates `root/rtc-edit-<pid>/`. `main` passes
      `std::env::temp_dir()`; tests pass their own folder, so parallel
      tests don't share one.
    - `pub fn path(&self) -> &Path` is `<that folder>/<scene_name>`.
    - `pub fn write(&self, desc: &SceneDescription, scene_dir: &Path) -> Result<(), String>`
      absolutizes a clone and `save`s it to `path()`.
    - `Drop` removes the folder, ignoring errors.
  - **`viewer`:**
    - `pub fn find(env: Option<&str>, crate_dir: &Path) -> Result<PathBuf, String>`.
      An env value that is set and non-empty wins, without an existence
      check, because the user asked for it. Otherwise the result is
      `crate_dir/../rtc-metal/target/release/rtc-metal`, or
      `Err("rtc-metal not found at <path>; build it or set RTC_METAL")`.
    - `pub struct Viewer`.
    - `pub fn start(exe: &Path, preview: &Path) -> Result<Viewer, String>`
      spawns `exe <preview> --window --watch`, with stdin null and
      stdout/stderr inherited. On error it returns
      `Err("starting <exe>: <error>")`.
    - `pub fn running(&mut self) -> bool` is `try_wait` returning
      `Ok(None)`.
    - `pub fn id(&self) -> u32`.
    - `Drop` kills and waits, ignoring errors.

- [ ] **Step 1: Write the failing tests**

`preview`:
- **`write_puts_a_loadable_scene_in_the_folder`:** the chapter15
  description, written, loads with `rtc_rust::scene::load(p.path())`.
  `p.path()` ends with `chapter15.json`, and its parent's name starts
  with `rtc-edit-`.
- **`dropping_the_preview_removes_its_folder`.**

`viewer`. A test script is a file in the test's temp folder with mode
0o755 (`std::os::unix::fs::PermissionsExt`). Tests:
- **`find_prefers_the_environment`:** `find(Some("/x/rtc-metal"), any)`
  is `/x/rtc-metal`.
- **`find_falls_back_to_the_sibling_crate`:** with a fake
  `root/rtc-metal/target/release/rtc-metal` file and
  `crate_dir = root/rtc-edit`, the result is that path. With `Some("")`,
  the same.
- **`find_reports_a_missing_viewer`:** `Err` contains `RTC_METAL`.
- **`running_is_true_while_the_child_runs`:** a script `exec sleep 30`;
  `running()` is true.
- **`running_is_false_after_the_child_exits`:** a script `exit 0`; poll
  `running()` for up to 2 s, and it becomes false.
- **`dropping_the_viewer_kills_it`:** after the drop,
  `kill -0 <id>` (via `std::process::Command`) fails.
- **`start_reports_a_missing_program`:** `Err` contains the exe path.

- [ ] **Step 2: Run them to verify they fail**

Run: `cargo test preview viewer`
Expected: compile errors, then assertion failures against stubs.

- [ ] **Step 3: Implement `preview.rs` and `viewer.rs`.**

- [ ] **Step 4: Run the tests**

Run: `cargo test && cargo clippy --all-targets`
Expected: all pass; clippy clean.

- [ ] **Step 5: Commit**
`ENHANCEMENT: rtc-edit writes a preview file and runs the viewer (milestone 17)`

---

### Task 5: the iced editor

**Files:**
- Modify: `RUST/rtc-edit/Cargo.toml` (add `iced = "0.14"`,
  `rfd = "0.15"`)
- Create: `RUST/rtc-edit/src/main.rs`
- Create: `RUST/rtc-edit/README.md`
- Modify: `Readme.adoc` (top level: one rtc-edit bullet next to
  rtc-metal's), `RUST/rtc-metal/README.md` (one line: rtc-edit drives
  `--watch`)

**Interfaces:**
- Consumes: `paths::{rebase, dir_of}`, `tree::{rows, json, Node, Row}`,
  `preview::Preview`, and `viewer::{find, Viewer}`.

**The application** is `iced::application(boot, update, view)` with:
- `.title(…)`, giving `rtc-edit — <file name>`;
- `.subscription(…)`;
- `.exit_on_close_request(false)`;
- `.window_size((900.0, 640.0))`.

State:
- `path: PathBuf`;
- `desc: SceneDescription`;
- `folded: HashSet<Node>`;
- `selected: Option<Node>`;
- `preview: Option<Preview>`;
- `viewer: Option<Viewer>`;
- `status: String`.

Messages:
- `Select(Node)` and `Toggle(Node)`;
- `Save`, `SaveAs` and `ShowViewer`;
- `CloseRequested(window::Id)`.

**Behaviour:**
- **`boot`** reads the scene path from `args().nth(1)`. Without one it
  prints `usage: rtc-edit <scene.json>` and exits with 2. If
  `SceneDescription::from_path` fails, it prints the error and exits
  with 1. Then it creates the preview, writes it, and runs the
  `ShowViewer` logic, which sets the status line.
- **`ShowViewer`:**
  - if a viewer is running, do nothing;
  - otherwise
    `find(env::var("RTC_METAL").ok().as_deref(), Path::new(env!("CARGO_MANIFEST_DIR")))`,
    then `Viewer::start(exe, preview.path())`;
  - the status is `viewer started` or the error;
  - with no preview (it failed to write), the status is
    `no preview file: <error>`.
- **`Save`:** `desc.save(&path)`; the status is `saved <file name>` or
  the error.
- **`SaveAs`:**
  - open `rfd::FileDialog::new()` with `set_directory(dir_of(&path))`,
    `set_file_name(<current name>)`, `add_filter("scene", &["json"])`,
    then `save_file()`;
  - `None` (cancelled) does nothing;
  - `Some(new)`: clone `desc`, `rebase(clone, dir_of(&path), dir_of(&new))`,
    then `save(&new)`. Only when both succeed, `desc = clone` and
    `path = new`;
  - the status is `saved <new name>` or the error.
- **The subscription** combines two:
  - `keyboard::listen()`, filtered to `KeyPressed` with
    `modifiers.command()`. The character, compared case-insensitively,
    decides: `s` gives Save, or SaveAs with shift; `r` gives ShowViewer;
  - `window::close_requests()`, giving `CloseRequested`.
- **`CloseRequested`:** drop the viewer and the preview (set them to
  `None`), then `iced::exit()`.

**`view`:**
- **The tree column** is 260 px wide inside a `scrollable`. Each row is a
  `button` showing `"  " × depth`, then `▾ `/`▸ ` when foldable, then the
  label. Pressing a row selects it; pressing a foldable row's marker
  toggles it. The selected row uses `button::primary`, the others
  `button::text`. Rows are 14 px text.
- **The inspector** is `text(json(&desc, node))` in
  `iced::font::Font::MONOSPACE`, inside a `scrollable`. When nothing is
  selected, it shows `select an item`.
- **The status line** is a 13 px `text` below both columns.

**README:**
- what rtc-edit is, and the milestone plan (17 read-only, 18 editing,
  19 add/delete and New/Open);
- the dependencies;
- how to build and run, including building rtc-metal first and
  `RTC_METAL`;
- the keys;
- the preview file's location.

- [ ] **Step 1: Build**

Run: `cargo build --release && cargo clippy --all-targets && cargo test`
Expected: builds; clippy clean; Tasks 2–4's tests still pass. `main.rs`
has no unit tests: its logic is the glue between the tested modules, and
it is checked by eye.

- [ ] **Step 2: Launch check**

Run: `(cd ../rtc-metal && cargo build --release)`, then
`cargo run --release -- ../rtc-metal/scenes/chapter15.json` in the
background for 8 s, then stop it.
Expected:
- no panic on stderr;
- a `rtc-edit-<pid>` folder exists in `$TMPDIR` while it runs, holding
  `chapter15.json` with absolute `.obj` paths;
- after the editor is stopped with SIGTERM, the viewer process is gone.
  If SIGTERM skips `Drop`, note it: the by-eye check is the real gate.

- [ ] **Step 3: By eye (the user)**
1. Open `../rtc-metal/scenes/chapter15.json`. The tree shows `camera`,
   `lights (1)`/`point`, and shapes with their comments. The viewer
   renders the teapots.
2. Fold and unfold `shapes`; select rows; the inspector shows each
   item's JSON.
3. Cmd-Shift-S into another folder (e.g. `/tmp/x/`); then
   `rtc-metal /tmp/x/chapter15.json --window` shows the teapots.
4. Cancel a Save As: nothing changes.
5. Close the viewer, press Cmd-R: it comes back. Press Cmd-R again: no
   second viewer.
6. Quit with Cmd-Q, and in a second run by closing the window. Each time
   the viewer closes and the `rtc-edit-<pid>` folder is gone.

   If Cmd-Q leaves them behind (spec Risks), make Cmd-Q a key binding in
   the subscription that sends the same close path, and record a ruling.

- [ ] **Step 4: Commit** (after the user confirms)
`ENHANCEMENT: rtc-edit shows a scene as a tree and saves it (milestone 17)`
