# Milestone 20: New and Open, adding and deleting items — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** rtc-edit starts a new scene, opens another one without restarting, and adds and deletes shapes and lights.

**Architecture:** A new `src/edit.rs` holds the scene changes as pure functions (new scene, add, delete, fold shifting), tested without a window. `src/main.rs` wires them to a button bar, keys, Cmd-N and Cmd-O, and replaces the preview and viewer when the scene switches.

**Tech Stack:** Rust 2024, iced 0.14, rfd 0.15 (file dialogs and the alert), rtc-rust's `SceneDescription`.

**Spec:** `RUST/rtc-edit/docs/superpowers/specs/2026-10-18-metal-20-new-open-add-delete-design.md`

## Global Constraints

- No new crates. rtc-edit already depends on `rtc-rust`, `serde_json`, `iced`, `rfd`.
- Tests in rtc-edit use descriptive snake_case names; run with `cargo test --lib` in `RUST/rtc-edit` (several filters: `cargo test --lib -- a b`).
- `cargo clippy --all-targets` clean; `cargo fmt`. The "future-incompat block v0.1.6" warning is pre-existing.
- Commits: prefix `ENHANCEMENT:` / `TEST:` / `DOCS:` / `FIX:`, subject ends ` (milestone 20)`, message ends with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`; commit with `git commit -F <file>`. Never push.
- UI work (Tasks 3 and 4) is checked by the user by eye before it is committed: stop after the tests and README, report, and commit only when told.
- New scene: camera from (0, 1.5, −5), to (0, 1, 0), up (0, 1, 0), 800×600, field of view π/3; `light`: point, position (−10, 10, −10), intensity (1, 1, 1); shapes: one plane, no transform, no material. Name `untitled.json`.
- New shapes: material `color` (0.8, 0.8, 0.8) only; sphere and cube `translate [0, 1, 0]`; plane no transform; cylinder and cone `min 0`, `max 1`, closed, no transform. New light: point at (−10, 10, −10), intensity (1, 1, 1).
- Status texts, verbatim: `a CSG needs both its shapes`, `a scene needs a light`, `add: s sphere  c cube  p plane  y cylinder  o cone  l light`. Alert: `Save changes to <name>?` with buttons `Save`, `Don't Save`, `Cancel`.

## Review Focus

1. **Adding into a folded group** (or with Shapes/Lights folded). Expected: the new row is visible and selected — its ancestors unfold. Task 3, by eye.
2. **Delete while a box of the selected item is shaded or being edited.** Expected: no stale `left`/`editing` pointing at the deleted item's boxes; the next item is selected cleanly. Task 3 uses `select`, which resets them.
3. **Cmd still marked as held after an rfd dialog or alert** (the key-up goes to the dialog). Expected: shortcuts and typing work normally afterwards. Task 4 clears `command` after every dialog.
4. **Save chosen in the alert for an untitled scene, then Save As cancelled.** Expected: New/Open is cancelled and the edits stay. Task 4: `save_as` returns whether it saved.
5. **Opening a scene with relative assets** (e.g. `rtc-metal/scenes/chapter15.json`'s OBJ). Expected: the viewer finds them — the preview is absolutized against the *new* scene's folder. Task 4, by eye.

---

### Task 1: `edit.rs` — the new scene and adding

**Files:**
- Create: `RUST/rtc-edit/src/edit.rs`
- Modify: `RUST/rtc-edit/src/lib.rs` (add `pub mod edit;`)

**Interfaces:**
- Consumes: `rtc_edit::tree::Node` (`Camera, Settings, Lights, Light(usize), Shapes, Shape(Vec<usize>)`; a `Shape` path indexes `desc.shapes`, then a group's children in order, a CSG's left = 0 and right = 1).
- Produces:
  - `#[derive(Debug, Clone, Copy, PartialEq, Eq)] pub enum Item { Sphere, Cube, Plane, Cylinder, Cone, Light }`
  - `impl Item { pub fn from_key(c: char) -> Option<Item> }` — `s c p y o l`.
  - `pub fn new_scene() -> SceneDescription`
  - `pub fn add(desc: &mut SceneDescription, selected: Option<&Node>, item: Item) -> Node`

**Placement rule (the spec's CSG wording made precise):** a CSG's two children are fixed, everything else can take more. With a shape selected: if it is a group, append to its children; otherwise insert right after it in its parent — and while that parent is a CSG, use the CSG itself as "the selected shape" and climb (a group inside a CSG still takes children; a shape directly in a CSG goes after the CSG, and after the outermost one for nested CSGs). Any other selection, or none: append to `desc.shapes`. A light always goes at the end of the lights; `light: Some(l)` becomes `lights: Some(vec![l, new])`, `light: None`.

- [ ] **Step 1: Write the failing tests** (in `edit.rs`'s `#[cfg(test)] mod tests`, with a fixture like `tree.rs`'s: shapes `[sphere, group[cube, plane], csg(left: group[sphere], right: sphere)]`, two point lights; a second fixture with a single `light`):
  - `the_new_scene_builds_and_reads_back_equal`: `new_scene().build(Path::new("."))` is `Ok`; `SceneDescription::from_json(&d.to_json())` equals `d`; camera, light and the one plane match the Global Constraints values.
  - `a_shape_goes_after_the_selected_top_level_shape`: select `Shape([0])`, add Cube → returns `Shape([1])`, `desc.shapes[1]` is a cube, 4 top-level shapes.
  - `a_shape_goes_after_the_selected_shape_inside_a_group`: select `Shape([1, 0])` → `Shape([1, 1])`, group has 3 children.
  - `a_shape_goes_at_the_end_of_a_selected_group`: select `Shape([1])` → `Shape([1, 2])`.
  - `a_shape_inside_a_csg_goes_after_the_csg`: select `Shape([2, 1])` → `Shape([3])`, the CSG unchanged.
  - `a_group_inside_a_csg_still_takes_a_shape`: select `Shape([2, 0])` → `Shape([2, 0, 1])`.
  - `with_no_shape_selected_a_shape_goes_at_the_end`: for `None`, `Camera`, `Lights`, `Light(0)`, `Shapes` → `Shape([3])` (fresh fixture each).
  - `added_shapes_have_their_defaults_and_build`: each of Sphere, Cube, Plane, Cylinder, Cone: the kind, transform (`[Translate([0.,1.,0.])]` for sphere/cube, empty otherwise), cylinder/cone `min Some(0.)`, `max Some(1.)`, `closed true`, material `color Some([0.8; 3])` and every other material field `None`; `desc.build(Path::new("."))` is `Ok`.
  - `adding_a_light_to_a_single_light_makes_a_list`: single-light fixture, add Light → returns `Light(1)`, `light` is `None`, `lights` has 2 with the old first.
  - `a_light_goes_at_the_end_of_the_lights`: two-light fixture → `Light(2)`.
  - `item_keys`: `from_key` maps `s c p y o l` and gives `None` for `x`.
- [ ] **Step 2: Run** `cargo test --lib edit` in `RUST/rtc-edit`. Expected: compile errors / FAIL (no `edit` items yet).
- [ ] **Step 3: Implement** `Item`, `from_key`, `new_scene`, `add` in `src/edit.rs`. Reach a shape mutably by walking the path (a small private `shape_at_mut`, or a helper returning the parent's child `Vec` — a CSG has no `Vec`, which is why the climb happens first).
- [ ] **Step 4: Run** `cargo test --lib` and `cargo clippy --all-targets`. Expected: all pass, no new warnings.
- [ ] **Step 5: Commit** `ENHANCEMENT: rtc-edit makes a new scene and adds shapes and lights to one (milestone 20)`.

### Task 2: `edit.rs` — deleting and folds

**Files:**
- Modify: `RUST/rtc-edit/src/edit.rs`

**Interfaces:**
- Consumes: Task 1's fixtures and helpers.
- Produces:
  - `pub fn can_delete(desc: &SceneDescription, node: &Node) -> Result<(), String>` — `Err("a CSG needs both its shapes")` for a shape whose parent is a CSG, `Err("a scene needs a light")` for a light when it is the only one, `Err` with a short reason (e.g. `nothing to delete`) for `Camera`, `Settings`, `Lights`, `Shapes`, a path that names nothing.
  - `pub fn delete(desc: &mut SceneDescription, node: &Node) -> Result<Node, String>` — checks `can_delete` first; returns the node to select: the next sibling (same path), else the previous one, else the parent (`Shape(prefix)`, or `Shapes`; for lights `Light(i)` or `Light(i - 1)`). When `lights` drops to one it stays a one-item `lights` list.
  - `pub fn shift(folded: &mut HashSet<Node>, at: &[usize], by: isize)` — `by = 1` after an insert at `at`; `by = -1` for a delete at `at`. A folded `Shape(p)` whose path has `at`'s parent as prefix and whose component at that depth is ≥ `at`'s last index moves by `by` there (its descendants with it); for `by = -1` the folds of `at` itself and below are dropped first. Other folds stay.

(The spec named `can_delete` as `-> bool`; returning the reason lets `delete` and the status line share it — the button uses `.is_ok()`.)

- [ ] **Step 1: Write the failing tests:**
  - `deleting_a_shape_in_a_group_selects_the_next_one`: delete `Shape([1, 0])` → `Ok(Shape([1, 0]))`, group has `[plane]`.
  - `deleting_the_last_shape_of_a_level_selects_the_previous_one`: delete `Shape([2])` → `Ok(Shape([1]))`.
  - `deleting_an_only_child_selects_its_group`: group with one child, delete it → `Ok(Shape([1]))`; delete the only top-level shape → `Ok(Shapes)`.
  - `deleting_a_group_removes_its_children`: delete `Shape([1])` → 2 top-level shapes left.
  - `a_csg_child_cant_be_deleted`: `Shape([2, 0])` → `Err("a CSG needs both its shapes")`, scene equal to before; `can_delete` agrees.
  - `the_last_light_cant_be_deleted`: single-light fixture, `Light(0)` → `Err("a scene needs a light")`, unchanged.
  - `deleting_a_light_selects_the_next_or_previous`: two lights, delete `Light(1)` → `Ok(Light(0))`, `lights` has one.
  - `heading_rows_cant_be_deleted`: `Camera`, `Settings`, `Lights`, `Shapes` → `Err`, unchanged.
  - `a_folded_group_stays_folded_after_an_insert_before_it`: folded `{Shape([1]), Shape([1, 0])}`, `shift(&mut f, &[0], 1)` → `{Shape([2]), Shape([2, 0])}`.
  - `a_folded_group_stays_folded_after_a_delete_before_it`: folded `{Shape([2])}`, `shift(&mut f, &[0], -1)` → `{Shape([1])}`.
  - `deleting_a_folded_group_drops_its_folds`: folded `{Shape([1]), Shape([1, 0]), Shape([2])}`, `shift(&mut f, &[1], -1)` → `{Shape([1])}` (the old `[2]`).
  - `folds_elsewhere_dont_move`: folded `{Lights, Shape([0]), Shape([1, 0])}`, `shift(&mut f, &[1, 1], 1)` → unchanged.
- [ ] **Step 2: Run** `cargo test --lib edit`. Expected: FAIL (functions missing).
- [ ] **Step 3: Implement** `can_delete`, `delete`, `shift`.
- [ ] **Step 4: Run** `cargo test --lib` and `cargo clippy --all-targets`. Expected: pass, clean.
- [ ] **Step 5: Commit** `ENHANCEMENT: rtc-edit deletes shapes and lights, and folds follow their items (milestone 20)`.

### Task 3: adding and deleting in the window

**Files:**
- Modify: `RUST/rtc-edit/src/main.rs`, `RUST/rtc-edit/README.md`

**Interfaces:**
- Consumes: `edit::{Item, add, can_delete, delete, shift}` (Tasks 1–2); existing `Editor::select(node)` (resets `left`, `editing`, texts, errors, remembered blur).
- Produces (for Task 4): `Editor.pending: Option<char>`, reset by Task 4's `switch`.

- [ ] **Step 1: Keys.** In `shortcut`, plain characters no box took — `h j k l i c x u a d s p y o` — become one new `Message::Key(char)`, replacing today's `Move` for h/j/k/l and `Letter` for i/c/x/u (arrow keys still send `Move`). `Editor::key(c)`:
  - `pending == Some('a')`: clear it; `Item::from_key(c)` → `Add(item)`, otherwise cancel (status cleared).
  - `pending == Some('d')`: clear it; `d` → `Delete`.
  - otherwise h/j/k/l → `arrow`, i/c/x/u → `letter` (as today); `a` or `d` with no box shaded or edited (`left` and `editing` both `None`) and an item selected or `a` → set `pending`; `a` sets the status to `add: s sphere  c cube  p plane  y cylinder  o cone  l light`.
  - Every other message except `Command` and `Focused` clears `pending`.
- [ ] **Step 2: `Editor::add_item(item)`**: `let node = edit::add(&mut desc, selected, item)`; for a new `Shape(p)`, `shift(&mut folded, &p, 1)` and unfold its ancestors (each `Shape(prefix)` and `Shapes`; `Lights` for a light); `select(node)`; `unsaved = unloaded = true`; status as after an edit (the `not in viewer — Cmd-L to load` wording `shown_status` already uses).
- [ ] **Step 3: `Editor::delete_selected()`**: on `Ok(next)`: for a deleted `Shape(p)`, `shift(&mut folded, &p, -1)`; `select(next)`; `unsaved = unloaded = true`. On `Err(reason)`: status = reason, scene unchanged.
- [ ] **Step 4: Button bar** above the tree: `+ sphere`, `+ cube`, `+ plane`, `+ cylinder`, `+ cone`, `+ light` sending `Add(item)`, and `Delete` sending `Delete` only while `selected` passes `can_delete` (no `on_press` otherwise, so iced greys it).
- [ ] **Step 5: README**: the button bar, `a` + letter, `d d`, what is refused and why, "no undo: open the file again and choose Don't Save" (the alert arrives in Task 4; word it so it reads right after Task 4 too).
- [ ] **Step 6: Run** `cargo test --lib`, `cargo clippy --all-targets`, `cargo build --release`. Expected: pass, clean.
- [ ] **Step 7: Stop for the by-eye check** (report; do not commit). The user checks: in `chapter10.json`, `a s` then Cmd-L shows a sphere at (0, 1, 0); `+ cylinder` with a group's child selected in `chapter16.json` (or `rtc-rust/scenes/dice-sentence.json`) lands in the group; `+ light` brightens; `d d` deletes and selects the next row; a CSG child and the last light are refused with the reason; adding into a folded group shows the new row; `a` then `x` cancels; h/j/k/l, i/c/x/u still work.
- [ ] **Step 8: Commit when told** `ENHANCEMENT: rtc-edit adds and deletes shapes and lights from a button bar and keys (milestone 20)`.

### Task 4: New, Open, and switching scenes

**Files:**
- Modify: `RUST/rtc-edit/src/main.rs`, `RUST/rtc-edit/README.md`

**Interfaces:**
- Consumes: `edit::new_scene` (Task 1); `Editor.pending` (Task 3); `Preview::new(root, scene_name)`, `Preview::write(desc, scene_dir)`, `Viewer` (dropping one sends SIGTERM if it still runs), `CHILDREN`.
- Produces: `Editor.path: Option<PathBuf>`; `fn switch(&mut self, desc: SceneDescription, path: Option<PathBuf>)`; `fn ask_unsaved(&mut self) -> bool`; `fn save_as(&mut self) -> bool` (true when saved); messages `New`, `Open`.

- [ ] **Step 1: `path: Option<PathBuf>`.** A helper `fn name(&self) -> String` gives the file name or `untitled.json` (title, alert, preview); `fn dir(&self) -> PathBuf` gives `dir_of(path)` or `std::env::current_dir()` (Save As, Open, preview, rebase).
- [ ] **Step 2: `switch`.** Lock `CHILDREN`; set `viewer = None` then `preview = None` (SIGTERM, folder removed); make and write a new `Preview` named `self.name()` against the new scene's dir; then set `desc`, `path`, and reset `folded`, `selected`, `texts`, `errors`, `left`, `editing`, `removed_blur`, `pending`, `unsaved`, `unloaded`, `preview_error`; release the lock and `show_viewer()`. `boot` builds the `Editor` once and calls `switch` — with `new_scene()` and `None` when there is no argument, else the file (a file that doesn't load still exits with the error, as today).
- [ ] **Step 3: Save on untitled.** `save()` with `path == None` calls `save_as()`. `save_as` defaults to `self.dir()` and `self.name()`, rebases from `self.dir()`, and returns whether it saved.
- [ ] **Step 4: `ask_unsaved`.** Not `unsaved` → `true`. Else `rfd::MessageDialog` (`MessageLevel::Warning`, description `Save changes to <name>?`, `MessageButtons::YesNoCancelCustom("Save".into(), "Don't Save".into(), "Cancel".into())`): `Custom("Save")` → `save()` then `!self.unsaved`; `Custom("Don't Save")` → `true`; anything else → `false`. Check by eye which `MessageDialogResult` variants rfd 0.15 returns for custom buttons on macOS and match those.
- [ ] **Step 5: New and Open.** Cmd-N → `New`: `ask_unsaved()` then `switch(new_scene(), None)`. Cmd-O → `Open`: `ask_unsaved()`, then `rfd::FileDialog` in `self.dir()` filtered to `json`, `pick_file`; `SceneDescription::from_path` → `switch(desc, Some(path))`, or on error the status line says `<file>: <error>` and nothing changes. After every rfd dialog or alert, set `command = false`.
- [ ] **Step 6: README**: `rtc-edit` with no argument, Cmd-N, Cmd-O, the alert, Cmd-S on an untitled scene, the old viewer offering its image on a switch.
- [ ] **Step 7: Run** `cargo test --lib`, `cargo clippy --all-targets`, `cargo build --release` (and `cargo build --release` in `RUST/rtc-metal`). Expected: pass, clean.
- [ ] **Step 8: Stop for the by-eye check** (report; do not commit). The user checks: `cargo run --release` with no argument shows the new scene and the viewer renders the floor; Cmd-S asks where to save; with an edit, Cmd-O shows the alert — Cancel keeps everything, Don't Save switches, the old viewer offers its image and a new viewer opens; Save on an untitled scene then cancelling Save As keeps the scene; opening `rtc-metal/scenes/chapter15.json` renders its teapot; Cmd-N gives `untitled.json`; Cmd-S/Cmd-L still work after a dialog.
- [ ] **Step 9: Commit when told** `ENHANCEMENT: rtc-edit starts new scenes and opens others, asking to save edits first (milestone 20)`.

### Task 5: a comment box for shapes

Spec: the amendment "a comment box (Task 5, 2026-10-19)".

**Files:**
- Modify: `RUST/rtc-edit/src/fields.rs`, `RUST/rtc-edit/src/main.rs`, `RUST/rtc-edit/README.md`

**Interfaces:**
- Consumes: `fields::{Field, fields, set}`, `main.rs`'s `describe` and `boxes`.
- Produces: `Field::Comment`.

- [ ] **Step 1: Write the failing tests** in `fields.rs`:
  - `a_shape_lists_its_comment_first`: a sphere with `"comment": "ball"` → `fields(..)[0] == (Field::Comment, "ball")`; one without → `(Field::Comment, "")`.
  - `a_group_and_a_csg_list_only_a_comment_box`: `fields` of a group and of a CSG → `vec![(Field::Comment, …)]`.
  - `a_non_string_comment_gets_no_box`: a sphere with `"comment": ["a", "b"]` → no `Field::Comment` in `fields`.
  - `typing_a_comment_sets_it`: `set(.., Field::Comment, "red ball")` → `comment == Some(json!("red ball"))`; on a group as well.
  - `an_empty_comment_removes_it`: `""` and `"   "` → `comment == None`, and the saved JSON has no `comment` key.
  - `lights_and_the_camera_have_no_comment_box`.
- [ ] **Step 2: Run** `cargo test --lib fields`. Expected: FAIL (`Field::Comment` missing).
- [ ] **Step 3: Implement** `Field::Comment`: `fields` puts it first for every shape whose comment is a string or absent (non-editable shapes get only it); `set` handles it before `parse`, for any shape, storing the text as typed, removing the key when the text trims empty.
- [ ] **Step 4: The box.** `describe(Comment)` → no heading, row label `comment`; `boxes` draws it about 248 wide (three 80-wide boxes and their spacing) and not monospace. Check that paste, Submit, Esc/u reset and `next_box` treat it like any box, and that nothing numeric (e.g. the paste filter, if any) rejects its text.
- [ ] **Step 5: README**: the comment box and that the tree shows it.
- [ ] **Step 6: Run** `cargo test --lib`, `cargo clippy --all-targets`, `cargo fmt`, `cargo build --release`. Expected: pass, clean.
- [ ] **Step 7: Stop for the by-eye check** (report; do not commit). The user checks: `a s`, then type a comment: the row shows `sphere — <comment>`; Cmd-S writes it; emptying it removes it; a group in `chapter16.json` shows a comment box alone; arrows and h/j/k/l reach the box.
- [ ] **Step 8: Commit when told** `ENHANCEMENT: rtc-edit edits a shape's comment, which labels its row (milestone 20)`.
