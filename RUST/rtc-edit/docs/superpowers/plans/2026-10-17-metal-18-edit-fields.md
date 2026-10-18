# rtc-edit, editing with fixed fields — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** the inspector gets text boxes for a simple shape's transform
and material, a cylinder's or cone's bounds, and a point light; edits
build up in the editor, and Cmd-L loads them into the viewer.

**Architecture:**
- **Two new library modules, tested without a window:** `trs.rs` reads a
  transform list as scale, rotation and translation and writes it back;
  `fields.rs` lists an item's boxes with their text and writes parsed
  text into the scene.
- **`main.rs`** keeps the typed text per box, marks unsaved and unloaded
  edits, and adds Cmd-L.
- rtc-rust and rtc-metal don't change.

**Tech Stack:** Rust 2024; rtc-rust's scene description types; iced 0.14
(`text_input`, `checkbox`); serde_json.

**Spec:** `RUST/rtc-edit/docs/superpowers/specs/2026-10-17-metal-18-edit-fields-design.md`

## Global Constraints

- **Branch** `metal/18-edit-fields`. Never push. Commits end with
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Commit subjects** end with `(milestone 18)`; prefixes
  `ENHANCEMENT:` / `TEST:` / `DOCS:` / `FIX:`.
- **Test names:** descriptive snake_case.
- **Dependencies:** none added. rtc-edit keeps `rtc-rust`, `serde_json`,
  `iced = "0.14"`, `rfd = "0.15"`.
- **Editable items:** sphere, cube, plane, cylinder and cone (transform
  and material; cylinder and cone also min, max, closed), and point
  lights (position, intensity); also inside groups and CSG. Nothing else.
- **Transform order:** scale, rotate x, rotate y, rotate z, translate,
  each at most once.
- **Rotation boxes** show degrees rounded to 6 decimals; radians are
  stored, and replaced only when that box is edited.
- **Keys:** Cmd-L loads, Cmd-S saves, Cmd-Shift-S Save As, Cmd-R shows
  the viewer.
- **Title:** `rtc-edit — <file name>`, followed by ` •` while there are
  unsaved edits.
- **Status line** shows `not in viewer — Cmd-L to load` while there are
  edits the viewer hasn't loaded.
- **Commands** (from `RUST/rtc-edit`): `cargo test --lib`,
  `cargo clippy --all-targets`; `cargo fmt` before committing.
- **YAGNI:** no camera editing, no undo, no quit prompt, no patterns.

## Review Focus

1. **`inf`, `NaN` and `1e400` typed into a box.** Rust's `f64` parser
   accepts them. Expected: rejected, box red, scene unchanged. Task 2's
   `set_rejects_text_that_isnt_a_finite_number` pins it.
2. **Editing one rotation box on a shape with unrounded radians.**
   Expected: the other two rotations keep their stored radians exactly.
   Task 2's `editing_one_rotation_keeps_the_others_exact` pins it.
3. **Clearing every material box.** Expected: the shape's `material`
   goes back to absent rather than saving `"material": {}`. Task 2's
   `clearing_the_last_material_value_removes_the_material` pins it.
4. **A scene with a single `"light"`** rather than a `"lights"` list.
   Expected: `Light(0)` is editable and the edit lands in `light`.
   Task 2's `a_single_light_is_editable` pins it.
5. **A Cmd shortcut while a box has focus.** Expected: Cmd-S, Cmd-L and
   Cmd-R still work. Checked by eye in Task 3.

---

### Task 1: transforms as scale, rotation and translation

**Files:**
- Create: `RUST/rtc-edit/src/trs.rs`; add `pub mod trs;` to `lib.rs`.

**Interfaces:**
- Produces, in `rtc_edit::trs`:
  - `#[derive(Debug, Clone, Copy, PartialEq)] pub struct Trs { pub scale: [f64; 3], pub rotate: [f64; 3], pub translate: [f64; 3] }`,
    rotations in radians. `Default` is scale `[1, 1, 1]`, the rest zero.
  - `pub fn from_ops(ops: &[TransformOp]) -> Option<Trs>` (an associated
    fn on `Trs`): `None` unless the ops come in the order scale,
    rotate x, rotate y, rotate z, translate, each at most once. Shear and
    Matrix never fit.
  - `pub fn to_ops(&self) -> Vec<TransformOp>`: that order, leaving out
    scale `[1, 1, 1]`, a rotation of exactly `0.0` and translate
    `[0, 0, 0]`.

- [ ] **Step 1: Write the failing tests**
  - **`an_empty_list_is_the_default`:** `from_ops(&[])` is
    `Some(Trs::default())`, and `Trs::default().to_ops()` is empty.
  - **`a_full_list_reads_and_writes_back`:** `[Scale([2,2,2]),
    RotateX(0.1), RotateY(0.2), RotateZ(0.3), Translate([1,2,3])]` reads
    as those values, and `to_ops` returns the same list.
  - **`a_partial_list_fills_in_the_rest`:** `[RotateY(0.5),
    Translate([0,1,0])]` reads as scale `[1,1,1]`, rotate `[0,0.5,0]`,
    translate `[0,1,0]`, and writes back unchanged.
  - **`no_op_steps_are_left_out_when_written`:** `Trs { scale:
    [1,1,1], rotate: [0,0,0], translate: [4,0,0] }.to_ops()` is
    `[Translate([4,0,0])]`.
  - **`lists_that_dont_fit_give_none`:** each of these gives `None`:
    `[Shear(..)]`, `[Matrix(..)]`, `[Translate(..), Scale(..)]`,
    `[RotateX(..), RotateX(..)]`, `[RotateZ(..), RotateX(..)]`.
  - **`every_scene_shape_reads_without_panicking`:** for every
    `scenes/*.json` in `../rtc-metal` and `../rtc-rust`, walk every shape
    including group and CSG children, and call `from_ops` on its
    transform. Whenever it gives `Some(t)` and the list holds no no-op
    step, `t.to_ops()` equals the list. Assert at least one shape in all
    gives `Some`.

- [ ] **Step 2: Run them to verify they fail**

Run: `cargo test --lib trs`
Expected: compile errors, then failures against stubs.

- [ ] **Step 3: Implement `trs.rs`.** `from_ops` walks the list with a
  stage index 0–4 that only moves forward: each op must belong to a later
  stage than the previous one.

- [ ] **Step 4: Run the tests**

Run: `cargo test --lib && cargo clippy --all-targets`
Expected: all pass; clippy clean.

- [ ] **Step 5: Commit**

`ENHANCEMENT: rtc-edit reads transforms as scale, rotation and translation (milestone 18)`

---

### Task 2: the fields of an item

**Files:**
- Create: `RUST/rtc-edit/src/fields.rs`; add `pub mod fields;` to
  `lib.rs`.

**Interfaces:**
- Consumes: `tree::Node`, `trs::Trs`.
- Produces, in `rtc_edit::fields`:
  - `#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)] pub enum Field { Position(usize), Rotation(usize), Scale(usize), Min, Max, Color(usize), Ambient, Diffuse, Specular, Shininess, Reflective, Transparency, RefractiveIndex, LightPosition(usize), LightIntensity(usize) }`;
    the `usize` is the axis, 0 to 2.
  - `pub fn fields(desc: &SceneDescription, node: &Node) -> Vec<(Field, String)>`.
  - `pub fn set(desc: &mut SceneDescription, node: &Node, field: Field, text: &str) -> Result<(), String>`.
  - `pub fn closed(desc: &SceneDescription, node: &Node) -> Option<bool>`:
    `Some` only for a cylinder or cone.
  - `pub fn set_closed(desc: &mut SceneDescription, node: &Node, closed: bool)`.

**Rules:**
- **Order of `fields`:** Position 0–2, Rotation 0–2, Scale 0–2 (only
  when `Trs::from_ops` gives `Some`); Min, Max (cylinder and cone);
  Color 0–2, Ambient, Diffuse, Specular, Shininess, Reflective,
  Transparency, RefractiveIndex. A point light gives LightPosition 0–2,
  LightIntensity 0–2. Anything else gives an empty list.
- **Text:** numbers use `{}` formatting of the `f64` (`1.0` shows as
  `1`). Rotations show `to_degrees()` rounded to 6 decimals. An unset
  value (material, min, max) shows as an empty string.
- **Parsing:** trim the text, parse `f64`, and reject anything not
  finite. An empty text is accepted only for Min, Max and the material
  fields, and means unset. A Color with one empty axis is an error.
- **Writing a transform field:** read `Trs::from_ops`, change the one
  component (Rotation: `deg.to_radians()`), write `to_ops`. A scale of 0
  on any axis is an error.
- **Writing a material field:** create `MaterialDescription::default()`
  when the shape has none. If the result equals the default (pattern and
  bump included), set `material` to `None`.
- **Lights** come from `light` when present, otherwise `lights[i]`, as in
  `tree`. Only `Point` is editable.
- **Shapes** are found by the node's path, as `tree::shape_at` does; add
  a mutable twin in `fields.rs`.
- Any `Err` leaves `desc` unchanged; a field the item doesn't have is an
  `Err`.

- [ ] **Step 1: Write the failing tests** on one `json!` scene: a sphere
  with `transform [scale [2,2,2], translate [0,1,0]]` and material color
  `[1, 0.2, 0.2]`; a cylinder with `min -1, max 1, closed true` and no
  material; a group holding a cube; a torus; an obj; a sphere with a
  shear; two lights, a point and an area.
  - **`a_sphere_lists_transform_and_material`:** the first nine entries
    are Position `0,1,0`, Rotation `0,0,0`, Scale `2,2,2`; Color texts
    are `1`, `0.2`, `0.2`; Ambient's text is empty.
  - **`a_cylinder_lists_its_bounds`:** includes Min `-1`, Max `1`, and
    `closed` is `Some(true)`.
  - **`a_point_light_lists_position_and_intensity`.**
  - **`a_shape_whose_transform_doesnt_fit_lists_material_only`:** the
    shear sphere's first field is `Color(0)`.
  - **`items_you_cant_edit_list_nothing`:** torus, obj, group, area
    light, `Camera`, `Settings`, `Lights`, `Shapes` all give `[]`.
  - **`setting_position_and_color_writes_the_scene`:** after
    `Position(0) = "3"` and `Color(1) = "0.5"`, the sphere's transform is
    `[Scale([2,2,2]), Translate([3,1,0])]`, its color `[1, 0.5, 0.2]`,
    and `desc.build(Path::new("."))` succeeds.
  - **`rotation_is_typed_in_degrees`:** `Rotation(1) = "90"` stores
    `RotateY(FRAC_PI_2)`, and `fields` shows `90`.
  - **`editing_one_rotation_keeps_the_others_exact`:** a shape with
    `[RotateX(0.785398), RotateZ(0.123456789)]`; after `Rotation(1) =
    "10"`, the RotateX and RotateZ values are bit-identical.
  - **`an_empty_material_box_unsets_the_value`:** `Color` set, then
    `Ambient = "0.3"`, then `Ambient = ""` gives `ambient: None`.
  - **`setting_a_value_creates_a_material`:** on the cylinder,
    `Reflective = "0.5"` gives a material with only `reflective`.
  - **`clearing_the_last_material_value_removes_the_material`:** then
    `Reflective = ""` gives `material: None`.
  - **`set_rejects_text_that_isnt_a_finite_number`:** `"abc"`, `"inf"`,
    `"NaN"`, `"1e400"` and `""` on Position(0) are `Err`, and `desc`
    equals its clone from before.
  - **`a_zero_scale_is_rejected`:** `Scale(2) = "0"` is `Err`, `desc`
    unchanged.
  - **`a_shape_in_a_group_is_editable`:** `Position(0) = "5"` on
    `Shape(vec![2, 0])` changes the cube.
  - **`a_single_light_is_editable`:** a scene with `"light"`;
    `LightPosition(1) = "7"` on `Light(0)` changes `desc.light`.
  - **`set_closed_toggles_the_caps`.**
  - **`whitespace_around_a_number_is_accepted`:** `" 2 "` sets 2.

- [ ] **Step 2: Run them to verify they fail**

Run: `cargo test --lib fields`
Expected: compile errors, then failures against stubs.

- [ ] **Step 3: Implement `fields.rs`.**

- [ ] **Step 4: Run the tests**

Run: `cargo test --lib && cargo clippy --all-targets`
Expected: all pass; clippy clean.

- [ ] **Step 5: Commit**

`ENHANCEMENT: rtc-edit lists an item's fields and writes them into the scene (milestone 18)`

---

### Task 3: the boxes in the editor, and Cmd-L

**Files:**
- Modify: `RUST/rtc-edit/src/main.rs`, `RUST/rtc-edit/README.md`.

**Interfaces:**
- Consumes: `fields::{Field, fields, set, closed, set_closed}`.

**Editor state and messages:**
- New fields on `Editor`: `texts: HashMap<Field, String>`,
  `errors: HashSet<Field>`, `unsaved: bool`, `unloaded: bool`.
- `Select(node)` also rebuilds `texts` from `fields(&desc, &node)` and
  clears `errors`.
- `Edit(Field, String)`: store the text; on `set` Ok remove the field
  from `errors` and set `unsaved` and `unloaded`; on `Err(e)` add it to
  `errors` and show `e` in the status line.
- `Closed(bool)`: `set_closed`, set `unsaved` and `unloaded`.
- `Load` (Cmd-L): write the preview through `CHILDREN`'s `Preview`
  (`write(&desc, dir_of(&path))`); on Ok clear `unloaded` and show
  `loaded into viewer`; on Err show the error. With no preview, show
  the same `no preview file: …` message as `show_viewer`.
- Save and Save As clear `unsaved` on success.
- When `errors` isn't empty, Load, Save and Save As still write, and
  their status message ends with ` (a field has an error)`.
- After Save As, rebuild `texts` (no paths in fields change, but keep
  the boxes in step with `desc`).

**View:**
- Inspector: for each `(Field, _)` in `fields`, a `text_input` with the
  stored text, `on_input(move |s| Message::Edit(field, s))`, monospace,
  about 80 px wide, with a red border while the field is in `errors`
  (a style closure that sets the border colour from
  `theme.palette().danger`).
- Headings `Transform`, `Cylinder`/`Cone` (Min, Max and the closed
  checkbox), `Material`, `Light`, as `text` rows; the three axes of a
  vector field in one `row!` with a short label (`position`, `rotation°`,
  `scale`, `color`, `position`, `intensity`); scalar fields one per row
  with their name.
- The JSON stays below the boxes, as now.
- Title adds ` •` while `unsaved`. The status line holds the status
  text, and at its right `not in viewer — Cmd-L to load` while
  `unloaded`.
- `shortcut` adds `"l" => Some(Message::Load)`.

**README:** the keys table gains Cmd-L; a paragraph on editing: which
items have boxes, red boxes, Cmd-L, the `•`, and that quitting discards
unsaved edits without asking.

- [ ] **Step 1: Implement the state, messages and view** as above.

- [ ] **Step 2: Build and lint**

Run: `cargo build --release && cargo test --lib && cargo clippy --all-targets`
Expected: builds; all tests pass; clippy clean.

- [ ] **Step 3: Check by eye** (with the user), running
  `cargo run --release -- ../rtc-metal/scenes/chapter10.json` (chapter15 is two teapots and a plane):
  1. Select a sphere; change its position and color; the title shows
     ` •` and the status line `not in viewer — Cmd-L to load`. Cmd-L:
     the viewer restarts and shows the change.
  2. Type `-`, `2.` and `abc` in a box: the text stays as typed; only
     `abc` is red.
  3. With a box focused, Cmd-L, Cmd-S and Cmd-R work. If iced's text
     box swallows them, switch the subscription to
     `iced::event::listen_with` and match key presses whatever their
     status, then check again.
  4. Cmd-S clears the ` •`.
  5. In a scene with a cylinder (`../rtc-metal/scenes/chapter13.json`),
     toggle closed and edit max; Cmd-L shows it.
  6. A shape with a shear shows material boxes only.

- [ ] **Step 4: Commit** after the user's approval:

`ENHANCEMENT: rtc-edit edits shapes and point lights in fixed fields (milestone 18)`

### Task 4: arrow keys and h/j/k/l for the tree and the boxes

**Files:**
- Modify: `RUST/rtc-edit/src/fields.rs`, `RUST/rtc-edit/src/tree.rs`,
  `RUST/rtc-edit/src/main.rs`, `RUST/rtc-edit/README.md`.

**Interfaces:**
- Consumes: `fields::fields`, `tree::{rows, Row, Node}`, Task 3's
  `left` (the shaded box) and `box_id`.
- Produces:
  - `fields::next_box(order: &[Field], from: Field, right: isize, down: isize) -> Field`
  - `tree::next_row(rows: &[Row], from: Option<&Node>, by: isize) -> Option<Node>`

- [ ] **Step 1: Write the failing tests.**

  In `fields.rs`, on the cylinder's order from `fields`:
  - `next_box_moves_along_a_vector_and_stops_at_its_end`:
    `Position(0)` right 1 → `Position(1)`; `Position(2)` right 1 →
    `Position(2)`; `Position(0)` right -1 → `Position(0)`.
  - `next_box_keeps_the_axis_going_down_and_clamps_it_on_a_single_box`:
    `Position(2)` down 1 → `Rotation(2)`; `Scale(2)` down 1 → `Min`;
    `Max` down 1 → `Color(0)`; `Ambient` up 1 → `Color(0)`.
  - `next_box_stops_at_the_first_and_last_row`: `Position(1)` up 1 →
    `Position(1)`; `RefractiveIndex` down 1 → `RefractiveIndex`.

  In `tree.rs`, on a scene with one light and a group holding two
  spheres:
  - `next_row_moves_to_the_next_and_previous_row`.
  - `next_row_skips_the_rows_a_fold_hides`: with the group folded, the
    row after the group is the next top-level shape (or `None` past the
    end, since rows stop).
  - `next_row_stops_at_the_first_and_last_row`: `Camera` by -1 →
    `Camera`; the last row by 1 → the last row.
  - `next_row_with_nothing_selected_picks_the_first_row`: `None` by 1
    and by -1 → `Camera`.

- [ ] **Step 2: Run them to verify they fail**

Run: `cargo test --lib next_`
Expected: FAIL to compile, `next_box` / `next_row` not found.

- [ ] **Step 3: Implement `next_box` and `next_row`.** A box's row is
  its variant (`std::mem::discriminant`), its column the axis (0 for a
  scalar). Going up or down, take the neighbouring row and the box at
  `min(axis, row width - 1)`. A `from` not in the list returns `from`.

- [ ] **Step 4: Run the tests**

Run: `cargo test --lib`
Expected: all pass.

- [ ] **Step 5: Wire the keys in `main.rs`.**
  - State: `editing: Option<Field>`, the box being edited, beside
    `left`.
  - `shortcut` emits, only for key presses with no modifiers:
    `Move(dx, dy)` for arrows and h/j/k/l when the event was
    **ignored** (a focused box takes ←/→ and letters); `Enter` when
    ignored; `Escape(captured: bool)` either way; `Letter(char)` for
    i, c and x when ignored.
  - `Move`: editing → nothing; shaded → `left = next_box(...)`; tree →
    ↑/↓ `select(next_row(...))`, ← fold / → unfold the selected
    foldable row.
  - `Enter`: shaded → focus it, `editing = left.take()`; tree with
    boxes → `left` = first box.
  - `Letter` on a shaded box (else nothing): `i` focuses it and
    `move_cursor_to_front`; `c` sets its text to `""` through `edit`,
    then focuses it; `x` calls `set(desc, node, field, "")` and only on
    Ok stores the empty text and marks unsaved and unloaded, leaving
    the box shaded.
  - `Escape(true)` (the box took it and lost focus): `left =
    editing.take()`. `Escape(false)`: `left = None`.
  - `Submit(f)`: `left = Some(f)`, `editing = None`, unfocus.
  - `Edit` / `Paste(f)`: `editing = Some(f)`.
  - `tab(step)`: always the editor's own order — from `editing`, else
    `left`, else first (Tab) or last (Shift-Tab) box; wrap; focus it and
    set `editing`, clear `left`. `focus_next` / `focus_previous` are no
    longer used.
  - `Clicked`: clear `left` and `editing`, then `Task::batch` of
    `operation::is_focused(box_id(f))` for each box, mapped to
    `Focused(f, bool)`; `true` sets `editing`.
  - `select` clears `editing` too.

- [ ] **Step 6: README.** The Editing section gets a short keys list:
  tree (↑↓/jk, ←→/hl, Enter), shaded box (arrows/hjkl, Enter, i, c,
  x, Esc),
  editing (Enter/Esc), Tab.

- [ ] **Step 7: Build and lint**

Run: `cargo build --release && cargo test --lib && cargo clippy --all-targets`
Expected: builds; all tests pass; clippy clean.

- [ ] **Step 8: Check by eye** (with the user), on
  `../rtc-metal/scenes/chapter10.json` and a scene with a group:
  1. j/k and ↑/↓ walk the tree; h folds, l unfolds a group.
  2. Enter shades the first box; h/j/k/l and arrows move the shading
     and stop at the edges.
  3. Enter edits; ↑/↓ change nothing; Enter shades it again; Esc goes
     back to the tree, and j moves the tree.
  4. Click a box, press ↓: nothing moves. Tab moves on from the clicked
     box.
  5. Cmd-L with the shading on still loads.
  6. On a shaded box: i edits with the cursor at the start; c empties
     it for typing; x empties ambient (default applies after Cmd-L)
     and leaves a position unchanged.

- [ ] **Step 9: Commit** after the user's approval:

`ENHANCEMENT: rtc-edit moves through the tree and boxes with arrows and h/j/k/l, and i, c, x (milestone 18)`
