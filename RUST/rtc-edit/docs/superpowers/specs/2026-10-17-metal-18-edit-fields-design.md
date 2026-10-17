# Milestone 18: rtc-edit edits with fixed fields

Branch `metal/18-edit-fields`, merged into develop with `--no-ff` when
done. The second of three milestones for the scene editor:

- 17: the program, a read-only tree, Save and Save As, and the viewer (done);
- 18: editing with fixed fields;
- 19: add and delete, New and Open, and the camera.

## Goal

The inspector gets text boxes for the values you change most when
composing a scene: a simple shape's transform and material, a cylinder's
or cone's bounds, and a point light. Edits build up in the editor;
Cmd-L loads them into the viewer, so several edits cost one restart.

Success:
- Select a sphere in `rtc-metal/scenes/chapter15.json`, change its
  position and color, press Cmd-L: the viewer shows the change.
- Typing `-`, `2.` or `abc` never fights you: the box keeps your text,
  and turns red while it doesn't parse.
- The title shows `•` after an edit; Cmd-S clears it.
- Opening and saving a scene without editing it changes no numbers.

## Out of scope

- The camera, and the milestone 16 minor where changing only the output
  size resets the orbit (19).
- Patterns and bump maps; area and spot lights; obj, group, CSG and
  torus settings. They stay visible in the read-only JSON.
- Adding and deleting items, New and Open (19).
- Undo. Reopen the file without saving.
- Asking before quitting. Cmd-Q in the app menu ends the process before
  the editor sees it (milestone 17), so unsaved edits are lost and the
  README says so.

## Decisions

- **Cmd-L loads the edits into the viewer.** Editing never writes the
  preview by itself. The status line says `not in viewer — Cmd-L to load`
  while the viewer is behind.
- **A box updates the scene as soon as its text parses.** The text you
  typed lives in the editor, not the scene, so partial numbers like `2.`
  stay as typed. Text that doesn't parse turns the box red and leaves the
  scene at its last good value.
- **Unsaved edits are marked, not guarded:** `•` in the title until
  Save or Save As. Quitting discards them without asking.
- **Editable items:**
  - sphere, cube, plane, cylinder and cone: transform and material;
  - cylinder and cone also: min, max, closed;
  - point light: position and intensity.

  This holds inside groups and CSG too. Everything else has no boxes.
- **Material boxes:** color (three numbers), ambient, diffuse, specular,
  shininess, reflective, transparency, refractive index. An empty box
  leaves the value unset, so the loader's default applies. Setting a
  value on a shape with no material creates a material holding only that
  value. A pattern or bump map is kept unchanged.
- **Transforms are position, rotation in degrees, and scale**, three
  numbers each (`src/trs.rs`):
  - A transform list fits when its operations come in the order scale,
    rotate x, rotate y, rotate z, translate, each at most once. Missing
    ones read as scale `1 1 1` and zero.
  - Shear, a matrix, a repeated operation or another order does not fit:
    the shape gets no transform boxes, and its material is still
    editable.
  - Writing rebuilds the list in that order and leaves out operations
    that do nothing (scale `1 1 1`, rotation `0`, translate `0 0 0`).
    Editing only the position of `[translate]` keeps it `[translate]`.
  - Rotations are stored in radians and shown in degrees, rounded to
    6 decimals. The stored radians are only replaced when that box is
    edited, so `0.785398` stays `0.785398` until you change it.
- **Validation:** numbers must be finite. A scale of 0 on any axis is
  rejected, because the matrix can't be inverted. Other values that
  parse but don't build (min above max, a refractive index of 0) are
  left to the viewer, whose overlay shows the loader's error and keeps
  the old image (milestone 16).

## rtc-edit

- **`src/trs.rs`**:
  - `Trs { scale: [f64; 3], rotate: [f64; 3], translate: [f64; 3] }`,
    rotations in radians; `Default` is scale `1 1 1`, the rest zero.
  - `Trs::from_ops(&[TransformOp]) -> Option<Trs>`: `None` when the list
    doesn't fit.
  - `Trs::to_ops(&self) -> Vec<TransformOp>`.
- **`src/fields.rs`**:
  - `Field` names one box: `Position(usize)`, `Rotation(usize)`,
    `Scale(usize)`, `Color(usize)`, `Ambient`, `Diffuse`, `Specular`,
    `Shininess`, `Reflective`, `Transparency`, `RefractiveIndex`, `Min`,
    `Max`, `LightPosition(usize)`, `LightIntensity(usize)`; the `usize`
    is the axis, 0 to 2.
  - `fields(desc, node) -> Vec<(Field, String)>`: the boxes for an item
    and their starting text, in display order. Empty for items you can't
    edit.
  - `set(desc, node, field, text) -> Result<(), String>`: parses `text`
    and writes it into the scene. On `Err` the scene is unchanged.
  - `closed(desc, node) -> Option<bool>` and
    `set_closed(desc, node, bool)`, for the cylinder's or cone's
    checkbox.
- **`src/main.rs`**:
  - New state: `texts: HashMap<Field, String>` for the selected item,
    rebuilt when the selection changes; `unsaved: bool`;
    `unloaded: bool`.
  - New messages: `Edit(Field, String)`, `Closed(bool)`, `Load`.
  - Selecting another item drops any red text; the scene keeps its last
    good value.
  - Cmd-L writes the preview and clears `unloaded`. If the preview can't
    be written, the status line says so.
  - Cmd-L, Cmd-S and Cmd-Shift-S with a red box write the scene with its
    last good values, and the status line warns that a field has an
    error.
- **The inspector:** the boxes come first, under *Transform*,
  *Material*, *Cylinder* (or *Cone*) and *Light*, with a vector's three
  boxes in one row. The item's JSON stays below them, read-only, and
  follows the edits.

## Error handling

- **Text that doesn't parse, or a zero scale:** the box turns red, and
  the scene is unchanged.
- **The preview can't be written on Cmd-L:** the status line shows the
  error, and `unloaded` stays set.
- **The edited scene doesn't build:** the viewer's overlay shows why.

## Testing

rtc-edit (unit tests, no window):
- `trs`:
  - fields to ops and back give the same `Trs`;
  - every fitting list survives `from_ops` then `to_ops` unchanged;
  - shear, a matrix, a repeated op and a wrong order give `None`;
  - every shape in every `scenes/*.json` file of rtc-metal and rtc-rust
    is read without panicking.
- `fields`:
  - a sphere, a cylinder and a point light list the expected fields
    with their starting text;
  - a shape whose transform doesn't fit lists material fields only;
  - obj, group, CSG, torus, area and spot lights, camera and settings
    list none;
  - `set` on position, rotation (degrees in, radians stored), color and
    min writes the scene, and the result still builds;
  - an empty material box unsets the value; setting one on a shape with
    no material creates it;
  - bad text and a zero scale return `Err` and leave the scene equal to
    before;
  - a shape inside a group can be edited.

By eye:
- Edit a sphere's position and color in chapter15, press Cmd-L: the
  viewer shows it.
- Type `-`, `2.` and `abc`: the box keeps the text, red for `abc`.
- The `•` appears after an edit, and Cmd-S clears it.
- Cmd-S, Cmd-L and Cmd-R work while a box has focus.
- Toggle a cylinder's closed box; the viewer shows the caps after Cmd-L.

## Risks

- **A focused text box may take the Cmd shortcuts.** If iced's
  `text_input` consumes them, the shortcuts listen to all key events
  rather than only ignored ones. Checked by eye.
- **Many boxes.** A cylinder has about 25; the inspector scrolls.
- **Rebuilding the transform list changes its text** on the first edit
  of a shape, for example dropping a `scale [1, 1, 1]`. The render is
  the same.
