# rtc-edit

A scene editor for the JSON scenes of [rtc-rust](../rtc-rust/README.md), written with
[iced](https://iced.rs). It has no image of its own: it writes a preview copy of the scene,
and [rtc-metal](../rtc-metal/README.md)'s `--window --watch` renders it live and reloads it
whenever it changes. macOS only, because the viewer is.

The editor grows in milestones:

- **17 (done):** the scene as a foldable tree, each item's JSON in the inspector, Save and
  Save As, and the viewer.
- **18 (done):** editing spheres, cubes, planes, cylinders and cones (position, rotation, scale and a
  basic material), and point lights.
- **19:** adding and deleting items, New and Open, and the camera.

## Dependencies

- [`rtc-rust`](../rtc-rust) (path dependency): the scene description types, loading and
  saving.
- `serde_json`: the inspector's JSON.
- `iced` 0.14: the window and widgets.
- `rfd`: the native Save panel, since iced has no file dialogs.

The editor never uses Metal itself; it starts `rtc-metal` as a separate program. No Cargo
workspace: the crate has its own `target/`, and the first build takes a minute or two.

## Building and running

```bash
(cd ../rtc-metal && cargo build --release)    # the viewer
cargo run --release -- ../rtc-metal/scenes/chapter15.json
cargo test                                    # no window needed
cargo clippy --all-targets
```

The editor looks for the viewer in `$RTC_METAL`, and otherwise at
`../rtc-metal/target/release/rtc-metal` next to this crate. Without it, the editor still
works and the status line says where it looked.

| key | does |
|---|---|
| Cmd-S | save the scene |
| Cmd-Shift-S | save it elsewhere; relative `.obj` and `.ppm` paths are rewritten so they still reach the same files |
| Cmd-L | load the edits into the viewer |
| Cmd-R | start the viewer again after closing it |

## Editing

Spheres, cubes, planes, cylinders and cones, also inside groups and CSG, get boxes for their
position, rotation (in degrees) and scale, and for their material: color, ambient, diffuse,
specular, shininess, reflective, transparency and refractive index. Cylinders and cones also
get min, max and closed; point lights get position and intensity; the camera gets from, to, up,
field of view (in degrees), width, height, antialias n and focal blur. Everything else is shown
read-only in the JSON below the boxes, and saved unchanged.

- The camera's n only matters to rtc-rust renders, as rtc-metal's progressive rendering
  antialiases anyway. An empty n removes antialias.
- Focal blur is set by an f-number, as on a camera: a bigger one is sharper. Typing one creates
  focal blur, focused on the point the camera looks at, with 32 samples; an empty f-number removes
  it. The file stores the lens radius as `aperture`, `focal distance / (2 × f-number)`, shown
  below the boxes as `lens radius`. Changing the focal distance keeps the f-number.
- A red box keeps its reason: the status line shows it while the box is shaded or being edited.
- A box updates the scene as you type. The shortcuts work while a box has focus.
- An empty material, min or max box leaves the value unset, so its default applies.
- A box whose text isn't a number turns red, and the scene keeps the last value that was.
- A transform the boxes can't show (shear, a matrix, or steps out of the order scale, rotate
  x, y, z, translate) has no boxes; the material can still be edited.
- Edits build up in the editor: Cmd-L loads them into the viewer, which then starts refining
  again, so several edits cost one restart. The status line says when the viewer is behind.
- The title shows ` •` while there are unsaved edits. Quitting doesn't ask: unsaved edits are
  lost.

The keyboard works in three places: the tree, a shaded box, and a box being edited. The
letters only act with no modifier held, so Cmd-L stays Load.

| where | key | does |
|---|---|---|
| tree | ↑/↓ or k/j | select the previous or next row |
| tree | ← or h, → or l | fold or unfold the selected item |
| tree | Enter | shade the item's first box |
| shaded box | arrows or h/j/k/l | move the shading; it stops at the edges |
| shaded box | Enter, i | edit it, with the cursor at the end or the start |
| shaded box | c | empty it and edit it |
| shaded box | x | empty it, if the box may be empty |
| shaded box | u | put the scene's value back, clearing a red box |
| shaded box | Esc | back to the tree |
| editing | Enter or Esc | leave the box, shaded; Esc on a red box also puts the scene's value back |
| anywhere | Tab, Shift-Tab | edit the next or previous box |

Clicking a box edits it. The tree doesn't scroll to follow the keys.

Saving writes the file in one step (a temporary file, then a rename), leaves out default
values, and keeps `comment` fields. It doesn't keep the file's exact text: `1` becomes
`1.0`, and keys follow a fixed order.

The preview is `rtc-edit-<pid>/<scene name>` in the system temp folder, with every asset path
absolute. Quitting the editor closes the viewer and removes the folder. Stopping it with
Ctrl-C or `kill` doesn't: the folder stays behind.

Quitting the editor or the viewer asks in the viewer whether to save the image; Save...
starts in the folder you started `rtc-edit` from.

## Layout

| path | what it holds |
|---|---|
| `src/main.rs` | the iced application: state, messages, keys, view |
| `src/paths.rs` | asset paths in a scene: rebasing for Save As, absolute for the preview |
| `src/tree.rs` | the scene as tree rows, and an item's JSON |
| `src/trs.rs` | a transform list as scale, rotation and translation |
| `src/fields.rs` | the boxes an item has, and writing typed text into the scene |
| `src/preview.rs` | the preview file and its folder |
| `src/viewer.rs` | finding and running `rtc-metal` |
| `docs/superpowers/` | design specs and implementation plans |
