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
- **19 (done):** the viewer's image on quit, and the camera.
- **20 (done):** New and Open, adding and deleting items, comments.
- **21 (done):** the editor and the viewer side by side.
- **23 (done):** the viewer's view into the scene's camera.

## Dependencies

- [`rtc-rust`](../rtc-rust) (path dependency): the scene description types, loading and
  saving.
- `serde_json`: the inspector's JSON.
- `iced` 0.14: the window and widgets.
- `rfd`: the native Save and Open panels and the unsaved-changes alert, since iced has no
  file dialogs.

The editor never uses Metal itself; it starts `rtc-metal` as a separate program. No Cargo
workspace: the crate has its own `target/`, and the first build takes a minute or two.

## Building and running

```bash
(cd ../rtc-metal && cargo build --release)    # the viewer
cargo run --release -- ../rtc-metal/scenes/chapter15.json
cargo run --release                           # no argument: a new scene
cargo test                                    # no window needed
cargo clippy --all-targets
```

The editor looks for the viewer in `$RTC_METAL`, and otherwise at
`../rtc-metal/target/release/rtc-metal` next to this crate. Without it, the editor still
works and the status line says where it looked.

The editor opens at the top-left of the screen, and each viewer it starts (on opening, New,
Open and Cmd-R) opens just right of it, tops aligned, fitted to the rest of that screen
(`rtc-metal --at`). A new viewer goes beside the editor where it is now; one already open
stays where it is. With less than 200 points left right of the editor, the viewer opens
centred.

| key | does |
|---|---|
| Cmd-N | a new scene, `untitled.json`: the camera, a light and a floor |
| Cmd-O | open another scene (a `.json` file) |
| Cmd-S | save the scene; on an untitled scene it asks where, like Save As |
| Cmd-Shift-S | save it elsewhere; relative `.obj` and `.ppm` paths are rewritten so they still reach the same files |
| Cmd-L | load the edits into the viewer |
| Cmd-R | start the viewer again after closing it (with one open, the status line says `viewer already running`) |
| Cmd-K | keep the followed view (see below) |

`rtc-edit` with no argument starts with a new scene. Cmd-N and Cmd-O replace the scene; with
unsaved edits (the ` •`) an alert first asks `Save changes to <name>?` with Save, Don't Save and
Cancel. Cancel, or a Save that fails or whose Save As is cancelled, keeps the scene and its
edits. A file that doesn't load leaves the scene as it is, and the status line says why. On a
switch the old viewer gets SIGTERM, so it offers to save its image, and a new viewer opens on the
new scene.

## The viewer's view

With the camera selected, a `viewer` column beside `from`, `to` and `up` shows the viewer's
camera, greyed, as you orbit and zoom in it (`no viewer` with none running, `waiting` until it
has written one). The two share it through a small file in the preview folder
([rtc-share](../rtc-share/README.md)); a viewer built from an older rtc-share prints
`rebuild rtc-metal` and the column stays at `waiting`.

- Ticking **follow the viewer**, above the `from` row, copies the viewer's `from` into the scene as you orbit, rounded to 4
  decimals, and with focal blur also its focal distance (zooming moves the focus with the
  camera; the f-number stays). The `from` boxes are read-only meanwhile, and each change is an
  unsaved edit.
- **Keep** (beside it while following, or Cmd-K) stops following and keeps the view; Cmd-S saves it.
- Unticking puts back the `from` and focal distance the camera had when you ticked.
- New, Open, or the viewer ending stop following without putting anything back; selecting
  another item doesn't stop it.
- Cmd-L isn't needed: the viewer already shows the view, and the next Cmd-L reloads it without
  a jump.

## Editing

Spheres, cubes, planes, cylinders and cones, also inside groups and CSG, get boxes for their
position, rotation (in degrees) and scale, and for their material: color, ambient, diffuse,
specular, shininess, reflective, transparency and refractive index. Cylinders and cones also
get min, max and closed; point lights get position and intensity; the camera gets from, to, up,
field of view (in degrees), width, height, antialias n and focal blur. Everything else is shown
read-only in the JSON below the boxes, and saved unchanged.

Every shape, groups, CSGs and OBJ files included, also gets a `comment` box above the other
boxes (alone, for those that have no others). The tree row shows it as `<type> — <comment>`. An
empty box removes the comment; a comment that isn't text gets no box and stays as it is.

- The camera's n only matters to rtc-rust renders, as rtc-metal's progressive rendering
  antialiases anyway. An empty n removes antialias.
- Focal blur is set by an f-number, as on a camera: a bigger one is sharper. Typing one creates
  focal blur, focused on the point the camera looks at, with 32 samples; an empty f-number removes
  it. The file stores the lens radius as `aperture`, `focal distance / (2 × f-number)`, shown
  below the boxes as `lens radius`. Changing the focal distance keeps the f-number. Emptying the f-number to retype it keeps the
  focal distance and samples until you select something else.
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
| tree | a, then s c p y o l | add a sphere, cube, plane, cylinder, cone or light; any other key cancels |
| tree | d d or Cmd-Delete | delete the selected shape or light |

## Adding and deleting

The button bar above the tree adds a sphere, cube, plane, cylinder, cone or light, and deletes
the selected item (the Delete button). `a` then a letter does the same from the keyboard (`a` shows the letters in
the status line), and `d d` or Cmd-Delete deletes. They work with no box shaded or edited.

A new shape goes at the end of a selected group, else right after the selected shape; inside a
CSG it goes after the CSG, whose two shapes are fixed. With nothing selected it goes at the end of
the scene. A new light goes at the end of the lights. The new item is selected, and the groups
that hide it unfold. Deleting selects the next row.

Two deletes are refused, with the reason in the status line: a shape directly inside a CSG (a
CSG needs both its shapes) and the last light (a scene needs a light). The camera, settings and
headings can't be deleted.

There is no undo: to go back, open the file again and choose Don't Save.

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
