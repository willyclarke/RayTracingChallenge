# rtc-edit

A scene editor for the JSON scenes of [rtc-rust](../rtc-rust/README.md), written with
[iced](https://iced.rs). It has no image of its own: it writes a preview copy of the scene,
and [rtc-metal](../rtc-metal/README.md)'s `--window` renders it live, reloading it when the
editor sends it a scene through a shared file. macOS only, because the viewer is.

The editor grows in milestones:

- **17 (done):** the scene as a foldable tree, each item's JSON in the inspector, Save and
  Save As, and the viewer.
- **18 (done):** editing spheres, cubes, planes, cylinders and cones (position, rotation, scale and a
  basic material), and point lights.
- **19 (done):** the viewer's image on quit, and the camera.
- **20 (done):** New and Open, adding and deleting items, comments.
- **21 (done):** the editor and the viewer side by side.
- **23 (done):** the viewer's view into the scene's camera.
- **24 (done):** undo and redo.
- **25 (done):** a native menu bar.
- **26 (done):** a command channel in the shared file; one viewer lives through New and Open.
- **27 (done):** GCD sources and callbacks instead of timers; an idle editor and viewer run no code.
- **28 (done):** a help window: the menus, keys, viewer and command line on one page.
- **29 (done):** About with a ray-traced icon and the build; Report an Issue….
- **30 (done):** a settings file: the scene and image folders and the editor's window.
- **31 (done):** an app bundle, RTC Edit with RTC Viewer inside.
- **32 (done):** File ▸ Open Recent, kept in the settings file.

Planned: a scene file type, opened from Finder; the viewer's samples and traversal;
importing a scene into the current one as a group (File ▸ Import Scene…).

## Dependencies

- [`rtc-rust`](../rtc-rust) (path dependency): the scene description types, loading and
  saving.
- `serde_json`: the inspector's JSON.
- `iced` 0.14: the window and widgets.
- `rfd`: the native Save and Open panels and the unsaved-changes alert, since iced has no
  file dialogs.
- `muda` 0.21: the native menu bar, since iced has none.
- `dispatch2`, `block2`: a GCD proc source that says when the viewer has exited.
- [`rtc-config`](../rtc-config) (path dependency): the settings file, shared with the
  viewer.

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

The editor looks for the viewer in `$RTC_METAL`, then inside the app bundle it runs from (see
below), and otherwise at `../rtc-metal/target/release/rtc-metal` next to this crate. Without it, the editor still
works and the status line says where it looked.

### Building the app

From `RUST/` (the alias in `RUST/.cargo/config.toml` works only there):

```bash
cargo xtask bundle            # builds both crates into "build/RTC Edit.app"
cargo xtask bundle --install  # …and replaces "~/Applications/RTC Edit.app"
```

`RTC Edit.app` holds the viewer as `Contents/Helpers/RTC Viewer.app`, so the Dock shows two
apps with their own names and icons, and the app finds its viewer wherever you put it.
Finder, the Dock and the menus say RTC Edit and RTC Viewer; the commands stay `rtc-edit` and
`rtc-metal`. Quit
the installed app before `--install`; a running copy keeps its old files until it quits.

The bundle is signed ad hoc (`codesign --sign -`): fine on the Mac that built it, with no Apple
account. Anyone else builds and signs their own the same way; a copy from another Mac would
need a Developer ID and notarisation to get past Gatekeeper.

Started from Finder, the Dock or Spotlight, rtc-edit starts with a new scene in `/`, which is
why the panels start in the [Settings](#settings) folders. Opening a scene from Finder (a
double-click, a drop on the Dock icon, `open -a rtc-edit x.json`) doesn't work yet: the app
opens with a new scene.

The editor opens at the top-left of the screen, and each viewer it starts (on opening and
Cmd-R) opens just right of it, without taking the focus, tops aligned, fitted to the rest of that screen
(`rtc-metal --at`). A new viewer goes beside the editor where it is now; one already open
stays where it is. With less than 200 points left right of the editor, the viewer opens
centred.

The menu bar holds every command, with its key; an item that can't act is greyed.

| menu | item | key | does |
|---|---|---|---|
| RTC Edit | About RTC Edit, Services, Hide RTC Edit, Hide Others, Show All, Quit RTC Edit | Cmd-H, Cmd-Opt-H, Cmd-Q | as in any Mac app; About shows the icon, the build (`Version 0.1.0 (<commit> <date>)`, ` dirty` for uncommitted changes) and the credits |
| File | New | Cmd-N | a new scene, `untitled.json`: the camera, a light and a floor |
| File | Open… | Cmd-O | open another scene (a `.json` file) |
| File | Open Recent ▸ | | the last 10 scenes opened, saved as or launched with, newest first; a name two of them share gets its folder; a moved or deleted scene says "not found" and leaves the list; Clear Menu empties it |
| File | Save | Cmd-S | save the scene; on an untitled scene it asks where, like Save As |
| File | Save As… | Cmd-Shift-S | save it elsewhere; relative `.obj` and `.ppm` paths are rewritten so they still reach the same files |
| Edit | Undo, Redo | Cmd-Z, Cmd-Shift-Z | undo, redo (see Undo); greyed with nothing to undo or redo |
| Edit | Delete | Cmd-Delete | delete the selected shape or light; greyed when it can't go, or while a box is shaded or edited |
| Scene | Add Sphere … Add Light | | add a sphere, cube, plane, cylinder, cone or light (see Adding and deleting) |
| Viewer | Load Edits | Cmd-L | load the edits into the viewer; greyed with no viewer |
| Viewer | Start Viewer | Cmd-R | start the viewer again after closing it; greyed while one runs |
| Viewer | Follow the Viewer | | tick or untick follow (see below); greyed with no viewer |
| Viewer | Keep Followed View | Cmd-K | keep the followed view (see below); greyed unless following |
| Help | RTC Edit Help | | open the help window beside the editor (centred if it doesn't fit), or bring it to the front; Esc closes it. Cmd-? opens the Help menu with macOS's search box, where this item is first |
| Help | Report an Issue… | | open a new GitHub issue with rtc-edit's and macOS's versions filled in |

The Help menu's Send RTC Edit Feedback isn't RTC Edit's: macOS adds it for users in Apple's beta
program, and it reports to Apple.

The icon, in About and in the Dock while rtc-edit runs, is rendered by rtc-metal from
`assets/about.json`; the app bundle's icon is the same scene at 1024 px, `assets/icon.json`.
After changing either scene, render it again from `../rtc-metal`:
`cargo run --release -- ../rtc-edit/assets/about.json --samples 1024 -o ../rtc-edit/assets/about.ppm`
(and the same for `icon.json` / `icon.ppm`).

The Edit menu has no Cut, Copy, Paste or Select All: Cmd-X, Cmd-C, Cmd-V and Cmd-A go to the
box being edited.

`rtc-edit` with no argument starts with a new scene. Cmd-N and Cmd-O replace the scene; with
unsaved edits (the ` •`) an alert first asks `Save changes to <name>?` with Save, Don't Save and
Cancel. Cancel, or a Save that fails or whose Save As is cancelled, keeps the scene and its
edits. A file that doesn't load leaves the scene as it is, and the status line says why. On a
switch the running viewer shows the new scene from its own camera, and the editor keeps the
focus; with no viewer running, none starts (Cmd-R starts one). To keep the old scene's image,
use the viewer's File ▸ Save Image… (Cmd-S) first.

## The viewer's view

With the camera selected, a `viewer` column beside `from`, `to` and `up` shows the viewer's
camera, greyed, as you orbit and zoom in it (`no viewer` with none running, `waiting` until it
has written one). The two share it through a small file, `$TMPDIR/rtc-edit-<pid>.share`
([rtc-share](../rtc-share/README.md)), made by the first viewer and removed when the editor
quits; the editor sleeps on it until the viewer writes, and sends the scenes to load through it
(Load Edits, New, Open); a viewer built from an older rtc-share prints
`rebuild rtc-metal` and the column stays at `waiting`.

- Ticking **follow the viewer**, above the `from` row, copies the viewer's `from` into the scene as you orbit, rounded to 4
  decimals, and with focal blur also its focal distance (zooming moves the focus with the
  camera and keeps the lens, so the f-number changes as the viewer's does). The `from` boxes are read-only meanwhile, and each change is an
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

Help ▸ RTC Edit Help (Cmd-?, then Enter on the item) shows these keys, the menus, the viewer's mouse and keys and the
command line in a window of its own; Esc closes it.

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
| tree | d d | delete the selected shape or light |

## Adding and deleting

The button bar above the tree and the Scene menu add a sphere, cube, plane, cylinder, cone or
light, and the Delete button and Edit ▸ Delete delete the selected item. `a` then a letter does
the same from the keyboard (`a` shows the letters in the status line), and `d d` or Cmd-Delete
deletes. The keys work with no box shaded or edited.

A new shape goes at the end of a selected group, else right after the selected shape; inside a
CSG it goes after the CSG, whose two shapes are fixed. With nothing selected it goes at the end of
the scene. A new light goes at the end of the lights. The new item is selected, and the groups
that hide it unfold. Deleting selects the next row.

Two deletes are refused, with the reason in the status line: a shape directly inside a CSG (a
CSG needs both its shapes) and the last light (a scene needs a light). The camera, settings and
headings can't be deleted.

## Undo

Cmd-Z takes back the last change to the scene and Cmd-Shift-Z puts it back; the status line
names the step (`undid: delete sphere`). One step is:

- everything typed into one box, from entering it to leaving it (`c` and the typing after it
  included), `x` on a shaded box, a click on `closed`;
- adding or deleting an item;
- a whole run of following the viewer (Keep ends it), and unticking follow.

Folding, selection and the viewer's orbit aren't steps. The last 100 steps are kept; a new
change clears what Cmd-Shift-Z would redo, and New and Open clear it all (Save doesn't). Undoing
back to the scene as saved removes the ` •`. After an undo the viewer is behind: Cmd-L loads it.
Cmd-Z while following stops following first.

Clicking a box edits it. The tree doesn't scroll to follow the keys.

Saving writes the file in one step (a temporary file, then a rename), leaves out default
values, and keeps `comment` fields. It doesn't keep the file's exact text: `1` becomes
`1.0`, and keys follow a fixed order.

The preview is `rtc-edit-<pid>/<scene name>` in the system temp folder, with every asset path
absolute. Quitting the editor closes the viewer and removes the folder. Stopping it with
Ctrl-C or `kill` doesn't: the folder stays behind.

Quitting the editor or the viewer asks in the viewer whether to save the image; Save...
starts in the folder you started `rtc-edit` from.

## Settings

rtc-edit and the viewer remember a few things in
`~/Library/Application Support/rtc-edit/config.json`:

```json
{
  "scene_folder": "/Users/me/rtc/scenes",
  "image_folder": "/Users/me/rtc/renders",
  "editor_frame": { "x": 0, "y": 25, "width": 900, "height": 640 },
  "recent_scenes": ["/Users/me/rtc/scenes/pencil.json"]
}
```

- `scene_folder`: where Open… and Save As… start for an untitled scene, set by the last
  Open… or Save As…. A scene with a file starts in its own folder. At first, `~/rtc/scenes`.
- `image_folder`: where the viewer's Save Image… starts, set by the last image saved. At
  first, `~/rtc/renders`.
- `editor_frame`: the editor's window when it last closed. A frame that no longer fits its
  screen goes back to the top-left at 900×640.
- `recent_scenes`: File ▸ Open Recent, newest first, at most 10 absolute paths.

The `~/rtc` folders are made when a panel is about to open in them. Every key is optional,
a file that isn't a JSON object counts as empty, and keys the programs don't know are kept.
Delete the file to forget everything.

## Layout

| path | what it holds |
|---|---|
| `src/main.rs` | the iced application: state, messages, keys, view |
| `src/paths.rs` | asset paths in a scene: rebasing for Save As, absolute for the preview |
| `src/tree.rs` | the scene as tree rows, and an item's JSON |
| `src/trs.rs` | a transform list as scale, rotation and translation |
| `src/fields.rs` | the boxes an item has, and writing typed text into the scene |
| `src/preview.rs` | the preview file and its folder |
| `src/recent.rs` | the Open Recent list: reading and writing it, adding, removing, and its menu labels |
| `src/viewer.rs` | finding and running `rtc-metal` |
| `src/menu.rs` | the menu bar: its commands, their labels and keys, when each is greyed, and the `muda` menus |
| `src/help.rs` | the help window's page: sections, chips, view |
| `src/windows.rs` | which window a key or a close belongs to |
| `src/about.rs` | About's build, credits and icon, and the GitHub issue URL |
| `build.rs` | the commit rtc-edit is built from, for About |
| `assets/` | the icons' scenes and their renders: `about.*` for About, `icon.*` for the app |
| `../xtask` | `cargo xtask bundle`: the app bundle, its `Info.plist` files and icons |
| `../rtc-config` | the settings file both programs share |
| `docs/superpowers/` | design specs and implementation plans |
