# Milestone 20: New and Open, adding and deleting items

Branch `metal/20-new-open-add-delete`, merged into develop with `--no-ff`
when done. The scene editor's milestones:

- 17: the program, a read-only tree, Save and Save As, and the viewer (done);
- 18: editing with fixed fields (done);
- 19: the viewer's image on quit, and the camera (done);
- 20: New and Open, adding and deleting items.

## Goal

The editor can start a scene from nothing, open another scene without
restarting, and add and delete shapes and lights, so a scene can be
composed in the editor from start to finish.

Success:
- `rtc-edit` with no argument shows a new scene, and the viewer renders
  its floor.
- `a s` adds a sphere; Cmd-L shows it in the viewer. `d d` deletes it.
- Cmd-O on another scene, with unsaved edits, asks to save them first;
  the old viewer offers its image and a new viewer opens on the new
  scene.
- Cmd-S on an untitled scene asks where to save it.

## Out of scope

- Undo. Open the same file again and choose Don't Save.
- Asking before the editor quits with unsaved edits. Cmd-Q in the app
  menu ends the process before the editor sees it (milestone 17).
- Moving items, copying, or adding groups, CSGs, OBJ files, tori, area
  and spot lights. They stay editable in the JSON file by hand.
- Remembering recent files or the last folder.

## Decisions

### New and Open

- **`rtc-edit` with no argument** starts with a new scene;
  `rtc-edit scene.json` works as before.
- **Cmd-N** replaces the scene with a new, untitled one:
  - camera: from (0, 1.5, −5), to (0, 1, 0), up (0, 1, 0), 800×600,
    field of view 60° (π/3);
  - lights: one white point light, position (−10, 10, −10), intensity
    (1, 1, 1), stored as `light`;
  - shapes: one plane, no transform, no material.
  - The title shows `untitled.json`. Cmd-S on an untitled scene acts as
    Save As.
- **Cmd-O** opens a file dialog in the current scene's folder (the
  working directory for an untitled scene), filtered to `.json`. A file
  that doesn't parse leaves the current scene, and the status line says
  why.
- **Unsaved edits** (the `•`): before New or Open, an alert asks "Save
  changes to `<name>`?" with Save, Don't Save and Cancel.
  - Save saves, through Save As for an untitled scene; cancelling that
    dialog, or a save that fails, cancels the New or Open.
  - Don't Save goes ahead. Cancel does nothing.
- **Switching scenes** (New, Open):
  - The old viewer gets SIGTERM, so it offers to save its image, named
    after the old scene (milestone 19). The old preview folder is
    removed.
  - A new preview is made, named after the new scene (`untitled.json`
    for New), and a new viewer starts on it.
  - The editor drops its selection, folds, boxes, red text, remembered
    focal blur, shading and editing, and clears `•` and `unloaded`.

### Adding

- **A button bar above the tree**: `+ sphere`, `+ cube`, `+ plane`,
  `+ cylinder`, `+ cone`, `+ light`, `Delete`. Delete is disabled while
  the selection can't be deleted.
- **Keys**, while the tree has the keys (no box shaded or edited):
  `a` then `s` sphere, `c` cube, `p` plane, `y` cylinder, `o` cone,
  `l` light. After `a`, the status line lists the letters; any other key
  or Esc cancels.
- **Where a shape goes:**
  - a shape selected: right after it, at its level (inside its group);
  - a group selected: at the end of that group;
  - a shape inside a CSG selected: right after the outermost CSG that
    holds it, since a CSG holds exactly two shapes;
  - Shapes, the camera, Settings, Lights, a light or nothing selected:
    at the end of the top-level shapes.
- **New shapes**, each with material `color` (0.8, 0.8, 0.8) and
  nothing else:
  - sphere: `translate [0, 1, 0]`;
  - cube: `translate [0, 1, 0]`;
  - plane: no transform;
  - cylinder and cone: `min 0`, `max 1`, closed, no transform.
- **A new light** is a white point light at (−10, 10, −10), added at the
  end of the lights. A scene holding a single `light` keeps it as the
  first entry of a new `lights` list.
- **After adding**, the new item is selected (its row shaded in the tree,
  no box shaded), the scene is marked `•` and `unloaded`, and the status
  line asks for Cmd-L as after any edit.

### Deleting

- **`Delete` or `d d`** deletes the selected shape (a group or CSG with
  everything in it) or light.
- **Refused**, with the reason in the status line and the scene
  unchanged:
  - a shape directly inside a CSG: "a CSG needs both its shapes";
  - the last light: "a scene needs a light" (the loader rejects a scene
    without one);
  - the camera, Settings, Lights and Shapes rows.
- **No confirmation.** There is no undo; open the file again and choose
  Don't Save.
- **The selection** moves to the next row at the same level, else the
  previous one, else the parent (the group, or the Shapes or Lights
  row).
- **Folds follow their items:** a folded group stays folded when a shape
  is added or deleted before it.
- Deleting marks the scene `•` and `unloaded`.

## rtc-edit

- **`src/edit.rs`** (new), changes to the scene with no UI:
  - `pub enum Item { Sphere, Cube, Plane, Cylinder, Cone, Light }`;
  - `pub fn new_scene() -> SceneDescription`;
  - `pub fn add(desc: &mut SceneDescription, selected: Option<&Node>,
    item: Item) -> Node`: inserts the item, returns its node;
  - `pub fn delete(desc: &mut SceneDescription, node: &Node)
    -> Result<Node, String>`: removes the item, returns the node to
    select; on `Err` the scene is unchanged;
  - `pub fn can_delete(desc: &SceneDescription, node: &Node) -> bool`,
    for the Delete button;
  - `pub fn shift(folded: &mut HashSet<Node>, at: &[usize], by: isize)`:
    moves the folds of the shapes at and after `at` (same parent) by
    `by` places, and drops the folds of a deleted shape and its
    children.
- **`src/main.rs`**:
  - `path: Option<PathBuf>`, `None` for an untitled scene.
  - `fn switch(&mut self, desc: SceneDescription, path: Option<PathBuf>)`,
    used at start-up, by New and by Open: replaces the preview and the
    viewer and resets the editor's state.
  - `fn ask_unsaved(&mut self) -> bool`: true when New or Open may go
    ahead; the alert is an `rfd::MessageDialog` with
    `MessageButtons::YesNoCancelCustom("Save", "Don't Save", "Cancel")`.
  - Messages `New`, `Open`, `Add(Item)`, `Delete`; shortcuts Cmd-N and
    Cmd-O; the pending key after `a` or the first `d`, cleared by any
    other key.
  - The button bar, a `row` of buttons above the tree.
- **`README.md`**: New, Open, the alert, the button bar, the keys, and
  how to undo.

## Error handling

- **Open of a file that doesn't parse:** status line, the scene unchanged.
- **The new preview can't be written:** as at start-up, the status line
  says so and Cmd-R reports it.
- **Save in the unsaved alert fails or is cancelled:** New or Open is
  cancelled, the edits kept.
- **A refused delete:** the reason in the status line.
- **An added shape the viewer can't build:** the viewer's overlay
  (milestone 16); the defaults above all build.

## Testing

rtc-edit (unit tests, `src/edit.rs`):
- `new_scene` builds with rtc-rust, and its JSON reads back equal;
- `add` places a shape after a top-level shape, after a shape inside a
  group, at the end of a selected group, after the CSG when a CSG's
  child is selected, and at the end of the top level with the camera,
  a light or nothing selected; it returns the new node;
- each added shape kind has its defaults, and the scene still builds;
- adding a light to a scene with `light` gives a two-item `lights`;
- `delete` removes a shape inside a group and returns the next row,
  the previous row, or the parent; it refuses a CSG child, the last
  light and the heading rows, leaving the scene equal to before;
  `can_delete` agrees;
- `shift` keeps a folded group folded after an add or a delete before
  it, and drops the folds of a deleted group.

By eye:
- `rtc-edit` with no argument: the new scene, the floor in the viewer.
- `a s`, Cmd-L: a sphere. A cylinder added inside a group in
  `chapter14.json` renders. `+ light` brightens the scene.
- `d d` deletes; a CSG child and the last light are refused with the
  reason.
- Cmd-O with edits: the alert; Don't Save switches, the old viewer
  offers its image, a new viewer opens. Cancel keeps everything.
- Cmd-S on an untitled scene opens Save As.

## Risks

- **An `rfd` alert while iced has focus** may leave a key or Cmd held in
  the editor's state. Checked by eye.
- **Two viewers on screen** during a switch: the old one waits for its
  image answer while the new one starts.

## Amendment: a comment box (Task 5, 2026-10-19)

The user found added shapes hard to tell apart in the tree: rows from a
file show `<type> — <comment>`, an added shape only its type, and the
editor had no way to set `comment`.

- **Every shape gets a comment box**, groups, CSGs, OBJ files and tori
  included; for those it is their only box. It sits above the first
  heading, labelled `comment`, with no heading of its own, about as wide
  as a row of three number boxes.
- **The text is stored as typed** as the shape's `comment` string; a box
  that is empty or only spaces removes the key. The tree row shows the
  new label at once, as other edits show at once.
- **A comment that isn't a string** (an array, say) gets no box, so it is
  never overwritten; it stays editable in the JSON file by hand.
- Lights and the camera get no comment box; their rows don't show one.
- The arrows, h/j/k/l, Tab, Enter, i/c/x/u, Esc and the red-box rules
  work on it as on any box; a comment is never red, since any text is
  valid.
