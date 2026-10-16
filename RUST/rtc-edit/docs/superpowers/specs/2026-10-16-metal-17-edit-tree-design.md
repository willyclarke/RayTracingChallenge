# Milestone 17: the rtc-edit scene editor, read-only

Branch `metal/17-edit-tree`, merged into develop with `--no-ff` when
done. This is the first of three milestones for the scene editor:

- 17: the program, a read-only tree, Save and Save As, and the viewer;
- 18: editing with fixed fields;
- 19: add and delete, New and Open.

The editor has no image of its own. It writes a preview file, and
`rtc-metal --window --watch` (milestone 16) shows it.

## Goal

`rtc-edit scene.json` opens a window with the scene as a tree, and starts
the Metal viewer on it. Save writes the scene back. Save As writes it to
another folder, and its `.obj` and `.ppm` files still load from there.

Success:
- Opening `rtc-metal/scenes/chapter15.json` shows its tree, and the viewer
  renders the teapots.
- Save As into another folder gives a file that loads its teapots.
- Closing the viewer and pressing Cmd-R opens it again.
- Quitting the editor closes the viewer and removes the preview file.
- Every scene file keeps its shape comments through a save.

## Out of scope

- Editing any value (18).
- Adding or deleting items, New and Open (19). The scene comes from the
  command line.
- An unsaved-changes prompt: nothing can change yet.
- A native menu bar. iced has none; the commands are keyboard shortcuts.
- Keeping the exact text of a file. Milestone 16 already writes `1` as
  `1.0` and orders keys by the structs.

## Decisions

- **A new crate, `RUST/rtc-edit`**, next to rtc-rust and rtc-metal, with
  its own `target/` and no Cargo workspace. Dependencies:
  - `rtc-rust` (path), for the scene description types and `save`;
  - `iced = "0.14"`;
  - `rfd`, for the native Save panel, since iced has no file dialogs.

  The editor never uses Metal. It talks to the viewer only through the
  preview file.
- **The viewer starts when the scene opens**, and Cmd-R starts it again if
  it was closed. Cmd-R does nothing while the viewer is still running.
- **Finding `rtc-metal`:** `$RTC_METAL` if it is set. Otherwise
  `../rtc-metal/target/release/rtc-metal` relative to this crate's
  folder (`env!("CARGO_MANIFEST_DIR")`).
- **The preview file** is `<temp dir>/rtc-edit-<pid>/<scene file name>`,
  so the viewer's title shows the scene's own name. It holds the scene
  with every asset path made absolute, written with `save`. Quitting the
  editor removes the folder.
- **Save As keeps asset paths relative.** A relative path is rewritten to
  reach the same file from the new folder (`models/teapot.obj` saved one
  folder down becomes `../models/teapot.obj`). An absolute path is left
  alone. Folders are compared after `canonicalize`, so `/tmp` and
  `/private/tmp` are the same folder. The relative path is computed by
  hand, so there is no extra dependency.
- **The inspector shows the selected item's JSON**, pretty-printed and
  read-only. A group's or CSG's children are left out, because they have
  their own rows. Milestone 18 replaces it with fields for the types you
  can edit.
- **Shapes keep their comments.** rtc-rust's `ShapeDescription` gains a
  `comment: Option<serde_json::Value>`, like `SceneDescription`'s.
  Without it, Save drops every shape comment: 17 scene files have them.
  `ShapeDescription` uses `#[serde(flatten)]`, so it cannot deny unknown
  keys and accepted the comment without storing it.

## rtc-rust

`src/scene.rs`: `ShapeDescription` gets

```rust
/// Free-form notes; ignored by the loader. The editor labels the shape
/// with it.
#[serde(default)]
#[serde(skip_serializing_if = "Option::is_none")]
pub comment: Option<serde_json::Value>,
```

`build` ignores it. The README's shape keys list it.

## rtc-edit

- **`src/paths.rs`**:
  - `for_each_asset(desc: &mut SceneDescription, f: impl FnMut(&mut String))`
    calls `f` on the `file` of every `obj` shape and every `uv_image`.
    It looks inside groups, CSG operands, nested, blended and perturbed
    patterns, texture maps and the six cube map faces.
  - `rebase(desc, from_dir, to_dir) -> Result<(), String>` rewrites
    relative paths for Save As.
  - `absolutize(desc, from_dir)` makes every relative path absolute, for
    the preview.
- **`src/tree.rs`**:
  - `Node` names one item: `Camera`, `Settings`, `Light(i)`, or
    `Shape(Vec<usize>)`. The shape path goes top-level index first, then
    the child index in a group (CSG: 0 left, 1 right).
  - `rows(desc, folded: &HashSet<Node>) -> Vec<Row>`, where
    `Row { node, depth, label, foldable }`. The rows are, in order:
    - `camera`;
    - `settings`, only when it is not the default;
    - `lights (n)` with one row per light;
    - `shapes (n)` with one row per shape, nested.

    The two sections fold, using `Node`s of their own. A shape's label is
    its type, then its comment when the comment is a string, e.g.
    `obj — smooth-shaded teapot`. Groups and CSG add `(n children)` and
    the operation respectively.
  - `json(desc, node) -> String` is the inspector's text.
- **`src/preview.rs`**: `Preview::new(scene_name)` makes the folder;
  `write(desc, scene_dir)` absolutizes a copy and saves it. Removing the
  folder happens on `Drop`.
- **`src/viewer.rs`**:
  - `find(env: Option<&str>, crate_dir: &Path) -> Result<PathBuf, String>`
    is pure, so it is tested without touching the environment.
  - `Viewer::start(exe, preview)` spawns
    `rtc-metal <preview> --window --watch`.
  - `running(&mut self)` uses `try_wait`.
  - Killing the child happens on `Drop`.
- **`src/main.rs`**: the iced application.
  - State: the scene path, the description, the tree's folded set, the
    selection, the preview, the viewer, and the status line.
  - Messages: Select, Toggle, Save, SaveAs, ShowViewer, and the window's
    close request.
  - A close request drops the viewer and the preview, then exits.

## The window

- **Layout:** the tree on the left, about 260 px wide and scrollable; the
  inspector on the right. Under them is a status line holding the last
  message, such as `saved chapter15.json`, `viewer started` or an error.
  The title is `rtc-edit — <file name>`.
- **The tree:** click a row to select it, and ▸/▾ to fold. The selected
  row is highlighted.
- **Keys:**
  - Cmd-S saves;
  - Cmd-Shift-S opens Save As;
  - Cmd-R shows the viewer;
  - Cmd-Q and closing the window quit.

## Error handling

- **The scene doesn't load at startup:** the message on stderr, exit
  code 1.
- **`rtc-metal` not found, or it fails to start:** the status line says
  so, naming the path tried. The editor keeps working.
- **Save or Save As fails:** the status line shows the error. `save` is
  atomic, so the old file is untouched.
- **Save As can't rebase a path** (its file doesn't exist, so its folder
  can't be canonicalized): it saves nothing, and the status line names
  the path.
- **The preview can't be written:** the status line says so, and no
  viewer is started.

## Testing

rtc-rust:
- `test_scene_14`: a shape's string comment survives a round trip, and
  `build` ignores it.
- `test_scene_10` still round-trips every scene file. It now also checks
  that the number of `"comment"` keys is the same after `to_json`.

rtc-edit (unit tests, no window):
- `paths`:
  - rebasing into a child, a parent and a sibling folder;
  - an absolute path is left alone;
  - `obj` paths inside a group and a CSG, and `uv_image` paths inside a
    texture map, a cube map and a perturbed pattern;
  - rebasing every `scenes/*.json` file of rtc-metal and rtc-rust into
    its own folder leaves it unchanged;
  - an `absolutize`d chapter15 saved in a temp folder loads and builds.
- `tree`:
  - rows, depths and labels for a scene with a group, a CSG, two lights
    and settings;
  - folding `shapes` hides all shape rows;
  - default settings give no `settings` row;
  - `json` leaves out a group's children.
- `viewer`: `find` prefers `$RTC_METAL`; without it it gives the path
  next to the crate; it is an error when that file doesn't exist.
- `preview`: `write` puts a loadable file in the folder; dropping the
  `Preview` removes the folder.

By eye:
- Open `rtc-metal/scenes/chapter15.json`: the tree shows camera, light
  and shapes with their comments, and the viewer starts.
- Fold and select; the inspector shows each item's JSON.
- Save As into another folder, then load that file with
  `rtc-metal <file> --window`: the teapots are there.
- Close the viewer, press Cmd-R: it comes back.
- Quit with Cmd-Q and by closing the window: the viewer closes and the
  `rtc-edit-<pid>` folder is gone.

## Risks

- **Quitting through the system menu might skip `Drop`.** If winit ends
  the process without unwinding, the viewer and preview survive the
  editor. The close-request handler covers the window's close button,
  and the by-eye check covers Cmd-Q. If Cmd-Q skips it, Cmd-Q becomes a
  key binding that goes through the same handler.
- **Unknown keys on shapes and patterns are still dropped on save.** They
  use `flatten`, so serde accepts any key there. `comment` was the only
  one in use; a typo like `"materail"` is lost silently, as it was
  ignored silently before.
- **iced's build time.** The first build of the crate takes minutes;
  rtc-rust and rtc-metal are not affected.
