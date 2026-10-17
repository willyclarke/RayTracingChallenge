# rtc-edit

A scene editor for the JSON scenes of [rtc-rust](../rtc-rust/README.md), written with
[iced](https://iced.rs). It has no image of its own: it writes a preview copy of the scene,
and [rtc-metal](../rtc-metal/README.md)'s `--window --watch` renders it live and reloads it
whenever it changes. macOS only, because the viewer is.

The editor grows in milestones:

- **17 (done):** the scene as a foldable tree, each item's JSON in the inspector, Save and
  Save As, and the viewer.
- **18:** editing spheres, cubes, planes, cylinders and cones (position, rotation, scale and a
  basic material), and point lights.
- **19:** adding and deleting items, New and Open.

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
| Cmd-R | start the viewer again after closing it |

Saving writes the file in one step (a temporary file, then a rename), leaves out default
values, and keeps `comment` fields. It doesn't keep the file's exact text: `1` becomes
`1.0`, and keys follow a fixed order.

The preview is `rtc-edit-<pid>/<scene name>` in the system temp folder, with every asset path
absolute. Quitting the editor closes the viewer and removes the folder. Stopping it with
Ctrl-C or `kill` doesn't: the folder stays behind.

Quitting the viewer itself writes its image, as `rtc-metal --window` always does, into the
folder you started `rtc-edit` from.

## Layout

| path | what it holds |
|---|---|
| `src/main.rs` | the iced application: state, messages, keys, view |
| `src/paths.rs` | asset paths in a scene: rebasing for Save As, absolute for the preview |
| `src/tree.rs` | the scene as tree rows, and an item's JSON |
| `src/preview.rs` | the preview file and its folder |
| `src/viewer.rs` | finding and running `rtc-metal` |
| `docs/superpowers/` | design specs and implementation plans |
