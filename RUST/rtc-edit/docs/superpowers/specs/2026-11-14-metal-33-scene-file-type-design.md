# Milestone 33: a scene file type, opened from Finder

Branch `metal/33-scene-file-type`, merged into develop with `--no-ff` when
done. Milestone 32 added Open Recent; after 33 come the viewer's settings,
+ csg / + group, and importing scenes (see Out of scope).

## Goal

Scenes are plain `.json` files, so Finder hands them to whatever app owns
JSON, and double-clicking one never reaches RTC Edit. Scenes should have
their own extension that RTC Edit owns, so a double-click, Open With or a
drop on the Dock icon opens the scene in the editor.

Success:
- Double-clicking a `.rtc` file opens it in RTC Edit, whether the editor
  is running or not, with the editor in front and the viewer behind.
- `.json` scenes still open everywhere they did, and Open With ▸ RTC Edit
  works for them without RTC Edit taking over `.json` double-clicks.
- New scenes and Save As default to `.rtc`.
- Opening from Finder behaves like Open…: the "Save changes?" alert,
  errors in the status line, the scene folder remembered, Open Recent.

## Out of scope

- **Merging several scenes into one** (raised while designing this
  milestone): when several files are dropped or opened at once, combine
  them into one scene, e.g. the largest canvas and every file's objects.
  It needs its own decisions: which camera wins, what to do with every
  file's lights (several sets overexpose), how to keep scenes built around
  the origin from overlapping, and name clashes. Proposed for later as
  **File ▸ Import Scene…**, which adds a file's objects to the current
  scene as one group you can move, built on + group. Here, only the first
  file opens.
- A document icon of its own for `.rtc` (Finder shows RTC Edit's icon).
- Renaming the repo's example scenes from `.json` to `.rtc`.
- The viewer (`rtc-metal`) registering for any file type; it loads by
  path, and either extension works.
- AppKit's own recent-documents list (`NSDocumentController`).
- Linux and Windows; rtc-edit is macOS-only (muda, AppKit).

## Decisions

- Extension `rtc`, type `no.sb21.rtc.scene`, described as "RTC Scene",
  conforming to `public.json`.
- RTC Edit is the **Owner** of `.rtc` and an **Alternate** for
  `public.json`.
- New scenes are `untitled.rtc`; Save As suggests `<stem>.rtc`; Open… and
  Save As accept both `rtc` and `json`. A plain Save writes back to the
  scene's own file, whatever its extension.
- Finder opens arrive through `application:openURLs:`, added with
  `class_addMethod` to winit's application delegate class in `boot()`.
  The spike (see below) ruled out the other ways.
- Several files at once: the first opens, the status line says how many
  were not opened.
- Window order after a Finder open: the editor in front, the viewer
  behind it.

### The spike

Done before this spec, code thrown away:
1. Replacing winit's delegate with `setDelegate:` crashes: winit's
   `app::send_event` looks up its own delegate and panics.
2. Installing an Apple Event handler for `odoc` in `boot()` never fires:
   AppKit installs its own afterwards, and macOS reports "RTC Edit cannot
   open files in the RTC Scene format".
3. `class_addMethod(application:openURLs:)` on `WinitApplicationDelegate`
   works on a cold launch and while running, and winit is unaffected.

## Design

### Bundle (`xtask/src/bundle.rs`)

`App` gains `documents: bool`: true for RTC Edit, false for RTC Viewer.
When true, `info_plist` adds, built from data and escaped like the other
values:

```xml
<key>UTExportedTypeDeclarations</key>
<array><dict>
  <key>UTTypeIdentifier</key><string>no.sb21.rtc.scene</string>
  <key>UTTypeDescription</key><string>RTC Scene</string>
  <key>UTTypeConformsTo</key><array><string>public.json</string></array>
  <key>UTTypeTagSpecification</key><dict>
    <key>public.filename-extension</key><array><string>rtc</string></array>
  </dict>
</dict></array>
<key>CFBundleDocumentTypes</key>
<array>
  <dict>
    <key>CFBundleTypeName</key><string>RTC Scene</string>
    <key>CFBundleTypeRole</key><string>Editor</string>
    <key>LSHandlerRank</key><string>Owner</string>
    <key>LSItemContentTypes</key><array><string>no.sb21.rtc.scene</string></array>
  </dict>
  <dict>
    <key>CFBundleTypeName</key><string>JSON</string>
    <key>CFBundleTypeRole</key><string>Editor</string>
    <key>LSHandlerRank</key><string>Alternate</string>
    <key>LSItemContentTypes</key><array><string>public.json</string></array>
  </dict>
</array>
```

`--install` already re-registers the bundle with Launch Services.

### Extensions (`rtc-edit/src/paths.rs`)

```rust
/// The extensions Open… and Save As accept, the scene type's first.
pub const SCENE_EXTENSIONS: [&str; 2] = ["rtc", "json"];

/// A scene's file name, or `untitled.rtc` for a scene with no file yet.
pub fn scene_name(path: Option<&Path>) -> String;   // default changes

/// The name Save As suggests: the scene's stem with `.rtc`
/// (`pencil.json` → `pencil.rtc`, no file → `untitled.rtc`).
pub fn save_as_name(path: Option<&Path>) -> String;
```

`main.rs`: both `add_filter("scene", …)` calls take `&SCENE_EXTENSIONS`;
Save As's `set_file_name` takes `save_as_name(self.path.as_deref())`.
The window title and the preview file keep using `scene_name`.

### Finder opens (`rtc-edit/src/finder.rs`, new)

```rust
/// Add `application:openURLs:` to winit's application delegate class.
/// On failure (no delegate, `class_addMethod` false) print
/// `rtc-edit: opening from Finder is off: <why>` and carry on.
pub fn install(mtm: MainThreadMarker);

/// The stream of opened files for a subscription: anything queued before
/// it started first, then each open as it comes. A restarted
/// subscription takes over from the old one.
pub fn opens() -> UnboundedReceiver<Vec<PathBuf>>;

/// The file paths among `urls` (as `absoluteString`s), in order;
/// non-file URLs skipped, percent-encoding decoded.
fn file_paths(urls: &[String]) -> Vec<PathBuf>;
```

Inside, a static `Mutex` holds an `Inbox`: a queue and an optional
sender. `Inbox::push(paths)` sends when the sender is there and still
open, else queues; `Inbox::take_over(sender)` sends the queue, then keeps
the sender. Both are plain code, tested without AppKit.

The added method only reads its arguments, calls `file_paths` and pushes;
it never touches winit. The `transmute` to `Imp`, `class_addMethod` and
the type encoding `v@:@@` carry SAFETY comments.

### Editor (`main.rs`)

- `boot()` calls `finder::install` before the scene is loaded.
- `subscription()` gains `Subscription::run(finder_opens)`, which maps
  `finder::opens()` to `Message::OpenFiles(Vec<PathBuf>)`.
- `open_files(paths)`:
  1. No paths: nothing.
  2. `ask_unsaved()`; cancelled → stop.
  3. `SceneDescription::from_path(first)`: ok → `switch`,
     `remember(SCENE_FOLDER)`, `note_recent`, then the status from
     `more_not_opened(paths.len())`; error → status `<name>: <error>`,
     the scene stays.
  4. `window::gain_focus(main)`.
- `more_not_opened(n: usize) -> String`, pure: empty for 0 or 1,
  `(N-1 more not opened)` otherwise, shown as
  `<name> (2 more not opened)`.
- A cold launch from Finder: `boot()` starts the untitled scene with
  nothing unsaved, so the Finder file replaces it without an alert. The
  command line is unchanged: `rtc-edit x.rtc` loads in `boot()`, and a bad
  file still exits.
- Opening the scene that is already open reloads it from disk, as Open…
  does.

## Testing

Unit tests:
- `bundle.rs`: RTC Edit's plist has both declarations, the type, the
  extension and both ranks; RTC Viewer's has neither key.
- `paths.rs`: `scene_name(None)` is `untitled.rtc`; `save_as_name` for
  `None`, `pencil.json`, `pencil.rtc`, a name without an extension and a
  name with dots (`a.b.json` → `a.b.rtc`).
- `finder.rs`: `file_paths` decodes `%20`, skips `https://…`, keeps
  order; `Inbox` hands queued paths over first, in order, sends later
  paths directly, and queues again after the receiver is dropped.
- `main.rs`: `more_not_opened` for 0, 1, 2 and 3.

Checked by eye after `cargo xtask bundle --install`:
1. Double-click a `.rtc` with the editor closed: it opens, editor in
   front, viewer behind.
2. With an unsaved edit, double-click another `.rtc`: "Save changes?";
   Cancel keeps the scene.
3. Open two `.rtc` files at once: the first opens, the status says one
   more was not opened.
4. Open With ▸ RTC Edit on a `.json`: it opens; double-clicking `.json`
   still goes to its usual app.
5. A broken `.rtc`: the error in the status line, the editor keeps going.
6. New, then Save As: `untitled.rtc` is suggested; Save As on
   `pencil.json` suggests `pencil.rtc`.
7. A Finder-opened file is first in Open Recent.

## Docs

rtc-edit's README: `.rtc` files open from Finder, Open With for `.json`;
milestone 33 marked done; a `src/finder.rs` row in the layout table; the
scene file type leaves the Planned line. That line already names
"importing a scene into the current one as a group (File ▸ Import
Scene…)", added with this spec.
