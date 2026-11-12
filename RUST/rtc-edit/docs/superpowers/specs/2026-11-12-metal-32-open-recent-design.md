# Milestone 32: Open Recent

Branch `metal/32-open-recent`, merged into develop with `--no-ff` when done.
Milestone 31 built the app bundle; after 32 come a scene file type with
opening from Finder, the viewer's settings, and + csg / + group.

## Goal

Getting back to a scene means Cmd-O and walking the Open panel to it. RTC
Edit should keep the scenes used last and offer them in File ▸ Open Recent,
the Mac way, across restarts.

Success:
- After a relaunch, the scenes opened or saved last are one click away in
  File ▸ Open Recent, newest first.
- Picking one behaves like Open…: the "Save changes?" alert when there are
  unsaved edits, and the scene folder is remembered.
- A moved or deleted scene says so in the status line and leaves the
  menu; nothing breaks.

## Out of scope

- Hiding or greying missing files when the menu is drawn (no disk checks
  for the menu).
- AppKit's own recent-documents list (`NSDocumentController`): its items
  arrive as `application:openFile:`, the Apple Event iced/winit don't pass
  on. That belongs to the Finder-open milestone, which can call
  `recent::add` too.
- A key for Open Recent or its entries; a configurable list length.
- Open Recent in the help window.

## Decisions

### What goes on the list

A scene goes to the front of the list when it is
- opened with Open…,
- saved with Save As…, or
- given on the command line at launch.

A plain Save doesn't touch the list: the scene is already at the front
from opening or Save As. Picking an entry moves it to the front.

At most 10 entries (the macOS default). A path already on the list moves
rather than appearing twice.

### Where it lives

The settings file (`rtc-config`) gets one key:

```rust
/// Scenes opened or saved last, newest first, as absolute paths.
pub const RECENT_SCENES: &str = "recent_scenes";
```

```json
"recent_scenes": ["/Users/me/rtc/scenes/pencil.json", "/Users/me/x/untitled.json"]
```

Nothing else in rtc-config changes; rtc-metal ignores the key.

Every change re-reads the list from the file, changes it and writes it
back with `set`, as `scene_folder` does. Two editors running at once keep
each other's entries; a write at the same instant loses at most one change.

### Labels

Each entry shows its file name, `pencil.json`. When two or more entries
share a file name, each of those gets ` — <parent folder name>`:
`untitled.json — examples`, `untitled.json — rtc-edit`. When the parent
names match too, the labels stay alike (rare; accepted). A path with no
file name shows the whole path.

### Missing and broken files

Picking an entry whose file doesn't exist puts `<name>: not found` in the
status line and removes the entry. No "Save changes?" alert: nothing would
open.

A file that exists but doesn't load (bad JSON, say) shows
`<name>: <error>` as Open… does, and **stays** on the list: it is likely
being fixed.

## Design

### `rtc-edit/src/recent.rs` (new)

Plain functions, no menu or AppKit, so all of it is unit-tested:

```rust
pub const MAX: usize = 10;

/// `map[RECENT_SCENES]` as paths: non-strings and empty strings skipped,
/// at most MAX; a missing or non-array key is an empty list.
pub fn from_config(map: &Map<String, Value>) -> Vec<PathBuf>;

/// The list as the settings file stores it.
pub fn to_value(list: &[PathBuf]) -> Value;

/// `path` made absolute (std::path::absolute: the working directory in
/// front of a relative path, symlinks left alone), put first, any equal
/// entry dropped, at most MAX.
pub fn add(list: &[PathBuf], path: &Path) -> Vec<PathBuf>;

/// The list without `path`.
pub fn remove(list: &[PathBuf], path: &Path) -> Vec<PathBuf>;

/// One menu label per entry, in order (see Labels).
pub fn labels(list: &[PathBuf]) -> Vec<String>;
```

### Menu (`menu.rs`)

- `Command` gains `OpenRecent(usize)` and `ClearRecent`.
- Ids: `recent-0` … `recent-9`, `clear-recent`. `command(id)` parses
  `recent-N` (N < `recent::MAX`) besides the fixed ids.
- `ALL` stays the fixed list, so the help page doesn't change. Neither new
  command has a key; `enabled()` is true for both (their items are not in
  `Menus::items`, so `refresh` never touches them; `set_recent` greys
  Clear Menu when the list is empty); `menu_of` is `"File"`;
  `label(ClearRecent)` is `"Clear Menu"`, and `label(OpenRecent(_))` is
  never used for the menu (entries take their label from `recent::labels`).
- File becomes New, Open…, **Open Recent ▸**, Save, Save As….
- `Menus` keeps `recent: Submenu` and gains
  `set_recent(&mut self, labels: &[String])`: remove the submenu's items,
  append one `MenuItem` per label with id `recent-<i>`, a separator, and
  Clear Menu. With no labels: only Clear Menu, greyed.
- `set_recent` runs only when the list changes; `refresh` is unchanged.

### Editor (`main.rs`)

- `Editor` gains `recent: Vec<PathBuf>`, filled in `boot()` from
  `recent::from_config(&config())`.
- Where the menus are made on the window's first `Opened`, `set_recent`
  runs straight after `Menus::new()`.
- `note_recent(path)`: `recent::add` on the list re-read from the file,
  `remember(RECENT_SCENES, …)`, store in `self.recent`, `set_recent` when
  the menus exist.
- `forget_recent(path)`: the same with `recent::remove`.
- Callers of `note_recent`: `open()` and `save_as()` on success, beside
  `remember(SCENE_FOLDER)`; `boot()` for the launch scene once loaded (a
  bad launch scene still exits as now).
- `menu_message`: `OpenRecent(n)` → `Message::OpenRecent(n)`,
  `ClearRecent` → `Message::ClearRecent`.
- `open_recent(n)`, for `Message::OpenRecent(n)`:
  1. No entry `n` (list shorter): nothing.
  2. File missing: status `<name>: not found`, `forget_recent`.
  3. `ask_unsaved()`; cancelled → stop.
  4. `SceneDescription::from_path`: ok → `switch`, `remember(SCENE_FOLDER)`,
     `note_recent`; error → status `<name>: <error>`, entry kept.
  Picking the open scene reloads it from disk, as Open… would.
- `Message::ClearRecent`: empty list, `remember`, `set_recent(&[])`.
- A failed settings write shows `couldn't save settings: …` (the existing
  `remember`); the list in memory and the menu still change.

## Testing

`recent.rs` unit tests:
- `add`: new path goes first; a duplicate moves to the front; the 11th
  entry drops the oldest; a relative path becomes absolute.
- `remove`: drops the entry, keeps the order; absent path is a no-op.
- `from_config`: missing key, non-array, non-string and empty items
  skipped, longer than MAX cut.
- `to_value` / `from_config` round trip.
- `labels`: all unique; a duplicate pair gets folders; three entries with
  one shared name (only the pair gets folders); no file name.

`menu.rs` tests:
- `recent-0`, `recent-9` and `clear-recent` map back; `recent-10`,
  `recent-x` and `recent-` don't.
- `ALL` and the help page's entries are unchanged (the existing tests).
- Both new commands are enabled.

`main.rs` glue, checked by eye:
1. Open a scene: it is first in File ▸ Open Recent.
2. Save As a new name: the new name is first.
3. Two `untitled.json` from different folders: both show their folder.
4. Rename an entry's file in Finder, pick it: "not found", entry gone.
5. Edit a scene, pick a recent one: "Save changes?" alert; Cancel keeps
   the scene.
6. Clear Menu: only a greyed Clear Menu.
7. Quit and relaunch (also `rtc-edit examples/x.json`): the list is there,
   the launch scene first.

## Docs

rtc-edit's README: an Open Recent line among the File commands, and
milestone 32 marked done.
