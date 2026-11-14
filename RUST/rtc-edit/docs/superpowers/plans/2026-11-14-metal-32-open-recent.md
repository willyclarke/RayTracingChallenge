# Open Recent Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** File ▸ Open Recent in RTC Edit: the last 10 scenes opened, saved as or launched with, newest first, kept in the settings file across restarts.

**Architecture:** A pure module `rtc-edit/src/recent.rs` owns the list (read from / written to the `recent_scenes` key, add, remove, labels). `menu.rs` gains `OpenRecent(n)` / `ClearRecent` commands and a submenu rebuilt by `Menus::set_recent`. `main.rs` glues them: `note_recent` / `forget_recent` keep file, `Editor::recent` and menu in step; `open_recent(n)` mirrors `open()`.

**Tech Stack:** Rust 2024, serde_json, iced 0.14, muda 0.21 (`Submenu::items`/`remove`/`append`).

**Spec:** `rtc-edit/docs/superpowers/specs/2026-11-12-metal-32-open-recent-design.md`

## Global Constraints

- Settings key `recent_scenes`: JSON array of absolute path strings, newest first; at most 10 (`recent::MAX`).
- Added by Open…, Save As… and the launch scene; a plain Save never touches the list; picking an entry moves it first.
- Labels: file name; ` — <parent folder name>` only on entries whose file name is shared. Clear Menu label: `Clear Menu`.
- Ids: `recent-0` … `recent-9`, `clear-recent`. Neither command is in `ALL`, has a key, or appears in the help page.
- Missing file on pick: status `<name>: not found`, entry removed, no "Save changes?" alert. Load error: status `<name>: <error>`, entry kept.
- Settings write failure: existing `couldn't save settings: <error>`; list in memory and menu still change.
- Tests never touch the real `$HOME` settings file; no test sets environment variables.
- Run tests from `RUST/rtc-edit` (or `RUST/rtc-config`); `cargo clippy --all-targets -- -D warnings` clean.
- Commit prefixes ENHANCEMENT/FIX/DOC, `(milestone 32)` at the end of the subject; every message ends with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. Do not push.

## Review Focus

1. The same scene reached two ways — `rtc-edit examples/x.json` at launch (relative) and later Open… (absolute) — must be one entry, not two. Task 1's `add_makes_relative_paths_absolute_and_dedupes` pins it.
2. A hand-edited or old settings file (`recent_scenes` of 15 entries, numbers, `""`, an object) must load without panicking and show at most 10 entries. Task 1's `from_config_*` tests.
3. Cancel in the "Save changes?" alert after picking an entry: the scene, the list and the menu stay as they were. By eye, Task 4 item 5.
4. Save As onto a name already on the list moves it to the top rather than duplicating it. Task 1's `add_moves_a_duplicate_to_the_front`.
5. A menu id outside the range (`recent-10`, `recent-`, `recent-x`) maps to no command, and `OpenRecent(n)` with `n` past the list does nothing. Task 2's `recent_ids_map_back_only_in_range`; Task 3's `open_recent` step 1.

---

### Task 1: the list (`rtc-config` key + `recent.rs`)

**Files:**
- Modify: `RUST/rtc-config/src/lib.rs` (after `EDITOR_FRAME`)
- Create: `RUST/rtc-edit/src/recent.rs`
- Modify: `RUST/rtc-edit/src/lib.rs` (add `pub mod recent;`, alphabetical)

**Interfaces:**
- Produces:
  - `rtc_config::RECENT_SCENES: &str = "recent_scenes"` with doc `/// Scenes opened or saved last, newest first, as absolute paths.`
  - `recent::MAX: usize = 10`
  - `recent::from_config(map: &rtc_config::Map<String, rtc_config::Value>) -> Vec<PathBuf>`
  - `recent::to_value(list: &[PathBuf]) -> rtc_config::Value`
  - `recent::add(list: &[PathBuf], path: &Path) -> Vec<PathBuf>`
  - `recent::remove(list: &[PathBuf], path: &Path) -> Vec<PathBuf>`
  - `recent::labels(list: &[PathBuf]) -> Vec<String>`

- [ ] **Step 1: Write the failing tests** in `recent.rs`'s `#[cfg(test)] mod tests` (build maps with `rtc_config::json!(...).as_object().unwrap().clone()`; paths with `PathBuf::from`):

  - `add_puts_a_new_path_first`: `add(&[/a/x.json], /b/y.json) == [/b/y.json, /a/x.json]`.
  - `add_moves_a_duplicate_to_the_front`: `add(&[/a, /b, /c], /c) == [/c, /a, /b]`.
  - `add_keeps_at_most_max`: 10 entries `/s/0.json`…`/s/9.json`, add `/s/new.json` → length 10, first `/s/new.json`, `/s/9.json` gone.
  - `add_makes_relative_paths_absolute_and_dedupes`: `let a = add(&[], Path::new("examples/x.json"));` → `a[0].is_absolute()` and `a[0].ends_with("examples/x.json")`; then `add(&a, &a[0])` has length 1.
  - `remove_drops_the_entry_and_keeps_order`: `remove(&[/a, /b, /c], /b) == [/a, /c]`; `remove(&[/a], /z) == [/a]`.
  - `from_config_reads_the_array`: `{"recent_scenes": ["/a.json", "/b.json"]}` → `[/a.json, /b.json]`.
  - `from_config_skips_bad_items_and_caps`: `{}` → empty; `{"recent_scenes": 3}` and `{"recent_scenes": {"a": 1}}` → empty; `["/a", 7, "", null, "/b"]` → `[/a, /b]`; 15 strings → first 10.
  - `to_value_round_trips`: `from_config` of `{"recent_scenes": to_value(&list)}` equals `list` for `[/a/x.json, /b/y.json]`.
  - `labels_are_file_names_when_unique`: `[/a/x.json, /b/y.json]` → `["x.json", "y.json"]`.
  - `labels_add_the_folder_to_shared_names`: `[/p/examples/untitled.json, /p/rtc-edit/untitled.json, /p/pencil.json]` → `["untitled.json — examples", "untitled.json — rtc-edit", "pencil.json"]`.
  - `labels_show_a_nameless_path_whole`: `[PathBuf::from("/")]` → `["/"]`.

- [ ] **Step 2: Run to see them fail**

  Run: `cd RUST/rtc-edit && cargo test --lib recent`
  Expected: compile error, `recent` functions not found.

- [ ] **Step 3: Add `RECENT_SCENES` to rtc-config and implement `recent.rs`**

  Module doc: the recent-scenes list behind File ▸ Open Recent, kept in the settings file. `add` makes the path absolute with `std::path::absolute` (falling back to the path as given if it errors), drops equal entries, inserts first, truncates to `MAX`. `labels`: name = `file_name().to_string_lossy()`, else `path.display()`; an entry gets ` — <parent's file_name>` (U+2014 with a space either side) when another entry has the same name; no parent name → no suffix.

- [ ] **Step 4: Run the tests**

  Run: `cd RUST/rtc-edit && cargo test --lib recent && cd ../rtc-config && cargo test && cd ../rtc-edit && cargo clippy --all-targets -- -D warnings`
  Expected: all pass, no warnings.

- [ ] **Step 5: Commit**

  ```bash
  git add RUST/rtc-config/src/lib.rs RUST/rtc-edit/src/recent.rs RUST/rtc-edit/src/lib.rs
  git commit -m "ENHANCEMENT: rtc-edit keeps a recent-scenes list in the settings file (milestone 32)"
  ```

### Task 2: the menu (`menu.rs`)

**Files:**
- Modify: `RUST/rtc-edit/src/menu.rs`

**Interfaces:**
- Consumes: `recent::MAX`.
- Produces:
  - `Command::OpenRecent(usize)`, `Command::ClearRecent`
  - `menu::command(id: &str) -> Option<Command>` also returns `OpenRecent(n)` for `recent-<n>`, `n < recent::MAX`, and `ClearRecent` for `clear-recent`.
  - `Menus::set_recent(&mut self, labels: &[String])`

- [ ] **Step 1: Write the failing tests** in `menu.rs`'s tests:

  - `recent_ids_map_back_only_in_range`: `command("recent-0") == Some(OpenRecent(0))`, `command("recent-9") == Some(OpenRecent(9))`, `command("clear-recent") == Some(ClearRecent)`; `command("recent-10")`, `command("recent-")`, `command("recent-x")`, `command("recent--1")` are `None`; `id(OpenRecent(3)) == "recent-3"`.
  - `recent_commands_are_file_items_outside_all`: `!ALL.contains(&OpenRecent(0)) && !ALL.contains(&ClearRecent)`; `menu_of` of both is `"File"`; `key` of both is `None`; `label(ClearRecent) == "Clear Menu"`; `enabled` of both with `State::default()` is true.

  The existing `every_command_id_maps_back` and help-page tests stay as they are and must keep passing.

- [ ] **Step 2: Run to see them fail**

  Run: `cd RUST/rtc-edit && cargo test --lib menu`
  Expected: compile error, no variant `OpenRecent`.

- [ ] **Step 3: Implement**

  - `id()` returns `&'static str`, so add `const RECENT_IDS: [&str; recent::MAX] = ["recent-0", …, "recent-9"]` and return `RECENT_IDS[n]` (`OpenRecent(n)` with `n >= MAX` never exists: `command` won't make one). `command(id)` checks `ALL`, then `clear-recent`, then `RECENT_IDS.iter().position`.
  - `label(OpenRecent(_)) == "Open Recent"` (the submenu's title; entries take labels from `recent::labels`).
  - `Menus` gains `recent: Submenu`, made empty (`Submenu::new("Open Recent", true)`) and inserted in File between Open… and Save. The File submenu can no longer use the `submenu(name, &[MenuItem])` helper as-is; build it with `Submenu::with_items("File", true, &[&new, &open, &recent, &save, &save_as])`.
  - `Menus::new()` calls `set_recent(&[])` so the submenu starts with a greyed Clear Menu.
  - `set_recent`: remove every `self.recent.items()` entry; append `MenuItem::with_id(RECENT_IDS[i], label, true, None)` per label; append `PredefinedMenuItem::separator()` (only when there are labels); append `MenuItem::with_id("clear-recent", "Clear Menu", !labels.is_empty(), None)`. These items are not pushed to `self.items`, so `refresh` never re-enables Clear Menu.

- [ ] **Step 4: Run the tests**

  Run: `cd RUST/rtc-edit && cargo test --lib && cargo clippy --lib --tests -- -D warnings`
  Expected: all pass, no warnings. The binary (`main.rs`'s exhaustive `menu_message`) won't build until Task 3; that is expected, so this task checks the library only.

- [ ] **Step 5: Commit**

  ```bash
  git add RUST/rtc-edit/src/menu.rs
  git commit -m "ENHANCEMENT: File ▸ Open Recent submenu and its commands (milestone 32)"
  ```

### Task 3: the editor (`main.rs`)

**Files:**
- Modify: `RUST/rtc-edit/src/main.rs` — `Editor` struct (~line 137), `boot()` (~236–288), `open()` (~430), `save_as()` (~963), `menu_message` (~1130), `handle` (`MenuEvent` arm ~1272), window `Opened` arm (~1334), `Message` enum (~150–200).

**Interfaces:**
- Consumes: `recent::{from_config, to_value, add, remove, labels}`, `rtc_config::RECENT_SCENES`, `Menus::set_recent`, `Command::{OpenRecent, ClearRecent}`.
- Produces: `Message::OpenRecent(usize)`, `Message::ClearRecent`; `Editor::recent: Vec<PathBuf>`; `Editor::note_recent(&mut self, path: &Path)`, `Editor::forget_recent(&mut self, path: &Path)`, `Editor::open_recent(&mut self, n: usize)`, `Editor::clear_recent(&mut self)`.

No unit tests: this is window glue (the spec checks it by eye in Task 4). Keep the logic in `recent.rs`; these methods only sequence calls.

- [ ] **Step 1: State and helpers**

  - `Editor.recent` initialised in `boot()` with `recent::from_config(&config())`.
  - A private `fn set_recent_list(&mut self, list: Vec<PathBuf>)`: `self.remember(RECENT_SCENES, recent::to_value(&list))`, `self.recent = list`, and when `self.menus` is `Some`, `set_recent(&recent::labels(&self.recent))`.
  - `note_recent(path)` = `set_recent_list(recent::add(&recent::from_config(&config()), path))`; `forget_recent(path)` the same with `remove`; `clear_recent()` = `set_recent_list(Vec::new())`.

- [ ] **Step 2: Callers**

  - `open()` success: `self.note_recent(&path)` after `remember(SCENE_FOLDER, …)` (clone the path before `switch` takes it).
  - `save_as()` success: `self.note_recent(&new)` beside `remember(SCENE_FOLDER, …)`.
  - `boot()`: in the `Some((desc, path))` arm, `note_recent(&path)` after `switch`.
  - Window `Opened`, where menus are made: after `Menus::new()`, `menus.set_recent(&recent::labels(&e.recent))` before storing them.

- [ ] **Step 3: Messages**

  - `menu_message`: `Command::OpenRecent(n) => Message::OpenRecent(n)`, `Command::ClearRecent => Message::ClearRecent`.
  - `handle`: `Message::OpenRecent(n) => e.open_recent(n)`, `Message::ClearRecent => e.clear_recent()`.
  - `open_recent(n)`:
    1. `let Some(path) = self.recent.get(n).cloned() else { return };`
    2. `!path.exists()` → `self.status = format!("{}: not found", file_name(&path))`, `forget_recent(&path)`, return.
    3. `if !self.ask_unsaved() { return; }`
    4. `SceneDescription::from_path(&path)`: `Ok` → `switch(desc, Some(path.clone()))`, `remember(SCENE_FOLDER, folder_value(&path))`, `note_recent(&path)`; `Err(e)` → `self.status = format!("{}: {e}", file_name(&path))`.

- [ ] **Step 4: Build, test, lint**

  Run: `cd RUST/rtc-edit && cargo build && cargo test && cargo clippy --all-targets -- -D warnings`
  Expected: builds; all tests pass; no warnings.

- [ ] **Step 5: Commit**

  ```bash
  git add RUST/rtc-edit/src/main.rs
  git commit -m "ENHANCEMENT: Open Recent opens, remembers and forgets scenes; Clear Menu empties it (milestone 32)"
  ```

### Task 4: docs and the eye check

**Files:**
- Modify: `RUST/rtc-edit/README.md` (File commands list; milestone list marks 32 done)

- [ ] **Step 1: README** — add a line among the File commands: `File ▸ Open Recent — the last 10 scenes opened, saved as or launched with; Clear Menu empties it.` Mark milestone 32 done in the same style as 31.

- [ ] **Step 2: Eye check, with the user** (`cd RUST/rtc-edit && cargo run -- examples/<scene>.json`), each item reported pass/fail:
  1. Open a scene: it is first in File ▸ Open Recent.
  2. Save As a new name: the new name is first.
  3. Two `untitled.json` from different folders: both show their folder.
  4. Rename an entry's file in Finder, pick it: "not found", entry gone.
  5. Edit a scene, pick a recent one: "Save changes?" alert; Cancel keeps the scene and the list.
  6. Clear Menu: only a greyed Clear Menu.
  7. Quit and relaunch with `cargo run -- examples/x.json`: the list is there, the launch scene first, and only once.

- [ ] **Step 3: Commit**

  ```bash
  git add RUST/rtc-edit/README.md
  git commit -m "DOC: rtc-edit's README has Open Recent and marks milestone 32 done"
  ```
