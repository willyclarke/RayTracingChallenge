# Milestone 25: a native menu bar — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** rtc-edit gets a macOS menu bar (rtc-edit, File, Edit, Scene, Viewer) that owns every Cmd key and greys what can't act.

**Architecture:** A pure part of `src/menu.rs` (commands, their ids, the editor state the greying needs, `enabled`) is unit-tested; a `muda` part builds, installs and refreshes the bar. `main.rs` installs it when the window opens, polls `MenuEvent::receiver()` every 50 ms and turns each event into the editor's existing message, and refreshes the bar after every message.

**Tech Stack:** Rust 2024, iced 0.14, muda 0.21 (macOS).

**Spec:** `RUST/rtc-edit/docs/superpowers/specs/2026-10-21-metal-25-menu-bar-design.md`

## Global Constraints

- Branch `metal/25-menu-bar`; commits end ` (milestone 25)`, prefixes `ENHANCEMENT:` / `TEST:` / `DOCS:` / `FIX:`; never push; never add `RUST/rtc-edit/untitled.json`.
- Menus and keys exactly as the spec's table: rtc-edit (About rtc-edit, Quit Cmd-Q); File (New Cmd-N, Open… Cmd-O, Save Cmd-S, Save As… Cmd-Shift-S); Edit (Undo Cmd-Z, Redo Cmd-Shift-Z, Delete Cmd-Delete); Scene (Add Sphere, Add Cube, Add Plane, Add Cylinder, Add Cone, Add Light; no keys); Viewer (Load Edits Cmd-L, Start Viewer Cmd-R, Follow the Viewer (check item, no key), Keep Followed View Cmd-K).
- No predefined Cut, Copy, Paste or Select All.
- Tests: `cargo test --lib` in `RUST/rtc-edit`; `cargo clippy --all-targets` clean; `cargo fmt`. The "future-incompat … block v0.1.6" warning is pre-existing.
- Test names: descriptive snake_case.

## Rulings made in this plan

- **`edit::can_delete` already exists** (`Result<(), String>`, tested in `edit.rs`), so no new function: the state uses `can_delete(..).is_ok()`.
- **Delete is also greyed while a box is shaded or edited** (`left` or `editing` set), as Cmd-Delete works today (README: "They work with no box shaded or edited"). Otherwise Cmd-Backspace in a box, which the menu now takes first, would delete the item.
- **No `Message::Menu(Command)`**: a `Message::Poll` drains the menu events, and each command becomes the existing message (`Command::StartViewer` → `ShowViewer`, `Follow` → `Follow(!following)`, …) handled by `handle`.

## Review Focus

1. **The 50 ms poll and a two-key wait** (`a` then `s`, `d d`): the poll must not cancel the wait. Task 2: `Poll` joins the messages that don't cancel `pending`.
2. **Cmd-Backspace or Cmd-L while typing in a box**: no item deleted; no `l` typed. Task 2: Delete greyed while `editing`; the `e.command` guard in `Message::Edit` stays.
3. **Clicking Follow the Viewer**: AppKit flips the tick itself on a click; when following doesn't start, the tick must go back. Task 2: `refresh` sets the tick whenever it differs from `following`.
4. **Each key acts once**: the menu sends it, and iced's copy of the key does nothing. Task 2: the Cmd table and `DELETE` are gone from `shortcut`; by eye.
5. **A viewer that ends**: Load Edits, Follow and Keep grey, Start Viewer ungreys, within a tick. Task 2: the state reads `e.live`; by eye.

---

### Task 1: the pure menu model and `History::can_undo/can_redo`

**Files:**
- Create: `RUST/rtc-edit/src/menu.rs` (`pub mod menu;` in `lib.rs`).
- Modify: `RUST/rtc-edit/src/history.rs`.

**Interfaces:**
- Produces:
  - `#[derive(Debug, Clone, Copy, PartialEq, Eq)] pub enum Command { New, Open, Save, SaveAs, Undo, Redo, Delete, Add(Item), Load, StartViewer, Follow, Keep }` (`Item` from `edit`).
  - `pub const ALL: [Command; 17]`, in menu order.
  - `pub fn id(command: Command) -> &'static str`: `new`, `open`, `save`, `save-as`, `undo`, `redo`, `delete`, `add-sphere`, `add-cube`, `add-plane`, `add-cylinder`, `add-cone`, `add-light`, `load`, `start-viewer`, `follow`, `keep`.
  - `pub fn command(id: &str) -> Option<Command>`.
  - `#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)] pub struct State { pub can_undo: bool, pub can_redo: bool, pub can_delete: bool, pub viewer_running: bool, pub following: bool }`.
  - `pub fn enabled(command: Command, state: &State) -> bool` (the spec's greying table).
  - `History::can_undo(&self) -> bool`, `History::can_redo(&self) -> bool`.

- [ ] **Step 1: Write the failing tests.**
  - `menu.rs`:
    - `every_command_id_maps_back`: for each `c` in `ALL`, `command(id(c)) == Some(c)`; the ids are distinct; `command("nonsense") == None`.
    - `undo_and_redo_need_a_step`: default state → both disabled; `can_undo` → Undo only; `can_redo` → Redo only.
    - `delete_needs_a_deletable_item`: `can_delete` false/true → Delete disabled/enabled.
    - `viewer_items_follow_the_viewer`: `viewer_running` false → Load, Follow disabled, StartViewer enabled; true → the opposite.
    - `keep_needs_following`: `following` false/true → Keep disabled/enabled.
    - `the_rest_are_always_enabled`: with `State::default()`, New, Open, Save, SaveAs and every `Add(_)` enabled.
  - `history.rs`: `can_undo_and_can_redo_follow_the_steps`: fresh → both false; after a `record` → `can_undo` only; after `undo` → `can_redo` only (and `can_undo` false).
- [ ] **Step 2: Run** `cargo test --lib menu history`. Expected: FAIL (missing items).
- [ ] **Step 3: Implement** `menu.rs`'s pure part and the two `History` methods.
- [ ] **Step 4: Run** `cargo test --lib`, `cargo clippy --all-targets`, `cargo fmt`. Expected: pass, clean.
- [ ] **Step 5: Commit** `ENHANCEMENT: rtc-edit models its menu commands and when each can act (milestone 25)`.

### Task 2: the bar, wiring, README

**Files:**
- Modify: `RUST/rtc-edit/Cargo.toml`, `RUST/rtc-edit/src/menu.rs`, `RUST/rtc-edit/src/main.rs`, `RUST/rtc-edit/README.md`.

**Interfaces:**
- Consumes: Task 1's `Command`, `ALL`, `id`, `command`, `State`, `enabled`, `can_undo`, `can_redo`.
- Produces: `pub struct Menus`; `Menus::new() -> Menus`; `Menus::install(&self)`; `Menus::refresh(&mut self, state: &State)`; `Message::Poll`.

- [ ] **Step 1: `Menus` in `menu.rs`** (`muda = "0.21"` in `Cargo.toml`). `Menus` keeps the `muda::Menu`, every `(Command, MenuItem)` but Follow, the Follow `CheckMenuItem`, and `last: Option<State>`.
  - `new()`: the five submenus of the Global Constraints, items made with `MenuItem::with_id(id(c), label, true, accel)` (`Accelerator::new(Modifiers::SUPER [| SHIFT], Code::…)`; Delete is `Code::Backspace`); rtc-edit has `PredefinedMenuItem::about(Some("About rtc-edit"), Some(AboutMetadata { name: Some("rtc-edit".into()), version: Some(env!("CARGO_PKG_VERSION").into()), ..Default::default() }))`, a separator and `PredefinedMenuItem::quit(None)`; `expect("<menu name> menu")` on each `Result`.
  - `install(&self)`: `init_for_nsapp()`.
  - `refresh(&mut self, state)`: when `Some(*state) != last`, `set_enabled(enabled(c, state))` on every item (Follow too) and store it; then, every call, `set_checked(state.following)` when `is_checked()` differs.
- [ ] **Step 2: Wire `main.rs`.**
  - `Editor.menus: Option<Menus>` (`None` at boot); on `window::Event::Opened`, when `None`, `Menus::new()`, `install()`, store.
  - `Message::Poll` (doc: "Time to read the menu bar's clicks and keys"), from `iced::time::every(Duration::from_millis(50))` in `subscription`; added to the messages that don't cancel `pending` in `handle`.
  - `Poll`: drain `MenuEvent::receiver().try_recv()`; for each id with `menu::command(&id.0) == Some(c)`, `handle(e, message(c, e.follower.following()))`; batch the tasks. Unknown ids are ignored.
  - `fn message(command: Command, following: bool) -> Message`: the matching existing message; `StartViewer` → `ShowViewer`, `Follow` → `Follow(!following)`.
  - `Editor::menu_state(&self) -> State`: `history.can_undo/can_redo()`; `can_delete` = selected item with `edit::can_delete(..).is_ok()`, `left` and `editing` both `None`; `viewer_running: live`; `following: follower.following()`.
  - `update`: after `handle`, `if let Some(m) = &mut e.menus { m.refresh(&state) }`.
  - `shortcut`: drop the Cmd-character table and the Cmd-Backspace/Delete arm; remove `DELETE` and its branch in `key`; `Key(OTHER)` stays. Update `shortcut`'s doc and the module doc's key mentions to name the menus.
- [ ] **Step 3: README**: the Cmd key table becomes a menus table (menu, items, keys, when greyed); `d d or Cmd-Delete` rows and the "Adding and deleting" text say Edit ▸ Delete (Cmd-Delete, the Delete key); Cut/Copy/Paste/Select All left to the boxes; milestone 25 in the list.
- [ ] **Step 4: Run** `cargo test --lib`, `cargo clippy --all-targets`, `cargo fmt`, `cargo build --release`. Expected: pass, clean.
- [ ] **Step 5: Stop for the by-eye check** (report; do not commit): each item by click and by key, once each; greying as you undo, select the camera, shade or edit a box, close and Cmd-R the viewer, follow; the Follow tick by menu and by checkbox; Cmd-Backspace and Cmd-L while typing in a box; `a` then `s` and `d d` still work; copy and paste in a box; About; Quit.
- [ ] **Step 6: Commit when told** `ENHANCEMENT: rtc-edit has a native menu bar that owns the Cmd keys (milestone 25)`.
