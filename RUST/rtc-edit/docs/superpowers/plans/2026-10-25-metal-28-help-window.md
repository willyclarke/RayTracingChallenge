# Milestone 28: help window Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Help ▸ rtc-edit Help (Cmd-?) opens a second iced window with one scrolling cheat-sheet page: menu commands (from the menu's own data), editor keys, viewer, command line.

**Architecture:** menu.rs gains pure `label`/`key`/`menu_of` that both the menu bar and the help page read. A new pure `help.rs` holds the page as data (`Section` → `Group` → `Entry` of `Chip`s) plus `help::view`. A new pure `windows.rs` decides what a key or a close does per window. main.rs switches from `iced::application` to `iced::daemon` with two window ids.

**Tech Stack:** Rust 2024, iced 0.14 (`iced::daemon`, `window::open`/`close`/`gain_focus`), muda 0.21.

**Spec:** `RUST/rtc-edit/docs/superpowers/specs/2026-10-25-metal-28-help-window-design.md`

All paths below are relative to `RUST/rtc-edit/`. Run tests with `cargo test` from there (baseline: 138 lib tests pass).

## Global Constraints

- Help menu item label `rtc-edit Help`, key Cmd-Shift-Slash (shown ⇧ ⌘ ?), never greyed; the Help submenu is `set_as_help_menu_for_nsapp`.
- Help window title `rtc-edit Help`, size 560×720, opened at the editor's top-right (`frame.x + frame.width, frame.y`).
- Editor window as today: 900×640 at `Point::ORIGIN`, title `rtc-edit — <name>[ •]`.
- Closing the editor → `iced::exit()` (`clean_up` still runs via `atexit`); closing help only clears `Editor.help`.
- In help, Esc closes it and other keys/clicks are ignored; Cmd (modifier) tracking passes from any window.
- Modifier chip order ⌃ ⌥ ⇧ ⌘; Shift+Slash → `?`, Backspace → `⌫`.
- No search, clickable entries, greying in help, remembered help frame, help book, localisation.
- Every commit message ends with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`; do not push.

## Refinement of the spec's data model

The spec's `Section { title, entries }` becomes `Section { title, groups, monospace }` with `Group { title, entries }`: the menu section groups by menu, editor keys by place, the viewer by Mouse / Keys / Menus, the command line by program. `monospace` is true only for the command line. `Chip` and `Entry` are as in the spec.

## Review Focus

1. Cmd-? pressed while help is in front: muda takes it (menu key), so it must raise help, not open a second one — `Editor.help` is `Some` until Closed arrives.
2. Closing help with the red button, then Cmd-?: Closed must clear `Editor.help` so a new window opens (test in Task 3: `closed(help) == Closed::Help`).
3. A pending `a`/`d` in the editor when a key is typed in help: the dropped key must not cancel or complete it (routing happens before `handle`, Task 4).
4. Esc in help while a box is edited in the editor: only help closes; the editor's box stays edited (route returns `CloseHelp`, never forwards).
5. Window events of help (Opened/Moved/Resized) must not overwrite the editor's tracked frame, else the viewer and next help open in the wrong place (Task 4 matches on `id == e.main`).

---

### Task 1: menu.rs label / key / menu_of and the Help command

**Files:**
- Modify: `src/menu.rs`
- Modify: `src/main.rs:984-1000` (`menu_message`), the `Message` enum, `handle`

**Interfaces:**
- Produces:
  - `Command::Help`; `ALL: [Command; 18]` with `Help` last; `id(Help) == "help"`.
  - `pub fn label(c: Command) -> &'static str`
  - `pub fn key(c: Command) -> Option<(Modifiers, Code)>`
  - `pub fn menu_of(c: Command) -> &'static str` — `"File" | "Edit" | "Scene" | "Viewer" | "Help"`
  - `Message::Help` in main.rs (handled in Task 4; here it does nothing).

- [ ] **Step 1: Write the failing tests** in `src/menu.rs` `mod tests`:

```rust
#[test]
fn labels_and_keys_are_the_menus() {
    let cmd = Modifiers::META;
    let shift = Modifiers::META | Modifiers::SHIFT;
    assert_eq!(label(Command::New), "New");
    assert_eq!(key(Command::New), Some((cmd, Code::KeyN)));
    assert_eq!(label(Command::Open), "Open…");
    assert_eq!(key(Command::SaveAs), Some((shift, Code::KeyS)));
    assert_eq!(key(Command::Redo), Some((shift, Code::KeyZ)));
    assert_eq!(key(Command::Delete), Some((cmd, Code::Backspace)));
    assert_eq!(label(Command::Add(Item::Cone)), "Add Cone");
    assert_eq!(key(Command::Add(Item::Cone)), None);
    assert_eq!(label(Command::Load), "Load Edits");
    assert_eq!(key(Command::StartViewer), Some((cmd, Code::KeyR)));
    assert_eq!(label(Command::Follow), "Follow the Viewer");
    assert_eq!(key(Command::Follow), None);
    assert_eq!(key(Command::Keep), Some((cmd, Code::KeyK)));
    assert_eq!(label(Command::Help), "rtc-edit Help");
    assert_eq!(key(Command::Help), Some((shift, Code::Slash)));
}

#[test]
fn every_command_has_a_label_and_a_menu() {
    for c in ALL {
        assert!(!label(c).is_empty(), "{c:?}");
        assert!(["File", "Edit", "Scene", "Viewer", "Help"].contains(&menu_of(c)), "{c:?}");
    }
    assert_eq!(menu_of(Command::SaveAs), "File");
    assert_eq!(menu_of(Command::Delete), "Edit");
    assert_eq!(menu_of(Command::Add(Item::Light)), "Scene");
    assert_eq!(menu_of(Command::Keep), "Viewer");
    assert_eq!(menu_of(Command::Help), "Help");
    assert!(enabled(Command::Help, &State::default()));
}
```

- [ ] **Step 2: Run** `cargo test --lib menu::` — Expected: FAIL to compile (`label`, `key`, `menu_of`, `Command::Help` missing).

- [ ] **Step 3: Implement** in `src/menu.rs`:
  - `Command::Help`, appended to `ALL` (length 18), `id` → `"help"`, `enabled` → `true` (add to the always arm).
  - `label`, `key`, `menu_of` as matches holding today's literals from `Menus::new`.
  - `Menus::new`: `item` becomes `|c: Command| -> MenuItem`, reading `label(c)` and `key(c)`; Follow's `CheckMenuItem` uses `label(Command::Follow)`. Add a `"Help"` submenu with `item(Command::Help)`, appended after Viewer, kept in a new field `help: Submenu`.
  - `install`: after `init_for_nsapp`, `self.help.set_as_help_menu_for_nsapp()`.
  - `main.rs`: `Message::Help` (doc: `/// Help ▸ rtc-edit Help (Cmd-?).`), `Command::Help => Message::Help` in `menu_message`, `Message::Help => {}` in `handle` (replaced in Task 4).

- [ ] **Step 4: Run** `cargo test` — Expected: all pass (138 + 2). `cargo clippy --all-targets` clean.

- [ ] **Step 5: Commit**

```bash
git add src/menu.rs src/main.rs
git commit -m "ENHANCEMENT: rtc-edit's menu labels and keys as data, and a Help menu (milestone 28)"
```

### Task 2: help.rs — the page as data

**Files:**
- Create: `src/help.rs`
- Modify: `src/lib.rs` (`pub mod help;`, alphabetical)

**Interfaces:**
- Consumes: `menu::{ALL, label, key, menu_of}` (Task 1).
- Produces:

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Chip { Modifier(&'static str), Key(&'static str), Mouse(&'static str), Word(&'static str) }
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Entry { pub chips: Vec<Chip>, pub text: &'static str }
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Group { pub title: &'static str, pub entries: Vec<Entry> }
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Section { pub title: &'static str, pub groups: Vec<Group>, pub monospace: bool }
pub fn chips(modifiers: Modifiers, code: Code) -> Vec<Chip>;
pub fn sections() -> Vec<Section>;
```

- [ ] **Step 1: Write the failing tests** in `src/help.rs` `mod tests`:

```rust
#[test]
fn chips_put_modifiers_in_mac_order() {
    use Chip::*;
    let cmd = Modifiers::META;
    let shift = Modifiers::META | Modifiers::SHIFT;
    assert_eq!(chips(shift, Code::KeyS), [Modifier("⇧"), Modifier("⌘"), Key("S")]);
    assert_eq!(chips(cmd, Code::Backspace), [Modifier("⌘"), Key("⌫")]);
    assert_eq!(chips(shift, Code::Slash), [Modifier("⇧"), Modifier("⌘"), Key("?")]);
    let all = Modifiers::CONTROL | Modifiers::ALT | Modifiers::SHIFT | Modifiers::META;
    assert_eq!(chips(all, Code::KeyZ), [Modifier("⌃"), Modifier("⌥"), Modifier("⇧"), Modifier("⌘"), Key("Z")]);
}

#[test]
fn every_menu_command_appears_once() {
    let menus = &sections()[0];
    let entries: Vec<&Entry> = menus.groups.iter().flat_map(|g| &g.entries).collect();
    assert_eq!(entries.len(), ALL.len());
    for c in ALL {
        let group = menus.groups.iter().find(|g| g.title == menu_of(c)).expect(menu_of(c));
        let found: Vec<_> = group.entries.iter().filter(|e| e.text == label(c)).collect();
        assert_eq!(found.len(), 1, "{c:?}");
        let want = key(c).map_or(vec![], |(m, k)| chips(m, k));
        assert_eq!(found[0].chips, want, "{c:?}");
    }
    let titles: Vec<_> = menus.groups.iter().map(|g| g.title).collect();
    assert_eq!(titles, ["File", "Edit", "Scene", "Viewer", "Help"]);
}

#[test]
fn every_menu_key_has_a_name() {
    for c in ALL {
        if let Some((m, k)) = key(c) {
            assert!(!chips(m, k).contains(&Chip::Key("…")), "{c:?}");
        }
    }
}

#[test]
fn sections_come_in_order_and_none_is_empty() {
    let s = sections();
    let titles: Vec<_> = s.iter().map(|s| s.title).collect();
    assert_eq!(titles, ["Menu commands", "Editor keys", "Viewer", "Command line"]);
    for section in &s {
        assert!(!section.groups.is_empty(), "{}", section.title);
        for g in &section.groups {
            assert!(!g.entries.is_empty(), "{} / {}", section.title, g.title);
        }
    }
    assert_eq!(s.iter().map(|s| s.monospace).collect::<Vec<_>>(), [false, false, false, true]);
}
```

- [ ] **Step 2: Run** `cargo test --lib help::` — Expected: FAIL to compile.

- [ ] **Step 3: Implement** `src/help.rs` (module doc: the help window's page, as data built from the menu and fixed tables):
  - `chips`: modifier chips in ⌃ ⌥ ⇧ ⌘ order, then one `Key`. Key names by match: `KeyK/L/N/O/R/S/Z` → the letter, `Backspace` → `"⌫"`, `Slash` → `"?"` with Shift else `"/"`, anything else → `"…"`.
  - Menu section: one `Group` per menu in `ALL`'s order (titles from `menu_of`), entries `label` + `chips(key)` (empty when `None`).
  - Editor keys: groups `Tree`, `Shaded box`, `Editing`, `Anywhere`, with the spec's rows (README key table order); `Word("then")`/`Word("or")` between keys, `Key("↑")`, `Key("Enter")`, `Key("Esc")`, `Key("Tab")`, `Key("Space")`; the "letters act only with no modifier held" line as an `Anywhere` entry with no chips.
  - Viewer (comment: hand-written from rtc-metal `src/window.rs`; keep in step with it): groups `Mouse` (`Mouse("drag")` orbits, `Mouse("scroll")` zooms), `Keys` (Space returns to the scene's view, `i` shows the overlay), `Menus` (⌘S Save Image…, ⌘1 BVH traversal, ⌘2 Metal traversal, ⌘] double samples, ⌘[ halve samples, ⌘H hide).
  - Command line, `monospace: true`, no chips: group `rtc-edit` (`rtc-edit [scene.json]`, `RTC_METAL=/path/to/rtc-metal rtc-edit scene.json`), group `rtc-metal` (`rtc-metal scene.json [-o out.ppm] [--samples N] [--window [--watch]] [--traversal bvh|metal]`, examples `rtc-metal scene.json -o out.ppm`, `rtc-metal scene.json --window --samples 1024`, `rtc-metal scene.json --window --watch`).

- [ ] **Step 4: Run** `cargo test` — Expected: all pass. `cargo clippy --all-targets` clean.

- [ ] **Step 5: Commit**

```bash
git add src/help.rs src/lib.rs
git commit -m "ENHANCEMENT: rtc-edit's help page as data, built from the menu (milestone 28)"
```

### Task 3: windows.rs — what a key or a close does per window, and help::view

**Files:**
- Create: `src/windows.rs`
- Modify: `src/lib.rs` (`pub mod windows;`), `src/help.rs` (add `view`)

**Interfaces:**
- Consumes: `help::{Section, Group, Entry, Chip}` (Task 2).
- Produces:

```rust
// windows.rs
pub const HELP_SIZE: Size = Size::new(560.0, 720.0);
pub fn help_at(editor: Rectangle) -> Point;            // editor's top-right
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Input { Escape, Modifiers, Other }
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Route { Editor, CloseHelp, Drop }
pub fn route(window: Id, main: Id, help: Option<Id>, input: Input) -> Route;
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Closed { Quit, Help, Other }
pub fn closed(window: Id, main: Id, help: Option<Id>) -> Closed;
// help.rs
pub fn view<M: 'static>(sections: Vec<Section>) -> iced::Element<'static, M>;
```
(`Id` = `iced::window::Id`; `Size`, `Point`, `Rectangle` from `iced`.)

- [ ] **Step 1: Write the failing tests** in `src/windows.rs` `mod tests`:

```rust
#[test]
fn help_opens_at_the_editors_top_right() {
    let frame = Rectangle::new(Point::new(10.0, 20.0), Size::new(900.0, 640.0));
    assert_eq!(help_at(frame), Point::new(910.0, 20.0));
}

#[test]
fn keys_act_in_the_editor_and_esc_closes_help() {
    let (main, help, other) = (Id::unique(), Id::unique(), Id::unique());
    for input in [Input::Escape, Input::Modifiers, Input::Other] {
        assert_eq!(route(main, main, Some(help), input), Route::Editor);
        assert_eq!(route(main, main, None, input), Route::Editor);
    }
    assert_eq!(route(help, main, Some(help), Input::Escape), Route::CloseHelp);
    assert_eq!(route(help, main, Some(help), Input::Other), Route::Drop);
    assert_eq!(route(help, main, Some(help), Input::Modifiers), Route::Editor);
    assert_eq!(route(other, main, Some(help), Input::Escape), Route::Drop);
}

#[test]
fn closing_the_editor_quits_and_closing_help_forgets_it() {
    let (main, help, other) = (Id::unique(), Id::unique(), Id::unique());
    assert_eq!(closed(main, main, Some(help)), Closed::Quit);
    assert_eq!(closed(help, main, Some(help)), Closed::Help);
    assert_eq!(closed(other, main, Some(help)), Closed::Other);
    assert_eq!(closed(other, main, None), Closed::Other);
}
```

And in `src/help.rs` `mod tests`:

```rust
#[test]
fn the_page_draws() {
    let _: iced::Element<'static, ()> = view(sections());
}
```

- [ ] **Step 2: Run** `cargo test --lib windows:: help::` — Expected: FAIL to compile.

- [ ] **Step 3: Implement** `windows.rs` (module doc: which window an input or a close belongs to; the rules in Global Constraints) and `help::view`:
  - `view`: `scrollable` of a padded column (padding 16); per section a bold title (size 18), per group a title (size 14, skipped when empty), per entry a `row![chips..., text]` with spacing 6. `Modifier`/`Key`/`Mouse` chips are `container(text(s).size(12))` with padding `[2, 6]`, border radius 4 and a background from `theme.extended_palette()`: modifiers `secondary.weak`, keys `primary.weak`, mouse `success.weak` (text colour the same pair's `.text`). `Word` is plain text. Monospace sections use `text(...).font(iced::Font::MONOSPACE)`.

- [ ] **Step 4: Run** `cargo test` — Expected: all pass. `cargo clippy --all-targets` clean.

- [ ] **Step 5: Commit**

```bash
git add src/windows.rs src/help.rs src/lib.rs
git commit -m "ENHANCEMENT: rtc-edit's help page view and its window rules (milestone 28)"
```

### Task 4: main.rs — iced::daemon with an editor and a help window

**Files:**
- Modify: `src/main.rs` (`main`, `boot`, `Editor`, `Message`, `update`, `handle`'s `Window` arm, `shortcut`, `subscription`, `view`)

**Interfaces:**
- Consumes: `Message::Help` (Task 1); `help::{sections, view}`; `windows::{HELP_SIZE, help_at, route, Input, Route, closed, Closed}` (Tasks 2–3).
- Produces: nothing for later tasks.

- [ ] **Step 1: Switch to a daemon.**
  - `main`: `iced::daemon(boot, update, view).title(title).subscription(subscription).run()` (no `window_size`/`position`).
  - `boot() -> (Editor, Task<Message>)`: `let (main, open) = window::open(window::Settings { size: Size::new(WIDTH, HEIGHT), position: window::Position::Specific(Point::ORIGIN), ..Default::default() })`; `Editor.main: window::Id = main`, `Editor.help: Option<window::Id> = None`; return `open.discard()` batched with whatever boot already returns (it returns `Editor` today, so just `open.discard()`).
  - `title(e: &Editor, id: window::Id) -> String`: `"rtc-edit Help"` when `Some(id) == e.help`, else today's closure body.
  - `view(e, id)`: `help::view(help::sections())` when `Some(id) == e.help`, else today's `view`, renamed `editor_view`.

- [ ] **Step 2: Window events by id.**
  - `Message::Window(window::Id, window::Event)`; `subscription` maps `window::events()` to `Message::Window(id, event)`; `handle`'s exclusion list matches `Message::Window(..)`.
  - The arm: `Closed` → `match closed(id, e.main, e.help) { Quit => return iced::exit(), Help => e.help = None, Other => {} }`. `Opened`/`Moved`/`Resized` only when `id == e.main` (menus install there, as today).

- [ ] **Step 3: Route keys by window.**
  - `Message::Shortcut(window::Id, Box<Message>)`; `shortcut` keeps its body but takes the id and wraps every `Some(m)` as `Message::Shortcut(id, Box::new(m))`.
  - In `update`, before `handle`: unwrap `Shortcut(id, m)` with `route(id, e.main, e.help, input)`, where `input` is `Escape` for `Message::Escape(_)`, `Modifiers` for `Message::Command(_)`, else `Other`. `Editor` → handle `*m`; `CloseHelp` → `return window::close(id)`; `Drop` → `return Task::none()`. Unwrapping first keeps a dropped key from cancelling a pending `a`/`d`.

- [ ] **Step 4: Open or raise help.** `Message::Help` in `handle`: with `Some(id) = e.help` → `return window::gain_focus(id)`; else `let (id, open) = window::open(window::Settings { size: HELP_SIZE, position: window::Position::Specific(help_at(e.frame)), ..Default::default() })`, `e.help = Some(id)`, `return open.discard()`.

- [ ] **Step 5: Run** `cargo test` and `cargo clippy --all-targets` — Expected: all pass, clean. `cargo build` succeeds.

- [ ] **Step 6: Check by eye** (`cargo run -- examples/<a scene>.json`), each item:
  - Cmd-? opens "rtc-edit Help" right of the editor, 560×720; Cmd-? again raises it, no second window.
  - Help ▸ has macOS's search box; the item is never greyed.
  - Esc in help closes it; the red button closes it; Cmd-? then opens a new one.
  - With help in front: Cmd-S, Cmd-Z, Cmd-R act on the editor; letters typed in help do nothing in the editor; a pending `a` in the editor survives a key typed in help.
  - Moving help doesn't move where the viewer opens (Cmd-R after closing the viewer).
  - Closing the editor quits: help closes, the viewer is asked to quit, `$TMPDIR/rtc-edit-<pid>` is gone.
  - Chips read in light and dark (System Settings ▸ Appearance).

- [ ] **Step 7: Commit**

```bash
git add src/main.rs
git commit -m "ENHANCEMENT: rtc-edit opens a help window from Help ▸ rtc-edit Help (milestone 28)"
```

### Task 5: README

**Files:**
- Modify: `README.md` (the keyboard section's intro and the source-files table)

- [ ] **Step 1:** Add before the key table: "Help ▸ rtc-edit Help (Cmd-?) shows these keys, the menus, the viewer's mouse and keys and the command line in a window of its own; Esc closes it." Add rows `src/help.rs | the help window's page: sections, chips, view` and `src/windows.rs | which window a key or a close belongs to` to the files table; `src/menu.rs`'s row mentions labels and keys.
- [ ] **Step 2: Commit**

```bash
git add README.md
git commit -m "DOC: rtc-edit's README points at the help window (milestone 28)"
```
