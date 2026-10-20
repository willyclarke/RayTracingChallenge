# Milestone 21: the editor and the viewer side by side — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The editor opens at the screen's top-left and every viewer it starts opens just right of it, fitted to the rest of the screen.

**Architecture:** rtc-metal learns `--at X,Y` (top-left of its window, top-left-origin points) and places itself with pure, tested geometry in `window.rs`; a size reload keeps its top-left. rtc-edit opens at (0, 0), tracks its frame from iced window events, and passes its top-right corner to each viewer it starts.

**Tech Stack:** Rust; rtc-metal on objc2 AppKit (`NSWindow`, `NSScreen`); rtc-edit on iced 0.14.

**Spec:** `RUST/rtc-edit/docs/superpowers/specs/2026-10-19-metal-21-window-layout-design.md`

## Global Constraints

- Branch `metal/21-window-layout`; commits end ` (milestone 21)` with prefixes `ENHANCEMENT:` / `TEST:` / `DOCS:` / `FIX:`; never push.
- `--at X,Y`: the window frame's top-left, in points, right and down from the top-left of the main screen; `y_appkit = main_screen_height − Y`.
- Less than **200** points of width right of X, or no screen holding the point: centred, as today.
- The viewer's top never goes above its screen's visible top (under the menu bar).
- A size reload keeps the window's top-left corner; room is measured from its left edge and top to the visible right and bottom edges.
- The image is never grown past the camera's width and height (`window_size`, unchanged).
- The editor's size stays 900×640; it opens at `Position::Specific(0, 0)`.
- `--at` without `--window` is a usage error, as `--watch` is.
- rtc-metal tests: `cargo test --release --lib`; rtc-edit: `cargo test --lib`; both `cargo clippy --all-targets` clean and `cargo fmt`. The "future-incompat … block v0.1.6" warning is pre-existing.
- Test names: descriptive snake_case.

## Review Focus

1. **The editor reports top 0 before macOS moves it below the menu bar.** Expected: the viewer's title bar is below the menu bar. Task 1: `place_moves_a_point_under_the_menu_bar_down`.
2. **The editor on a second screen left of or above the main one** (negative X, or Y outside the main screen). Expected: the viewer lands on that screen. Task 1: `place_finds_the_screen_holding_the_point` with a screen at negative X.
3. **The editor pushed to the right edge.** Expected: viewer centred, not a sliver. Task 1: `place_gives_none_without_200_points_of_room`.
4. **A size reload after the user dragged the viewer elsewhere.** Expected: it stays where the user put it, top-left fixed. Task 1: `room_from_measures_to_the_visible_right_and_bottom` + by eye.
5. **Dragging the editor while `a` is pending.** Expected: the pending key survives (window events are not keys). Task 2: the update preamble exempts `Message::Window`; by eye.

---

### Task 1: rtc-metal places its window with `--at`

**Files:**
- Modify: `RUST/rtc-metal/src/window.rs` (geometry fns, `run`, the size reload ~line 588), `RUST/rtc-metal/src/main.rs` (flag, `usage`, `open_window`), `RUST/rtc-metal/README.md`

**Interfaces:**
- Produces: `rtc-metal <scene> --window [--watch] --at X,Y`, X and Y any finite numbers (negative allowed), written by Rust's `{}` for `f64`.
- In `window.rs` (pub, so the lib tests reach them):
  - `pub fn parse_at(s: &str) -> Option<(f64, f64)>`
  - `pub fn place(at: (f64, f64), main_height: f64, screens: &[(NSRect, NSRect)]) -> Option<(NSPoint, NSSize)>` — `screens` are `(frame, visibleFrame)` pairs; returns the window frame's top-left in AppKit coordinates and the frame room.
  - `pub fn room_from(top_left: NSPoint, visible: NSRect) -> NSSize` — from a top-left (AppKit) to the visible right and bottom edges.
  - `run(…, at: Option<(f64, f64)>, …)` gains the parameter.

- [ ] **Step 1: Write the failing tests** in `window.rs`'s test module. A main screen `frame (0,0,1440,900)`, `visible (0,0,1440,875)` (25-point menu bar):
  - `parse_at_reads_two_numbers`: `"900,0"` → `Some((900.0, 0.0))`; `"-1440,12.5"` → `Some((-1440.0, 12.5))`.
  - `parse_at_rejects_bad_values`: `"5"`, `"a,b"`, `"1,2,3"`, `"inf,0"`, `"NaN,0"`, `""` → `None`.
  - `place_puts_the_window_right_of_the_point`: `at (900, 30)` → top-left `(900, 870)`, room `(540, 870)`.
  - `place_moves_a_point_under_the_menu_bar_down`: `at (900, 0)` → top-left `(900, 875)`, room `(540, 875)`.
  - `place_finds_the_screen_holding_the_point`: add a second screen `frame (-1280,0,1280,1024)`, `visible` the same; `at (-380, 0)` (its Y measured from the main screen's top, so AppKit y = 900) → top-left `(-380, 900)`, room `(380, 900)`.
  - `place_gives_none_without_200_points_of_room`: `at (1300, 30)` → `None`; `at (1240, 30)` (exactly 200) → `Some`.
  - `place_gives_none_off_every_screen`: `at (5000, 30)` → `None`.
  - `room_from_measures_to_the_visible_right_and_bottom`: top-left `(100, 800)`, visible `(0,0,1440,875)` → `(1340, 800)`.
- [ ] **Step 2: Run** `cargo test --release --lib` in `RUST/rtc-metal`. Expected: FAIL (functions missing).
- [ ] **Step 3: Implement** `parse_at`, `room_from`, `place` (the holding screen is the first whose `frame` contains the converted point, right and top edges inclusive; clamp the y to the visible top; room via `room_from`; `None` below 200 width).
- [ ] **Step 4: Run** the tests. Expected: PASS.
- [ ] **Step 5: Wire it.**
  - `main.rs`: `--at` takes the next argument through `parse_at`, else `usage()`; `--at` without `--window` → `usage()`; the usage text lists `--at X,Y`; pass it through `open_window` to `run`.
  - `run`: with `Some(at)`, collect `NSScreen::screens` as `(frame, visibleFrame)`, take the main screen's height from `NSScreen::mainScreen`'s frame; on `Some((top_left, room))` size the window with `window_size` against `contentRectForFrameRect_styleMask` of a rect of size `room`, and call `setFrameTopLeftPoint(top_left)` instead of `center()`. Otherwise `center()` as now.
  - The size reload: remember the frame's top-left (`x`, `y + height`) before `setContentSize`; take the room from `room_from(top_left, screen.visibleFrame())` turned into content room as today; after `setContentSize`, `setFrameTopLeftPoint(top_left)`.
- [ ] **Step 6: README** (`RUST/rtc-metal/README.md`): `--at X,Y` with the coordinate convention, the 200-point fallback, and that a size reload keeps the top-left.
- [ ] **Step 7: Run** `cargo test --release --lib`, `cargo clippy --all-targets`, `cargo fmt`, `cargo build --release`. Expected: pass, clean.
- [ ] **Step 8: Stop for the by-eye check** (report; do not commit). The user checks from `RUST/rtc-metal`: `cargo run --release -- scenes/chapter10.json --window --at 900,0` opens right of x = 900 under the menu bar; `--at 1400,0` opens centred; without `--at` centred; `--at 5` prints usage; with `--watch`, editing the scene's camera width keeps the title bar in place.
- [ ] **Step 9: Commit when told** `ENHANCEMENT: rtc-metal --at places the window, and a size reload keeps its top-left (milestone 21)`.

### Task 2: rtc-edit opens at the left and starts viewers beside itself

**Files:**
- Modify: `RUST/rtc-edit/src/viewer.rs`, `RUST/rtc-edit/src/main.rs`, `RUST/rtc-edit/README.md`

**Interfaces:**
- Consumes: rtc-metal's `--at X,Y` (Task 1).
- Produces: `pub fn viewer_args(preview: &Path, at: (f64, f64)) -> Vec<OsString>` in `viewer.rs`; `Viewer::start(exe: &Path, preview: &Path, at: (f64, f64))`.

- [ ] **Step 1: Write the failing test** in `viewer.rs`: `viewer_args_put_the_window_at_the_point`: `viewer_args(Path::new("/p/s.json"), (900.0, 0.0))` == `["/p/s.json", "--window", "--watch", "--at", "900,0"]`; `(-1440.0, 12.5)` gives `"-1440,12.5"`.
- [ ] **Step 2: Run** `cargo test --lib viewer`. Expected: FAIL.
- [ ] **Step 3: Implement** `viewer_args` and use it in `Viewer::start`, which gains `at`; update the existing viewer tests' calls (any point, e.g. `(0.0, 0.0)`).
- [ ] **Step 4: Editor frame.**
  - `main()`: `.position(iced::window::Position::Specific(iced::Point::ORIGIN))` beside the existing `.window_size((900.0, 640.0))`.
  - `Editor.frame: iced::Rectangle`, starting `(0, 0, 900, 640)`.
  - `subscription` batches the existing `listen_with(shortcut)` with `iced::window::events().map(|(_, e)| Message::Window(e))`.
  - `Message::Window(iced::window::Event)`: `Opened { position, size }` sets the size and, when `Some`, the position; `Moved(p)` the position; `Resized(s)` the size; other events nothing.
  - The update preamble that clears `pending` exempts `Message::Window` too, like `Command` and `Focused`.
  - `show_viewer` passes `(frame.x + frame.width, frame.y)`.
- [ ] **Step 5: README**: the editor opens at the top-left, each viewer beside it (where the editor is now), centred when there's no room.
- [ ] **Step 6: Run** `cargo test --lib`, `cargo clippy --all-targets`, `cargo fmt`, `cargo build --release` (and `cargo build --release` in `RUST/rtc-metal`, which the editor runs). Expected: pass, clean.
- [ ] **Step 7: Stop for the by-eye check** (report; do not commit). The user checks: `cargo run --release -- ../rtc-metal/scenes/chapter10.json` — editor top-left, viewer beside it, tops aligned; move the editor, Cmd-O another scene → the new viewer beside the editor's new spot; close the viewer, Cmd-R → beside; `a` then drag the editor, then `s` still adds a sphere; a second screen if available.
- [ ] **Step 8: Commit when told** `ENHANCEMENT: rtc-edit opens at the left and starts each viewer beside itself (milestone 21)`.
