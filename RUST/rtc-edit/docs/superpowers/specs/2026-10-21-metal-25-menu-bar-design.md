# Milestone 25: a native menu bar for the scene editor

Branch `metal/25-menu-bar`, merged into develop with `--no-ff` when done.
Milestones 17–21, 23 and 24 built the scene editor; 26 is planned as a
help window built on the same commands.

## Goal

The editor's commands are only keys, which have to be known. A native
macOS menu bar makes them discoverable, shows each key beside its command,
and greys out what can't act, as other Mac apps do.

Success:
- The menu bar shows rtc-edit, File, Edit, Scene and Viewer; each item works
  by click and by its key.
- Items that can't act are greyed, and Follow the Viewer shows a tick while
  following.
- Copy, paste, cut and select all still work in the boxes.

## The spike (2026-10-21, throwaway)

An iced 0.14 app with a `muda` 0.21 menu bar, checked by eye:
- the menu installs (`Menu::init_for_nsapp`) when the window opens, and an
  iced subscription polling `MenuEvent::receiver()` gets its clicks and
  keys;
- a menu item's key arrives twice: as a menu event and as an iced key event;
- `muda`'s predefined Cut, Copy, Paste and Select All take Cmd-X/C/V/A from
  iced's text boxes, which then do nothing; without them, the boxes work;
- the predefined Quit quits.

## Out of scope

- A Help menu: milestone 26, with the help window.
- Menus in the viewer (rtc-metal has its own).
- The two-letter keys (`a s`, `d d`, hjkl, …): they stay keys.

## Decisions

### The menus

| Menu | Items (key) |
|---|---|
| rtc-edit | About rtc-edit, Quit (Cmd-Q) |
| File | New (Cmd-N), Open… (Cmd-O), Save (Cmd-S), Save As… (Cmd-Shift-S) |
| Edit | Undo (Cmd-Z), Redo (Cmd-Shift-Z), Delete (Cmd-Delete) |
| Scene | Add Sphere, Add Cube, Add Plane, Add Cylinder, Add Cone, Add Light |
| Viewer | Load Edits (Cmd-L), Start Viewer (Cmd-R), Follow the Viewer (tick), Keep Followed View (Cmd-K) |

- No predefined Cut, Copy, Paste or Select All (see the spike).
- Each item sends a command the editor already has: `Save`, `Undo`,
  `Add(Item::Sphere)`, `Follow(bool)` (the opposite of the tick), `Keep`, …

### One path per command

- The menu owns every Cmd key: `shortcut` loses its Cmd table (S, Shift-S,
  N, O, R, L, K, Z, Shift-Z, Cmd-Delete). Plain keys, Tab, Esc and the
  arrows stay in `shortcut`.
- Cmd-Delete is Cmd-Backspace (the Mac's Delete key); forward-delete with
  Cmd no longer deletes.

### Greying

| Item | Enabled when |
|---|---|
| Undo / Redo | there is a step to undo / redo |
| Delete | an item is selected and `edit::delete` would accept it |
| Load Edits, Follow the Viewer | a viewer is running |
| Start Viewer | no viewer is running |
| Keep Followed View | following |
| everything else | always |

- Follow the Viewer is a check item, ticked while following.
- The menus are refreshed after every message; only items whose state
  changed are touched.

## Code

- `Cargo.toml`: `muda = "0.21"`.
- `src/menu.rs`:
  - pure: `pub enum Command { New, Open, Save, SaveAs, Undo, Redo, Delete,
    Add(Item), Load, StartViewer, Follow, Keep }`; `fn id(Command) ->
    &'static str` and `pub fn command(id: &str) -> Option<Command>`; `pub
    struct State { can_undo, can_redo, can_delete, viewer_running,
    following: bool }`; `pub fn enabled(Command, &State) -> bool`;
  - `muda`: `pub struct Menus` (the bar and the items that change);
    `Menus::new() -> Menus`; `install(&self)` (once, when the window
    opens); `refresh(&mut self, &State)`.
- `src/history.rs`: `can_undo()`, `can_redo()`.
- `src/edit.rs`: `pub fn can_delete(desc: &SceneDescription, node: &Node)
  -> bool` (a delete on a copy).
- `src/main.rs`: `Editor.menus: Option<Menus>` (built and installed on the
  window's `Opened`); a 50 ms subscription turns menu events into
  `Message::Menu(Command)`, handled as the existing messages; after every
  message `update` builds the `State` and refreshes; the Cmd table goes.
- `README.md`: the menus, and the key table without the Cmd rows that moved
  (their keys are listed with the menus).

## Error handling

- A menu id the editor doesn't know (none expected) is ignored.
- Building the menus can't fail at run time with fixed items; `muda`'s
  `Result`s are unwrapped with a message naming the menu.

## Testing

`menu.rs` (unit tests, no window):
- every `Command`'s id maps back to it;
- `enabled` follows the greying table;
`edit.rs`: `can_delete` is true for a sphere and false for a CSG's child,
the last light, the camera and a heading.
`history.rs`: `can_undo` / `can_redo` before and after a step and an undo.

By eye: each item by click and by key, once each (no double action);
greying as the state changes; the Follow tick; copy and paste in a box;
About and Quit.
