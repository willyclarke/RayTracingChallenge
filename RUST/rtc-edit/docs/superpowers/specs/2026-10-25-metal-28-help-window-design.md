# Milestone 28: a help window for the scene editor

Branch `metal/28-help-window`, merged into develop with `--no-ff` when done.
Milestones 17–21 and 23–27 built the scene editor and its menu bar.

## Goal

The editor's keys, the viewer's mouse and keys, and the command line are
only written down in the READMEs. A help window gives one short page that
shows them all, opened from the menu bar as other Mac apps do.

Success:
- Help ▸ rtc-edit Help opens a second window to the right of the editor
  (centred when it wouldn't fit); choosing it again brings it to the front.
  Cmd-? is macOS's own Help-menu search, which lists the item (see the
  2026-10-25 eye check below).
- The window lists, in one scrolling page: the menu commands with their
  keys, the editor's keys by place, the viewer's mouse, keys and menus, and
  the command line.
- The menu commands are built from the menu's own data, so they can't
  drift; a test checks that every command appears.
- Keys are drawn as coloured chips: modifiers, keys and mouse actions each
  have their own colour, readable in light and dark.
- Closing the editor quits the app, as today; closing help only closes help.

## Out of scope

- Search, clickable entries, greying commands that can't act.
- Remembering the help window's size or position.
- A help book (needs an app bundle).
- Building the viewer table from rtc-metal's code (another crate).
- Localisation.
- The README stays the long reference; the help page is a cheat sheet.

## Decisions

### A second iced window

- `iced::daemon(boot, update, view)` replaces `iced::application`. A
  daemon opens no window of its own and doesn't quit when a window closes.
- `boot` opens the editor window with `window::open` (`WIDTH`×`HEIGHT` at
  `Position::Specific(Point::ORIGIN)`, as today) and keeps its id in
  `Editor.main`; the help window's id, while open, is `Editor.help:
  Option<window::Id>`.
- `title(e, id)` and `view(e, id)` choose by window id: the editor as today,
  and `"rtc-edit Help"` / `help::view(...)` for help.
- `Message::Window(window::Id, window::Event)`: frame tracking (Opened,
  Moved, Resized) only for `main`; the menus install on main's Opened, as
  today.
- `window::Event::Closed`: for `main`, `iced::exit()` (`clean_up` still runs
  from `atexit`, so the preview folder and viewer go as today); for help,
  `Editor.help = None`. A Closed for any other id is ignored.
- Keys: `shortcut` already gets the window id. In the main window it works
  as today; in help, Esc closes help and every other key is ignored, so
  nothing typed in help reaches the tree or the boxes.
- App-wide menu commands (Save, Undo, Add, Load…) keep acting on the editor
  while help is in front.

### The Help menu

- New `Command::Help`, in a new Help menu after Viewer, labelled
  `rtc-edit Help`, no key of its own. The menu is
  `set_as_help_menu_for_nsapp`, so macOS adds its search box.
- Never greyed.
- With no help window: `window::open` at the editor's top-right
  (`frame.x + frame.width, frame.y`), about 560×720. With one open:
  `window::gain_focus(id)`.

### The menu's data (menu.rs)

Three pure functions taken out of `Menus::new`:
- `label(Command) -> &'static str`: "New", "Open…", …, "Add Sphere", …,
  "rtc-edit Help".
- `key(Command) -> Option<(Modifiers, Code)>`: today's accelerators, plus
  `None` for the Add items, Follow and Help.
- `menu_of(Command) -> &'static str`: File, Edit, Scene, Viewer or Help.

`Menus::new`'s `item` takes only the command and reads `label` and `key`, so
the menu bar and the help page use the same data. `Command::Help` joins
`ALL` (17 → 18). Ids and accelerators are otherwise unchanged.

### The help page's data (help.rs)

Pure, with no iced types in the data:

```rust
pub enum Chip {
    Modifier(&'static str), // ⌃ ⌥ ⇧ ⌘
    Key(&'static str),      // N, ⌫, ?, ↑, Esc, Tab, Space
    Mouse(&'static str),    // drag, scroll
    Word(&'static str),     // "then", "or": drawn without a box
}
pub struct Entry { pub chips: Vec<Chip>, pub text: &'static str }
pub struct Section { pub title: &'static str, pub entries: Vec<Entry> }

pub fn sections() -> Vec<Section>;
pub fn chips(modifiers: Modifiers, code: Code) -> Vec<Chip>;
```

- `chips` gives the modifiers in macOS order (⌃ ⌥ ⇧ ⌘) and then the key;
  Backspace is shown as `⌫`.
- A command with no key gets no chips; its keyboard route (`a` then `s`, …)
  is listed under Editor keys.

The four sections, in order:

1. **Menu commands**, grouped by `menu_of` in `ALL`'s order: each command's
   `label` and `chips(key)`.
2. **Editor keys**, by place, following the README's key table (rows in
   its order, its wording shortened):
   - Tree: ↑/↓ or `k`/`j` select; ←/`h`, →/`l` fold, unfold; Enter shades
     the first box; `a` then `s` `c` `p` `y` `o` `l` adds a sphere, cube,
     plane, cylinder, cone or light; `d` `d` deletes.
   - Shaded box: arrows or `h` `j` `k` `l` move the shading; Enter / `i`
     edit with the cursor at the end / start; `c` empties and edits; `x`
     empties; `u` puts the scene's value back; Esc back to the tree.
   - Editing: Enter or Esc leave the box (Esc on a red box puts the scene's
     value back).
   - Anywhere: Tab / ⇧Tab edit the next / previous box.
   - A line: the letters act only with no modifier held.
3. **Viewer** (hand-written; source: rtc-metal `window.rs`): drag orbits,
   scroll zooms, Space returns to the scene's view, `i` shows the overlay;
   menus: ⌘S Save Image…, ⌘1 / ⌘2 BVH / Metal traversal, ⌘] / ⌘[ double /
   halve samples, ⌘H hide.
4. **Command line**, monospace lines: `rtc-edit [scene.json]`;
   `RTC_METAL=/path/to/rtc-metal` picks the viewer; `rtc-metal scene.json
   [-o out.ppm] [--samples N] [--window [--watch]] [--traversal bvh|metal]`;
   a one-line example of each.

The README also gets a line pointing at Help ▸ rtc-edit Help.

### Drawing (help.rs)

`help::view(&[Section]) -> Element<'_, Message>`: a scrollable column with a
bold title per section and one row per entry (the chips, then the text).
Each chip is a small rounded container with a background colour per kind,
taken from the theme's palette so it reads in light and dark. Word chips
are plain text. The command-line section is monospace text with no chips.

## Tests

- menu.rs: `label` and `key` give what `Menus::new` used before;
  `key(Help)` is `None`; every command in `ALL` has a non-empty
  label; `menu_of(Help)` is "Help".
- help.rs: every command in `ALL` appears exactly once in the menu section,
  with its label and `chips(key)` when it has one; `chips` gives ⇧ ⌘ S for
  Cmd-Shift-S, ⌘ ⌫ for Cmd-Backspace; the
  sections are titled in order and none is empty.
- main.rs, as small plain functions where cheap: what a Closed does for
  main / help / another id; whether a key acts in a window (main yes, help
  only Esc); where help opens (right of the tracked frame).
- By eye: Help ▸ rtc-edit Help opens help beside the editor and raises it when open; Esc
  and the close button close it; menu keys act on the editor while help is
  in front; closing the editor quits, help included, and cleans up the
  preview folder and the viewer; chips read in light and dark.

## Errors

- `window::open` reports no failure we can act on; nothing to handle.
- Help open when the editor closes: `iced::exit()` takes it with the app.

## Eye check, 2026-10-25

Cmd-Shift-/ never reached the menu item: macOS's App Shortcut "Show Help
menu" (⇧⌘/) drops the Help menu open with its search box first. The user
chose to follow the Mac way: the item has no key, and Cmd-? then Enter
opens help. A review fix also centres help when it wouldn't fit beside the
editor on its screen (`windows::help_position`).
