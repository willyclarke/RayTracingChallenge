# Milestone 21: the editor and the viewer side by side

Branch `metal/21-window-layout`, merged into develop with `--no-ff` when
done. The scene editor's milestones:

- 17: the program, a read-only tree, Save and Save As, and the viewer (done);
- 18: editing with fixed fields (done);
- 19: the viewer's image on quit, and the camera (done);
- 20: New and Open, adding and deleting items, comments (done);
- 21: the editor and the viewer side by side.

## Goal

Today both windows open centred, so the viewer covers the editor. The
editor should open at the left of the screen and every viewer it starts
should open just right of it, so both can be seen without moving
anything.

Success:
- `rtc-edit ../rtc-metal/scenes/chapter10.json`: the editor at the
  screen's top-left, the viewer beside it, tops aligned, nothing
  overlapping.
- Move the editor, then Cmd-O another scene: the new viewer opens beside
  the editor where it is now.
- A width change in the camera boxes and Cmd-L: the viewer resizes and
  its title bar stays put.
- `rtc-metal scene.json --window` on its own still opens centred.

## Out of scope

- Moving an open viewer when the editor moves; only new viewers follow.
- Remembering window positions between runs (a config file is a later
  idea).
- Changing the editor's size (900×640) or letting the viewer grow past
  the camera's width and height.
- Putting the viewer left of the editor when there is no room on the
  right.

## Decisions

### Where the windows go

- **The editor** opens at the top-left of the main screen,
  `Position::Specific(0, 0)`; macOS moves it below the menu bar.
- **A viewer the editor starts** — at start-up, on New, on Open, and on
  Cmd-R — opens with its window's top-left corner at the editor's
  current top-right corner. Its image is shrunk uniformly to fit the
  room between there and the right and bottom edges of that screen's
  visible area, and never grown past the camera's width and height (as
  now, `window_size`).
- **The editor tracks its frame**: its position and size as it is moved
  or resized, so each new viewer goes beside the editor where it is, on
  the screen the editor is on. An open viewer stays where it is.
- **No room:** when less than 200 points of width remain right of the
  editor, or no screen holds the point, the viewer opens centred, as
  today.
- **The title bar stays below the menu bar:** the editor's frame may
  still say top 0 before macOS has moved it, so the viewer moves its top
  down into the screen's visible area.

### Size reloads

- A reload that changes the camera's width or height (milestone 19)
  resizes the viewer keeping its window's **top-left** corner; today
  AppKit keeps the bottom-left, so the title bar jumps.
- The room for that resize is measured from the window's current left
  edge to the right edge of its screen's visible area, and from its top
  to the bottom, so a viewer beside the editor stays beside it.

### How the position gets across

- **`rtc-metal --at X,Y`**: the top-left corner of the window's frame,
  in points, measured right and down from the top-left of the main
  screen (the coordinates iced and winit use). Only meaningful with
  `--window`; without `--at` the window is centred.
- rtc-metal converts to AppKit's coordinates, which measure up from the
  bottom of the main screen: `y_appkit = main_screen_height − Y`.
- A malformed value (`--at 5`, `--at a,b`, a non-finite number) is an
  error at start-up, like any other bad flag.

## rtc-metal

- `src/main.rs`: parses `--at X,Y` into `Option<(f64, f64)>` and passes
  it to `window::run`.
- `src/window.rs`:
  - `fn place(at: (f64, f64), main_height: f64, visible: &[NSRect])
    -> Option<(NSPoint, NSSize)>`: from the point and the screens'
    visible frames, the window's top-left in AppKit coordinates (moved
    down to the visible top if above it) and the room for the window's
    frame (to the visible right and bottom edges); `None` when no screen
    holds the point or less than 200 points of width remain. The caller
    turns the frame room into content room with
    `contentRectForFrameRect:styleMask:`. Pure, so it is tested without
    a window.
  - At start-up, with `--at` and a `Some` from `place`, the window is
    sized with `window_size` against that room and placed with
    `setFrameTopLeftPoint`; otherwise `center()` as now.
  - The size reload keeps the top-left corner and measures the room from
    the window's left edge and top.

## rtc-edit

- `src/main.rs`:
  - the window opens at `Position::Specific(0, 0)`;
  - `Editor.frame`: the editor's outer top-left and inner size, starting
    at (0, 0) and 900×640, updated from `iced::window::events()`
    (`Opened`, `Moved`, `Resized`);
  - `show_viewer` passes the frame's top-right corner to the viewer.
- `src/viewer.rs`: `Viewer::start` takes `at: (f64, f64)` and adds
  `--at X,Y`; a small pure helper builds the arguments.
- `README.md`: the layout, and `--at` in rtc-metal's README.

## Error handling

- **No room, or no screen at the point:** the viewer opens centred.
- **A bad `--at`:** rtc-metal exits with the error, as for other flags;
  the editor only ever passes numbers it formatted itself.
- **The editor's frame is stale** (moved and no event yet): the viewer
  lands where the editor was; harmless.

## Testing

rtc-metal (unit tests):
- `place` on one screen: the top-left and room for a point right of the
  editor;
- a point on a second screen lands on that screen;
- a point with less than 200 points of room, and a point on no screen,
  give `None`;
- a point above the visible top (under the menu bar) is moved down;
- parsing `--at 900,0`, and rejecting `--at 5` and `--at a,b`;
- a size reload keeps the window's top-left (the function computing the
  new frame, tested without a window).

rtc-edit (unit tests):
- the viewer's arguments include `--at` with the frame's top-right
  corner.

By eye:
- start-up layout with `chapter10.json`;
- move the editor, Cmd-O: the new viewer beside it; Cmd-R likewise;
- a second screen, if there is one;
- a width change keeps the viewer's title bar in place;
- `rtc-metal …/chapter10.json --window` alone opens centred.

## Risks

- **winit's reported position under the menu bar**: macOS moves a window
  placed at y = 0 down, and a `Moved` event may or may not follow; the
  viewer's clamp covers either case.
- **Full screen**: a viewer in full screen on a size reload; the
  milestone 19 deferred minor already notes `setContentSize` there. Not
  made worse; checked by eye if convenient.
