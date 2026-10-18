# Milestone 19: the viewer asks to save its image, and the camera gets boxes

Branch `metal/19-quit-image-camera`, merged into develop with `--no-ff`
when done. The scene editor's milestones:

- 17: the program, a read-only tree, Save and Save As, and the viewer (done);
- 18: editing with fixed fields (done);
- 19: the viewer's image on quit, and the camera;
- 20: New and Open, adding and deleting items.

## Goal

Quitting the viewer, however it happens, asks whether to save the
rendered image and where. The camera gets boxes in the editor like the
shapes have.

Success:
- Quitting `rtc-metal scene.json --window` with Cmd-Q, by closing the
  window, or with `kill <pid>` shows "Save the rendered image?". Save…
  writes the image where you choose; Quit Without Saving writes nothing.
- Quitting `rtc-edit` shows the same alert in the viewer; before, the
  viewer was killed and its image lost.
- `-o out.ppm` still writes there without asking.
- Selecting the camera in the editor shows its boxes; an edit and Cmd-L
  move the viewer's view. Changing only width or height keeps the orbit.

## Out of scope

- New and Open, adding and deleting items (20).
- Remembering the last save folder in a config file, and a key that
  copies the viewer's orbited view into the scene's camera: later.
- Image formats other than `.ppm`.
- Asking before the editor quits with unsaved edits.
- Motion blur settings; they stay read-only in the JSON.

## Decisions

### rtc-metal: asking on quit

- **The question is asked in `applicationShouldTerminate:`**, which every
  way of quitting goes through: Cmd-Q, closing the last window, and
  SIGTERM. The image is no longer written in `applicationWillTerminate:`.
- **With `-o`**, the viewer quits at once and writes the image there.
- **Without `-o`**, the viewer activates itself and runs an `NSAlert`:
  message "Save the rendered image?", buttons **Save…** and **Quit
  Without Saving**.
  - Save… runs an `NSSavePanel` in the current directory, named
    `<scene stem>.ppm`, allowing `.ppm`. Save writes the image at the
    size last shown, as now, then quits. Cancel quits without saving.
  - Quit Without Saving quits.
- **A write that fails** shows a second alert with the error, then
  quits.
- **SIGTERM** sets an atomic flag in a signal handler (the only safe
  work there). A 0.25 s `NSTimer`, running with or without `--watch`,
  sees the flag and calls `NSApp.terminate`. A second SIGTERM while the
  alert is up does nothing more.
- Without `--window`, nothing changes.
- **Where the image goes is a value, not a callback detail:**
  `enum Save { To(String), Ask(String) }` — `To` is the `-o` path,
  `Ask` the default name. `run` takes a `Save`, and `on_close` gets the
  path chosen: `FnOnce(pixels, width, height, &Path)`.

### rtc-metal: keeping the orbit

- A reload resets the view only when the camera changed in something
  other than width and height. The image keeps the window's size, as
  now; width and height only set the window's first size.

### rtc-edit: stopping the viewer

- `Viewer`'s `Drop` sends SIGTERM instead of `Child::kill` (SIGKILL), and
  doesn't wait: the editor exits, and removes its preview folder, while
  the viewer waits for your answer. The image is already in the viewer's
  memory. The editor never reaps the child before sending, so its pid
  can't have been reused.

### rtc-edit: the camera's boxes

- **Boxes**, under headings:
  - *View*: `from`, `to`, `up`, three numbers each;
  - *Lens*: field of view in degrees (rounded to 6 decimals, stored in
    radians, replaced only when that box is edited), width, height;
  - *Antialias*: n;
  - *Focal blur*: aperture, focal distance, samples.
- **Validation:** width, height, n and samples are whole numbers of at
  least 1; field of view is above 0 and below 180°; numbers are finite.
  A zero `up`, or `from` equal to `to`, is left to the viewer's overlay.
- **Optional parts:**
  - n empty removes `antialias`; n typed creates it.
  - aperture empty removes `focal_blur`, and empties the other two
    boxes. Aperture typed with no focal
    blur creates it with focal distance = the distance from `from` to
    `to`, and samples = 32; their boxes then show those values.
  - focal distance and samples need a number while focal blur exists.
- rtc-metal ignores `antialias` (progressive rendering antialiases
  already); n matters to rtc-rust renders of the saved scene. The README
  says so.
- The arrow keys and h/j/k/l move over these boxes like any other.

## rtc-metal

- `src/main.rs`: builds `Save::To(out)` when `-o` was given, else
  `Save::Ask(<scene stem>.ppm)`, and writes the PPM to the path
  `on_close` receives. `fn image_name(scene: &Path) -> String` gives the
  default name.
- `src/window.rs`:
  - `Save`, and `run(…, save: Save, on_close)`.
  - `applicationShouldTerminate:` decides through
    `fn ask(save: &Save) -> Option<PathBuf>` (alert and panel for `Ask`,
    the path for `To`), resolves the image, calls `on_close` when there
    is a path, and answers `NSTerminateNow`.
  - The SIGTERM flag, its handler, and the timer.
  - `Renderer::reload` compares cameras with `same_view(a, b) -> bool`,
    which ignores width and height.

## rtc-edit

- `src/viewer.rs`: `Drop` sends SIGTERM (`libc`-free: an
  `extern "C" fn kill(pid, sig)` declaration, like the `atexit` one).
- `src/fields.rs`: `Field` gains `CameraFrom(usize)`, `CameraTo(usize)`,
  `CameraUp(usize)`, `FieldOfView`, `Width`, `Height`, `Antialias`,
  `Aperture`, `FocalDistance`, `FocalSamples`; `fields` and `set` handle
  `Node::Camera`.
- `src/main.rs`: `describe` gets the four headings.
- `README.md`: the camera's boxes and the alert on quit.

## Error handling

- **The image can't be written:** an alert with the error, then quit.
- **SIGTERM to a viewer that already quit:** nothing happens.
- **A camera that doesn't build:** the viewer's overlay shows why and
  keeps the old image (milestone 16).
- **Bad text in a camera box:** red, the scene unchanged (milestone 18).

## Testing

rtc-metal (unit tests):
- `image_name` gives `<stem>.ppm` for a scene path.
- `same_view` ignores width and height and notices `from`, field of
  view and focal blur.
- a reload changing only width keeps the orbit; one changing `from`
  resets it (the milestone 16 reload tests' style).

rtc-edit (unit tests):
- the camera lists its boxes with their starting text, and n and the
  focal blur boxes are empty when those parts are absent;
- field of view typed in degrees is stored in radians;
- 0 and 180 degrees, width 0, width 2.5 and samples -1 are rejected,
  the scene unchanged;
- typing n creates `antialias`, emptying it removes it;
- typing an aperture creates focal blur with the derived distance and
  32 samples; emptying it removes focal blur; emptying focal distance
  is rejected;
- dropping a `Viewer` stops a child that exits on SIGTERM (the existing
  test, now with SIGTERM).

By eye:
- `rtc-metal …/chapter10.json --window`: Cmd-Q, the close button and
  `kill <pid>` each show the alert; Save… writes where chosen; Quit
  Without Saving writes nothing; Cancel in the panel quits.
- With `-o x.ppm`, no alert and `x.ppm` is written.
- Quitting `rtc-edit` brings up the alert in the viewer.
- In the editor, change the camera's `from` and Cmd-L: the view moves.
  Change only width: the orbit stays.
- Type an aperture on a scene without focal blur, Cmd-L: the blur
  shows.

## Risks

- **A modal alert inside `applicationShouldTerminate:`** may misbehave
  in AppKit. Checked by eye first; the fallback is answering
  `NSTerminateLater`, showing the alert as a sheet, and replying with
  `replyToApplicationShouldTerminate:` when it closes.
- **The editor's terminal**: after the editor exits, the viewer keeps
  running until answered; a shell prompt may come back first.
