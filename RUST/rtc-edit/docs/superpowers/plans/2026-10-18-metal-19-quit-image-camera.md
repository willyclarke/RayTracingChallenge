# rtc-edit, the viewer's image on quit and the camera — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** quitting the viewer, however it happens, asks whether to save the
image and where; the editor's camera gets boxes.

**Architecture:**
- **rtc-metal** decides where the image goes from a `Save` value (`-o`
  path, or ask with a default name). `applicationShouldTerminate:` runs
  the alert and save panel. SIGTERM sets a flag that a timer turns into
  a normal quit. A reload keeps the orbit unless the camera changed in
  more than its size.
- **rtc-edit** stops the viewer with SIGTERM, and `fields.rs` lists and
  sets the camera's boxes.
- rtc-rust doesn't change.

**Tech Stack:** Rust 2024; objc2 0.6 / objc2-app-kit 0.3 (`NSAlert`,
`NSSavePanel`); iced 0.14; rtc-rust's scene types.

**Spec:** `RUST/rtc-edit/docs/superpowers/specs/2026-10-18-metal-19-quit-image-camera-design.md`

## Global Constraints

- **Branch** `metal/19-quit-image-camera`. Never push. Commits end with
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Commit subjects** end with `(milestone 19)`; prefixes
  `ENHANCEMENT:` / `TEST:` / `DOCS:` / `FIX:`.
- **Test names:** descriptive snake_case.
- **No new crates.** rtc-metal may enable more `objc2-app-kit` features
  (`NSAlert`, `NSSavePanel`, `NSPanel`, and what they require). Signals
  and `kill` go through `extern "C"` declarations, like rtc-edit's
  `atexit`. SIGTERM is 15 on macOS.
- **Alert copy:** message `Save the rendered image?`; buttons `Save…`
  (first, default) and `Quit Without Saving`.
- **Save panel:** starts in the current directory, name
  `<scene stem>.ppm`, only `.ppm`.
- **Focal blur created from an aperture:** focal distance = |to − from|,
  samples = 32.
- **Field of view** shown in degrees rounded to 6 decimals, stored in
  radians; valid above 0 and below 180.
- **Whole-number boxes** (width, height, n, samples): at least 1.
- **Commands:**
  - rtc-metal: `cargo test --release --lib` (needs the GPU),
    `cargo clippy --all-targets`.
  - rtc-edit: `cargo test --lib`, `cargo clippy --all-targets`.
  - `cargo fmt` before committing.

## Review Focus

1. **Quitting again while the alert is up** (Cmd-Q, a second SIGTERM,
   or the editor quitting twice). Expected: one alert, no crash, no
   second write. Checked by eye in Task 3.
2. **Quitting before any frame was drawn.** Expected: no alert, nothing
   written, the viewer quits. Task 2's
   `nothing_is_saved_before_a_frame_is_drawn` pins it.
3. **The editor's preview folder disappearing while the viewer's alert
   is up.** Expected: at most an overlay error, the alert still works
   and saves. Checked by eye in Task 3.
4. **An aperture typed on a camera whose `from` equals `to`.**
   Expected: focal blur is created with focal distance 0, and the
   viewer's overlay reports anything that doesn't build. Task 4's
   `an_aperture_with_from_equal_to_to_gives_distance_zero` pins it.
5. **A width-only reload after the user orbited.** Expected: the orbit
   and the stopwatch restart as for any reload, but the view stays
   where the user left it. Task 1's
   `reload_with_only_a_new_size_keeps_the_view` pins it.

---

### Task 1: rtc-metal keeps the orbit when only the size changes

**Files:**
- Modify: `RUST/rtc-metal/src/window.rs` (`try_reload`, tests).

**Interfaces:**
- Produces: `fn same_view(a: &CameraDescription, b: &CameraDescription) -> bool`
  (private to `window.rs`): equal in everything but `width` and
  `height`.

- [ ] **Step 1: Write the failing tests** in `window.rs`'s tests:
  - `same_view_ignores_the_size`: a chapter7 camera vs. the same with
    `width = 999, height = 3` → true; with `from` changed, with
    `field_of_view` changed, and with `focal_blur` set → false each.
  - `reload_with_only_a_new_size_keeps_the_view`: as
    `reload_with_the_same_camera_keeps_the_view`, but changing
    `desc.camera.width = 77` instead of a material; after `reload`,
    `r.camera.transform()` approx-equals the dragged view, `r.moved` is
    still true, `(hsize, vsize) == (40, 20)`, and `r.camera_desc ==
    desc.camera`.

- [ ] **Step 2: Run them to verify they fail**

Run (from `RUST/rtc-metal`): `cargo test --release --lib same_view reload_with_only`
Expected: `same_view` not found; once it exists as a stub, the reload
test fails because the view resets.

- [ ] **Step 3: Implement `same_view`** (compare clones with width and
  height zeroed) and use it in `try_reload`'s `new_view` condition in
  place of `!=`.

- [ ] **Step 4: Run the tests**

Run: `cargo test --release --lib`
Expected: all pass.

- [ ] **Step 5: Commit**

`ENHANCEMENT: rtc-metal keeps the orbit when only the scene's image size changes (milestone 19)`

---

### Task 2: rtc-metal asks where to save the image on quit

**Files:**
- Modify: `RUST/rtc-metal/Cargo.toml` (app-kit features),
  `RUST/rtc-metal/src/window.rs`, `RUST/rtc-metal/src/main.rs`,
  `RUST/rtc-metal/README.md` (the `--window` paragraph on quitting).

**Interfaces:**
- Produces (in `window.rs`, public):
  - `pub enum Save { To(String), Ask(String) }`: `To` is the `-o`
    path, `Ask` the default name.
  - `pub fn image_name(scene: &Path) -> String`: `<stem>.ppm`.
  - `pub type OnClose = Box<dyn FnOnce(Vec<[f32; 4]>, usize, usize, &Path) -> Result<(), String>>`.
  - `run(…, save: Save, on_close: OnClose)`: `save` added before
    `on_close`.
- Private: `fn choose(save: &Save, frames: u32, ask: impl FnOnce(&str) -> Option<PathBuf>) -> Option<PathBuf>`:
  `None` when `frames == 0`; the `To` path without calling `ask`; for
  `Ask(name)`, whatever `ask(name)` returns.

- [ ] **Step 1: Write the failing tests** in `window.rs`'s tests:
  - `image_name_is_the_scene_stem_with_ppm`: `scenes/chapter7.json` →
    `chapter7.ppm`; `/tmp/rtc-edit-1/a.b.json` → `a.b.ppm`.
  - `a_given_path_is_used_without_asking`: `choose(&Save::To("x.ppm"),
    3, |_| panic!())` → `Some("x.ppm")`.
  - `asking_returns_the_answer`: `Ask("a.ppm")` with an `ask` that
    asserts it got `"a.ppm"` and returns `Some("/d/b.ppm")` → that path;
    one returning `None` → `None`.
  - `nothing_is_saved_before_a_frame_is_drawn`: frames 0 with `To` and
    with `Ask` → `None`, and `ask` not called.

- [ ] **Step 2: Run them to verify they fail**

Run: `cargo test --release --lib image_name a_given asking nothing_is`
Expected: compile errors, `Save` / `image_name` / `choose` not found.

- [ ] **Step 3: Implement the pure parts:** `Save`, `image_name`,
  `choose`. `main.rs` uses `image_name` for its default `-o` too.

- [ ] **Step 4: Run the tests**

Run: `cargo test --release --lib`
Expected: all pass.

- [ ] **Step 5: Wire the quit.**
  - `Delegate` keeps `save: Save`.
  - `applicationShouldTerminate:` (replacing the work in
    `applicationWillTerminate:`): a `Cell<bool>` guard makes a second
    call while asking answer `NSTerminateCancel`. Otherwise
    `choose(&save, acc.frames(), ask)`; with a path, resolve the image
    and call `on_close`; on `Err(e)` run an alert showing `e`. Answer
    `NSTerminateNow`.
  - `ask(name)`: activate the app (`activateIgnoringOtherApps(true)`),
    run the alert (Global Constraints' copy); on `Save…` run an
    `NSSavePanel` with `directoryURL` = the current directory,
    `nameFieldStringValue` = `name`, allowed type `ppm`, and return its
    URL's path on OK; `None` otherwise.
  - SIGTERM: a `static QUIT: AtomicBool`; an `extern "C"` handler that
    only stores `true`; installed in `run` with `signal(15, handler)`.
  - The timer: always scheduled (not only with `--watch`), interval
    `WATCH_INTERVAL`; its method checks `QUIT` (swap to false, then
    `NSApplication::terminate`), then the watched file if any.
  - `main.rs`: `Save::To(o)` when `-o` was given, else
    `Save::Ask(image_name(path))`; `on_close` writes the PPM to the
    path and returns the write error as `Err` instead of exiting;
    without `--window`, `-o`'s default stays `image_name` and the
    write stays as now.
  - README: `--window` asks on quit; `-o` writes without asking.

- [ ] **Step 6: Build and lint**

Run: `cargo build --release && cargo test --release --lib && cargo clippy --all-targets`
Expected: builds; all pass; clippy clean.

- [ ] **Step 7: Commit** (the by-eye check is with Task 3)

`ENHANCEMENT: rtc-metal asks whether to save the image when the window quits (milestone 19)`

---

### Task 3: rtc-edit asks the viewer to quit

**Files:**
- Modify: `RUST/rtc-edit/src/viewer.rs`, `RUST/rtc-edit/README.md`.

**Interfaces:**
- Consumes: Task 2's SIGTERM handling in the viewer.

- [ ] **Step 1: Replace the failing test.** `dropping_the_viewer_kills_it`
  becomes `dropping_the_viewer_asks_it_to_quit`: the script is
  `trap 'echo term > <root>/got; exit 0' TERM; sleep 30 & wait`; after
  `drop(v)`, `<root>/got` appears within 2 s. (Without waiting, a
  SIGTERMed child stays a zombie, so `kill -0` can't tell.)

- [ ] **Step 2: Run it to verify it fails**

Run (from `RUST/rtc-edit`): `cargo test --lib dropping_the_viewer`
Expected: FAIL — SIGKILL can't be trapped, `got` never appears.

- [ ] **Step 3: Implement.** `Drop` calls `kill(self.child.id() as i32, 15)`
  through an `extern "C" fn kill(pid: c_int, sig: c_int) -> c_int`, and
  doesn't wait.

- [ ] **Step 4: Run the tests**

Run: `cargo test --lib && cargo clippy --all-targets`
Expected: all pass; clippy clean.

- [ ] **Step 5: README.** Replace "Quitting the viewer itself writes its
  image … into the folder you started `rtc-edit` from" with: quitting
  the editor or the viewer asks in the viewer whether to save the
  image; Save… starts in the folder you started `rtc-edit` from.

- [ ] **Step 6: Check by eye** (with the user; Tasks 2 and 3):
  1. `rtc-metal scenes/chapter10.json --window` (from `RUST/rtc-metal`):
     Cmd-Q shows the alert; Save… → the panel, in the current folder,
     named `chapter10.ppm`; Save writes it.
  2. Again with the close button and Quit Without Saving: nothing
     written. Again with `kill <pid>` and Cancel in the panel: nothing
     written, the viewer quits.
  3. With `-o /tmp/x.ppm`: no alert, `/tmp/x.ppm` written.
  4. Quit while the alert is up (Cmd-Q in the alert's app menu, or a
     second `kill`): still one alert.
  5. `rtc-edit ../rtc-metal/scenes/chapter10.json`, then Cmd-Q in the
     editor: the viewer comes to the front with the alert; saving
     works although the editor's preview folder is gone.

- [ ] **Step 7: Commit** after the user's approval:

`ENHANCEMENT: rtc-edit asks the viewer to quit, so it offers to save its image (milestone 19)`

---

### Task 4: the camera's boxes

**Files:**
- Modify: `RUST/rtc-edit/src/fields.rs`, `RUST/rtc-edit/src/main.rs`
  (`describe`), `RUST/rtc-edit/README.md`.

**Interfaces:**
- Consumes: milestone 18's `fields`, `set`, `next_box`.
- Produces: `Field` variants `CameraFrom(usize)`, `CameraTo(usize)`,
  `CameraUp(usize)`, `FieldOfView`, `Width`, `Height`, `Antialias`,
  `Aperture`, `FocalDistance`, `FocalSamples`, listed in that order for
  `Node::Camera`.

- [ ] **Step 1: Write the failing tests** in `fields.rs`, on the test
  scene's camera (`width 10, height 5, field_of_view 1.0, from [0,0,-5],
  to [0,0,0], up [0,1,0]`):
  - `the_camera_lists_its_boxes`: `CameraFrom(2)` = `"-5"`, `Width` =
    `"10"`, `FieldOfView` = `"57.29578"` (1 rad, rounded to 6 decimals
    as milestone 18's rotation text), `Antialias`, `Aperture`,
    `FocalDistance`, `FocalSamples` = `""`.
  - `field_of_view_is_typed_in_degrees`: set `"90"` → stored
    `FRAC_PI_2`.
  - `camera_values_out_of_range_are_rejected`: `FieldOfView` `"0"` and
    `"180"`, `Width` `"0"` and `"2.5"`, `Height` `"-1"`, `Antialias`
    `"0"` → `Err`, scene equal to before.
  - `typing_n_adds_antialias_and_emptying_it_removes_it`.
  - `an_aperture_creates_focal_blur`: set `Aperture` `"0.25"` →
    `focal_blur == Some(FocalBlur { aperture: 0.25, focal_distance:
    5.0, samples: 32 })`, and `fields` shows `"5"` and `"32"`.
  - `emptying_the_aperture_removes_focal_blur`: after the above, `""`
    → `None`, and the three focal boxes show `""`.
  - `focal_distance_needs_a_number_while_focal_blur_exists`: `""` →
    `Err`.
  - `an_aperture_with_from_equal_to_to_gives_distance_zero`.
  - `items_you_cant_edit_list_nothing` drops the camera from its list.

- [ ] **Step 2: Run them to verify they fail**

Run: `cargo test --lib camera aperture antialias field_of_view focal`
Expected: compile errors, the variants don't exist.

- [ ] **Step 3: Implement** the variants, `fields` and `set` for
  `Node::Camera`. Whole numbers parse as `usize` from the trimmed text
  and must be ≥ 1. `main.rs`'s `describe` gives the headings *View*
  (`from`, `to`, `up`), *Lens* (`field of view°`, `width`, `height`),
  *Antialias* (`n`), *Focal blur* (`aperture`, `focal distance`,
  `samples`).

- [ ] **Step 4: Run the tests**

Run: `cargo test --lib && cargo clippy --all-targets`
Expected: all pass; clippy clean.

- [ ] **Step 5: README.** The camera in the Editing section's list of
  items with boxes; n only matters to rtc-rust renders, as rtc-metal's
  progressive rendering antialiases anyway; an empty aperture removes
  focal blur.

- [ ] **Step 6: Check by eye** (with the user), in
  `rtc-edit ../rtc-metal/scenes/chapter10.json`:
  1. Select the camera: the four headings and their boxes.
  2. Change `from`, Cmd-L: the view moves.
  3. Orbit in the viewer, change only width, Cmd-L: the orbit stays.
  4. Type an aperture of 0.3, Cmd-L: the blur shows; empty it, Cmd-L:
     sharp again.
  5. h/j/k/l move over the camera's boxes.

- [ ] **Step 7: Commit** after the user's approval:

`ENHANCEMENT: rtc-edit edits the camera in fixed fields (milestone 19)`

### Task 4 amendment: f-number, red-box reasons, reset

Spec: the "Amendment after the Task 4 check by eye" section. Replaces
the aperture parts of Task 4; the rest of Task 4 stands.

**Interfaces:**
- `Field::Aperture` becomes `Field::FNumber`.
- `pub fn lens_radius(desc: &SceneDescription) -> Option<f64>` in
  `fields.rs`: the camera's `aperture` when focal blur is set.
- `main.rs`: `errors` becomes `HashMap<Field, String>` (the reason).

- [ ] **Step 1: Write the failing tests** in `fields.rs`, replacing the
  aperture tests (test camera: from `[0,0,-5]`, to `[0,0,0]`):
  - `an_f_number_creates_focal_blur`: set `FNumber` `"8"` →
    `FocalBlur { aperture: 0.3125, focal_distance: 5.0, samples: 32 }`
    (5 / 16); `fields` shows `"8"`, `"5"`, `"32"`;
    `lens_radius == Some(0.3125)`.
  - `emptying_the_f_number_removes_focal_blur`: then `""` → `None`, the
    three boxes `""`, `lens_radius == None`.
  - `changing_the_focal_distance_keeps_the_f_number`: with f-number 8,
    set `FocalDistance` `"10"` → aperture `0.625`, `FNumber` still
    shows `"8"`.
  - `the_f_number_must_be_above_zero`: `"0"` and `"-2"` → `Err`, scene
    unchanged.
  - `an_aperture_of_zero_shows_an_infinite_f_number`: a camera with
    `focal_blur {aperture: 0, focal_distance: 5, samples: 8}` shows
    `FNumber` `"∞"`.
  - `an_f_number_with_from_equal_to_to_gives_distance_zero` (replaces
    the aperture version): focal distance 0, aperture 0.

- [ ] **Step 2: Run them to verify they fail**

Run: `cargo test --lib -- f_number focal lens`
Expected: compile errors, `FNumber` / `lens_radius` not found.

- [ ] **Step 3: Implement** `FNumber` and `lens_radius` in `fields.rs`.
  In `main.rs`:
  - `describe`: *Focal blur* row labels `f-number`, `focal distance`,
    `samples`; below the focal blur boxes a `text` line
    `lens radius <r>` from `lens_radius` when it is `Some`.
  - `errors: HashMap<Field, String>`; `edit` stores the reason. The
    status line shows `<label>: <reason>` (the row label from
    `describe`, with the axis as `x`/`y`/`z` for vectors) while the
    shaded or edited box is in `errors`; otherwise the usual status.
  - `fn reset(&mut self, field)`: the text from `fields` for that
    field, remove it from `errors`. `Letter('u')` on a shaded box calls
    it; `Escape(true)` resets first when the box it leaves is in
    `errors`. `shortcut` adds `"u"` to the letters.
  - README: the f-number and lens radius; `u` and Esc-on-red in the
    keys table; a red box's reason shows in the status line.

- [ ] **Step 4: Run the tests**

Run: `cargo test --lib && cargo clippy --all-targets && cargo build --release`
Expected: all pass; clippy clean; builds.

- [ ] **Step 5: Check by eye** (with the user), on chapter10:
  1. Type an f-number of 4, Cmd-L: blur; the lens radius line shows;
     16 is sharper than 4.
  2. Type `abc` in width, shade it: the status line says why; `u`
     restores the value. Type `abc` again, Esc: restored.
  3. Set focal distance with f-number kept: the lens radius changes.

- [ ] **Step 6: Commit** Task 4 and this amendment together, after the
  user's approval, with Task 4's subject.

### Task 5: the viewer resizes its window to the scene's width and height

Spec: the "Amendment: the viewer follows width and height" section.

**Files:**
- Modify: `RUST/rtc-metal/src/window.rs`, `RUST/rtc-metal/README.md`.

**Interfaces:**
- `Renderer` gains `resize: Option<(usize, usize)>`: set by
  `try_reload` to the file's `(width, height)` when they differ from
  the previous `camera_desc`'s; taken (and cleared) by the view after a
  reload.
- `fn window_size(width: usize, height: usize, scale: f64, room: Option<NSSize>) -> NSSize`:
  the content size `run` computes today (pixels / backing scale, shrunk
  uniformly to fit `room`, never grown), extracted so `run` and the
  resize share it.

- [ ] **Step 1: Write the failing tests** in `window.rs`'s tests:
  - `reload_with_a_new_size_asks_for_that_window_size`: chapter7, set
    `desc.camera.width = 77`, reload → `r.resize == Some((77, height))`
    with the file's height; the orbit kept as in Task 1's test.
  - `reload_with_the_same_size_asks_for_no_resize`: a material change
    → `r.resize == None`.
  - `window_size_fits_the_room_and_never_grows`: `(800, 600)` at scale
    2 with no room → `400 × 300`; with room `200 × 300` → `200 × 150`;
    with room `4000 × 4000` → `400 × 300`.

- [ ] **Step 2: Run them to verify they fail**

Run (from `RUST/rtc-metal`): `cargo test --release --lib -- resize window_size`
Expected: compile errors, `resize` / `window_size` not found.

- [ ] **Step 3: Implement.** `try_reload` sets `resize`; `run` uses
  `window_size`; `reload_if_changed`, after a good reload, takes
  `resize` and calls `setContentSize` on the view's window with
  `window_size(w, h, backing scale, room from the window's screen)`.
  The drawable then follows through `drawableSizeWillChange`, as for a
  hand resize. README: width and height in the file resize the window
  on reload.

- [ ] **Step 4: Run the tests**

Run: `cargo test --release --lib && cargo clippy --all-targets && cargo build --release`
Expected: all pass; clippy clean; builds.

- [ ] **Step 5: Check by eye** (with the user): in rtc-edit on
  chapter10, orbit, change width to half and Cmd-L: the window gets
  narrower, the orbit stays; a size larger than the screen is fitted.

- [ ] **Step 6: Commit** after the user's approval:

`ENHANCEMENT: rtc-metal resizes its window when the scene's width or height changes (milestone 19)`
