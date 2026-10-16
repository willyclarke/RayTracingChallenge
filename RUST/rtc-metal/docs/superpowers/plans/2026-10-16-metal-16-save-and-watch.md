# Saving scenes and `--window --watch` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** rtc-rust saves scene descriptions back to JSON, and rtc-metal's
`--window --watch` reloads the scene file whenever it changes.

**Architecture:**
- **rtc-rust:** the serde description types in `rtc-rust/src/scene.rs`
  gain `Serialize`/`PartialEq`, plus `to_json` and an atomic `save`.
- **rtc-metal:** a small `Watch` (mtime polling) decides when the file
  changed; `Renderer::reload` swaps in the new scene. An `NSTimer` on the
  view polls every 0.25 s, so a paused window still reloads.

**Tech Stack:** Rust 2024, serde/serde_json, objc2 0.6 /
objc2-foundation 0.3 (`NSTimer`, in its default features), Metal.

**Spec:** `RUST/rtc-metal/docs/superpowers/specs/2026-10-16-metal-16-save-and-watch-design.md`

## Global Constraints

- **Branch** `metal/16-save-and-watch`. Never push. Commits end with
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Commit subjects:**
  - rtc-rust commits: `ENHANCEMENT: rtc-rust …`;
  - rtc-metal commits end with `(milestone 16)`;
  - prefixes `ENHANCEMENT:` / `TEST:` / `DOCS:` / `FIX:`.
- **Test names:**
  - rtc-rust loader tests continue `test_scene_N` from 10;
  - rtc-metal tests use descriptive snake_case names, like `window.rs`'s.
- **Existing behaviour:** all existing tests pass, and every existing
  render is unchanged.
- **Reload restarts** accumulation from one sample, and the stopwatch
  from zero.
- **Poll interval:** 0.25 s.
- **Status on a bad reload:** `reload failed: <message>`.
- **Commands:**
  - rtc-rust (from `RUST/rtc-rust`): `cargo test --lib`,
    `cargo clippy --lib`, `cargo test --doc`;
  - rtc-metal (from `RUST/rtc-metal`): `cargo test --release --lib`,
    `cargo clippy --release`.
- **Formatting:** run `cargo fmt` in each crate touched before committing.
- **YAGNI:** no asset-path rewriting (milestone 17), and no file-system
  event APIs.

## Review Focus

1. **A watched file rewritten twice within the same mtime tick.**
   Expected: the second change is still seen. APFS has nanosecond mtimes;
   `Watch` compares exact `SystemTime`s. Task 3's test sets the mtimes
   explicitly.
2. **A reload while the window is paused at its target.** Expected: it
   redraws and refines again. Task 5 resumes drawing after a successful
   reload; checked by eye.
3. **Changing only a material.** Expected: the orbit/zoom is kept,
   because the camera description is equal. Task 4's "same camera" test
   covers it.
4. **A reload after the window was resized.** Expected: the new scene
   renders at the window's current size, not the file's
   `width`/`height`. Task 4's tests assert the size.
5. **The scene file deleted and recreated** (some editors save that way).
   Expected: no error while it is missing, and a reload once it is back.
   Task 3's test covers the missing file.

---

### Task 1: rtc-rust saves scenes

**Files:**
- Modify: `RUST/rtc-rust/src/scene.rs` (derives, attributes,
  `to_json`/`save`, module doc, tests)
- Modify: `RUST/rtc-rust/README.md` (lines ~32–33 and ~114, "load-only")
- Modify: `RUST/rtc-rust/CLAUDE.md` (`scene` paragraph)

**Interfaces:**
- Produces:
  - `#[derive(Serialize, PartialEq, Clone)]` on every description type
    (`SceneDescription`, `CameraDescription`, `Antialias`, `FocalBlur`,
    `MotionBlur`, `LightDescription`, `Settings`, `PathTracing`,
    `TransformOp`, `ShapeDescription`, `ShapeKind`, `CsgOp`,
    `MaterialDescription`, `BumpDescription`, `PatternDescription`,
    `PatternKind`, `MappingKind`, `UvDescription`), keeping existing
    derives;
  - `SceneDescription::to_json(&self) -> String`;
  - `SceneDescription::save<P: AsRef<Path>>(&self, path: P) -> Result<(), SceneError>`.
  Task 4 relies on `CameraDescription: Clone + PartialEq`.

- [ ] **Step 1: Write the failing tests** in `scene.rs`'s test module:
  - `test_scene_10`: for every `scenes/*.json` (read with
    `std::fs::read_dir("scenes")`), `from_json(&text)? ==
    from_json(&desc.to_json())?`. Fail with the file name. Assert that at
    least 5 files were checked.
  - `test_scene_11`: this input
    ```json
    {"camera":{"width":10,"height":5,"field_of_view":1.0,"from":[0,0,-5],"to":[0,0,0],"up":[0,1,0]},
     "light":{"type":"point","position":[-10,10,-10],"intensity":[1,1,1]},
     "shapes":[{"type":"sphere"}]}
    ```
    serializes so that `serde_json::from_str::<serde_json::Value>(&desc.to_json())?`
    equals a `serde_json::json!` of the same scene, written the way serde
    writes it: `usize` fields as integers (`"width": 10`), `f64` fields
    as floats (`"from": [0.0, 0.0, -5.0]`). Also assert the text contains
    none of `"transform"`, `"settings"`, `"closed"`, `null`.
  - `test_scene_12`: a scene with `"comment": {"why": "test"}` and
    `"lights": [two point lights]` round-trips into an equal description.
    Its `to_json()` contains `"lights"`, not `"light":`, and the comment
    object.
  - `test_scene_13`: `save` to `std::env::temp_dir().join("test_scene_13.json")`
    the description from `scenes/cover.json`. Then
    `load(&path)` builds a world with the same shape count, light count
    and camera `hsize`/`vsize` as `load("scenes/cover.json")`.
    `SceneDescription::from_path(&path)? == desc`. Assert no
    `test_scene_13.json.tmp` exists afterwards.

- [ ] **Step 2: Run them to see them fail**

Run: `cargo test --lib test_scene_1` (from `RUST/rtc-rust`)
Expected: compile errors — no `to_json`, no `save`, no `PartialEq` on
`SceneDescription`.

- [ ] **Step 3: Implement**
  - Add `Serialize, PartialEq, Clone` to every description type.
  - Add `#[serde(skip_serializing_if = "Option::is_none")]` to every
    `Option` field (including `comment`, `light`, `lights`, the pattern
    and bump seeds, cylinder/cone `min`/`max`, torus radii, and every
    `MaterialDescription` field);
  - `skip_serializing_if = "Vec::is_empty"` on `ShapeDescription.transform`;
  - a `fn is_false(b: &bool) -> bool` for `closed`;
  - `fn is_default_settings(s: &Settings) -> bool` (`*s == Settings::default()`)
    for `SceneDescription.settings`.
  - `to_json`: `serde_json::to_string_pretty(self)`; serializing these
    types cannot fail, so `expect` with a message.
  - `save`: write to `<path>.tmp` (append `.tmp` to the file name), then
    `std::fs::rename` over `path`. Map I/O errors to
    `SceneError(format!("writing {}: {e}", path.display()))`.
  - Doc comments as in the spec's code block. Module doc line 1 becomes
    `//! JSON scene description — load and save.`

- [ ] **Step 4: Run the tests, then the whole suite**

Run: `cargo test --lib test_scene && cargo test --lib && cargo test --doc && cargo clippy --lib`
Expected:
- all `test_scene_*` pass, including 10–13;
- the full suite passes;
- the doctests pass;
- clippy is clean.

If `test_scene_10` fails on a flattened pattern or shape, the
`flatten` + internally tagged enum is serializing differently from what it
reads. Fix the serde attributes, not the test.

- [ ] **Step 5: Docs.** Update the README's two "load-only" passages and
  CLAUDE.md's `scene` paragraph:
  - drop "Load-only" / "never writes scenes back";
  - say `to_json`/`save` write a description back as pretty JSON that loads
    into an equal description;
  - defaults are left out, and `save` replaces the file atomically.

- [ ] **Step 6: Commit**

```bash
cargo fmt
git add src/scene.rs README.md CLAUDE.md
git commit -m "ENHANCEMENT: rtc-rust saves scene descriptions as JSON" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 2: rtc-metal's scene files round-trip

**Files:**
- Modify: `RUST/rtc-metal/src/scene.rs` (test module only)

**Interfaces:**
- Consumes: Task 1's `to_json`, `PartialEq`.

- [ ] **Step 1: Write the test.** `every_scene_file_round_trips`: for
  every `scenes/*.json` in rtc-metal, `from_json(&text) ==
  from_json(&desc.to_json())`, failing with the file name. Assert at
  least 10 files were checked.

- [ ] **Step 2: Run it**

Run: `cargo test --release --lib every_scene_file_round_trips` (from `RUST/rtc-metal`)
Expected: PASS. Task 1 already provides the behaviour; this test pins
rtc-metal's own scene files.

- [ ] **Step 3: Commit**

```bash
git add src/scene.rs
git commit -m "TEST: rtc-metal scene files round-trip through to_json (milestone 16)" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 3: `Watch` — has the file changed?

**Files:**
- Create: `RUST/rtc-metal/src/watch.rs`
- Modify: `RUST/rtc-metal/src/lib.rs` (`pub mod watch;`)

**Interfaces:**
- Produces:
  - `pub struct Watch`;
  - `pub fn new(path: impl Into<PathBuf>) -> Watch`;
  - `pub fn path(&self) -> &Path`;
  - `pub fn changed(&mut self) -> bool`.

- [ ] **Step 1: Write the failing tests** in `watch.rs`. Use a file in
  `std::env::temp_dir()` with a test-specific name, and set mtimes with
  `std::fs::File::options().write(true).open(..)?.set_modified(t)`.
  - `changed_reports_each_new_modification_once`:
    - write the file and make the `Watch`; `changed()` is false;
    - set the mtime to `t0 + 1 s`; `changed()` is true, then false;
    - set it to `t0 + 1 s + 1 ns`; true again.
  - `a_missing_file_is_unchanged_until_it_reappears`:
    - make a `Watch`, then delete the file; `changed()` is false, twice;
    - recreate it with a new mtime; `changed()` is true.
  - `watching_a_file_that_does_not_exist_yet`: `Watch::new` on a missing
    path; `changed()` is false. Create the file; `changed()` is true.

- [ ] **Step 2: Run them to see them fail**

Run: `cargo test --release --lib watch`
Expected: compile error — no `watch` module.

- [ ] **Step 3: Implement.** `Watch { path: PathBuf, modified:
  Option<SystemTime> }`.
  - `new` records `fs::metadata(path).and_then(|m| m.modified()).ok()`.
  - `changed`: read the same way. A missing file returns false and
    leaves `modified` as is, so a file that comes back with the old mtime
    is not a change. A different `Some(t)` stores it and returns true.

- [ ] **Step 4: Run the tests**

Run: `cargo test --release --lib watch && cargo clippy --release`
Expected: 3 pass; clippy clean.

- [ ] **Step 5: Commit**

```bash
cargo fmt
git add src/watch.rs src/lib.rs
git commit -m "ENHANCEMENT: rtc-metal Watch notices a changed scene file (milestone 16)" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 4: `Renderer::reload`

**Files:**
- Modify: `RUST/rtc-metal/src/window.rs`
- Modify: `RUST/rtc-metal/src/orbit.rs` (only if the helper below lands there)
- Modify: `RUST/rtc-metal/src/main.rs` (use the helper)

**Interfaces:**
- Consumes: `CameraDescription: Clone + PartialEq` (Task 1).
- Produces, all in `window.rs`:
  - `Renderer.camera_desc: CameraDescription`, the description the
    current scene was loaded with;
  - `Renderer.reload_error: Option<String>`;
  - `fn reload(&mut self, desc: &SceneDescription, base_dir: &Path) -> Result<(), String>`,
    which also sets `reload_error` (`None` on success; `Some(msg)` on
    error, and returns `Err(msg)`);
  - `pub fn orbit_of(c: &CameraDescription) -> Orbit`, the orbit
    `main.rs::open_window` builds today, moved into `orbit.rs` as
    `Orbit::of(c)` or kept as a free function. Pick one; `main.rs` uses it.
  - `run` gains a `camera_desc: CameraDescription` parameter (Task 5
    adds the watch).

- [ ] **Step 1: Write the failing tests** in `window.rs`'s tests. Add a
  helper that builds a `Renderer` from a `SceneDescription` parsed from
  `scenes/chapter7.json`, at 40×20 via `resized`. Each test adds 2 frames
  first, then asserts:
  - `reload_with_the_same_camera_keeps_the_view`:
    - move the orbit (`r.orbit.drag(0.3, 0.0)`, then `r.restart()`), and
      remember `r.camera.transform()`;
    - change a sphere's material colour in a cloned description, then
      `reload`;
    - `Ok`; `acc.frames() == 0`; `elapsed == ZERO`; `camera` size
      still 40×20; `camera.transform()` unchanged; `moved` still true.
  - `reload_with_a_new_camera_takes_it_as_home`:
    - set `from` to `[0, 3, -8]` in the clone, then `reload`;
    - `Ok`; size 40×20; `camera.transform()` equals the transform built
      from the new from/to/up; `home.0` equals it too; `moved` false;
      `camera_desc` is the new one.
  - `a_broken_reload_keeps_the_old_scene`:
    - set `lights = Some(vec![])` and `light = None` in the clone (which
      `build` rejects);
    - `reload` is `Err`; `reload_error` is `Some` containing `"lights"`;
      `acc.frames() == 2`; `camera_desc` unchanged;
    - a following good reload clears `reload_error`.
  - `reload_keeps_the_traversal`: `set_traversal(TRAVERSAL_METAL)`, then
    a good reload; `scene.info.counts.traversal == TRAVERSAL_METAL`, and
    `acc.add_frame` succeeds.
  - `an_unsupported_reload_keeps_the_old_scene`: replace the light with a
    `LightDescription::Spot`; `Err` containing `"spot"`; frames still 2.

- [ ] **Step 2: Run them to see them fail**

Run: `cargo test --release --lib reload`
Expected: compile errors — no `reload`, `camera_desc` or `reload_error`.

- [ ] **Step 3: Implement `reload`.** Order matters: everything fallible
  happens before any field changes.
  1. `desc.build(base_dir)` → `(world, file_camera)`.
  2. The camera to use:
     - `desc.camera == self.camera_desc` → `self.camera`;
     - otherwise `resized(&file_camera, self.camera.hsize, self.camera.vsize)`,
       with the new `orbit_of(&desc.camera)`.
  3. `GpuScene::from_world(&world, &camera)`, set its traversal to the
     current one, then `self.gpu.accumulator(&scene)`.
  4. Only now assign:
     - `scene`, `acc` and `camera_desc`;
     - for a new camera, also `camera`, `orbit`, `home` and `moved = false`;
     - `elapsed = ZERO`, `last_frame = now`, `reload_error = None`.
  Errors become `Strings` (`SceneError.0`) and set `reload_error`.
  - `show_status`: when `reload_error` is `Some(e)`, append
    ` · reload failed: {e}` to the status text.
  - `main.rs::open_window` passes `desc.camera.clone()` to `run` and uses
    `orbit_of`.

- [ ] **Step 4: Run the tests and the suite**

Run: `cargo test --release --lib && cargo clippy --release`
Expected: all pass, including the 5 new tests; clippy clean.

- [ ] **Step 5: Commit**

```bash
cargo fmt
git add src/window.rs src/main.rs src/orbit.rs
git commit -m "ENHANCEMENT: rtc-metal --window can reload its scene (milestone 16)" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 5: `--watch` — the timer and the flag

**Files:**
- Modify: `RUST/rtc-metal/src/window.rs` (`run` takes the watch; timer;
  `checkWatchedFile:`; module doc)
- Modify: `RUST/rtc-metal/src/main.rs` (flag, usage, help, module doc)
- Modify: `RUST/rtc-metal/README.md` (`--window` bullets)

**Interfaces:**
- Consumes:
  - `Watch` (Task 3);
  - `Renderer::reload`, `reload_error` (Task 4);
  - the existing `RtcView::resume()`.
- Produces: `run(..., camera_desc: CameraDescription, watch: Option<Watch>, ...)`.

- [ ] **Step 1: Wire it up**
  - **Where the watch lives:** `Renderer` gets `watch: Option<Watch>`.
  - **Timer:** in `run`, when the watch is `Some`, schedule
    `NSTimer::scheduledTimerWithTimeInterval_target_selector_userInfo_repeats(0.25, &view, sel!(checkWatchedFile:), None, true)`,
    with a `// SAFETY:` comment. The view outlives the timer, which
    repeats until the app exits.
  - **`RtcView` method** `#[unsafe(method(checkWatchedFile:))] fn check_watched_file(&self, _timer: &NSTimer)`:
    - borrow the renderer; if the watch reports `changed()`, read
      `SceneDescription::from_path` and call `reload` with the file's
      parent directory;
    - on success, `resume()`;
    - on error, `eprintln!("rtc-metal: reload failed: {e}")`, set
      `reload_error` (a parse error, too), and `setNeedsDisplay(true)`
      so a paused window shows the message;
    - don't hold the `RefCell` borrow across `resume()`.
  - **`main.rs`:**
    - parse `--watch` into a `bool`;
    - `--watch` without `--window` → `usage()`;
    - `open_window` passes `Some(Watch::new(path))` when set.
  - **Docs:**
    - usage line: `[--window [--watch]]`;
    - help line: `--watch      with --window: reload the scene whenever its file changes`;
    - the `main.rs` and `window.rs` module docs;
    - the README `--window` bullets.

- [ ] **Step 2: Build, test, lint**

Run: `cargo build --release && cargo test --release --lib && cargo clippy --release && ./target/release/rtc-metal scenes/chapter7.json --watch; echo "exit $?"`
Expected:
- all tests pass and clippy is clean;
- the last command prints the usage and `exit 2`.

- [ ] **Step 3: Launch check.** Copy `scenes/chapter7.json` to the
  scratchpad (the base dir doesn't matter, since it has no assets). Run
  `--window --watch` on the copy in the background for 6 s, logging to the
  scratchpad. Rewrite the copy's first sphere colour with a script after
  2 s.
  Expected: an empty log (no `reload failed`).

- [ ] **Step 4: By eye (the user).** On a copy of `scenes/chapter7.json`,
  in a text editor:
  - move a sphere and save, while refining and while paused;
  - orbit first, then change only a colour: the view stays;
  - break the JSON: the overlay shows `reload failed: …`;
  - fix it: the message goes away.

- [ ] **Step 5: Commit** once the user confirms.

```bash
cargo fmt
git add src/window.rs src/main.rs README.md
git commit -m "ENHANCEMENT: rtc-metal --window --watch reloads a changed scene (milestone 16)" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```
