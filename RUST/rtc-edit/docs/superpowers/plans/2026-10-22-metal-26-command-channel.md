# Milestone 26: a command channel, and one viewer through New and Open — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The shared file carries load commands as well as the camera, each side sleeps in `os_sync_wait_on_address` until the other writes, and one viewer lives through New and Open without taking the focus.

**Architecture:** `rtc-share` replaces `flock` with per-block seqlocks and futex-style waits. rtc-metal gains a listener thread that hands each `Load` to the main queue (`dispatch2`), a fresh/non-fresh load on `Renderer`, Hide items and File ▸ Save Image…. rtc-edit keeps one share file for its life, sends loads instead of restarting the viewer, and follows the camera through a `Subscription::run_with` stream.

**Tech Stack:** Rust 2024; macOS 14.4+ `os_sync_wait_on_address` / `os_sync_wake_by_address_all` (libSystem); objc2 0.6, objc2-app-kit 0.3, `dispatch2 = "0.3"`; iced 0.14.

**Spec:** `RUST/rtc-edit/docs/superpowers/specs/2026-10-22-metal-26-command-channel-design.md`

## Global Constraints

- Branch `metal/26-command-channel`; commits end ` (milestone 26)`, prefixes `ENHANCEMENT:` / `TEST:` / `DOCS:` / `FIX:`; never push; never add `RUST/rtc-edit/untitled.json`.
- Each crate is its own Cargo project: run `cargo test`, `cargo clippy --all-targets`, `cargo fmt` in the crate's folder (`rtc-share`, `rtc-metal`, `rtc-edit`). The "future-incompat … block v0.1.6" warning is pre-existing.
- `rtc-share` stays one file (`src/lib.rs`); `SIGNATURE` still hashes it.
- Share file: `$TMPDIR/rtc-edit-<pid>.share`. Load path limit: 1024 bytes.
- Viewer menus: rtc-metal (Hide rtc-metal Cmd-H, Hide Others Cmd-Opt-H, Show All, separator, Quit Cmd-Q); File (Save Image… Cmd-S); Traversal and Samples unchanged.
- The viewer's 100 ms timer, the editor's Tick (viewer death only) and the 50 ms menu poll stay (milestone 27).
- Test names: descriptive snake_case.

## Rulings made in this plan

- **`wait_camera` returns `Result<Option<Snapshot>, String>`** (spec: `Option<Snapshot>`): the spec's error handling needs the wait's failure to reach the editor's status line.
- **`camera_seq` and `generation` share one 8-byte word**, and `wait_camera` sleeps on all 8 bytes: a bump of either wakes it, so `end_viewer` can't lose a wake between the waiter's check and its sleep.
- **Each side opens its own mapping by path** for its waiting thread (the viewer's listener, the editor's stream), rather than sharing one `Shared` across threads.
- **The quit writer becomes reusable** (`OnClose` → `Writer`, `Fn` not `FnOnce`), so Save Image… and quit use the same `write_ppm`.
- **`Snapshot.seq` becomes `u32`.**

## Review Focus

1. **A load arriving while the quit alert or a save panel is up:** the main queue also drains in modal run loops, so `deliver` can run then; it must apply cleanly, never panic on a `RefCell` borrow. Tasks 2–3: no `Renderer` borrow is held across `runModal` (Save Image… and quit); by eye with Save Image…'s panel open, then Cmd-L in the editor.
2. **The editor quits while a viewer runs:** the viewer gets SIGTERM and asks as today; the share file is removed; no thread keeps the editor alive. Task 4: `clean_up` removes `share_path`; by eye.
3. **A viewer closed by the user, then Cmd-R:** the old camera-stream thread ends, the new one follows the new viewer, and the generation in the subscription key changes. Task 4: `end_viewer` on death; Task 1's `wait_camera_ends_when_the_viewer_ends` test.
4. **New/Open to an unsaved scene (`untitled.json`) twice in a row:** the preview path is the same; `fresh` still resets the orbit. Task 3: `a_fresh_load_with_the_same_camera_resets_the_view`.
5. **A reader racing a writer:** never a torn camera or path. Task 1: `a_reader_never_sees_a_half_written_camera`.

---

### Task 1: `rtc-share` — seqlocks and waits

**Files:**
- Modify: `RUST/rtc-share/src/lib.rs`, `RUST/rtc-share/README.md`.
- Modify (callers, rename only): `RUST/rtc-metal/src/window.rs` (`publish`: `write` → `write_camera`), `RUST/rtc-edit/src/main.rs` (`tick`: `Shared::read` → `Shared::read_camera`), `RUST/rtc-edit/src/follow.rs` (test helper `shot(seq: u32)`).

**Interfaces:**
- Produces:
  - Layout `#[repr(C)] struct Layout { signature: u64, size: u64, camera_seq: u32, generation: u32, from: [f64; 3], to: [f64; 3], up: [f64; 3], focal_distance: f64, load_seq: u32, fresh: u32, len: u32, path: [u8; 1024] }` (`camera_seq` and `generation` at an 8-aligned offset, together).
  - `Shared::create(path: &Path) -> Result<Shared, String>`, `Shared::open(path: &Path) -> Result<Shared, String>` (as today, no `flock`).
  - `write_camera(&self, from: [f64; 3], to: [f64; 3], up: [f64; 3], focal_distance: f64)`; `read_camera(&self) -> Option<Snapshot>` (`None` while `camera_seq == 0`).
  - `wait_camera(&self, last: u32, generation: u32) -> Result<Option<Snapshot>, String>`: `Ok(None)` once `generation()` differs from `generation`; `Ok(Some(s))` with `s.seq != last`.
  - `send_load(&self, path: &Path, fresh: bool) -> Result<(), String>`; `load_seq(&self) -> u32`; `wait_load(&self, last: u32) -> Result<Load, String>`.
  - `#[derive(Clone, Debug, PartialEq)] pub struct Load { pub seq: u32, pub path: PathBuf, pub fresh: bool }`.
  - `generation(&self) -> u32`; `end_viewer(&self)` (bump `generation`, wake the camera word).
  - `Snapshot { seq: u32, from, to, up, focal_distance }`.
  - `unsafe impl Send` stays; add `unsafe impl Sync` (every access is atomic or seqlocked).

- [ ] **Step 1: Write the failing tests** (rename the existing ones' `write`/`read` to `write_camera`/`read_camera`; keep them):
  - `a_sent_load_reads_back_with_its_path_and_fresh`: `a = create`, `b = open`; `a.send_load("/x/y.json", true)`; `b.wait_load(0)? == Load { seq: 2, path: "/x/y.json".into(), fresh: true }`; then `send_load("/z.json", false)`; `b.wait_load(2)?.fresh == false`, path `/z.json`, `b.load_seq() == 4`.
  - `wait_load_sleeps_until_a_send`: `b.wait_load(b.load_seq())` on a thread; main thread sleeps 50 ms, asserts the thread hasn't finished (`JoinHandle::is_finished`), sends, joins, gets the path.
  - `wait_camera_returns_a_new_write`: `b.wait_camera(0, b.generation())` on a thread; `a.write_camera(...)`; join → `Ok(Some(s))` with `s.seq == 2` and the values.
  - `wait_camera_ends_when_the_viewer_ends`: thread waits as above with no write; after 50 ms `a.end_viewer()`; join → `Ok(None)`; also `b.wait_camera(0, g_before)` after the bump returns `Ok(None)` at once.
  - `a_reader_never_sees_a_half_written_camera`: writer thread writes `[k; 3]` for `from`, `to`, `up` and `focal_distance = k` for `k` in `1..=20_000`; reader loops `read_camera` until `focal_distance == 20_000.`, asserting every snapshot has all ten values equal.
  - `a_path_over_1024_bytes_is_refused`: `send_load` of a 1025-byte path is `Err` containing `"1024"`; `load_seq()` stays 0; a 1024-byte path is `Ok`.
- [ ] **Step 2: Run** `cargo test` in `rtc-share`. Expected: FAIL to compile (missing items).
- [ ] **Step 3: Implement.**
  - FFI (libSystem): `os_sync_wait_on_address(addr: *mut c_void, value: u64, size: usize, flags: u32) -> i32`, `os_sync_wake_by_address_all(addr: *mut c_void, size: usize, flags: u32) -> i32`; `OS_SYNC_WAIT_ON_ADDRESS_SHARED = 1`, `OS_SYNC_WAKE_BY_ADDRESS_SHARED = 1`. A wait returning -1 with `EINTR` loops; any other errno is `Err`. A wake's -1 (`ENOENT`, no waiters) is ignored.
  - Words are reached as `&AtomicU32` / `&AtomicU64` on the mapping; data copied with `read_volatile` / `write_volatile`.
  - Seqlock write: `s = seq.load(Relaxed)`; `seq.store(s.wrapping_add(1), Relaxed)`; `fence(Release)`; write data; `seq.store(s.wrapping_add(2), Release)`; wake (`camera`: the 8-byte word, size 8; `load`: `load_seq`, size 4).
  - Seqlock read: loop { `s1 = seq.load(Acquire)`; odd → spin; copy; `fence(Acquire)`; `seq.load(Relaxed) == s1` → done }.
  - `wait_camera`: loop { `w = word.load(Acquire)` (u64 of both halves); generation half ≠ `generation` → `Ok(None)`; camera half even and ≠ `last` → `Ok(read_camera())` with the read's seq; else `os_sync_wait_on_address(word, w, 8, SHARED)` }.
  - `wait_load`: same on `load_seq` (size 4), returning the seqlocked copy.
  - `create` writes the header with every word 0; `open` checks `signature`/`size` as today.
  - Module doc: replace the `flock` paragraph with the seqlock/wait description and the macOS 14.4 floor. README: the same in a few lines.
- [ ] **Step 4: Run** `cargo test`, `cargo clippy --all-targets`, `cargo fmt` in `rtc-share`; then `cargo build` and `cargo test` in `rtc-metal` and `rtc-edit` (callers renamed). Expected: all pass, clean.
- [ ] **Step 5: Commit** `ENHANCEMENT: rtc-share wakes the other side with os_sync_wait_on_address and carries load commands (milestone 26)`.

### Task 2: rtc-metal — Hide items and File ▸ Save Image…

**Files:**
- Modify: `RUST/rtc-metal/src/window.rs`, `RUST/rtc-metal/src/main.rs`, `RUST/rtc-metal/README.md`.

**Interfaces:**
- Produces:
  - `pub type Writer = Box<dyn Fn(Vec<[f32; 4]>, usize, usize, &Path) -> Result<(), String>>` (replaces `OnClose`; `run`'s `on_close: OnClose` becomes `write: Writer`; `Delegate.on_close: Cell<Option<OnClose>>` becomes `write: Writer`).
  - `fn offered_name(save: &Save) -> String`: `Save::Ask(name)` → `name`; `Save::To(path)` → its file name.
  - Action `saveImage:` on the delegate's responder (`RtcView`, like `useBvhTraversal:`).

- [ ] **Step 1: Write the failing test** `offered_name_is_the_ask_name_or_the_out_file_name`: `offered_name(&Save::Ask("a.ppm".into())) == "a.ppm"`; `offered_name(&Save::To("dir/b.ppm".into())) == "b.ppm"`.
- [ ] **Step 2: Run** `cargo test offered_name` in `rtc-metal`. Expected: FAIL (not defined).
- [ ] **Step 3: Implement.**
  - `offered_name`; `Writer`, with the quit path calling `(ivars.write)(pixels, w, h, &path)`.
  - `saveImage:`: `save_path(mtm, &offered_name(&save))`; on `Some(path)`, resolve the image (as the quit path does, no borrow held across the panel) and write; `Err` → `eprintln!` and `show_error`. Rendering is not paused.
  - App menu: `hide:` "Hide rtc-metal" `h`; `hideOtherApplications:` "Hide Others" `h` with `setKeyEquivalentModifierMask(Option | Command)`; `unhideAllApplications:` "Show All" no key; separator; Quit. A new File menu with "Save Image…" `s` → `saveImage:`, placed between the app menu and Traversal.
  - README: the menus and Save Image….
- [ ] **Step 4: Run** `cargo test`, `cargo clippy --all-targets`, `cargo fmt`, `cargo build --release` in `rtc-metal`. Expected: pass, clean.
- [ ] **Step 5: Commit** `ENHANCEMENT: rtc-metal has Hide items and File ▸ Save Image… (milestone 26)`.

### Task 3: rtc-metal — taking loads

**Files:**
- Create: `RUST/rtc-metal/src/loads.rs` (`pub mod loads;` in `lib.rs`).
- Modify: `RUST/rtc-metal/Cargo.toml` (`dispatch2 = "0.3"`), `RUST/rtc-metal/src/window.rs`, `RUST/rtc-metal/src/main.rs`, `RUST/rtc-metal/README.md`.

**Interfaces:**
- Consumes: Task 1's `Shared::open`, `load_seq`, `wait_load`, `Load`; Task 2's `offered_name`.
- Produces:
  - `Renderer::load(&mut self, desc: &SceneDescription, base_dir: &Path, fresh: bool) -> Result<(), String>`; today's `reload(desc, base_dir)` becomes `load(desc, base_dir, false)` (rename callers and tests). `fresh` makes the scene's camera the view and home even when `same_view` holds.
  - `impl Save { fn for_scene(&self, scene: &Path) -> Save }`: `Ask(_)` → `Ask(image_name(scene))`; `To(p)` → `To(p)` unchanged.
  - `fn activates_at_start(shared: bool) -> bool` (`!shared`); `run` gains `activate: bool` and calls `app.activate()` only when true.
  - `pub fn deliver(load: Load)` in `window.rs`: main thread only; finds the view through a `thread_local!` set by `run`, and loads `load.path` with `load.fresh` through the same path as the `--watch` reload (refactor `reload_if_changed` into `reload_from(&self, path: &Path, fresh: bool)`), then sets the title to the file name and `Delegate.save` to `save.for_scene(path)` (`title: RefCell<String>`, `save: RefCell<Save>`). Never activates or raises the window.
  - `pub fn listen(share: &Path) -> Result<(), String>` in `loads.rs`: opens its own `Shared`, takes `load_seq()` as the baseline, spawns a thread looping on `wait_load`, sending each `Load` to `DispatchQueue::main().exec_async(move || window::deliver(load))`; an `Err` prints `rtc-metal: loads stopped: {e}` and ends the thread.

- [ ] **Step 1: Write the failing tests** (`window.rs` tests, using the existing `chapter7()` fixture):
  - `a_fresh_load_with_the_same_camera_resets_the_view`: drag the orbit, `moved = true`, `restart()`; `load(&desc, "scenes", true)`; `!r.moved`, `r.camera.transform()` approx equals the scene camera's, `r.acc.frames() == 0`.
  - `a_load_that_is_not_fresh_keeps_the_view`: as `reload_with_the_same_camera_keeps_the_view`, through `load(.., false)` (the old test renamed).
  - `save_for_a_new_scene_offers_its_name`: `Save::Ask("old.ppm".into()).for_scene("/t/new.json")` is `Ask("new.ppm")`; `Save::To("o.ppm".into()).for_scene(..)` is `To("o.ppm")` (derive `PartialEq, Debug` on `Save`).
  - `the_viewer_activates_only_when_run_alone`: `activates_at_start(false)` and `!activates_at_start(true)`.
- [ ] **Step 2: Run** `cargo test` in `rtc-metal`. Expected: FAIL (missing items).
- [ ] **Step 3: Implement** the Produces block. `main.rs`: with `--share`, call `loads::listen(path)` before `run` (an `Err` prints and the viewer runs without loads) and pass `activate: share.is_none()`; usage text: `--share PATH` "write the camera to PATH and take load commands from it; don't take the focus at start". README: the same, and that the editor no longer passes `--watch`.
- [ ] **Step 4: Run** `cargo test`, `cargo clippy --all-targets`, `cargo fmt`, `cargo build --release` in `rtc-metal`. Expected: pass, clean.
- [ ] **Step 5: Commit** `ENHANCEMENT: rtc-metal takes load commands from the shared file (milestone 26)`.

### Task 4: rtc-edit — one viewer, loads, and the camera stream

**Files:**
- Modify: `RUST/rtc-edit/src/viewer.rs`, `RUST/rtc-edit/src/main.rs`, `RUST/rtc-edit/README.md`.

**Interfaces:**
- Consumes: Task 1's `Shared`, `send_load`, `wait_camera`, `generation`, `end_viewer`, `Snapshot`; Task 3's viewer behaviour.
- Produces:
  - `pub fn share_path(root: &Path, pid: u32) -> PathBuf` in `viewer.rs`: `root/rtc-edit-<pid>.share`.
  - `viewer_args(preview: &Path, at: (f64, f64), share: &Path) -> Vec<OsString>`: `<preview> --window --at X,Y --share <share>` (no `--watch`; `share` no longer optional). `Viewer::start` takes `share: &Path` too.
  - `Message::Camera(Snapshot)`, `Message::CameraLost(String)`.
  - `fn camera_stream(key: &(PathBuf, u32)) -> impl Stream<Item = Message>`: opens its own `Shared` at `key.0`; a `std::thread` loops `wait_camera(last, key.1)`, sending through `iced::futures::channel::mpsc::unbounded`; `Ok(Some(s))` → `Camera(s)`, `last = s.seq`; `Ok(None)` or a closed channel → end; `Err(e)` → `CameraLost(e)` and end. An `open` error yields `CameraLost`.

- [ ] **Step 1: Write the failing tests** (`viewer.rs`):
  - `share_path_is_named_for_the_pid`: `share_path("/t".as_ref(), 42) == PathBuf::from("/t/rtc-edit-42.share")`.
  - Update `viewer_args`'s existing test(s) to: `viewer_args("/p/s.json", (10.0, 20.0), "/t/x.share")` == `["/p/s.json", "--window", "--at", "10,20", "--share", "/t/x.share"]`, and assert no `"--watch"`.
- [ ] **Step 2: Run** `cargo test --lib viewer` in `rtc-edit`. Expected: FAIL.
- [ ] **Step 3: Implement.**
  - `viewer.rs`: `share_path`, the new `viewer_args`/`start`; module doc without `--watch`.
  - `Editor.share`: made by `show_viewer` when `None` (`Shared::create(&share_path(&temp_dir(), pid))`; error → status `no shared file: {e}`, no viewer), then kept. `clean_up` also removes `share_path(&temp_dir(), pid)` (errors ignored).
  - `show_viewer`: as today minus the per-viewer `Shared`; start with the editor's share path.
  - `switch`: no longer drops the viewer or the share, nor calls `show_viewer`; stops following, clears `seen`, drops then remakes the preview as today; when `live`, writes the new preview and `send_load(preview.path(), true)` (`Err` → status). Doc: "the running viewer shows `desc` from its own camera".
  - `load`: write the preview, then `send_load(preview.path(), false)`; no viewer → status `no viewer running`.
  - `tick`: only `live`; when it turns false (dead viewer), `share.end_viewer()`; no camera read.
  - `subscription`: add `Subscription::run_with((path, generation), camera_stream)` while `live` and `share` is `Some` (`generation` read from `share`).
  - `Message::Camera(s)`: if `seen` differs in `seq`, set `seen`, `follow(s)`. `Message::CameraLost(e)`: stop following, status `lost the viewer's camera: {e}`. Both join the messages that don't cancel `pending`.
  - Module docs and README: no `--watch`; New/Open keep the viewer; Cmd-R doesn't take the focus; milestone 26.
- [ ] **Step 4: Run** `cargo test --lib`, `cargo clippy --all-targets`, `cargo fmt`, `cargo build --release` in `rtc-edit` (and `cargo build --release` in `rtc-metal`). Expected: pass, clean.
- [ ] **Step 5: Stop for the by-eye check** (report; do not commit): Cmd-N and Cmd-O keep the focus and the viewer shows the new scene with a reset orbit (twice with Untitled); Cmd-L keeps the orbit; following tracks; close the viewer, Cmd-R, following works again; Cmd-R opens the viewer behind the editor; Cmd-H / Hide Others / Show All in the viewer; Save Image… writes the image, also with Cmd-L sent while its panel is up; quitting the editor with a viewer running asks and removes `$TMPDIR/rtc-edit-<pid>.share`; `rtc-metal --window --watch scenes/chapter7.json` alone activates and reloads.
- [ ] **Step 6: Commit when told** `ENHANCEMENT: rtc-edit keeps one viewer through New and Open and follows its camera without a timer (milestone 26)`.
