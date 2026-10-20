# Milestone 23: the viewer's view into the scene's camera — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The editor shows the viewer's camera live and, with **follow** ticked, copies the viewer's `from` (and focal distance) into the scene; **Keep** keeps it, unticking restores.

**Architecture:** A one-file crate `rtc-share` maps a small file `MAP_SHARED` in both programs, guarded by `flock` and a compile-time signature of its own source. rtc-metal (`--share <path>`) writes its camera on every change; rtc-edit creates the file per viewer, reads it on a 10 Hz timer, and a pure `Follower` applies, keeps or restores.

**Tech Stack:** Rust 2024; `extern "C"` `mmap`/`munmap`/`flock` (no new crate); rtc-metal on objc2 AppKit; rtc-edit on iced 0.14 (gains its `tokio` feature for `iced::time::every`).

**Spec:** `RUST/rtc-edit/docs/superpowers/specs/2026-10-20-metal-23-view-to-camera-design.md`

## Global Constraints

- Branch `metal/23-view-to-camera`; commits end ` (milestone 23)`, prefixes `ENHANCEMENT:` / `TEST:` / `DOCS:` / `FIX:`; never push.
- `rtc-share` is one file, `RUST/rtc-share/src/lib.rs` (tests included), so `SIGNATURE = fnv1a(include_str!("lib.rs"))` covers the whole layout.
- `View` is `#[repr(C)]`: `signature: u64, size: u64, seq: u64, from: [f64; 3], to: [f64; 3], up: [f64; 3], focal_distance: f64` (104 bytes).
- FNV-1a 64: offset basis `0xcbf29ce484222325`, prime `0x100000001b3`.
- macOS constants: `PROT_READ = 1`, `PROT_WRITE = 2`, `MAP_SHARED = 1`, `MAP_FAILED = !0usize as *mut c_void`, `LOCK_SH = 1`, `LOCK_EX = 2`, `LOCK_UN = 8`.
- Shared reads and writes copy the whole `View` with `ptr::read_volatile` / `ptr::write_volatile` while holding the lock.
- The editor's file: `<preview folder>/view`, created anew for every viewer it starts; `--share` only with `--window`.
- Status words in the viewer column: `no viewer`, `waiting`. Status line on Cmd-R with a viewer open: `viewer already running`.
- Followed values go into the scene rounded to 4 decimals (readable boxes and files); the viewer column shows 3 decimals.
- Timer: `iced::time::every(Duration::from_millis(100))`.
- Tests: rtc-share `cargo test`; rtc-metal `cargo test --release --lib`; rtc-edit `cargo test --lib`; all `cargo clippy --all-targets` clean and `cargo fmt`. The "future-incompat … block v0.1.6" warning is pre-existing.
- Test names: descriptive snake_case.

## Review Focus

1. **A viewer change that leaves `from` as it was** (traversal switch, window resize, Cmd-L reload) bumps `seq`. Expected: no unsaved mark. Task 4: `an_unchanged_from_is_no_edit`.
2. **The viewer quits while following.** Expected: following stops, the boxes keep the last value and become editable, the column says `no viewer`. Task 4: `stopping_keeps_the_values`; Task 5 by eye.
3. **A stale rtc-metal** (built before rtc-share changed). Expected: it prints why and runs without sharing; the editor shows `waiting`. Task 1: `open_refuses_another_signature`.
4. **New or Open while following.** Expected: the new scene isn't touched by the old cache. Task 4: `stopping_keeps_the_values`; Task 5 calls `stop` in `switch`.
5. **Unticking after Cmd-S while following.** Expected: the pre-tick `from` comes back as an unsaved edit. Task 4: `unticking_restores_from_and_focal_distance` asserts `untick` returns `true`.

---

### Task 1: the `rtc-share` crate

**Files:**
- Create: `RUST/rtc-share/Cargo.toml` (`name = "rtc-share"`, `edition = "2024"`, no dependencies), `RUST/rtc-share/src/lib.rs`, `RUST/rtc-share/README.md` (what it is, the layout, the signature, `flock`).

**Interfaces:**
- Produces:
  - `pub const fn fnv1a(s: &str) -> u64`; `pub const SIGNATURE: u64`.
  - `#[derive(Clone, Copy, Debug, PartialEq)] pub struct Snapshot { pub seq: u64, pub from: [f64; 3], pub to: [f64; 3], pub up: [f64; 3], pub focal_distance: f64 }` (the spec's `Camera`, renamed: rtc-metal already has a `Camera`).
  - `pub struct Shared` (owns the `File` and the mapping; `unsafe impl Send`; unmaps on `Drop`).
  - `Shared::create(path: &Path) -> Result<Shared, String>`: removes any existing file first (a viewer still quitting keeps the old one), creates it, `set_len(104)`, maps, writes the header, `seq` 0.
  - `Shared::open(path: &Path) -> Result<Shared, String>`: maps; `Err` if missing, shorter than `View`, or `signature`/`size` differ — the message names the path and says `rebuild rtc-metal`.
  - `Shared::write(&self, from: [f64; 3], to: [f64; 3], up: [f64; 3], focal_distance: f64)`: `seq + 1`.
  - `Shared::read(&self) -> Option<Snapshot>`: `None` while `seq == 0`.

- [ ] **Step 1: Write the failing tests** in `lib.rs` (files under `std::env::temp_dir().join(format!("rtc-share-test-{}", std::process::id()))`, one name per test, removed at the end):
  - `fnv1a_matches_known_answers`: `fnv1a("") == 0xcbf29ce484222325`, `fnv1a("a") == 0xaf63dc4c8601ec8c`.
  - `signature_is_the_hash_of_this_file`: `SIGNATURE == fnv1a(include_str!("lib.rs"))`.
  - `read_is_none_until_the_first_write`: `create(p)?.read() == None`.
  - `a_write_through_one_mapping_is_read_through_another`: `let a = create(p)?; let b = open(p)?; a.write([1.,2.,3.],[0.,1.,0.],[0.,1.,0.],5.)`; `b.read() == Some(Snapshot { seq: 1, from: [1.,2.,3.], to: [0.,1.,0.], up: [0.,1.,0.], focal_distance: 5. })`; a second write gives `seq` 2.
  - `open_refuses_another_signature`: after `create`, overwrite bytes 0..8 with `!SIGNATURE` via `std::fs` (write at offset 0); `open(p)` is `Err` containing `rebuild rtc-metal`.
  - `open_refuses_a_short_or_missing_file`: a 10-byte file → `Err`; a missing path → `Err`.
  - `create_replaces_an_old_file`: `create`, `write`, `create` again → `read() == None`.
- [ ] **Step 2: Run** `cargo test` in `RUST/rtc-share`. Expected: FAIL to compile (items missing).
- [ ] **Step 3: Implement** the items above. `fnv1a` loops over `s.as_bytes()` with a `while` (const fn). `read`/`write` take `flock(fd, LOCK_SH|LOCK_EX)` and `LOCK_UN` around the volatile copy.
- [ ] **Step 4: Run** `cargo test`, `cargo clippy --all-targets`, `cargo fmt`. Expected: 7 pass, clean.
- [ ] **Step 5: Commit** `ENHANCEMENT: rtc-share maps a camera both programs can read, signed by its source (milestone 23)`.

### Task 2: rtc-metal writes its camera with `--share`

**Files:**
- Modify: `RUST/rtc-metal/Cargo.toml` (`rtc-share = { path = "../rtc-share" }`), `RUST/rtc-metal/src/orbit.rs`, `RUST/rtc-metal/src/window.rs`, `RUST/rtc-metal/src/main.rs`, `RUST/rtc-metal/README.md`.

**Interfaces:**
- Consumes: `Shared::open`, `Shared::write` (Task 1).
- Produces: `rtc-metal <scene> --window [--watch] [--at X,Y] [--share PATH]`.
  - `Orbit::to(&self) -> Tuple`, `Orbit::up(&self) -> Tuple`.
  - `pub fn shot(orbit: &Orbit, focal_distance: f64) -> ([f64; 3], [f64; 3], [f64; 3], f64)` in `window.rs`: from, to, up as arrays (x, y, z).
  - `run(…, at, share: Option<Shared>, save, on_close)`.

- [ ] **Step 1: Write the failing tests:**
  - `orbit.rs` `orbit_keeps_its_to_and_up`: `Orbit::new(point(0,1,-5), point(0,1,0), vector(0,1,0))`, after `drag(30.0, 10.0)`, `to()` is `point(0,1,0)` and `up()` `vector(0,1,0)`.
  - `window.rs` `shot_gives_the_orbits_camera`: the same orbit, `shot(&o, 4.5)` == `([0.,1.,-5.], [0.,1.,0.], [0.,1.,0.], 4.5)` within 1e-9 per component; after `o.zoom(…)` its `from` equals `o.from()`.
- [ ] **Step 2: Run** `cargo test --release --lib`. Expected: FAIL (missing).
- [ ] **Step 3: Implement** `to`, `up`, `shot`.
- [ ] **Step 4: Wire it.**
  - `main.rs`: `--share PATH` → `Shared::open`; on `Err(e)` print `rtc-metal: {e}` to stderr and carry `None`; `--share` without `--window` → `usage()` (as `--watch`, `--at`); usage text and doc comment list `[--share PATH]`; pass through `open_window` to `run`.
  - `window.rs`: `Renderer` gains `share: Option<Shared>` and `fn publish(&self)` writing `shot(&self.orbit, self.camera.focal_distance)`; call it at the end of `Renderer::restart` (every camera change restarts) and once in `run` after the renderer is built. A publish with an unchanged camera is harmless (Review Focus 1 is the editor's).
- [ ] **Step 5: README**: `--share PATH`, what is written and when, the signature check and its message.
- [ ] **Step 6: Run** `cargo test --release --lib`, `cargo clippy --all-targets`, `cargo fmt`, `cargo build --release`. Expected: pass, clean.
- [ ] **Step 7: Commit** `ENHANCEMENT: rtc-metal --share writes its camera on every change (milestone 23)`.

### Task 3: rtc-edit starts each viewer with a shared file

**Files:**
- Modify: `RUST/rtc-edit/Cargo.toml` (`rtc-share = { path = "../rtc-share" }`; `iced = { version = "0.14", features = ["tokio"] }`), `RUST/rtc-edit/src/viewer.rs`, `RUST/rtc-edit/src/preview.rs`, `RUST/rtc-edit/src/main.rs`.

**Interfaces:**
- Consumes: `Shared::create` (Task 1); rtc-metal's `--share PATH` (Task 2).
- Produces: `viewer_args(preview: &Path, at: (f64, f64), share: Option<&Path>) -> Vec<OsString>`; `Viewer::start(exe, preview, at, share: Option<&Path>)`; `Preview::dir(&self) -> &Path`; `Editor.share: Option<Shared>` (the current viewer's file).

- [ ] **Step 1: Write the failing test** `viewer_args_pass_the_shared_file`: `viewer_args(Path::new("/p/s.json"), (0.0, 0.0), Some(Path::new("/p/view")))` ends with `["--share", "/p/view"]`; with `None` it is the 5 arguments of today. Update the existing `viewer_args_put_the_window_at_the_point` call with `None`.
- [ ] **Step 2: Run** `cargo test --lib viewer`. Expected: FAIL.
- [ ] **Step 3: Implement** `viewer_args`, `Viewer::start`, `Preview::dir`.
- [ ] **Step 4: Wire it** in `main.rs`:
  - `show_viewer`: with a viewer running, set the status `viewer already running` and return. Otherwise `Shared::create(&preview.dir().join("view"))`; `Ok` → `self.share = Some(s)`, start with `Some(path)`; `Err(e)` → `self.share = None`, start without, and the status after start reads `viewer started; no shared view: {e}`.
  - `switch`: `self.share = None` before the preview is replaced.
- [ ] **Step 5: Run** `cargo test --lib`, `cargo clippy --all-targets`, `cargo fmt`, `cargo build --release`. Expected: pass, clean.
- [ ] **Step 6: Commit** `ENHANCEMENT: rtc-edit gives each viewer a shared camera file, and says when Cmd-R finds one running (milestone 23)`.

### Task 4: following, as pure state

**Files:**
- Create: `RUST/rtc-edit/src/follow.rs` (`pub mod follow;` in `lib.rs`).
- Modify: `RUST/rtc-edit/src/fields.rs` (extract the focal distance setter).

**Interfaces:**
- Consumes: `Snapshot` (Task 1).
- Produces:
  - `pub fn set_focal_distance(blur: &mut FocalBlur, distance: f64)` in `fields.rs`: keeps the f-number (the existing logic in `set_camera`'s `FocalDistance` arm, which now calls it).
  - `#[derive(Default)] pub struct Follower` with `cache: Option<([f64; 3], Option<f64>)>` (the `from` and, with focal blur, the focal distance at tick time).
  - `following(&self) -> bool`; `start(&mut self, camera: &CameraDescription)`; `update(&mut self, camera: &mut CameraDescription, shot: &Snapshot) -> bool` (only while following; rounds to 4 decimals; `true` when the camera changed); `keep(&mut self)`; `untick(&mut self, camera: &mut CameraDescription) -> bool` (restores, stops, `true` when changed); `stop(&mut self)` (stops without restoring).

- [ ] **Step 1: Write the failing tests** in `follow.rs` (a camera `from [0,1.5,-5]`, `to [0,1,0]`, `up [0,1,0]`; a snapshot `from [0.81234, 2.1, -4.87]`, `focal_distance 6.0`):
  - `nothing_follows_until_started`: `update` → `false`, camera unchanged.
  - `following_copies_from_rounded`: after `start`, `update` → `true`, `from == [0.8123, 2.1, -4.87]`.
  - `an_unchanged_from_is_no_edit`: a second `update` with the same snapshot (any `seq`) → `false`.
  - `focal_distance_follows_only_with_focal_blur`: without `focal_blur` it stays `None`; with `FocalBlur { focal_distance: 5.0, aperture: 0.25, samples: 32 }` it becomes `6.0` and `aperture == 0.3` (f-number 10 kept).
  - `unticking_restores_from_and_focal_distance`: start, update, `untick` → `true`, `from` and focal distance back, `following() == false`.
  - `keeping_keeps_the_values`: start, update, `keep` → `following() == false`, `from` stays the followed one; a later `update` → `false`.
  - `stopping_keeps_the_values`: as keep, via `stop`.
  - In `fields.rs` the existing focal distance tests stay green after the extraction.
- [ ] **Step 2: Run** `cargo test --lib follow`. Expected: FAIL.
- [ ] **Step 3: Implement** `set_focal_distance` (extract) and `Follower`.
- [ ] **Step 4: Run** `cargo test --lib`, `cargo clippy --all-targets`, `cargo fmt`. Expected: pass, clean.
- [ ] **Step 5: Commit** `ENHANCEMENT: rtc-edit follows, keeps and restores the viewer's camera as pure state (milestone 23)`.

### Task 5: the viewer column, follow, Keep

**Files:**
- Modify: `RUST/rtc-edit/src/main.rs`, `RUST/rtc-edit/README.md`.

**Interfaces:**
- Consumes: `Editor.share` (Task 3), `Follower`, `Snapshot` (Tasks 1, 4).
- Produces: `Message::Tick`, `Message::Follow(bool)`, `Message::Keep`; Cmd-K → `Keep`.

- [ ] **Step 1: State and messages.** `Editor` gains `seen: Option<Snapshot>` and `follower: Follower`. `subscription` batches `iced::time::every(Duration::from_millis(100)).map(|_| Message::Tick)`. The update preamble that clears `pending` exempts `Tick` too.
  - `Tick`: no viewer running → `seen = None`, `follower.stop()`. Else `share.read()`; when its `seq` differs from `seen`'s: `seen = it`, and `follower.update(&mut desc.camera, &it)` returning `true` sets `unsaved` and refreshes the camera's boxes (the existing `refresh_camera` path, so a box being edited isn't overwritten).
  - `Follow(true)` → `follower.start(&desc.camera)`, then apply `seen` at once if any. `Follow(false)` → `untick`; `true` sets `unsaved` and refreshes the boxes.
  - `Keep` (button or Cmd-K, `"k"` in `shortcut`'s Cmd table) → `follower.keep()`.
  - `switch` → `follower.stop()`, `seen = None`.
- [ ] **Step 2: View.** In `boxes`, for the camera only, after the third box of the `from`, `to` and `up` rows: the viewer's values as 3-decimal text, greyed (`size 13`, monospace), or `no viewer` / `waiting` on the `from` row; after `from`'s also `checkbox("follow")` (`on_toggle(Message::Follow)`, only while a viewer runs) and, while following, `button("Keep").on_press(Message::Keep)`. While following, the `from` boxes have no `on_input` (read-only).
- [ ] **Step 3: README**: the viewer column, follow, Keep and Cmd-K (key table), what unticking restores, when following stops, Cmd-R's `viewer already running`, and milestone 23 in the list.
- [ ] **Step 4: Run** `cargo test --lib`, `cargo clippy --all-targets`, `cargo fmt`, `cargo build --release` (and `cargo build --release` in `RUST/rtc-metal`). Expected: pass, clean.
- [ ] **Step 5: Stop for the by-eye check** (report; do not commit). From `RUST/rtc-edit`, `cargo run --release -- ../rtc-metal/scenes/chapter10.json`, select the camera: orbit and zoom the viewer, the column follows; tick follow, the `from` boxes follow and the title shows ` •`; Keep, Cmd-S, Cmd-O the same file: the kept view; tick, orbit, untick: the old `from` back; Cmd-L after Keep: no jump; Cmd-R with the viewer open: `viewer already running`; quit the viewer while following: `no viewer`, boxes editable; a focal blur scene (`../rtc-metal/scenes/` one with `focal_blur`, or type an f-number): zoom while following keeps the focus after Keep and Cmd-L.
- [ ] **Step 6: Commit when told** `ENHANCEMENT: rtc-edit shows the viewer's camera and follows it into the scene (milestone 23)`.
