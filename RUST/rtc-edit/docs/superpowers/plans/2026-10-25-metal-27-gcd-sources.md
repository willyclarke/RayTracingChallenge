# GCD Sources Instead of Timers Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the viewer's 250 ms timer and the editor's 100 ms Tick and 50 ms menu poll with GCD dispatch sources and callbacks, so an idle editor and viewer run no code.

**Architecture:** rtc-metal gets a signal source for SIGTERM on the main queue and a vnode-source file watcher (`FileEvents`, file ⇄ folder switching on a private serial queue) that hands reloads to the main queue. rtc-edit gets a proc source per viewer feeding an iced subscription, sends held loads on the new viewer's first camera write (the camera stream starts from the seq recorded at Start Viewer), and takes muda's menu events through `set_event_handler` into a subscription. rtc-share's `end_viewer` becomes one 64-bit add.

**Tech Stack:** Rust 2024, dispatch2 0.3 (`DispatchSource`, `DispatchQueue`), block2 `RcBlock`, objc2 AppKit, iced 0.14 `Subscription::run` / `run_with`, muda 0.21, futures unbounded channels.

**Spec:** `RUST/rtc-edit/docs/superpowers/specs/2026-10-25-metal-27-gcd-sources-design.md`

## Global Constraints

- Each crate (`RUST/rtc-share`, `RUST/rtc-metal`, `RUST/rtc-edit`) is its own Cargo project; run cargo with `--manifest-path` or from the crate.
- rtc-metal tests: `cargo test --release --lib -- --skip test_chap_17_gpu_perlin_noise_matches_cpu` (that probe fails before this branch); plus `cargo test --release --test share`.
- `dispatch2 = "0.3"` in both rtc-metal and rtc-edit; `block2` as a direct dependency where a handler block is built (same version the lock already has).
- Source types: `&raw const dispatch2::_dispatch_source_type_{signal,proc,vnode}` cast to `dispatch_source_type_t`; `DispatchSource::new` panics on NULL, so it is only called with a valid handle (an open descriptor, a pid we started, SIGTERM).
- Every source is created inactive, given its handlers, then `activate()`d; its owner keeps the `DispatchRetained<DispatchSource>`.
- Handler closures given to a non-main queue are `Send + Sync + 'static`; nothing non-`Send` crosses a queue.
- No polling timer may remain: no `NSTimer`, no `iced::time::every`.
- Commits: `ENHANCEMENT:`/`TEST:`/`DOCS:`/`FIX:` prefix, subject ends ` (milestone 27)`, body ends `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`; explicit `git add` paths; never add `RUST/rtc-edit/untitled.json` or untracked files in `RUST/rtc-edit/examples/`.
- UI behaviour is checked by the user's eye before Task 4's and Task 5's commits; stop and report there.

## Review Focus

1. **A real editor's save** (vim writes a backup and renames; VS Code and TextEdit save atomically): every save reloads once, and later saves keep working. Task 2: `a_new_file_renamed_over_the_scene_is_reported_and_watched` (two saves in a row); by eye with a real editor.
2. **The watched file deleted for good, or its folder removed:** no crash, no busy loop, the old image stays; a file made again later is picked up. Task 2: `a_deleted_file_made_again_is_reported`.
3. **SIGTERM twice, or while the save panel is up:** the second does nothing; the alert has no Cancel. Task 3: `QUIT` and the `asking` guard keep today's behaviour; by eye with `kill -TERM` twice.
4. **A viewer killed with SIGKILL, or crashing:** the proc source reports it; "viewer closed"; Cmd-R works. Task 4: `the_exit_stream_reports_a_child_that_exits`; by eye with `kill -9`.
5. **The menu after Cmd-R, New, Open or a long idle:** every menu command still arrives (the handler is set once per subscription start). Task 5: by eye after each.

---

### Task 1: rtc-share's `end_viewer` as one 64-bit add

**Files:**
- Modify: `RUST/rtc-share/src/lib.rs` (`end_viewer`, `generation_word` if unused afterwards)

**Interfaces:**
- Consumes: nothing new.
- Produces: unchanged API (`end_viewer`, `generation`).

- [ ] **Step 1: Make every access to the camera word 64-bit.** No new test: the mixed-size race is undefined behaviour in Rust's model, not observable; the existing tests pin the behaviour.
  - `end_viewer`: `camera_word().fetch_add(1 << 32, Release)`, then the wake.
  - `generation()`: `(camera_word().load(Acquire) >> 32) as u32`.
  - The camera seqlock: its odd/even steps become a `compare_exchange_weak` loop on the 64-bit word that changes only the low half (`(w & !0xffff_ffff) | low.wrapping_add(1) as u64`), since `end_viewer` may change the high half at the same time; readers take `load(Acquire) as u32`. The load block keeps its own 32-bit `seq_write`/`seq_read`.
  - Remove `camera_seq` and `generation_word`.

- [ ] **Step 2: Run the suite**

Run: `cargo test --lib --manifest-path RUST/rtc-share/Cargo.toml`
Expected: 14 passed (incl. `wait_camera_ends_when_the_viewer_ends`, `a_viewer_is_ready_from_saying_so_until_it_ends`)

- [ ] **Step 3: Commit**

`FIX: rtc-share moves the viewer generation with one 64-bit add (milestone 27)` — `git add RUST/rtc-share/src/lib.rs`

### Task 2: rtc-metal's `FileEvents`, a vnode-source file watcher

**Files:**
- Modify: `RUST/rtc-metal/src/watch.rs` (keep `Watch`; add `FileEvents`)
- Modify: `RUST/rtc-metal/src/window.rs` (`run` starts `FileEvents` when `watch` is `Some`; `pub fn reload_watched()` via `WINDOW`; `checkQuitAndWatchedFile:` stops calling `reload_if_changed`)
- Modify: `RUST/rtc-metal/Cargo.toml` (`block2`), `RUST/rtc-metal/README.md`

**Interfaces:**
- Produces: `pub struct FileEvents` and `impl FileEvents { pub fn start(path: &Path, notify: impl Fn() + Send + Sync + 'static) -> Result<FileEvents, String> }` — dropping it cancels its sources. `pub fn reload_watched()` in `window.rs` (main thread: the view's `reload_if_changed`).

- [ ] **Step 1: Write the failing tests in `watch.rs`**

Helper: `fn events(path) -> (FileEvents, std::sync::mpsc::Receiver<()>)` whose `notify` sends `()`; `fn next(rx) -> bool` = `rx.recv_timeout(Duration::from_secs(2)).is_ok()`; drain with `while rx.try_recv().is_ok() {}` between steps.

```rust
#[test]
fn a_write_in_place_is_reported() // write the file → next(rx)
#[test]
fn a_new_file_renamed_over_the_scene_is_reported_and_watched()
    // write tmp, rename over → next; drain; write tmp2, rename over again → next;
    // drain; write in place → next
#[test]
fn a_deleted_file_made_again_is_reported()
    // remove → (may or may not notify; drain after 200 ms); assert !next within 300 ms;
    // write the path again → next; drain; write in place → next
#[test]
fn start_fails_for_a_missing_folder() // FileEvents::start("/nonexistent/x.json") is Err
```

Each test uses its own folder under `std::env::temp_dir()` named with the test and the pid.

- [ ] **Step 2: Run to verify they fail**

Run: `cargo test --release --lib watch:: --manifest-path RUST/rtc-metal/Cargo.toml`
Expected: FAIL to compile, `FileEvents` not found

- [ ] **Step 3: Implement `FileEvents` in `watch.rs`**

State `Arc<Mutex<Inner>>` with `path`, `notify: Arc<dyn Fn() + Send + Sync>`, `source: Option<DispatchRetained<DispatchSource>>`, all driven on a private serial queue `DispatchQueue::new("rtc-metal.watch", None)`. `open(path)` = `libc::open(c_path, O_EVTONLY)`.

- `watch_file(inner)`: open the path; if it opens, a vnode source with mask `WRITE | EXTEND | ATTRIB | DELETE | RENAME`; its handler reads `data()`: on `DELETE | RENAME` → cancel, then `watch_file` again and `notify()` if the file opened, else `watch_folder`; otherwise `notify()`. If it doesn't open → `watch_folder`.
- `watch_folder(inner)`: open the parent folder, vnode source mask `WRITE`; handler: if the path opens → cancel this source, `watch_file`, `notify()`. A folder that doesn't open → `Err` from `start` (at start) or stop watching (later; stderr).
- Each source's cancel handler closes its descriptor. `Drop for FileEvents` cancels the current source.
- `start` returns `Err` only when neither the file nor its folder opens.

- [ ] **Step 4: Run the tests**

Run: same as Step 2. Expected: 4 new tests pass, plus the existing `Watch` tests.

- [ ] **Step 5: Wire it into the window**

In `run`, when `watch` is `Some`, `FileEvents::start(path, || DispatchQueue::main().exec_async(window::reload_watched))`; keep it alive in the delegate (`file_events: Option<FileEvents>`); an `Err` shows on the overlay as `reload_error` and stderr. `reload_watched` borrows nothing across `reload_if_changed`. The timer selector keeps only its SIGTERM check (Task 3 removes it).

- [ ] **Step 6: Full suite, clippy, release build**

Run: `cargo test --release --lib -- --skip test_chap_17_gpu_perlin_noise_matches_cpu` and `cargo clippy --release --all-targets` and `cargo build --release` in `RUST/rtc-metal`
Expected: all pass (177 + 4), no warnings

- [ ] **Step 7: Commit**

`ENHANCEMENT: rtc-metal's --watch sleeps on a GCD vnode source instead of polling (milestone 27)` — add `src/watch.rs src/window.rs Cargo.toml Cargo.lock README.md`

### Task 3: rtc-metal's SIGTERM signal source; the timer goes

**Files:**
- Modify: `RUST/rtc-metal/src/window.rs` (SIGTERM source; remove `NSTimer`, `checkQuitAndWatchedFile:`, `WATCH_INTERVAL`, `on_sigterm`, the `signal` extern's handler use)
- Modify: `RUST/rtc-metal/README.md`

**Interfaces:**
- Consumes: Task 2's `FileEvents` (the timer no longer reloads).
- Produces: `fn sigterm_source() -> DispatchRetained<DispatchSource>` (main queue), kept in the delegate.

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn a_sigterm_source_reports_the_signal_without_killing() {
    // ignore SIGTERM, a signal source on a test serial queue whose handler
    // sends on a channel, activate, libc::kill(getpid(), SIGTERM),
    // assert rx.recv_timeout(2 s).is_ok(); restore nothing (the test
    // process ignoring SIGTERM is harmless)
}
```

Make the source's creation a function taking the queue and the handler: `fn signal_source(sig: i32, queue: &DispatchQueue, handler: impl Fn() + Send + Sync + 'static) -> DispatchRetained<DispatchSource>`, which also sets `signal(sig, SIG_IGN)`.

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test --release --lib a_sigterm_source --manifest-path RUST/rtc-metal/Cargo.toml`
Expected: FAIL to compile, `signal_source` not found

- [ ] **Step 3: Implement `signal_source`; in `run`, `signal_source(SIGTERM, &DispatchQueue::main(), || { QUIT.store(true, SeqCst); terminate })`**

The handler gets `MainThreadMarker::new()` (it runs on the main queue) for `NSApplication::sharedApplication`. Remove the `NSTimer`, its selector, `WATCH_INTERVAL` and `on_sigterm`. `QUIT` stays for `offers_cancel`.

- [ ] **Step 4: Full suite, clippy, release build** (as Task 2 Step 6)

Expected: all pass (+1), no warnings; `grep -n NSTimer src/window.rs` finds nothing

- [ ] **Step 5: Commit**

`ENHANCEMENT: rtc-metal quits on SIGTERM through a GCD signal source; its timer goes (milestone 27)` — add `src/window.rs README.md`

### Task 4: rtc-edit's viewer-exit source and camera-start; Tick goes

**Files:**
- Modify: `RUST/rtc-edit/src/viewer.rs` (`Viewer::pid`, `exit_source`)
- Modify: `RUST/rtc-edit/src/main.rs` (`Message::ViewerExited`, `viewer_exit` stream, `camera_start`, camera stream key, held-load flush on `Camera`; remove `Message::Tick`, `tick()`, its `every`)
- Modify: `RUST/rtc-edit/Cargo.toml` (`dispatch2 = "0.3"`, `block2`), `RUST/rtc-edit/README.md`

**Interfaces:**
- Produces: `impl Viewer { pub fn pid(&self) -> u32 }`; `pub fn exit_source(pid: u32, notify: impl Fn() + Send + Sync + 'static) -> DispatchRetained<DispatchSource>` (proc source, `DISPATCH_PROC_EXIT`, global queue, activated); `fn viewer_exit(pid: &u32) -> impl Stream<Item = Message>` in main.rs; camera stream key `(PathBuf, u32, u32)` = (share path, generation, start seq).

- [ ] **Step 1: Write the failing tests in `viewer.rs`**

```rust
#[test]
fn the_exit_source_reports_a_child_that_exits() {
    // spawn `sleep 0.1`; exit_source(child.id(), send on channel);
    // assert rx.recv_timeout(3 s).is_ok(); child.wait()
}
#[test]
fn the_exit_source_is_created_for_a_child_that_already_exited() {
    // spawn `true`, sleep 200 ms (exited, not reaped); exit_source(...) must not
    // panic; whether it fires is not asserted (the stream's after-check covers it)
}
```

And in rtc-share's test style, a test in `main.rs` is not possible (binary); put `camera_stream`'s start rule in a pure helper in `viewer.rs`: `pub fn camera_start(share: &Shared) -> u32` = `share.read_camera().map_or(0, |s| s.seq)`, tested:

```rust
#[test]
fn camera_start_is_the_last_write_so_waiting_skips_it() {
    // create a share, write_camera once; s = camera_start(&share);
    // wait_camera(s, generation) in a thread does not finish within 100 ms;
    // write_camera again → it returns the new seq
}
```

- [ ] **Step 2: Run to verify they fail**

Run: `cargo test --lib viewer:: --manifest-path RUST/rtc-edit/Cargo.toml`
Expected: FAIL to compile, `exit_source` / `camera_start` / `pid` not found

- [ ] **Step 3: Implement `pid`, `exit_source`, `camera_start` in `viewer.rs`**

- [ ] **Step 4: Run the tests** — Expected: 3 new pass

- [ ] **Step 5: Wire it into `main.rs`**

- `show_viewer`: before `Viewer::start`, `self.camera_start = camera_start(share)`; on success `self.live = true`.
- `subscription`: while `live` and a viewer exists, `Subscription::run_with(pid, viewer_exit)` and the camera stream with key `(own_share_path(), generation, camera_start)` (it waits from `key.2`).
- `viewer_exit(pid)`: unbounded channel; `exit_source(pid, send ViewerExited)` kept alive inside the stream (`rx.map(move |m| { let _ = &source; m })`); after activating, if `CHILDREN`'s viewer is not `running`, send `ViewerExited` once.
- `Message::ViewerExited` → `viewer_ended()` (idempotent through `live`); add it to the list that doesn't cancel `pending`.
- `camera(shot)`: after the `is_ready` check, send a held load (`send_or_hold` with the preview path) before following.
- Remove `Message::Tick`, `tick()` and its `iced::time::every`.

- [ ] **Step 6: Full suite, clippy, release build**

Run: `cargo test --lib`, `cargo clippy --release`, `cargo build --release` in `RUST/rtc-edit`
Expected: 135 + 3 pass; no new warnings

- [ ] **Step 7: Stop for the by-eye check (report; do not commit)**

Close the viewer → "viewer closed"; `kill -9` the viewer → "viewer closed"; Cmd-R then Cmd-N with `RTC_METAL=/tmp/slow-viewer.sh`; following after Cmd-R; Cmd-Q with the viewer open.

- [ ] **Step 8: Commit after the user's go-ahead**

`ENHANCEMENT: rtc-edit learns of a closed viewer from a GCD proc source; Tick goes (milestone 27)` — add `src/viewer.rs src/main.rs Cargo.toml Cargo.lock README.md`

### Task 5: rtc-edit's menu events by callback; the poll goes

**Files:**
- Modify: `RUST/rtc-edit/src/main.rs` (`Message::MenuEvent(muda::MenuEvent)`, `menu_events` stream; remove `Message::Poll`, its `every`; drop iced's `tokio` feature if the build allows)
- Modify: `RUST/rtc-edit/Cargo.toml`, `RUST/rtc-edit/README.md`, the spec's by-eye list if anything changed

**Interfaces:**
- Produces: `fn menu_events() -> UnboundedReceiver<Message>`: `MenuEvent::set_event_handler(Some(move |e| { let _ = tx.unbounded_send(Message::MenuEvent(e)); }))`.

- [ ] **Step 1: No unit test** — muda's handler fires only from AppKit menu actions; the mapping is `menu::command` (already tested). Ruling-worthy if the executor finds a testable seam.

- [ ] **Step 2: Implement** — `subscription` adds `Subscription::run(menu_events)`; `Message::MenuEvent(e)` maps through `menu::command` and `menu_message` exactly as `Poll` did for one event; remove `Poll` and its `every`; remove `"tokio"` from iced's features and the `Duration` import if unused, keeping it only if the build fails without it.

- [ ] **Step 3: Full suite, clippy, release build** (as Task 4 Step 6)

Expected: 138 pass; `grep -n "time::every" src/main.rs` finds nothing

- [ ] **Step 4: Stop for the by-eye check (report; do not commit)**

Every menu command (File, Edit, View items, Follow the Viewer tick) works, also after Cmd-R, New and Open; Activity Monitor shows rtc-edit and rtc-metal at 0% CPU when idle; `rtc-metal --window --watch scenes/chapter7.json` alone reloads on a vim save, an in-place save, and delete-and-recreate; `kill -TERM` the viewer twice.

- [ ] **Step 5: Commit after the user's go-ahead**

`ENHANCEMENT: rtc-edit takes menu events by callback; no timer is left (milestone 27)` — add `src/main.rs Cargo.toml Cargo.lock README.md`
