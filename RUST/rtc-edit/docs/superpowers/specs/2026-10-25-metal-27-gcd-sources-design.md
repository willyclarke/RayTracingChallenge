# Milestone 27: GCD sources instead of timers

Branch `metal/27-gcd-sources`, merged into develop with `--no-ff` when
done. Milestone 26 gave the editor and viewer a command channel that
sleeps in `os_sync_wait_on_address`; 28 is the help window.

## Goal

Four timers still poll: the viewer's 250 ms `NSTimer` (the SIGTERM flag
and the `--watch` file's modification time), the editor's 100 ms Tick
(a dead viewer, and a held load waiting for the viewer to be ready) and
its 50 ms menu poll. Replace each with an event: a GCD dispatch source
(signal, proc, vnode) or a callback, so an idle editor and viewer run no
code at all.

Success:
- No polling timer is left in either app.
- An idle editor and viewer show 0% CPU in Activity Monitor.
- Every behaviour from milestone 26's by-eye checklist is unchanged.
- `--watch` follows an atomic save (write a temporary file, rename it
  over the scene), a save in place, and a file deleted and made again.

## Out of scope

- The help window: milestone 28.
- The remaining 26 minors: a load during the quit alert; a writer killed
  mid-write; the spec's stale lines; `hands_focus_back`'s always-true
  `quiet`; `load()` checking `live`; `end_viewer` only when live (a proc
  source per pid reports every exit, which mostly removes it).

## Decisions

### Dispatch sources

A dispatch source joins a kernel event (`kqueue`) to a queue: created
inactive with a type, a handle, a mask and a queue, given an event
handler block (and a cancel handler), then activated; `get_data` says
which mask bits fired, and events arriving while the handler runs are
merged. `dispatch2` 0.3 (`DispatchSource`, `block2` handlers) builds
them; the owner keeps the source alive, as dropping it releases it.

### The viewer (`rtc-metal`)

- **SIGTERM:** at start `signal(SIGTERM, SIG_IGN)` and a signal source on
  the main queue. Its handler sets `QUIT` (still read by
  `offers_cancel`) and calls `terminate`. The `signal()` handler
  `on_sigterm` goes.
- **`--watch`:** `watch.rs` becomes a vnode-source watcher on the main
  queue; `Watch::changed` (the modification-time filter) stays as the
  last check.
  - The file, opened `O_EVTONLY`, is watched for `WRITE | EXTEND |
    ATTRIB`: each calls the reload check (`reload_if_changed`), so an
    event that leaves the modification time alone does nothing.
  - On `DELETE | RENAME` the file's source is cancelled. If the path
    opens again (the usual case after an atomic save) the new file is
    watched and the reload check runs at once; if not, the scene's
    folder is watched for `WRITE` until the path opens, then the file
    is watched again and the check runs.
  - Each source's cancel handler closes its descriptor.
- The 250 ms `NSTimer`, `checkQuitAndWatchedFile:` and
  `WATCH_INTERVAL` go.
- The SIGTERM source lives in the delegate, the watcher in the renderer
  (as `Watch` does today), for the app's life.

### The editor (`rtc-edit`)

- Depends on `dispatch2 = "0.3"`.
- **A dead viewer:** `Viewer::pid()`. While a viewer runs,
  `Subscription::run_with(pid, viewer_exit)`: the stream makes a proc
  source for `EXIT` on a global queue whose handler sends
  `Message::ViewerExited` into an unbounded channel. A process that died
  before the source was active is never reported, so after activating
  the stream checks the child once and sends `ViewerExited` if it has
  exited. `ViewerExited` runs today's `viewer_ended` ("viewer closed",
  `end_viewer`, following stops, the held load is dropped). `live` is
  set by `show_viewer` when the viewer starts.
- **Held loads and the camera:** the viewer calls `viewer_ready` before
  its first `write_camera`, so the new viewer's first camera write means
  it is ready, and the held load is sent on that `Camera`. `show_viewer`
  records the camera seq before starting the viewer; the camera
  stream's key becomes `(path, generation, start_seq)` and it waits from
  `start_seq`, so it never yields the last viewer's camera. The
  `is_ready` checks stay.
- **The menu:** `Subscription::run(menu_events)` calls
  `MenuEvent::set_event_handler` with a closure that forwards each event
  (sent by muda on the main thread) into an unbounded channel; the
  stream yields `Message::Menu(event)`, mapped through `menu::command`
  as `Poll` mapped them. A restarted subscription sets the handler again
  for its new channel; a send into a closed channel is ignored.
- `Message::Poll`, `Message::Tick`, `tick()` and both
  `iced::time::every` subscriptions go, and iced's `tokio` feature if
  nothing else needs it.

### The shared file (`rtc-share`)

- `end_viewer` moves the generation with one 64-bit
  `fetch_add(1 << 32)` on the camera word, so no access to the word is
  of another size than its other accesses (the 26 review's
  mixed-size-atomics minor).

## Error handling

- The viewer's SIGTERM source can't be made: SIGTERM keeps its default
  action (the viewer dies without the save alert); stderr says why.
- The `--watch` file can't be opened: its folder is watched instead; the
  folder can't be opened either: the overlay says why and `--watch` is
  off.
- The editor's proc source can't be made: `ViewerExited` is sent at once
  if the viewer is gone; otherwise the status line says the viewer can't
  be watched, and `Viewer::running` is checked again on the next Start
  Viewer.
- A burst of saves: the source merges events into one or two handler
  calls; `Watch::changed` lets one reload through per new modification
  time.

## Testing

`rtc-share`:
- `wait_camera_ends_when_the_viewer_ends` and the ready tests pass with
  the 64-bit `end_viewer`.

`rtc-metal` (`Watch` driven by its source on a test queue, each change
awaited on a channel):
- a write in place is reported;
- a new file renamed over the scene is reported, and a later write to it
  too;
- a deleted file made again is reported once it is back.

`rtc-edit`:
- the viewer-exit stream yields `ViewerExited` for a child that exits
  (`sleep 0.1`);
- and for a child that exited before the stream started;
- a camera stream started from `start_seq` does not yield the camera
  written before it.

By eye: the menu bar's commands; closing the viewer ("viewer closed");
Cmd-R then Cmd-N with `RTC_METAL=/tmp/slow-viewer.sh`; Cmd-Q with the
viewer open (the alert, no Cancel); `rtc-metal --window --watch` alone
with an atomic save, an in-place save, and delete-and-recreate; Activity
Monitor showing both apps at 0% CPU when idle.
