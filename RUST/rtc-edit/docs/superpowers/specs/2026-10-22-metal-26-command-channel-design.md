# Milestone 26: a command channel, and one viewer through New and Open

Branch `metal/26-command-channel`, merged into develop with `--no-ff` when
done. Milestone 23 shared the viewer's camera through `rtc-share`; 25 gave
the editor a menu bar. 27 replaces the remaining timers with GCD sources;
28 is the help window.

## Goal

New and Open kill the viewer and start another, which takes the focus from
the editor. The editor and the viewer talk through a polled file (`--watch`)
and an `flock`ed mapping read on a timer. Make it macOS-native: the shared
file carries commands as well as the camera, each side sleeps in
`os_sync_wait_on_address` until the other writes, and one viewer lives
through New and Open.

Success:
- Cmd-N and Cmd-O leave the focus in the editor; the running viewer shows
  the new scene from its own camera.
- Cmd-L reloads the viewer without `--watch`, keeping the orbit.
- Following the viewer reads the camera when it changes, not on a timer.
- Cmd-R opens the viewer without taking the focus.
- The viewer has Hide (Cmd-H), Hide Others, Show All and File ▸ Save
  Image… (Cmd-S).

## Out of scope

- GCD sources for SIGTERM, viewer death, `--watch` and the menu poll:
  milestone 27. The viewer's 100 ms timer, the editor's Tick that notices a
  dead viewer, and the 50 ms menu poll stay.
- The help window: milestone 28.
- The quit alert knowing an image was just saved: it asks as today.

## Decisions

### The shared file (`rtc-share`)

- Lives at `$TMPDIR/rtc-edit-<pid>.share`, outside the preview folder
  (whose Drop removes the folder). Made by the editor's first Start Viewer,
  kept for the editor's life, removed by `clean_up`. Every viewer the
  editor starts uses it.
- Still one source file; the `SIGNATURE` hash of it guards the layout.
- Layout (`#[repr(C)]`):

| Block | Fields | Writer |
|---|---|---|
| header | `signature: u64`, `size: u64` | `create` |
| viewer | `generation: u32` | editor |
| camera | `seq: u32`, `from`, `to`, `up: [f64; 3]`, `focal_distance: f64` | viewer |
| load | `seq: u32`, `fresh: u32`, `len: u32`, `path: [u8; 1024]` | editor |

- **No `flock`.** Each block has one writer, so a seqlock: the writer sets
  `seq` odd, writes, sets it to the next even number, then
  `os_sync_wake_by_address_all(&seq, 4, OS_SYNC_WAKE_BY_ADDRESS_SHARED)`. A
  reader copies the block and retries while `seq` is odd or changed during
  the copy. Words are accessed as atomics on the mapping.
- **Waiting.** `os_sync_wait_on_address(&seq, last, 4,
  OS_SYNC_WAIT_ON_ADDRESS_SHARED)` sleeps until `seq` is no longer `last`
  (or a wake); `EINTR` and spurious wakes loop.
- **Generation.** The editor bumps `generation` when a viewer ends and
  wakes the camera `seq` address; a camera waiter started for another
  generation returns `None`.
- API:
  - `Shared::create(path)`, `Shared::open(path)` (as today);
  - `write_camera(&self, from, to, up, focal_distance)`,
    `read_camera(&self) -> Option<Snapshot>` (as today's `write`/`read`),
    `wait_camera(&self, last: u32, generation: u32) -> Option<Snapshot>`;
  - `send_load(&self, path: &Path, fresh: bool) -> Result<(), String>`
    (`Err` past 1024 bytes), `load_seq(&self) -> u32`,
    `wait_load(&self, last: u32) -> Result<Load, String>` with
    `pub struct Load { pub seq: u32, pub path: PathBuf, pub fresh: bool }`;
  - `generation(&self) -> u32`, `end_viewer(&self)` (bump and wake).
- Requires macOS 14.4 (`os_sync_wait_on_address`).

### The viewer (`rtc-metal`)

- The editor starts it as `rtc-metal <preview> --window --at X,Y --share
  <file>`, without `--watch`. `--share` now means: write the camera, take
  load commands, and don't activate at start. `--watch` stays, unchanged,
  for running it alone.
- At start it takes the load block's `seq` as the baseline, so an old
  command is not replayed. A thread waits in `wait_load`; each `Load` goes
  to the main thread with `dispatch_async` on the main queue (`dispatch2`).
- On a `Load`:
  - the scene at `path` is read and built through the existing `reload`;
    `path` becomes the scene's path (the title, and the name Save Image and
    the quit alert offer);
  - `fresh` (New, Open): the scene's own camera becomes the view and its
    home even when it equals the old one, so the orbit resets;
  - not `fresh` (Load Edits): today's `--watch` reload, keeping the orbit
    unless the file's camera changed;
  - a scene that fails keeps the old image and shows why, as today;
  - the window is neither activated nor raised.
- Menus:

| Menu | Items (key) |
|---|---|
| rtc-metal | Hide rtc-metal (Cmd-H), Hide Others (Cmd-Opt-H), Show All, Quit (Cmd-Q) |
| File | Save Image… (Cmd-S) |
| Traversal, Samples | unchanged |

- Save Image… runs the existing save panel, then writes the image as it
  stands (rendering carries on); a failed write shows the existing alert.

### The editor (`rtc-edit`)

- `Editor.share: Option<Shared>`, made by the first Start Viewer.
- Start Viewer (Cmd-R): as today, but the args drop `--watch` and pass the
  editor's share file.
- New, Open (`switch`): the viewer is not dropped; following stops and the
  history clears as today. The old preview is dropped before the new one is
  made (they share the `rtc-edit-<pid>` folder). With a viewer running, the
  new preview is written and `send_load(preview, true)`; with none, nothing
  starts (today a switch starts one).
- Load Edits (Cmd-L): write the preview, `send_load(preview, false)`.
- Following: a subscription `Subscription::run_with((share path,
  generation), camera_stream)` while a viewer runs. Its thread opens the
  file and loops on `wait_camera`, yielding `Message::Camera(Snapshot)`,
  handled by `Follower::update` as today's Tick reads were. Tick no longer
  reads the camera.
- When a viewer ends (dropped, or found dead by Tick): `end_viewer()`.

## Error handling

- The shared file can't be made: the status line says why; no viewer.
- A viewer from another build: refuses the file by `SIGNATURE` and exits
  with that message, as today.
- A load path over 1024 bytes: `send_load` errs; the status line says so.
- `os_sync_wait_on_address` fails other than `EINTR`: the waiting thread
  ends. The viewer stops taking loads (stderr says why) but still writes
  its camera; the editor stops following and the status line says so.

## Testing

`rtc-share` (unit tests; two mappings of one temp file, threads for
processes):
- a camera write reads back;
- a sent load reads back with its path and `fresh`;
- `wait_load` sleeps until a send and returns it;
- `wait_camera` returns `None` once `end_viewer` is called;
- a reader never sees a half-written camera (a writer thread hammering it);
- a path over 1024 bytes is refused;
- a file of another signature is refused.

`rtc-metal`:
- a fresh load resets the orbit and home even with the same camera;
- a load that isn't fresh keeps the orbit;
- a load makes its path the scene's path.

`rtc-edit`:
- `viewer_args` has no `--watch` and passes the share file;
- the share path is `$TMPDIR/rtc-edit-<pid>.share`.

By eye: Cmd-N and Cmd-O keep the focus and the viewer shows the new scene
with a reset orbit; Cmd-L keeps the orbit; following tracks; closing the
viewer then Cmd-R works; Cmd-R opens the viewer behind the editor; Cmd-H,
Hide Others and Show All in the viewer; Save Image… writes the image;
`rtc-metal --window --watch scene.json` alone still activates and reloads.
