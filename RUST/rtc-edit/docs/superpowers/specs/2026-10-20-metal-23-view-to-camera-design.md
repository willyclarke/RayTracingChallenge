# Milestone 23: the viewer's view into the scene's camera

Branch `metal/23-view-to-camera`, merged into develop with `--no-ff` when
done. Milestones 17–21 built the scene editor; 22 gave the viewer's quit
alert a Cancel button.

## Goal

Orbiting and zooming in the viewer finds a good shot faster than typing
numbers, but the shot stays in the viewer. The editor should show the
viewer's camera live and, on request, follow it into the scene's camera
boxes, so the shot can be kept and saved.

Success:
- With the camera selected, the inspector shows the viewer's from, to and
  up beside the scene's, updating while the viewer is orbited.
- Ticking **follow** makes the `from` boxes track the viewer; **Keep**
  (or Cmd-K) keeps the followed value as an unsaved edit; unticking puts
  back the `from` the boxes had when follow was ticked.
- Cmd-S then saves the followed view; the next Cmd-L shows the same view
  in the viewer, without a jump.

## Out of scope

- Undo. The cached `from` is a one-level snapshot; an undo would keep a
  history of whole scenes and needs nothing from this milestone.
- Sending the editor's camera to the viewer live: Cmd-L stays the way
  edits reach the viewer.
- Following `to` or `up`: the viewer never changes them (it orbits around
  `to` and keeps `up`), so they are shown but not followed.
- Field of view, size, antialias and focal blur samples: the viewer
  doesn't change them.

## Decisions

### The channel: a shared, memory-mapped file

- A new crate, `RUST/rtc-share`, a path dependency of both rtc-metal and
  rtc-edit. One file, `src/lib.rs`, so everything that defines the layout
  is in it.
- The shared struct, `#[repr(C)]`, in a file mapped `MAP_SHARED` by both
  programs:

  ```rust
  pub struct View {
      signature: u64,     // SIGNATURE of the rtc-share the editor was built from
      size: u64,          // size_of::<View>()
      seq: u64,           // +1 on every write by the viewer; 0 = none yet
      from: [f64; 3],
      to: [f64; 3],
      up: [f64; 3],
      focal_distance: f64,
  }
  ```

- `pub const SIGNATURE: u64 = fnv1a(include_str!("lib.rs"))`: a
  compile-time FNV-1a hash of the crate's own source, a `const fn`. Any
  change to `lib.rs` changes it, so a viewer built from an older
  rtc-share is caught. FNV-1a, not SHA-256: it only has to notice a
  change, and needs no dependency or `build.rs`.
- Locking: `flock` on the file around every read and write, exclusive
  for writes, shared for reads. The kernel drops it when a process dies,
  so a crashed viewer can't hang the editor. Not real-time; a lock held
  for a copy of 104 bytes.
- `libc`-style `extern "C"` declarations (`mmap`, `munmap`, `flock`), as
  rtc-metal and rtc-edit already declare `signal`, `kill` and `atexit`;
  no new dependency.
- API:
  - `Shared::create(path) -> Result<Shared, String>`: the editor; makes
    (or replaces) the file, sized to `View`, header written, `seq` 0.
  - `Shared::open(path) -> Result<Shared, String>`: the viewer; maps an
    existing file and checks signature and size, else an `Err` saying to
    rebuild.
  - `write(&self, from, to, up, focal_distance)`: bumps `seq`.
  - `read(&self) -> Option<Camera>` (from, to, up, focal distance, seq);
    `None` while `seq` is 0 (a stale viewer never writes, so the editor
    sees `None` too).
  - Unmapped on drop.

### rtc-metal

- `--share <path>`, only with `--window` (a usage error otherwise, like
  `--watch`). It opens the file with `Shared::open`; on `Err` it prints
  the message and runs without sharing.
- It writes its camera — `Orbit::from()`, `to`, `up` and the current
  focal distance — at start-up, after every orbit or zoom, on space, and
  after a reload that changes the camera.
- Zooming already scales the focal distance with the orbit's radius (so
  the focus stays on what it was on); it is shared for that reason.

### rtc-edit

- **The file:** `<preview folder>/view`, created by `Shared::create` each
  time the editor starts a viewer (start-up, New, Open, Cmd-R), as a
  switch replaces the preview folder. `Viewer::start` passes
  `--share <path>`. The editor keeps the `Shared` for the current viewer.
- **Reading:** a timer subscription reads it 10 times a second; a new
  `seq` updates the shown values (and, following, the boxes). An
  unchanged `seq` changes nothing.
- **The camera's inspector** shows a viewer column beside `from`, `to`
  and `up`: the viewer's values, read-only and greyed, or one word in
  their place: `no viewer` (none running), `waiting` (`seq` still 0;
  also what a stale viewer shows — it prints why on its own terminal).
- **follow**, a checkbox beside `from`:
  - ticking caches the boxes' current `from` and focal distance;
  - while ticked, each new `seq` writes the viewer's `from` into the
    scene and the `from` boxes, and its focal distance when the scene
    has focal blur (the f-number is kept, as when the focal distance box
    is edited). These are unsaved edits (` •`). The `from` boxes are
    read-only while following.
  - unticking restores the cached `from` and focal distance.
- **Keep** (and Cmd-K), shown while following: stop following, keep the
  current values.
- **Following stops without restoring** when the scene is replaced (New,
  Open) or the viewer ends; selecting another item doesn't stop it, it
  only hides the controls.
- **Cmd-L** is unchanged: while following, the viewer already shows the
  view; the next Cmd-L reloads it with the same `from`, so the orbit
  doesn't jump.
- **Cmd-R with a viewer open** says `viewer already running` in the
  status line instead of nothing.

## Error handling

- `Shared::create` fails (temp folder unwritable): the viewer starts
  without `--share`, the column shows `no viewer`, the status line says
  why.
- A stale viewer (signature or size mismatch): it prints `rebuild
  rtc-metal` and runs without sharing; the editor shows `waiting`.
- The viewer crashes: its lock goes with it; the editor shows the last
  values until it notices the viewer ended, then `no viewer`.

## Testing

rtc-share (unit tests, no window):
- `create` then `open` agree; `read` is `None` until the first `write`;
  a `write` through one mapping is read through another, with `seq` 1.
- `open` refuses a file with another signature, a short file, and a
  missing one.
- `SIGNATURE` is the FNV-1a of the source (a known-answer test of `fnv1a`
  on a short string).

rtc-metal:
- `--share` without `--window` is a usage error.
- The camera written for an orbit: `Orbit::from()`, `to`, `up` (the pure
  part, without a window).

rtc-edit:
- `viewer_args` includes `--share <path>`.
- Following, keeping and unticking as pure state changes: tick caches,
  a new `seq` writes `from` (and focal distance only with focal blur),
  the same `seq` writes nothing, untick restores, Keep keeps, a switch
  stops without restoring.

By eye:
- orbit and zoom: the column follows; tick follow: the boxes follow;
  Keep, Cmd-S, reopen: the saved view; untick: the old view back;
- Cmd-L after Keep: no jump;
- focal blur scene: zoom while following keeps the focus;
- quit the viewer while following; Cmd-R with a viewer open.

## Risks

- **Timer cost:** a 10 Hz iced subscription redraws the editor only on a
  new `seq`; to be checked by eye that the editor stays idle otherwise.
- **macOS `flock` on the temp folder:** local APFS, supported.
