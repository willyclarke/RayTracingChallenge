# rtc-share

The viewer's camera and the editor's load commands, shared between
[rtc-metal](../rtc-metal/README.md) and [rtc-edit](../rtc-edit/README.md) through a small
memory-mapped file, so the editor can follow the viewer's view into the scene's camera and
tell the viewer which scene to show.

- The editor makes the file (`Shared::create`) and starts the viewer with
  `--share <path>`; the viewer maps it (`Shared::open`).
- The viewer writes its camera (`write_camera`) whenever it changes; the editor sleeps in
  `wait_camera` until it does. The editor sends a scene (`send_load`, `fresh` for a new
  scene); the viewer sleeps in `wait_load` until one comes. `end_viewer` moves the viewer
  generation on, so the editor's camera waiter for an ended viewer returns `None`.
- The viewer calls `viewer_ready` once it takes loads; until `is_ready` says so for the
  current generation, the editor holds its loads, as one sent earlier may be missed.
- The layout, `#[repr(C)]`: `signature`, `size`, the camera word (seq in the low half,
  viewer generation in the high half), `from`, `to`, `up` and the focal distance (`f64`),
  then the load's `seq`, `fresh`, `len`, a 1024-byte path, and `ready` (one more than the
  generation of the viewer that last said it takes loads).
- Each block has one writer and is a seqlock: its `seq` is odd while it is written, and a
  reader retries while `seq` is odd or changed. After a write the writer calls
  `os_sync_wake_by_address_all`; the reader sleeps in `os_sync_wait_on_address` (both with
  the shared flag, macOS 14.4+) until the `seq` it saw changes. No locks, so a crashed
  viewer can't hang the editor.
- `SIGNATURE` is an FNV-1a hash of `src/lib.rs`, computed while compiling. `open` refuses a
  file with another signature or size, so a viewer built from another version of this crate
  says `rebuild rtc-metal` instead of misreading the file. Any change to `lib.rs`, a comment
  included, changes it; keep the crate to that one file so the hash covers the whole layout.

No dependencies: `mmap`, `munmap`, `os_sync_wait_on_address` and
`os_sync_wake_by_address_all` are declared directly. `cargo test` needs no window.
