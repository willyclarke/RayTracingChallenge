# rtc-share

The viewer's camera, shared between [rtc-metal](../rtc-metal/README.md) and
[rtc-edit](../rtc-edit/README.md) through a small memory-mapped file, so the editor can show
the viewer's view and follow it into the scene's camera.

- The editor makes the file (`Shared::create`) and starts the viewer with
  `--share <path>`; the viewer maps it (`Shared::open`) and writes its camera
  (`write`) whenever it changes; the editor reads it (`read`) when it likes.
- The layout, `#[repr(C)]`: `signature`, `size`, `seq` (+1 per write, 0 before the first),
  then `from`, `to`, `up` and the focal distance, all `f64`; 104 bytes.
- `SIGNATURE` is an FNV-1a hash of `src/lib.rs`, computed while compiling. `open` refuses a
  file with another signature or size, so a viewer built from another version of this crate
  says `rebuild rtc-metal` instead of misreading the file. Any change to `lib.rs`, a comment
  included, changes it; keep the crate to that one file so the hash covers the whole layout.
- Every read and write holds an `flock` on the file (shared for reads, exclusive for
  writes). The kernel releases it when a process dies, so a crashed viewer can't hang the
  editor.

No dependencies: `mmap`, `munmap` and `flock` are declared directly. `cargo test` needs no
window.
