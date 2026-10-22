# rtc-metal

A Metal GPU port of [rtc-rust](../rtc-rust/README.md), the Rust implementation of
_The Ray Tracer Challenge_ (Jamis Buck). macOS only. It renders the same JSON scene
files on the GPU, and its tests compare the images with rtc-rust's CPU renders, which
stay the reference: pixel for pixel for deterministic renders, statistically for
sampled and path-traced ones.

## Dependencies

- [`rtc-rust`](../rtc-rust) (path dependency): loads scenes into a `World` and `Camera`,
  and writes PPM files.
- `objc2`, `objc2-foundation`: calling Objective-C APIs from Rust.
- `objc2-metal`: devices, buffers, compute pipelines, and the ray tracing API
  (acceleration structures).
- `objc2-metal-kit`, `objc2-quartz-core`, `objc2-app-kit`: the `--window` view, its
  window, menus and status overlay.

The Metal shaders in `shaders/` are embedded with `include_str!` and compiled when the
program starts, so there is no separate shader build step. No Cargo workspace: the crate
sits next to rtc-rust and has its own `target/`.

## Building and running

```bash
cargo build --release
cargo test --release --lib       # unit, probe and image-compare tests (needs the GPU)
cargo clippy --release

cargo run --release -- scenes/chapter15.json                    # → chapter15.ppm
cargo run --release -- scenes/soft_shadows.json --samples 256   # average 256 sampled frames
cargo run --release -- scenes/cornell_three_lights.json --window
cargo run --release -- my_scene.json --window --watch            # reload when the file changes
cargo run --release -- scenes/chapter10.json --window --at 900,0 # window's top-left at x 900
cargo run --release -- scenes/chapter15.json --traversal metal  # Metal's ray tracing API
cargo run --release -- --help
```

- Without `--samples`, a scene renders deterministically, like rtc-rust.
- `--samples N` averages N frames with random light samples and jittered pixels, for soft
  shadows and anti-aliasing.
- `--window` shows the sampled frames refining live:
  - drag to orbit, scroll to zoom, space to return to the scene's view;
  - `i` toggles the status overlay (samples, a stopwatch to the target, frame time, size,
    traversal);
  - resize or go full screen to render at the new size;
  - the app menu has Hide rtc-metal (Cmd-H), Hide Others (Cmd-Opt-H), Show All and Quit (Cmd-Q);
  - File ▸ Save Image… (Cmd-S) writes the image as it stands where you pick, while it keeps
    refining;
  - the Traversal menu switches between BVH (Cmd-1) and Metal (Cmd-2);
  - the Samples menu sets the target (presets, Unlimited, Cmd-] doubles, Cmd-[ halves);
  - quitting (Cmd-Q, closing the window, or SIGTERM) asks whether to save the image at the
    window's size and, if so, where (a save panel offering `<scene>.ppm` in the working
    directory); with `-o` it writes there without asking. Cancel (or Esc) in the alert keeps
    the window open, and cancelling the save panel asks again; a SIGTERM (the scene editor
    closing its viewer) offers no Cancel;
  - with `--watch`, the window reloads the scene whenever its file changes and refines it
    again, keeping your view unless the camera changed in more than width and height. A new width or
    height in the file resizes the window, keeping its top-left corner and fitted to the room
    from there to the screen's right and bottom edges (the whole screen when under 200 points
    are left either way, moving the window back onto the screen if it no longer fits), and
    keeps your view. A file that doesn't load
    keeps the old image, and the overlay shows why. The scene editor
    [rtc-edit](../rtc-edit/README.md) uses this as its viewer.
  - with `--at X,Y`, the window's top-left goes at that point, in points right and down from
    the primary display's top-left (moved below the menu bar if needed), and the image is
    fitted to the room right of and below it on that screen; with less than 200 points of
    room either way, or off every screen, the window is centred as without `--at`. rtc-edit uses this to put its viewer
    beside itself.
  - with `--share PATH`, the window writes its camera (from, to, up and focal distance) to
    `PATH`, a file the scene editor made with [rtc-share](../rtc-share/README.md), at start-up
    and whenever it changes: orbiting, zooming, space, or a reload with a new camera. rtc-edit
    shows it and can follow it into the scene. A file made by another version of rtc-share is
    refused with `rebuild rtc-metal`, and the window runs without sharing.
    Through the same file the editor sends scenes to load (New and Open from the scene's own
    camera, Load Edits keeping the view), which a thread sleeping in
    `os_sync_wait_on_address` hands to the main queue; the window is never raised for them.
    With `--share` the viewer doesn't take the focus at start, so the editor keeps the
    keyboard. rtc-edit no longer passes `--watch`.
- `--traversal bvh|metal` picks how rays search the bounded shapes: our own flattened BVH
  walk (the default) or Metal's ray tracing API.

Without `--window`, images are written as PPM files to the working directory (`-o` names the
file). The scenes in `scenes/` are rtc-rust's chapter scenes, plus soft-shadow,
multiple-light and Cornell box scenes.

## What the GPU supports

- **Shapes:** spheres, planes, cubes, cylinders, cones, groups, OBJ triangle meshes (smooth
  or flat) and CSG.
- **Materials:** Phong shading, reflection and refraction, and stripe, gradient, ring,
  checkers and perturbed patterns.
- **Lights:** point and area lights, up to 8 per scene.
- **Rendering:** focal blur, and rtc-rust's path-traced indirect light (chapter 17).

Not yet supported: spot lights, motion blur, bump maps, texture maps and tori. Scenes that
use them are rejected with a message naming the feature.

## Layout

| path | what it holds |
|---|---|
| `src/gpu.rs` | the Metal device, pipelines, one-shot render, and the `Accumulator` for sampled frames |
| `src/scene.rs`, `src/types.rs` | rtc-rust's `World` flattened into GPU buffers; the `#[repr(C)]` types shared with the shaders |
| `src/bvh.rs`, `src/accel.rs` | the flattened BVH, and Metal acceleration structures |
| `src/window.rs`, `src/orbit.rs` | the `--window` view and the orbiting camera |
| `src/probes.rs`, `src/compare.rs` | tests: per-function probe kernels, and whole-image comparisons with the CPU |
| `shaders/` | the Metal shading language sources |
| `docs/superpowers/` | design specs and implementation plans, one per milestone |

The port was built in milestones, each on its own `metal/NN-name` branch merged into
`develop`; the specs in `docs/superpowers/specs/` record the decisions made along the way.
