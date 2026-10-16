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
  - the Traversal menu switches between BVH (Cmd-1) and Metal (Cmd-2);
  - the Samples menu sets the target (presets, Unlimited, Cmd-] doubles, Cmd-[ halves);
  - quitting writes the image at the window's size.
- `--traversal bvh|metal` picks how rays search the bounded shapes: our own flattened BVH
  walk (the default) or Metal's ray tracing API.

Images are written as PPM files to the working directory (`-o` names the file). The scenes in
`scenes/` are rtc-rust's chapter scenes, plus soft-shadow, multiple-light and Cornell box
scenes.

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
