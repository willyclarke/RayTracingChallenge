# rtc-rust

Rust implementation of [The Ray Tracer Challenge](http://raytracerchallenge.com/) (Jamis Buck),
built chapter by chapter with the book's test-driven approach — one of three implementations
in this repository, alongside C++ and Zig (see the repo root `Readme.adoc`). All 16 core chapters, the three
online bonus chapters (soft shadows, bounding boxes/BVH, texture mapping), and all eight
chapter 17 "next steps" (area lights, spotlights, focal blur, motion blur, anti-aliasing,
texture maps, normal perturbation, torus) are implemented, plus path-traced indirect lighting.

The crate is a library plus the `rtc` scene-rendering binary. External dependencies: `rayon`
(parallel rendering) and `serde`/`serde_json` (scene loading, confined to the `scene` module).

## Building and testing

```bash
cargo build --release
cargo test --lib                 # all unit tests (named test_chap_N after book chapters)
cargo test --doc                 # doctests
cargo test --release --lib cornell_box -- --ignored --nocapture --test-threads=1  # benchmarks
```

Renders are written as binary PPM files to the working directory.

## JSON scene format

The `rtc` binary reads a JSON scene description and renders it (load-only — the renderer
never writes scenes back):

```bash
cargo run --release -- scenes/cover.json            # → cover.ppm
cargo run --release -- scene.json -o render.ppm
cargo run --release -- scene.json --orbit 240 -o frames/scene.ppm   # film: 240 numbered frames
cargo run --release -- help                         # man-style scene-format reference
cargo run --release -- help shapes                  # one topic, with copyable JSON examples
```

`rtc help` covers the whole format (topics: scene, camera, lights, transforms, shapes,
materials, patterns, settings; `rtc help all` prints everything), so the sections below
are also available from the binary itself.

`scenes/cover.json` — the book's cover image (Appendix A1), translated from the appendix's
YAML — is a complete example. The dice scenes are larger ones: three marbled dice built from
nested CSG (rounded cube minus 21 pip spheres, materials on the CSG leaves) in a checkered
room, with perturbed-stripe patterns and adaptive anti-aliasing. Their Python generators
share `scenes/dicelib.py` (die/room/camera/light builders plus a 3×5 dot-matrix dice font):

```bash
python3 scenes/dice-light-area.py            # area light, soft shadows
python3 scenes/dice-light-spot.py            # spotlight aimed at the stack
python3 scenes/dice-sentence.py "HELLO"      # text spelled in small dice, colors per letter
cargo run --release -- scenes/dice-sentence.json
```

Each writes the JSON of the same name (optional trailing `width height` arguments;
`dice-sentence.py` auto-frames the camera to the text length and supports A–Z, 0–9 and
basic punctuation). The loader lives in `src/scene.rs`
(`scene::load(path) -> (World, Camera)`), with field-level errors and `test_scene_*` tests.

### Filming an orbit

`--orbit <frames>` renders a film instead of a still: the camera circles the look-at point
once, rotating about the camera's `up` axis at constant height and distance, and writes
`<stem>_0000.ppm` … `<stem>_<frames-1>.ppm` next to the `-o` path. The last frame stops one
step short of 360° so the sequence loops seamlessly. The scene is parsed and built once;
only the camera transform changes per frame. `rtc` prints the ffmpeg command that merges the
frames, and `film.sh` runs both steps:

```bash
./film.sh scenes/dice-light-spot.json            # 240 frames at 24 fps → 10 s dice-light-spot.mp4
./film.sh scenes/small.json frames=48 fps=12     # both optional; bare numbers work too
ffmpeg -framerate 24 -i frames/small/small_%04d.ppm -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" \
    -c:v libx264 -pix_fmt yuv420p small.mp4      # what film.sh runs after rendering
```

Render time scales with the frame count, so preview with a copy of the scene at a small
`width`/`height` (or a few frames) before the full run.

`screensaver/` turns a film into a macOS screen saver: `RtcSaverView.swift` is a
`ScreenSaverView` that loops the bundled MP4 with `AVPlayerLooper` (muted, aspect-fill, one
instance per display), and `build.sh` compiles it with `swiftc`, packs the film into
`build/RtcSaver.saver`, ad-hoc signs it, and optionally installs it:

```bash
./film.sh scenes/cover.json frames=480 fps=48     # render at the display's resolution first
./screensaver/build.sh cover.mp4 install          # -> ~/Library/Screen Savers/RtcSaver.saver
```

Then pick RtcSaver under "Other" in the screen saver section of System Settings > Wallpaper
(macOS 26 folded the Screen Saver pane into Wallpaper). The settings tile shows the film's
first frame, which `build.sh` extracts with ffmpeg as `thumbnail.png`/`thumbnail@2x.png`.
Needs Xcode's command-line tools.
Third-party savers run in Apple's sandboxed legacy host, so the film must live inside the
bundle; rebuild and reinstall after re-rendering. The System Settings preview thumbnail may
stay black for video savers; the full-screen saver is unaffected. Open scenes (`small.json`,
`cover.json`) orbit cleanly; the dice scenes are closed rooms, so most of the turn looks at
the walls from outside — move the walls out or drop them to film those.

Design rules:

- **Load-only.** JSON describes a scene to build; the renderer never writes scenes back.
- **A separate scene layer.** Plain serde-derived structs/enums (`scene` module) mirror the
  core types and `build()` a `World` + `Camera`. Core types stay serde-free; trait objects
  (`Shape`, `Pattern`) are dispatched by enums with a `"type"` tag.
- **Human-writable.** Transforms are op lists, not raw matrices. Bulk data stays in files:
  meshes are referenced as `.obj` paths, image textures as `.ppm` paths.
- **Computed state is rebuilt, not stored.** World ids, bounds, and the BVH are derived at
  load (`divide` runs when `bvh_threshold` is set). Perlin noise is seeded and deterministic,
  so a seed is its complete state.

### Top level

```json
{
  "camera": { ... },
  "light": { ... },
  "settings": { ... },
  "shapes": [ ... ]
}
```

One light per scene, matching `World`. `settings` is optional. An optional top-level
`"comment"` field (any JSON value) is ignored by the loader; unknown fields elsewhere are
rejected, so typos fail loudly.

### Conventions

| JSON | Meaning |
|---|---|
| `[x, y, z]` | point or vector (which one is clear from the field) |
| `[r, g, b]` | color, 0.0–1.0 floats |
| angles | radians |

### Camera

```json
{
  "width": 1000,
  "height": 1000,
  "field_of_view": 1.0472,
  "from": [0, 1.5, -5],
  "to": [0, 1, 0],
  "up": [0, 1, 0],
  "antialias": { "n": 3 },
  "focal_blur": { "aperture": 0.1, "focal_distance": 5.0, "samples": 32 },
  "motion_blur": { "samples": 16 }
}
```

`antialias`, `focal_blur`, and `motion_blur` are optional and map to the
`with_antialias` / `with_focal_blur` / `with_motion_blur` builders.

### Light

Tagged by `"type"` — the three `Light` constructors:

```json
{ "type": "point", "position": [-10, 10, -10], "intensity": [1, 1, 1] }
```

```json
{
  "type": "area",
  "corner": [-1, 2, 4], "uvec": [2, 0, 0], "usteps": 4,
  "vvec": [0, 2, 0], "vsteps": 4, "intensity": [1.5, 1.5, 1.5]
}
```

```json
{
  "type": "spot",
  "position": [0, 5, 0], "target": [0, 0, 0], "intensity": [1, 1, 1],
  "inner_angle": 0.3, "outer_angle": 0.5
}
```

### Settings

```json
{
  "path_tracing": { "samples": 64, "depth": 4 },
  "bvh_threshold": 4
}
```

Both optional: `path_tracing` calls `set_path_tracing`; `bvh_threshold` runs `divide` on
each top-level group after the scene is assembled.

### Transforms

An ordered list, applied left to right (i.e. listed order = the order they act on the object):

```json
"transform": [
  { "scale": [0.5, 0.5, 0.5] },
  { "rotate_y": 0.7854 },
  { "translate": [0, 1, 0] }
]
```

Ops: `translate`, `scale`, `rotate_x`, `rotate_y`, `rotate_z`, `shear` (6 numbers:
`[xy, xz, yx, yz, zx, zy]`), and `matrix` (4 rows of 4 numbers, row-major escape hatch).
Omitted = identity.

### Shapes

Tagged by `"type"`; all take optional `transform`, `material`, and `motion` (velocity vector
for motion blur, → `set_motion`). Primitive-specific fields:

```json
{ "type": "sphere" }
{ "type": "plane" }
{ "type": "cube" }
{ "type": "cylinder", "min": 0, "max": 1, "closed": true }
{ "type": "cone",     "min": -1, "max": 0, "closed": true }
{ "type": "torus",    "major_radius": 1.0, "minor_radius": 0.25 }
{ "type": "obj",      "file": "models/teapot.obj" }
{ "type": "group",    "children": [ ...shapes... ] }
{ "type": "csg",      "operation": "difference", "left": { ...shape... }, "right": { ...shape... } }
```

Cylinder/cone `min`/`max` default to ±infinity and `closed` to false. Groups and CSG nest
their children inline as JSON subtrees; world ids are assigned during `build()` in the usual
`add_shape` order. `obj` loads via `load_obj_file` and yields a group.

### Material

All fields optional, defaulting to `Material::default()`:

```json
{
  "color": [1, 0.9, 0.9],
  "ambient": 0.1, "diffuse": 0.9, "specular": 0.9, "shininess": 200.0,
  "reflective": 0.0, "transparency": 0.0, "refractive_index": 1.0,
  "pattern": { ... },
  "bump": { "seed": 7, "amplitude": 0.3, "frequency": 4.0 }
}
```

### Patterns

Tagged by `"type"`, each with an optional `transform`. The four book patterns take two colors:

```json
{ "type": "stripe",   "a": [1, 1, 1], "b": [0, 0, 0] }
{ "type": "gradient", "a": [1, 0, 0], "b": [0, 0, 1] }
{ "type": "ring",     "a": [1, 1, 1], "b": [0.5, 0.5, 0.5] }
{ "type": "checkers", "a": [1, 1, 1], "b": [0, 0, 0] }
```

Combinators nest patterns where colors would go:

```json
{ "type": "nested",    "a": { ... }, "b": { ... } }
{ "type": "blended",   "a": { ... }, "b": { ... } }
{ "type": "perturbed", "pattern": { ... }, "amplitude": 0.4, "frequency": 2.0, "seed": 3 }
```

Texture mapping (bonus chapter):

```json
{
  "type": "texture_map",
  "mapping": "spherical",
  "uv": { "type": "uv_checkers", "width": 16, "height": 8, "a": [0, 0.5, 0], "b": [1, 1, 1] }
}
```

`mapping`: `spherical` | `planar` | `cylindrical`. UV patterns: `uv_checkers`,
`uv_align_check` (five colors: `main`, `ul`, `ur`, `bl`, `br`), `uv_image` (`{"file": "x.ppm"}`).

```json
{
  "type": "cube_map",
  "left": { ...uv... }, "right": { ...uv... }, "front": { ...uv... },
  "back": { ...uv... }, "up": { ...uv... }, "down": { ...uv... }
}
```

### Example scene

```json
{
  "camera": {
    "width": 400, "height": 300, "field_of_view": 1.0472,
    "from": [0, 1.5, -5], "to": [0, 1, 0], "up": [0, 1, 0]
  },
  "light": { "type": "point", "position": [-10, 10, -10], "intensity": [1, 1, 1] },
  "shapes": [
    { "type": "plane", "material": { "pattern": { "type": "checkers", "a": [1, 1, 1], "b": [0.2, 0.2, 0.2] } } },
    {
      "type": "sphere",
      "transform": [ { "translate": [0, 1, 0] } ],
      "material": { "color": [0.1, 0.2, 0.8], "reflective": 0.3 }
    }
  ]
}
```

### Implementation notes

- serde stays confined to `src/scene.rs`: plain description structs/enums deserialize the
  JSON, then `build()` instantiates the real `World`/`Camera` — the core types are serde-free.
- `scene::load(path)` resolves relative asset paths (`.obj`, `.ppm`) against the scene
  file's directory.
- `src/main.rs` is the CLI: `rtc <scene.json> [-o out.ppm] [--orbit <frames>]`; without
  `-o` the output is `<scene-stem>.ppm` in the working directory. `--orbit` keeps the parsed
  `SceneDescription` (`from_path` + `build`) to read the camera's `from`/`to`/`up` and moves
  the eye with `world::orbit_from` per frame.
- Acceptance test: `test_scene_cover_putting_it_together` renders `scenes/cover.json`.
