# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Context

Rust implementation of "The Ray Tracer Challenge" book, built up chapter by chapter. Rust edition 2024. The one external dependency is `rayon`, used for parallel rendering.

## Commands

```bash
cargo build              # debug build
cargo build --release
cargo test --lib         # run all unit tests
cargo test --lib --ignored   # run only ignored tests (projectile demos, benchmarks)
cargo test --release --lib cornell_box -- --ignored --nocapture   # Cornell box benchmarks (point light, area light, area light + AA, focal blur), log render time
cargo test --lib test_chap_13 # run one chapter's tests by name prefix
cargo test --doc         # run documentation examples (doctests)
cargo clippy --lib       # lint
cargo fmt                # auto-format
```

## Architecture

The crate is a library (`src/lib.rs`); `src/main.rs` is a placeholder. Rendering flows: a `Camera` casts rays through a `World` of shapes lit by a `Light`, and `World::color_at` returns the shaded color (with reflection/refraction recursion).

**Foundation types:**

**`tuple`** — The core type. `Tuple { x, y, z, w: f64 }` is used for everything:
- Points: `w = 1.0`, Vectors: `w = 0.0`, Colors: RGB stored in x/y/z
- All arithmetic operators overloaded; float equality via `approx_eq()` (ε = 1e-9)

**`matrix`** — `Matrix2`, `Matrix3`, `Matrix4` as fixed `[[f64; N]; N]` arrays. `Matrix4` is the main type: `inverse()`, `transpose()`, `translation()`/`scaling()`/`rotation_{x,y,z}()`, identity, and `Mul<Tuple>` for transforming points/vectors.

**`canvas`** — 2D pixel buffer (`Vec<Tuple>`). Writes binary PPM via `to_ppm()` / `write_ppm(path)`; reads PPM via `from_ppm` (P3), `from_ppm_binary` (P6) and `read_ppm(path)`. Indexed by `canvas[(x, y)]`. PPM images are written to the working directory.

**`math`** — `approx_eq(a, b)` and `EPSILON`, used throughout for float comparison.

**`noise`** — `Perlin` improved noise (`noise`, `vector`, `octaves`), seeded and deterministic. Used by `patterns::perturbedpattern` (jitters any pattern's lookup point — marble) and `Material.bump` (tilts the object-space normal in `prepare_computations` — bumpy surfaces).

**Rendering pipeline:**

**`shape` / `shapes`** — `Shape` is the object-safe trait every primitive implements; shared state (id, transform, material) lives in `ShapeData`, exposed via `data()`/`data_mut()`. Each primitive implements `local_intersect` and `local_normal_at` in object space; the trait handles the world↔object transform. Primitives: `sphere`, `plane`, `cube`, `cylinder`, `cone`.

**`pattern` / `patterns`** — `Pattern` trait for material surface patterns: `stripe`, `gradient`, `ring`, `checkers`, plus `nested`/`blended` combinators and a `test` pattern for unit tests. Texture mapping (bonus chapter): `uvpattern` (`UvPattern` trait — `UvCheckers`, `UvAlignCheck`, `UvImage`) and `texturemap` (`TextureMap` with spherical/planar/cylindrical `UvMap`, `CubeMap` with one UV pattern per face).

**`world`** — Holds the shapes and light. `intersect`, `is_shadowed` (allocation-free any-hit walk via `Shape::local_occludes`), `color_at`, `shade_hit`, `reflected_color`, `refracted_color`, `prepare_computations` (builds `Computations`, including `n1`/`n2` for refraction), plus `render`/`render_parallel` and `view_transform`. Shapes get a world-assigned id via `add_shape`.

**`intersection`** — `Intersection { t, object_id }` and `Intersections`, a `t`-sorted collection (`push` inserts in order; `hit()` returns the first non-negative).

**`ray`, `camera`, `light`, `material`** — `Ray` (origin/direction); `Camera` (view rays via `ray_for_pixel`/`ray_for_subpixel`; `with_antialias(n)` enables edge-detected n×n supersampling and `with_focal_blur(aperture, focal_distance, samples)` depth of field via `ray_for_lens`, both in `render_parallel`); `Light` (jittered rectangular area light — `point_light` is the 1×1 case — with `intensity_at` for soft shadows and Phong `lighting` averaged over the sample points; `Sequence` is the jitter generator); `Material` (color, ambient/diffuse/specular/shininess, reflective, transparency, refractive_index, optional pattern, optional `Bump`).

**`log` / `color`** — `logi!()`, `logd!()`, `loge!()` macros with timestamps; `Color` enum for ANSI codes. Used in tests and demos.

## Tests

Tests live inline at the bottom of each module under `#[cfg(test)]`, named by book chapter (`test_chap_1_05`, `test_chap_13_9`, etc.) so a chapter's tests share a `test_chap_N` prefix. They return `Result<(), String>` (or `std::io::Result<()>` for I/O). `#[ignore]` marks the projectile trajectory demos and the Cornell box benchmarks. Bonus-chapter tests use a `test_bonus_*` prefix (`test_bonus_soft_shadows_N`, `test_bonus_texture_N`). Some `*_putting_it_all_together` tests render a scene to a PPM. Public helpers additionally carry doctests (run with `cargo test --doc`).
