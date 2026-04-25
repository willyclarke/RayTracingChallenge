# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Context

Rust implementation of "The Ray Tracer Challenge" book. Each chapter builds on the last — currently at Chapter 4 (matrix transformations). No external dependencies; Rust edition 2024.

## Commands

```bash
cargo build              # debug build
cargo build --release
cargo test --lib         # run all unit tests
cargo test --lib --ignored  # run only ignored tests (projectile demos)
cargo test --lib test_chap_3_20  # run a single test by name
cargo clippy --lib       # lint
cargo fmt                # auto-format
```

## Architecture

The crate is a library (`src/lib.rs`) with six modules. `src/main.rs` is a placeholder.

**Module dependency order:**
```
math    → tuple → matrix
                → canvas
log / color     (standalone utilities)
```

**`tuple`** — The core type. `Tuple { x, y, z, w: f64 }` is used for everything:
- Points: `w = 1.0`, Vectors: `w = 0.0`, Colors: RGB stored in x/y/z
- All arithmetic operators overloaded; float equality via `approx_eq()` (ε = 1e-9)

**`matrix`** — `Matrix2`, `Matrix3`, `Matrix4` as fixed `[[f64; N]; N]` arrays. `Matrix4` is the main type: supports `inverse()`, `transpose()`, `translation()`, identity, and `Mul<Tuple>` for transforming points/vectors. Display renders with ANSI color (0=yellow, 1=green, negatives=red).

**`canvas`** — 2D pixel buffer (`Vec<Tuple>`). Writes to binary PPM via `to_ppm()` / `write_ppm(path)`. Indexed by `canvas[(x, y)]`. PPM test images are written to the working directory.

**`math`** — Single `approx_eq(a, b)` helper used throughout for float comparison.

**`log` / `color`** — `logi!()`, `logd!()`, `loge!()` macros with timestamps; `Color` enum for ANSI codes. Used in tests and demos for visibility.

## Tests

All tests live inline at the bottom of each module under `#[cfg(test)]`. Named by chapter: `test_chap_1_05`, `test_chap_3_20`, etc. Return `Result<(), String>` (or `std::io::Result<()>` for I/O tests). `#[ignore]` marks the two projectile trajectory demos that generate PPM output.
