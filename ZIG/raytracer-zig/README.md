# raytracer-zig

Zig implementation of [The Ray Tracer Challenge](http://raytracerchallenge.com/)
(Jamis Buck), built chapter by chapter alongside the C++ and Rust versions in this
repository (see the repo root `Readme.adoc`).

**Status: work in progress** — currently around chapter 10: tuples, matrices, canvas,
rays, camera and world with point lights, sphere and plane shapes, and the beginnings
of patterns (stripes). The `chapNputtingtogether.ppm` files in this directory are the
book's "putting it together" renders produced so far.

## Building and testing

```bash
zig build run      # build and run the app (renders to PPM in this directory)
zig build test     # run the unit tests
```

Sources live in `src/` (one module per book concept: `matrix.zig`, `canvas.zig`,
`world.zig`, `shapes/`, `patterns/`, `lights/`).
