//! Man-page-style help topics for the `rtc` binary (`rtc help [topic]`).
//!
//! Binary-only module (declared in `main.rs`, not part of the library). The
//! text mirrors what `scene.rs` actually accepts — when the loader grows a
//! field, add it here and in `README.md`.

/// `(name, one-line summary, text)` per topic, in `rtc help all` print order.
pub const TOPICS: &[(&str, &str, &str)] = &[
    (
        "scene",
        "file structure, conventions, a minimal example",
        SCENE,
    ),
    (
        "camera",
        "size, view point, anti-aliasing, focal and motion blur",
        CAMERA,
    ),
    ("lights", "point, area (soft shadows), spot", LIGHTS),
    (
        "transforms",
        "translate / scale / rotate / shear op lists",
        TRANSFORMS,
    ),
    (
        "shapes",
        "sphere, plane, cube, cylinder, cone, torus, obj, group, csg",
        SHAPES,
    ),
    (
        "materials",
        "Phong fields, reflection, refraction, bump mapping",
        MATERIALS,
    ),
    (
        "patterns",
        "two-color patterns, combinators, texture maps",
        PATTERNS,
    ),
    ("settings", "path tracing, BVH threshold", SETTINGS),
];

/// Print one topic (or the topic list when `topic` is None, or everything for
/// "all"). Returns false when the topic is unknown.
pub fn print(topic: Option<&str>) -> bool {
    match topic {
        None => {
            println!("usage: rtc help <topic>\n");
            println!("Topics (also: 'all' for the full reference):\n");
            for (name, summary, _) in TOPICS {
                println!("    {name:<12} {summary}");
            }
            println!("\nComplete example scenes live in scenes/ (cover.json, dice-light-*.json).");
            true
        }
        Some("all") => {
            for (i, (_, _, text)) in TOPICS.iter().enumerate() {
                if i > 0 {
                    println!();
                }
                println!("{}", text.trim_end());
            }
            true
        }
        Some(name) => {
            // accept singular/plural sloppiness: "shape" finds "shapes"
            let hit = TOPICS.iter().find(|(n, _, _)| {
                *n == name || n.trim_end_matches('s') == name.trim_end_matches('s')
            });
            match hit {
                Some((_, _, text)) => {
                    println!("{}", text.trim_end());
                    true
                }
                None => {
                    eprintln!("rtc: no help topic '{name}'");
                    eprintln!(
                        "topics: {} (or 'all')",
                        TOPICS
                            .iter()
                            .map(|(n, _, _)| *n)
                            .collect::<Vec<_>>()
                            .join(", ")
                    );
                    false
                }
            }
        }
    }
}

const SCENE: &str = r#"SCENE
    A scene file is one JSON object with a camera, exactly one light,
    optional render settings, and a list of shapes. Rendered with:

        rtc scene.json [-o out.ppm]

    Unknown fields are rejected (typos fail loudly). A top-level "comment"
    field (any JSON value) is ignored. Relative file paths in the scene
    (.obj models, .ppm textures) resolve against the scene file's directory.

    Conventions: [x, y, z] arrays are points or vectors depending on the
    field; colors are [r, g, b] with 0.0-1.0 floats; angles are radians.

EXAMPLE
    {
      "comment": "smallest useful scene",
      "camera": {
        "width": 400, "height": 300, "field_of_view": 1.047,
        "from": [0, 1.5, -5], "to": [0, 1, 0], "up": [0, 1, 0]
      },
      "light": { "type": "point", "position": [-10, 10, -10], "intensity": [1, 1, 1] },
      "shapes": [
        { "type": "plane" },
        { "type": "sphere",
          "transform": [ { "translate": [0, 1, 0] } ],
          "material": { "color": [0.1, 0.2, 0.8] } }
      ]
    }

SEE ALSO
    rtc help camera | lights | transforms | shapes | materials | patterns | settings
"#;

const CAMERA: &str = r#"CAMERA
    View point, image size, and per-pixel effects.

    width, height        image size in pixels
    field_of_view        horizontal FOV in radians
    from, to, up         eye position, look-at point, up vector
    antialias            optional: { "n": 3 } - edge-detected n x n supersampling
    focal_blur           optional: { "aperture": 0.1, "focal_distance": 5.0, "samples": 32 }
                         depth of field; objects at focal_distance stay sharp
    motion_blur          optional: { "samples": 16 } - jittered shutter times;
                         pairs with a shape's "motion" velocity vector

EXAMPLE
    "camera": {
      "width": 1000, "height": 1000, "field_of_view": 1.047,
      "from": [-5, 4, -8], "to": [0, 1, 0], "up": [0, 1, 0],
      "antialias": { "n": 3 },
      "focal_blur": { "aperture": 0.1, "focal_distance": 8.0, "samples": 32 }
    }
"#;

const LIGHTS: &str = r#"LIGHTS
    Exactly one light per scene, tagged by "type".

    point       position, intensity - hard shadows
    area        jittered rectangle - soft shadows; corner is one corner,
                uvec/vvec its edges, split into usteps x vsteps sample cells
    spot        point light restricted to a cone: full intensity inside
                inner_angle, fading to nothing at outer_angle (radians,
                measured from the axis toward "target")

EXAMPLES
    "light": { "type": "point", "position": [-10, 10, -10], "intensity": [1, 1, 1] }

    "light": {
      "type": "area",
      "corner": [-8, 12, -10], "uvec": [4, 0, 0], "usteps": 4,
      "vvec": [0, 0, 4], "vsteps": 4, "intensity": [1, 1, 1]
    }

    "light": {
      "type": "spot",
      "position": [-3, 10, -4], "target": [0, 2, 0],
      "intensity": [1.2, 1.2, 1.2],
      "inner_angle": 0.22, "outer_angle": 0.42
    }
"#;

const TRANSFORMS: &str = r#"TRANSFORMS
    An ordered list of ops; listed order is the order they act on the object
    ([scale, translate] scales first, then moves). Omitted = identity.
    Shapes, patterns, and OBJ models all take a "transform".

    { "translate": [x, y, z] }
    { "scale":     [x, y, z] }
    { "rotate_x":  r }              also rotate_y, rotate_z (radians)
    { "shear":     [xy, xz, yx, yz, zx, zy] }
    { "matrix":    [[..4],[..4],[..4],[..4]] }   row-major escape hatch

EXAMPLE
    "transform": [
      { "scale": [0.5, 0.5, 0.5] },
      { "rotate_y": 0.7854 },
      { "translate": [0, 1, 0] }
    ]
"#;

const SHAPES: &str = r#"SHAPES
    Tagged by "type". Every shape also takes optional "transform",
    "material", and "motion" ([x, y, z] velocity for motion blur).

    sphere      unit sphere at the origin
    plane       infinite xz-plane through the origin
    cube        axis-aligned, -1..1 on each axis
    cylinder    y-axis; optional min, max (default +-infinity), closed (false)
    cone        double cone along y; same min / max / closed fields
    torus       xz-plane ring; major_radius (1.0), minor_radius (0.25)
    obj         { "file": "models/teapot.obj" } - Wavefront OBJ, loads as a group
    group       { "children": [ ...shapes... ] } - transforms apply to the subtree
    csg         { "operation": "union" | "intersection" | "difference",
                  "left": {shape}, "right": {shape} }
                a material on a group or csg node recolors its whole subtree

EXAMPLES
    { "type": "cylinder", "min": 0, "max": 2, "closed": true,
      "material": { "color": [0.9, 0.6, 0.1] } }

    { "type": "csg", "operation": "difference",
      "left":  { "type": "cube" },
      "right": { "type": "sphere", "transform": [ { "scale": [1.3, 1.3, 1.3] } ] } }

    { "type": "group",
      "transform": [ { "translate": [0, 1, 0] } ],
      "children": [ { "type": "sphere" }, { "type": "obj", "file": "models/teddy.obj" } ] }
"#;

const MATERIALS: &str = r#"MATERIALS
    All fields optional; defaults in parentheses.

    color              [r, g, b] (white); ignored where a pattern applies
    ambient            (0.1)   diffuse    (0.9)   specular   (0.9)
    shininess          (200)   reflective (0.0)
    transparency       (0.0)   refractive_index (1.0)  - glass is ~1.5
    pattern            optional surface pattern (rtc help patterns)
    bump               optional Perlin normal perturbation:
                       { "amplitude": 0.3, "frequency": 4.0, "seed": 7 }

EXAMPLE
    "material": {
      "color": [0.373, 0.404, 0.55],
      "diffuse": 0.2, "ambient": 0, "specular": 1, "shininess": 200,
      "reflective": 0.7, "transparency": 0.7, "refractive_index": 1.5
    }
"#;

const PATTERNS: &str = r#"PATTERNS
    Tagged by "type"; every pattern takes an optional "transform".

    Two-color patterns (a, b):
    stripe      alternates along x        gradient    blends a -> b along x
    ring        concentric in xz          checkers    3-D checkerboard

    Combinators (nest patterns where colors would go):
    nested      { "a": {pattern}, "b": {pattern} }
    blended     { "a": {pattern}, "b": {pattern} } - averages the two
    perturbed   { "pattern": {pattern}, "amplitude": 0.4, "frequency": 2.0,
                  "seed": 3 } - Perlin-jitters the lookup point (marble)

    Texture mapping:
    texture_map { "mapping": "spherical" | "planar" | "cylindrical", "uv": {uv} }
    cube_map    { "left": {uv}, "front": {uv}, "right": {uv},
                  "back": {uv}, "up": {uv}, "down": {uv} }

    UV patterns (used inside texture_map / cube_map):
    uv_checkers    { "width": 16, "height": 8, "a": [..], "b": [..] }
    uv_align_check { "main": [..], "ul": [..], "ur": [..], "bl": [..], "br": [..] }
    uv_image       { "file": "texture.ppm" }

EXAMPLES
    "pattern": { "type": "checkers", "a": [1, 1, 1], "b": [0.14, 0.13, 0.16],
                 "transform": [ { "scale": [2, 2, 2] } ] }

    "pattern": {
      "type": "perturbed", "amplitude": 1.1, "frequency": 0.9, "seed": 11,
      "pattern": { "type": "stripe", "a": [0.93, 0.08, 0.45], "b": [0.45, 0.02, 0.28],
                   "transform": [ { "scale": [0.35, 0.35, 0.35] } ] }
    }

    "pattern": {
      "type": "texture_map", "mapping": "spherical",
      "uv": { "type": "uv_image", "file": "textures/earth.ppm" }
    }
"#;

const SETTINGS: &str = r#"SETTINGS
    Optional top-level "settings" object.

    path_tracing     { "samples": 64, "depth": 4 } - path-traced indirect
                     lighting: hemisphere samples at the first hit replace
                     the ambient term; depth diffuse bounces per path
    bvh_threshold    integer; groups with at least this many children are
                     split into a bounding-volume hierarchy after load
                     (recommended for OBJ meshes)

EXAMPLE
    "settings": { "path_tracing": { "samples": 64, "depth": 4 }, "bvh_threshold": 4 }
"#;
