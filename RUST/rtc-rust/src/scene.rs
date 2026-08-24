//! JSON scene description — load-only.
//!
//! A scene file describes a camera, one light, optional render settings, and a
//! tree of shapes; [`load`] turns it into a ready-to-render `(World, Camera)`
//! pair. The serde-derived description types mirror the core types without
//! touching them: trait objects (`Shape`, `Pattern`) are dispatched by enums
//! with a `"type"` tag, transforms are op lists, and bulk data stays in
//! referenced files (`.obj` meshes, `.ppm` textures) resolved relative to the
//! scene file's directory. Computed state (world ids, bounds, the BVH) is
//! rebuilt at load, never stored. The format is documented in `README.md`.

use serde::Deserialize;
use std::fmt;
use std::path::Path;

use crate::camera::Camera;
use crate::canvas::Canvas;
use crate::light::Light;
use crate::material::{Bump, Material};
use crate::matrix::Matrix4;
use crate::pattern::Pattern;
use crate::patterns::blendedpattern::BlendedPattern;
use crate::patterns::checkerspattern::CheckersPattern;
use crate::patterns::gradientpattern::GradientPattern;
use crate::patterns::nestedpattern::NestedPattern;
use crate::patterns::perturbedpattern::PerturbedPattern;
use crate::patterns::ringpattern::RingPattern;
use crate::patterns::stripepattern::StripePattern;
use crate::patterns::texturemap::{CubeMap, TextureMap, UvMap};
use crate::patterns::uvpattern::{UvAlignCheck, UvCheckers, UvImage, UvPattern};
use crate::shape::Shape;
use crate::shapes::cone::Cone;
use crate::shapes::csg::{Csg, CsgOperation};
use crate::shapes::cube::Cube;
use crate::shapes::cylinder::Cylinder;
use crate::shapes::group::Group;
use crate::shapes::plane::Plane;
use crate::shapes::sphere::Sphere;
use crate::shapes::torus::Torus;
use crate::tuple::Tuple;
use crate::world::{World, view_transform};

/// Error loading or building a scene, with a human-readable message.
#[derive(Debug)]
pub struct SceneError(pub String);

impl fmt::Display for SceneError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for SceneError {}

/// Top level of a scene file: a camera, one light (matching `World`),
/// optional settings, and the shape tree.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SceneDescription {
    /// Free-form notes; ignored by the loader (JSON has no comment syntax).
    #[serde(default)]
    pub comment: Option<serde_json::Value>,
    pub camera: CameraDescription,
    pub light: LightDescription,
    #[serde(default)]
    pub settings: Settings,
    #[serde(default)]
    pub shapes: Vec<ShapeDescription>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CameraDescription {
    pub width: usize,
    pub height: usize,
    pub field_of_view: f64,
    pub from: [f64; 3],
    pub to: [f64; 3],
    pub up: [f64; 3],
    #[serde(default)]
    pub antialias: Option<Antialias>,
    #[serde(default)]
    pub focal_blur: Option<FocalBlur>,
    #[serde(default)]
    pub motion_blur: Option<MotionBlur>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Antialias {
    /// Supersampling grid size: `n`×`n` rays on edge pixels.
    pub n: usize,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FocalBlur {
    pub aperture: f64,
    pub focal_distance: f64,
    pub samples: usize,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MotionBlur {
    pub samples: usize,
}

/// The three `Light` constructors, tagged by `"type"`.
#[derive(Debug, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum LightDescription {
    Point {
        position: [f64; 3],
        intensity: [f64; 3],
    },
    Area {
        corner: [f64; 3],
        uvec: [f64; 3],
        usteps: usize,
        vvec: [f64; 3],
        vsteps: usize,
        intensity: [f64; 3],
    },
    Spot {
        position: [f64; 3],
        target: [f64; 3],
        intensity: [f64; 3],
        inner_angle: f64,
        outer_angle: f64,
    },
}

#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Settings {
    /// → `World::set_path_tracing(samples, depth)`.
    #[serde(default)]
    pub path_tracing: Option<PathTracing>,
    /// Run `World::divide(root, threshold)` on every root shape after assembly.
    #[serde(default)]
    pub bvh_threshold: Option<usize>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PathTracing {
    pub samples: usize,
    pub depth: i32,
}

/// One transform step; a shape's `"transform"` is an ordered list of these,
/// applied first-to-last (`[scale, translate]` scales, then translates).
#[derive(Debug, Clone, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TransformOp {
    Translate([f64; 3]),
    Scale([f64; 3]),
    RotateX(f64),
    RotateY(f64),
    RotateZ(f64),
    /// `[xy, xz, yx, yz, zx, zy]`
    Shear([f64; 6]),
    /// Raw row-major 4×4 escape hatch.
    Matrix([[f64; 4]; 4]),
}

impl TransformOp {
    fn matrix(&self) -> Matrix4 {
        match *self {
            TransformOp::Translate([x, y, z]) => Matrix4::translation(x, y, z),
            TransformOp::Scale([x, y, z]) => Matrix4::scaling(x, y, z),
            TransformOp::RotateX(r) => Matrix4::rotation_x(r),
            TransformOp::RotateY(r) => Matrix4::rotation_y(r),
            TransformOp::RotateZ(r) => Matrix4::rotation_z(r),
            TransformOp::Shear([xy, xz, yx, yz, zx, zy]) => {
                Matrix4::shearing(xy, xz, yx, yz, zx, zy)
            }
            TransformOp::Matrix(rows) => Matrix4::new(rows),
        }
    }
}

/// Fold an op list into one matrix. Listed order is the order the ops act on
/// the object, so each successive op is premultiplied.
/// # Examples
/// ```
/// use rtc_rust::matrix::Matrix4;
/// use rtc_rust::scene::{TransformOp, transform_from_ops};
/// use rtc_rust::tuple::Tuple;
///
/// let ops = [
///     TransformOp::Scale([2.0, 2.0, 2.0]),
///     TransformOp::Translate([1.0, 0.0, 0.0]),
/// ];
/// let got = transform_from_ops(&ops) * Tuple::point(1.0, 0.0, 0.0);
/// assert!(got.approx_eq(Tuple::point(3.0, 0.0, 0.0))); // scaled first, then moved
/// ```
pub fn transform_from_ops(ops: &[TransformOp]) -> Matrix4 {
    ops.iter()
        .fold(Matrix4::identity(), |acc, op| op.matrix() * acc)
}

/// A shape node: a tagged kind plus the fields every shape shares.
#[derive(Debug, Deserialize)]
pub struct ShapeDescription {
    #[serde(flatten)]
    pub shape: ShapeKind,
    #[serde(default)]
    pub transform: Vec<TransformOp>,
    #[serde(default)]
    pub material: Option<MaterialDescription>,
    /// Velocity vector for motion blur (→ `set_motion`).
    #[serde(default)]
    pub motion: Option<[f64; 3]>,
}

#[derive(Debug, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum ShapeKind {
    Sphere,
    Plane,
    Cube,
    Cylinder {
        #[serde(default)]
        min: Option<f64>,
        #[serde(default)]
        max: Option<f64>,
        #[serde(default)]
        closed: bool,
    },
    Cone {
        #[serde(default)]
        min: Option<f64>,
        #[serde(default)]
        max: Option<f64>,
        #[serde(default)]
        closed: bool,
    },
    Torus {
        #[serde(default)]
        major_radius: Option<f64>,
        #[serde(default)]
        minor_radius: Option<f64>,
    },
    /// Wavefront OBJ file, loaded via `World::load_obj_file` into a group.
    Obj {
        file: String,
    },
    Group {
        children: Vec<ShapeDescription>,
    },
    Csg {
        operation: CsgOp,
        left: Box<ShapeDescription>,
        right: Box<ShapeDescription>,
    },
}

#[derive(Debug, Clone, Copy, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CsgOp {
    Union,
    Intersection,
    Difference,
}

impl From<CsgOp> for CsgOperation {
    fn from(op: CsgOp) -> Self {
        match op {
            CsgOp::Union => CsgOperation::Union,
            CsgOp::Intersection => CsgOperation::Intersection,
            CsgOp::Difference => CsgOperation::Difference,
        }
    }
}

/// All fields optional, on top of `Material::new()` defaults.
#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MaterialDescription {
    pub color: Option<[f64; 3]>,
    pub ambient: Option<f64>,
    pub diffuse: Option<f64>,
    pub specular: Option<f64>,
    pub shininess: Option<f64>,
    pub reflective: Option<f64>,
    pub transparency: Option<f64>,
    pub refractive_index: Option<f64>,
    pub pattern: Option<PatternDescription>,
    pub bump: Option<BumpDescription>,
}

/// Perlin normal perturbation; the noise is deterministic, so the seed is its
/// complete state.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BumpDescription {
    pub amplitude: f64,
    pub frequency: f64,
    #[serde(default)]
    pub seed: Option<u64>,
}

#[derive(Debug, Deserialize)]
pub struct PatternDescription {
    #[serde(flatten)]
    pub pattern: PatternKind,
    #[serde(default)]
    pub transform: Vec<TransformOp>,
}

#[derive(Debug, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum PatternKind {
    Stripe {
        a: [f64; 3],
        b: [f64; 3],
    },
    Gradient {
        a: [f64; 3],
        b: [f64; 3],
    },
    Ring {
        a: [f64; 3],
        b: [f64; 3],
    },
    Checkers {
        a: [f64; 3],
        b: [f64; 3],
    },
    Nested {
        a: Box<PatternDescription>,
        b: Box<PatternDescription>,
    },
    Blended {
        a: Box<PatternDescription>,
        b: Box<PatternDescription>,
    },
    Perturbed {
        pattern: Box<PatternDescription>,
        amplitude: f64,
        frequency: f64,
        #[serde(default)]
        seed: Option<u64>,
    },
    TextureMap {
        mapping: MappingKind,
        uv: UvDescription,
    },
    CubeMap {
        left: Box<UvDescription>,
        front: Box<UvDescription>,
        right: Box<UvDescription>,
        back: Box<UvDescription>,
        up: Box<UvDescription>,
        down: Box<UvDescription>,
    },
}

#[derive(Debug, Clone, Copy, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MappingKind {
    Spherical,
    Planar,
    Cylindrical,
}

impl From<MappingKind> for UvMap {
    fn from(m: MappingKind) -> Self {
        match m {
            MappingKind::Spherical => UvMap::Spherical,
            MappingKind::Planar => UvMap::Planar,
            MappingKind::Cylindrical => UvMap::Cylindrical,
        }
    }
}

#[derive(Debug, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum UvDescription {
    UvCheckers {
        width: usize,
        height: usize,
        a: [f64; 3],
        b: [f64; 3],
    },
    UvAlignCheck {
        main: [f64; 3],
        ul: [f64; 3],
        ur: [f64; 3],
        bl: [f64; 3],
        br: [f64; 3],
    },
    /// PPM image texture, path resolved relative to the scene file.
    UvImage { file: String },
}

fn point(v: [f64; 3]) -> Tuple {
    Tuple::point(v[0], v[1], v[2])
}

fn vector(v: [f64; 3]) -> Tuple {
    Tuple::vector(v[0], v[1], v[2])
}

fn color(v: [f64; 3]) -> Tuple {
    Tuple::color(v[0], v[1], v[2])
}

fn build_light(desc: &LightDescription) -> Light {
    match *desc {
        LightDescription::Point {
            position,
            intensity,
        } => Light::point_light(point(position), color(intensity)),
        LightDescription::Area {
            corner,
            uvec,
            usteps,
            vvec,
            vsteps,
            intensity,
        } => Light::area_light(
            point(corner),
            vector(uvec),
            usteps,
            vector(vvec),
            vsteps,
            color(intensity),
        ),
        LightDescription::Spot {
            position,
            target,
            intensity,
            inner_angle,
            outer_angle,
        } => Light::spotlight(
            point(position),
            point(target),
            color(intensity),
            inner_angle,
            outer_angle,
        ),
    }
}

fn build_camera(desc: &CameraDescription) -> Camera {
    let mut camera = Camera::new(desc.width, desc.height, desc.field_of_view).with_transform(
        view_transform(point(desc.from), point(desc.to), vector(desc.up)),
    );
    if let Some(aa) = &desc.antialias {
        camera = camera.with_antialias(aa.n);
    }
    if let Some(fb) = &desc.focal_blur {
        camera = camera.with_focal_blur(fb.aperture, fb.focal_distance, fb.samples);
    }
    if let Some(mb) = &desc.motion_blur {
        camera = camera.with_motion_blur(mb.samples);
    }
    camera
}

fn build_uv(desc: &UvDescription, base_dir: &Path) -> Result<Box<dyn UvPattern>, SceneError> {
    Ok(match desc {
        UvDescription::UvCheckers {
            width,
            height,
            a,
            b,
        } => Box::new(UvCheckers::new(*width, *height, color(*a), color(*b))),
        UvDescription::UvAlignCheck {
            main,
            ul,
            ur,
            bl,
            br,
        } => Box::new(UvAlignCheck::new(
            color(*main),
            color(*ul),
            color(*ur),
            color(*bl),
            color(*br),
        )),
        UvDescription::UvImage { file } => {
            let path = base_dir.join(file);
            let canvas = Canvas::read_ppm(&path)
                .map_err(|e| SceneError(format!("uv_image {}: {e}", path.display())))?;
            Box::new(UvImage::new(canvas))
        }
    })
}

fn build_pattern(
    desc: &PatternDescription,
    base_dir: &Path,
) -> Result<Box<dyn Pattern>, SceneError> {
    let mut pattern: Box<dyn Pattern> = match &desc.pattern {
        PatternKind::Stripe { a, b } => Box::new(StripePattern::new(color(*a), color(*b))),
        PatternKind::Gradient { a, b } => Box::new(GradientPattern::new(color(*a), color(*b))),
        PatternKind::Ring { a, b } => Box::new(RingPattern::new(color(*a), color(*b))),
        PatternKind::Checkers { a, b } => Box::new(CheckersPattern::new(color(*a), color(*b))),
        PatternKind::Nested { a, b } => Box::new(NestedPattern::new(
            build_pattern(a, base_dir)?,
            build_pattern(b, base_dir)?,
        )),
        PatternKind::Blended { a, b } => Box::new(BlendedPattern::new(
            build_pattern(a, base_dir)?,
            build_pattern(b, base_dir)?,
        )),
        PatternKind::Perturbed {
            pattern,
            amplitude,
            frequency,
            seed,
        } => {
            let p =
                PerturbedPattern::new(build_pattern(pattern, base_dir)?, *amplitude, *frequency);
            Box::new(match seed {
                Some(s) => p.with_seed(*s),
                None => p,
            })
        }
        PatternKind::TextureMap { mapping, uv } => Box::new(TextureMap::new(
            build_uv(uv, base_dir)?,
            UvMap::from(*mapping),
        )),
        PatternKind::CubeMap {
            left,
            front,
            right,
            back,
            up,
            down,
        } => Box::new(CubeMap::new(
            build_uv(left, base_dir)?,
            build_uv(front, base_dir)?,
            build_uv(right, base_dir)?,
            build_uv(back, base_dir)?,
            build_uv(up, base_dir)?,
            build_uv(down, base_dir)?,
        )),
    };
    if !desc.transform.is_empty() {
        pattern.set_transform(transform_from_ops(&desc.transform));
    }
    Ok(pattern)
}

fn build_material(desc: &MaterialDescription, base_dir: &Path) -> Result<Material, SceneError> {
    let mut m = Material::new();
    if let Some(c) = desc.color {
        m.color = color(c);
    }
    if let Some(v) = desc.ambient {
        m.ambient = v;
    }
    if let Some(v) = desc.diffuse {
        m.diffuse = v;
    }
    if let Some(v) = desc.specular {
        m.specular = v;
    }
    if let Some(v) = desc.shininess {
        m.shininess = v;
    }
    if let Some(v) = desc.reflective {
        m.reflective = v;
    }
    if let Some(v) = desc.transparency {
        m.transparency = v;
    }
    if let Some(v) = desc.refractive_index {
        m.refractive_index = v;
    }
    if let Some(p) = &desc.pattern {
        m.pattern = Some(build_pattern(p, base_dir)?);
    }
    if let Some(b) = &desc.bump {
        let bump = Bump::new(b.amplitude, b.frequency);
        m.bump = Some(match b.seed {
            Some(s) => bump.with_seed(s),
            None => bump,
        });
    }
    Ok(m)
}

/// Add `desc`'s subtree to the world (under `parent` if given) and return its
/// world id.
fn add_shape_to_world(
    world: &mut World,
    parent: Option<usize>,
    desc: &ShapeDescription,
    base_dir: &Path,
) -> Result<usize, SceneError> {
    let transform = transform_from_ops(&desc.transform);
    let material = match &desc.material {
        Some(m) => Some(build_material(m, base_dir)?),
        None => None,
    };

    // Primitives: build the box, configure it, add it.
    let primitive: Option<Box<dyn Shape>> = match &desc.shape {
        ShapeKind::Sphere => Some(Box::new(Sphere::new())),
        ShapeKind::Plane => Some(Box::new(Plane::new())),
        ShapeKind::Cube => Some(Box::new(Cube::new())),
        ShapeKind::Cylinder { min, max, closed } => {
            let mut c = Cylinder::new();
            c.minimum = min.unwrap_or(f64::NEG_INFINITY);
            c.maximum = max.unwrap_or(f64::INFINITY);
            c.closed = *closed;
            Some(Box::new(c))
        }
        ShapeKind::Cone { min, max, closed } => {
            let mut c = Cone::new();
            c.minimum = min.unwrap_or(f64::NEG_INFINITY);
            c.maximum = max.unwrap_or(f64::INFINITY);
            c.closed = *closed;
            Some(Box::new(c))
        }
        ShapeKind::Torus {
            major_radius,
            minor_radius,
        } => Some(Box::new(Torus::with_radii(
            major_radius.unwrap_or(1.0),
            minor_radius.unwrap_or(0.25),
        ))),
        _ => None,
    };
    if let Some(mut shape) = primitive {
        shape.set_transform(transform);
        if let Some(m) = material {
            shape.set_material(m);
        }
        if let Some(v) = desc.motion {
            shape.set_motion(vector(v));
        }
        return Ok(match parent {
            Some(p) => world.add_child(p, shape),
            None => world.add_shape(shape),
        });
    }

    // Composites: add the node first, then recurse into its children.
    let id = match &desc.shape {
        ShapeKind::Group { children } => {
            let mut g = Group::new();
            g.set_transform(transform);
            let g_id = match parent {
                Some(p) => world.add_child(p, Box::new(g)),
                None => world.add_shape(Box::new(g)),
            };
            for child in children {
                add_shape_to_world(world, Some(g_id), child, base_dir)?;
            }
            g_id
        }
        ShapeKind::Csg {
            operation,
            left,
            right,
        } => {
            let mut c = Csg::new((*operation).into());
            c.set_transform(transform);
            let c_id = match parent {
                Some(p) => world.add_child(p, Box::new(c)),
                None => world.add_shape(Box::new(c)),
            };
            // children[0] = left operand, children[1] = right operand
            add_shape_to_world(world, Some(c_id), left, base_dir)?;
            add_shape_to_world(world, Some(c_id), right, base_dir)?;
            c_id
        }
        ShapeKind::Obj { file } => {
            let path = base_dir.join(file);
            let top_id = world
                .load_obj_file(&path, transform)
                .map_err(|e| SceneError(format!("obj {}: {e}", path.display())))?;
            // load_obj_file adds the group at the root; re-parent when nested.
            if let Some(p) = parent {
                world.shapes[top_id - 1].data_mut().parent = Some(p);
                world.shapes[p - 1].add_child_id(top_id);
            }
            top_id
        }
        _ => unreachable!("primitives handled above"),
    };
    if let Some(m) = &desc.material {
        world.set_material_recursive(id, &build_material(m, base_dir)?);
    }
    if let Some(v) = desc.motion {
        world.shapes[id - 1].set_motion(vector(v));
    }
    Ok(id)
}

impl SceneDescription {
    /// Parse a scene from JSON text.
    /// # Examples
    /// ```
    /// use rtc_rust::scene::SceneDescription;
    ///
    /// let json = r#"{
    ///   "camera": { "width": 10, "height": 10, "field_of_view": 1.047,
    ///               "from": [0, 1.5, -5], "to": [0, 1, 0], "up": [0, 1, 0] },
    ///   "light": { "type": "point", "position": [-10, 10, -10], "intensity": [1, 1, 1] },
    ///   "shapes": [ { "type": "sphere", "transform": [ { "translate": [0, 1, 0] } ] } ]
    /// }"#;
    /// let scene = SceneDescription::from_json(json).unwrap();
    /// let (world, camera) = scene.build(std::path::Path::new(".")).unwrap();
    /// assert_eq!(world.shapes.len(), 1);
    /// assert_eq!(camera.hsize, 10);
    /// ```
    pub fn from_json(json: &str) -> Result<Self, SceneError> {
        serde_json::from_str(json).map_err(|e| SceneError(format!("scene JSON: {e}")))
    }

    /// Build the world and camera. `base_dir` anchors relative asset paths
    /// (`.obj`, `.ppm`); pass the scene file's directory.
    pub fn build(&self, base_dir: &Path) -> Result<(World, Camera), SceneError> {
        let mut world = World::new();
        world.set_light(build_light(&self.light));
        for shape in &self.shapes {
            add_shape_to_world(&mut world, None, shape, base_dir)?;
        }
        if let Some(pt) = &self.settings.path_tracing {
            world.set_path_tracing(pt.samples, pt.depth);
        }
        world.build_bounds();
        if let Some(threshold) = self.settings.bvh_threshold {
            let roots: Vec<usize> = world
                .shapes
                .iter()
                .filter(|s| s.data().parent.is_none())
                .map(|s| s.id())
                .collect();
            for root in roots {
                world.divide(root, threshold);
            }
            world.build_bounds();
        }
        Ok((world, build_camera(&self.camera)))
    }
}

/// Load a scene file and build it; relative asset paths resolve against the
/// file's directory.
pub fn load<P: AsRef<Path>>(path: P) -> Result<(World, Camera), SceneError> {
    let path = path.as_ref();
    let json = std::fs::read_to_string(path)
        .map_err(|e| SceneError(format!("reading {}: {e}", path.display())))?;
    let scene = SceneDescription::from_json(&json)?;
    scene.build(path.parent().unwrap_or(Path::new(".")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::approx_eq;

    /// Scene loader - a minimal scene builds a world and camera
    #[test]
    fn test_scene_1() -> Result<(), String> {
        let json = r#"{
          "camera": { "width": 160, "height": 120, "field_of_view": 1.047,
                      "from": [0, 1.5, -5], "to": [0, 1, 0], "up": [0, 1, 0] },
          "light": { "type": "point", "position": [-10, 10, -10], "intensity": [1, 1, 1] },
          "shapes": [
            { "type": "sphere",
              "transform": [ { "translate": [0, 1, 0] } ],
              "material": { "color": [0.1, 0.2, 0.8], "reflective": 0.3 } }
          ]
        }"#;
        let scene = SceneDescription::from_json(json).map_err(|e| e.to_string())?;
        let (world, camera) = scene.build(Path::new(".")).map_err(|e| e.to_string())?;

        let s = world.shapes[0].as_ref();
        let chk = camera.hsize == 160
            && camera.vsize == 120
            && world.light.is_some()
            && world.shapes.len() == 1
            && s.material().color.approx_eq(Tuple::color(0.1, 0.2, 0.8))
            && approx_eq(s.material().reflective, 0.3)
            && (*s.transform() * Tuple::point(0.0, 0.0, 0.0))
                .approx_eq(Tuple::point(0.0, 1.0, 0.0));
        if chk {
            Ok(())
        } else {
            Err("minimal scene must build a world and camera".into())
        }
    }

    /// Scene loader - transform ops apply in listed order
    #[test]
    fn test_scene_2() -> Result<(), String> {
        let ops = [
            TransformOp::Scale([2.0, 2.0, 2.0]),
            TransformOp::RotateY(std::f64::consts::FRAC_PI_2),
            TransformOp::Translate([5.0, 0.0, 0.0]),
        ];
        let m = transform_from_ops(&ops);
        // (1,0,0) --scale--> (2,0,0) --rot_y--> (0,0,-2) --translate--> (5,0,-2)
        let got = m * Tuple::point(1.0, 0.0, 0.0);
        if got.approx_eq(Tuple::point(5.0, 0.0, -2.0)) {
            Ok(())
        } else {
            Err(format!("ops must apply in listed order, got {got:?}"))
        }
    }

    /// Scene loader - all three light types build their constructors' shapes
    #[test]
    fn test_scene_3() -> Result<(), String> {
        let point_l = build_light(
            &serde_json::from_str(
                r#"{ "type": "point", "position": [1, 2, 3], "intensity": [1, 1, 1] }"#,
            )
            .map_err(|e| e.to_string())?,
        );
        let area_l = build_light(
            &serde_json::from_str(
                r#"{ "type": "area", "corner": [-1, 2, 4], "uvec": [2, 0, 0], "usteps": 4,
                     "vvec": [0, 2, 0], "vsteps": 2, "intensity": [1.5, 1.5, 1.5] }"#,
            )
            .map_err(|e| e.to_string())?,
        );
        let spot_l = build_light(
            &serde_json::from_str(
                r#"{ "type": "spot", "position": [0, 5, 0], "target": [0, 0, 0],
                     "intensity": [1, 1, 1], "inner_angle": 0.3, "outer_angle": 0.5 }"#,
            )
            .map_err(|e| e.to_string())?,
        );
        let chk = point_l.position.approx_eq(Tuple::point(1.0, 2.0, 3.0))
            && point_l.samples == 1
            && area_l.usteps == 4
            && area_l.vsteps == 2
            && area_l.samples == 8
            && area_l.position.approx_eq(Tuple::point(0.0, 3.0, 4.0))
            && spot_l.spot.is_some()
            && approx_eq(spot_l.spot.as_ref().unwrap().cos_inner, 0.3f64.cos());
        if chk {
            Ok(())
        } else {
            Err("light types must map onto the Light constructors".into())
        }
    }

    /// Scene loader - patterns build, including combinators and texture maps
    #[test]
    fn test_scene_4() -> Result<(), String> {
        let desc: PatternDescription = serde_json::from_str(
            r#"{ "type": "checkers", "a": [1, 1, 1], "b": [0, 0, 0],
                 "transform": [ { "scale": [2, 2, 2] } ] }"#,
        )
        .map_err(|e| e.to_string())?;
        let checkers = build_pattern(&desc, Path::new(".")).map_err(|e| e.to_string())?;

        let desc: PatternDescription = serde_json::from_str(
            r#"{ "type": "perturbed", "amplitude": 0.4, "frequency": 2.0, "seed": 7,
                 "pattern": { "type": "stripe", "a": [1, 0, 0], "b": [0, 0, 1] } }"#,
        )
        .map_err(|e| e.to_string())?;
        let perturbed = build_pattern(&desc, Path::new(".")).map_err(|e| e.to_string())?;

        let desc: PatternDescription = serde_json::from_str(
            r#"{ "type": "texture_map", "mapping": "spherical",
                 "uv": { "type": "uv_checkers", "width": 16, "height": 8,
                         "a": [0, 0.5, 0], "b": [1, 1, 1] } }"#,
        )
        .map_err(|e| e.to_string())?;
        let texture = build_pattern(&desc, Path::new(".")).map_err(|e| e.to_string())?;

        // checkers at origin returns `a`; the others just have to build.
        let chk = checkers
            .color_at(Tuple::point(0.0, 0.0, 0.0))
            .approx_eq(Tuple::color(1.0, 1.0, 1.0))
            && perturbed.color_at(Tuple::point(0.0, 0.0, 0.0)).w == 0.0
            && texture.color_at(Tuple::point(0.0, 1.0, 0.0)).w == 0.0;
        if chk {
            Ok(())
        } else {
            Err("patterns must build from their descriptions".into())
        }
    }

    /// Scene loader - groups and CSG wire up parent/child ids
    #[test]
    fn test_scene_5() -> Result<(), String> {
        let json = r#"{
          "camera": { "width": 10, "height": 10, "field_of_view": 1.047,
                      "from": [0, 0, -5], "to": [0, 0, 0], "up": [0, 1, 0] },
          "light": { "type": "point", "position": [-10, 10, -10], "intensity": [1, 1, 1] },
          "shapes": [
            { "type": "group",
              "transform": [ { "translate": [0, 1, 0] } ],
              "children": [
                { "type": "sphere" },
                { "type": "csg", "operation": "difference",
                  "left":  { "type": "cube" },
                  "right": { "type": "sphere", "transform": [ { "scale": [1.2, 1.2, 1.2] } ] } }
              ] }
          ]
        }"#;
        let scene = SceneDescription::from_json(json).map_err(|e| e.to_string())?;
        let (world, _) = scene.build(Path::new(".")).map_err(|e| e.to_string())?;

        // ids in add order: 1 group, 2 sphere, 3 csg, 4 cube (left), 5 sphere (right)
        let group = world.shapes[0].as_ref();
        let csg = world.shapes[2].as_ref();
        let chk = world.shapes.len() == 5
            && group.children() == Some(&[2, 3][..])
            && csg.csg_operation() == Some(CsgOperation::Difference)
            && csg.children() == Some(&[4, 5][..])
            && world.shapes[3].data().parent == Some(3)
            && world.shapes[1].data().parent == Some(1);
        if chk {
            Ok(())
        } else {
            Err("group/CSG trees must wire parent and child ids".into())
        }
    }

    /// Scene loader - settings map to path tracing and BVH subdivision
    #[test]
    fn test_scene_6() -> Result<(), String> {
        let json = r#"{
          "camera": { "width": 10, "height": 10, "field_of_view": 1.047,
                      "from": [0, 0, -5], "to": [0, 0, 0], "up": [0, 1, 0] },
          "light": { "type": "point", "position": [-10, 10, -10], "intensity": [1, 1, 1] },
          "settings": { "path_tracing": { "samples": 16, "depth": 3 }, "bvh_threshold": 2 },
          "shapes": [
            { "type": "group", "children": [
              { "type": "sphere", "transform": [ { "translate": [-2, 0, 0] } ] },
              { "type": "sphere", "transform": [ { "translate": [-1, 0, 0] } ] },
              { "type": "sphere", "transform": [ { "translate": [1, 0, 0] } ] },
              { "type": "sphere", "transform": [ { "translate": [2, 0, 0] } ] }
            ] }
          ]
        }"#;
        let scene = SceneDescription::from_json(json).map_err(|e| e.to_string())?;
        let (world, _) = scene.build(Path::new(".")).map_err(|e| e.to_string())?;

        // divide(2) must have reorganized the 4 spheres into sub-groups.
        let top_children = world.shapes[0].children().unwrap();
        let has_subgroup = top_children
            .iter()
            .any(|&id| world.shapes[id - 1].children().is_some());
        let chk = world.gi_samples == 16 && world.gi_depth == 3 && has_subgroup;
        if chk {
            Ok(())
        } else {
            Err("settings must enable path tracing and BVH subdivision".into())
        }
    }

    /// Scene loader - malformed scenes are rejected with an error
    #[test]
    fn test_scene_7() -> Result<(), String> {
        let unknown_field = r#"{
          "camera": { "width": 10, "height": 10, "field_of_view": 1.047,
                      "from": [0, 0, -5], "to": [0, 0, 0], "up": [0, 1, 0],
                      "focal_lenght": 5 },
          "light": { "type": "point", "position": [0, 0, 0], "intensity": [1, 1, 1] }
        }"#;
        let bad_shape = r#"{
          "camera": { "width": 10, "height": 10, "field_of_view": 1.047,
                      "from": [0, 0, -5], "to": [0, 0, 0], "up": [0, 1, 0] },
          "light": { "type": "point", "position": [0, 0, 0], "intensity": [1, 1, 1] },
          "shapes": [ { "type": "teapot" } ]
        }"#;
        let missing_obj = r#"{
          "camera": { "width": 10, "height": 10, "field_of_view": 1.047,
                      "from": [0, 0, -5], "to": [0, 0, 0], "up": [0, 1, 0] },
          "light": { "type": "point", "position": [0, 0, 0], "intensity": [1, 1, 1] },
          "shapes": [ { "type": "obj", "file": "no_such_file.obj" } ]
        }"#;
        let chk = SceneDescription::from_json(unknown_field).is_err()
            && SceneDescription::from_json(bad_shape).is_err()
            && SceneDescription::from_json(missing_obj)
                .map_err(|e| e.to_string())?
                .build(Path::new("."))
                .is_err();
        if chk {
            Ok(())
        } else {
            Err("malformed scenes must be rejected".into())
        }
    }

    /// Scene loader - OBJ shapes load the referenced model into a group
    #[test]
    fn test_scene_8() -> Result<(), String> {
        let json = r#"{
          "camera": { "width": 10, "height": 10, "field_of_view": 1.047,
                      "from": [0, 0, -5], "to": [0, 0, 0], "up": [0, 1, 0] },
          "light": { "type": "point", "position": [-10, 10, -10], "intensity": [1, 1, 1] },
          "shapes": [
            { "type": "obj", "file": "models/icosahedron.obj",
              "transform": [ { "scale": [0.5, 0.5, 0.5] } ],
              "material": { "color": [1, 0, 0] } }
          ]
        }"#;
        let scene = SceneDescription::from_json(json).map_err(|e| e.to_string())?;
        let (world, _) = scene
            .build(Path::new(env!("CARGO_MANIFEST_DIR")))
            .map_err(|e| e.to_string())?;

        let top = world.shapes[0].as_ref();
        let n_children = top.children().map_or(0, |c| c.len());
        let chk = n_children > 0
            && world.shapes[1]
                .material()
                .color
                .approx_eq(Tuple::color(1.0, 0.0, 0.0));
        if chk {
            Ok(())
        } else {
            Err(format!(
                "OBJ scene must load triangles, got {n_children} children"
            ))
        }
    }

    /// Scene loader - the book's cover image (Appendix A1) renders from JSON.
    /// The scene file's camera is wallpaper-sized, so the test renders it
    /// downscaled; `cargo run --release -- scenes/cover.json` does the full one.
    #[test]
    fn test_scene_cover_putting_it_together() -> std::io::Result<()> {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("scenes/cover.json");
        let json = std::fs::read_to_string(&path)?;
        let mut scene = SceneDescription::from_json(&json).map_err(std::io::Error::other)?;
        scene.camera.width = 216; // full scene, 1/16 of the file's resolution
        scene.camera.height = 140;
        let (world, camera) = scene
            .build(path.parent().unwrap())
            .map_err(std::io::Error::other)?;
        let canvas = world.render_parallel(camera);
        canvas.write_ppm("test_scene_cover_putting_it_together.ppm")?;
        Ok(())
    }
}
