//! Flatten an rtc-rust `World` + `Camera` into GPU buffers.
//!
//! Only what the shaders implement so far is accepted; anything else is an
//! error rather than being silently dropped from the image.

use crate::bvh::{self, Leaf};
use crate::types::{
    GpuBvhNode, GpuCamera, GpuLight, GpuMaterial, GpuSceneCounts, GpuSceneInfo, GpuShape,
    PATTERN_CHECKERS, PATTERN_GRADIENT, PATTERN_NONE, PATTERN_RING, PATTERN_STRIPE, SHAPE_CONE,
    SHAPE_CUBE, SHAPE_CYLINDER, SHAPE_PLANE, SHAPE_SPHERE, SHAPE_TRIANGLE, float4,
};
use rtc_rust::camera::Camera;
use rtc_rust::light::Light;
use rtc_rust::material::Material;
use rtc_rust::matrix::Matrix4;
use rtc_rust::pattern::Pattern;
use rtc_rust::patterns::checkerspattern::CheckersPattern;
use rtc_rust::patterns::gradientpattern::GradientPattern;
use rtc_rust::patterns::ringpattern::RingPattern;
use rtc_rust::patterns::stripepattern::StripePattern;
use rtc_rust::shape::Shape;
use rtc_rust::shapes::cone::Cone;
use rtc_rust::shapes::cube::Cube;
use rtc_rust::shapes::cylinder::Cylinder;
use rtc_rust::shapes::plane::Plane;
use rtc_rust::shapes::sphere::Sphere;
use rtc_rust::shapes::triangle::Triangle;
use rtc_rust::shapes::triangleuv::TriangleUV;
use rtc_rust::tuple::Tuple;
use rtc_rust::world::World;
use std::any::Any;

pub struct GpuScene {
    pub info: GpuSceneInfo,
    /// Unbounded shapes first, then the bounded ones in BVH leaf order.
    pub shapes: Vec<GpuShape>,
    /// Book chapter 15's bounding volume hierarchy over the bounded shapes;
    /// empty when there are none.
    pub nodes: Vec<GpuBvhNode>,
    /// rtc-rust shape id per entry of `shapes`.
    pub ids: Vec<usize>,
    /// Lookup tables the shapes and materials refer to by word offset.
    pub tables: Vec<u32>,
}

impl GpuScene {
    pub fn from_world(world: &World, camera: &Camera) -> Result<Self, String> {
        if world.gi_samples > 0 {
            return Err("path tracing is not supported so far".into());
        }
        let light = world.light.as_ref().ok_or("scene has no light")?;
        let mut leaves = Vec::new();
        let identity = Matrix4::identity();
        for root in world.shapes.iter().filter(|s| s.data().parent.is_none()) {
            flatten(world, root.as_ref(), &identity, &identity, &mut leaves)?;
        }
        let bvh = bvh::build(leaves)?;
        Ok(Self {
            info: GpuSceneInfo {
                camera: gpu_camera(camera)?,
                light: gpu_light(light)?,
                counts: GpuSceneCounts {
                    shape_count: bvh.shapes.len() as u32,
                    unbounded_count: bvh.unbounded_count,
                    node_count: bvh.nodes.len() as u32,
                },
                frame: 0,
                sampling: 0,
            },
            shapes: bvh.shapes,
            nodes: bvh.nodes,
            ids: bvh.ids,
            tables: Vec::new(),
        })
    }

    /// Point the scene at `camera`, keeping the uploaded shapes and light.
    pub fn set_camera(&mut self, camera: &Camera) -> Result<(), String> {
        self.info.camera = gpu_camera(camera)?;
        Ok(())
    }

    pub fn width(&self) -> usize {
        self.info.camera.width as usize
    }

    pub fn height(&self) -> usize {
        self.info.camera.height as usize
    }
}

/// The camera's anti-aliasing sample count is ignored: a deterministic frame
/// has one ray through each pixel centre, and sampled frames jitter that ray
/// within the pixel instead, anti-aliasing as they accumulate.
pub fn gpu_camera(camera: &Camera) -> Result<GpuCamera, String> {
    if camera.aperture > 0.0 || camera.motion_samples > 1 {
        return Err("focal blur and motion blur are not supported so far".into());
    }
    let inverse = camera
        .transform()
        .inverse()
        .ok_or("camera transform is not invertible")?;
    Ok(GpuCamera {
        inverse: (&inverse).into(),
        half_width: camera.half_width as f32,
        half_height: camera.half_height as f32,
        pixel_size: camera.pixel_size as f32,
        width: camera.hsize as u32,
        height: camera.vsize as u32,
    })
}

/// rtc-rust's point light is already the 1x1 area light: `corner` is its
/// position and the edges are zero.
pub fn gpu_light(light: &Light) -> Result<GpuLight, String> {
    if light.spot.is_some() {
        return Err("spot lights are not supported so far".into());
    }
    Ok(GpuLight {
        corner: float4(light.corner),
        uvec: float4(light.uvec),
        vvec: float4(light.vvec),
        intensity: float4(light.intensity),
        usteps: light.usteps as u32,
        vsteps: light.vsteps as u32,
        samples: light.samples as u32,
        _pad: 0,
    })
}

pub fn gpu_material(m: &Material) -> Result<GpuMaterial, String> {
    if m.bump.is_some() {
        return Err("bump maps are not supported yet".into());
    }
    let (pattern, pattern_inverse, pattern_a, pattern_b) = match &m.pattern {
        None => (
            PATTERN_NONE,
            Matrix4::identity(),
            Tuple::color(0.0, 0.0, 0.0),
            Tuple::color(0.0, 0.0, 0.0),
        ),
        Some(p) => {
            let (kind, a, b) = gpu_pattern(p.as_ref())?;
            (kind, *p.data().transform_inv(), a, b)
        }
    };
    Ok(GpuMaterial {
        color: float4(m.color),
        ambient: m.ambient as f32,
        diffuse: m.diffuse as f32,
        specular: m.specular as f32,
        shininess: m.shininess as f32,
        reflective: m.reflective as f32,
        transparency: m.transparency as f32,
        refractive_index: m.refractive_index as f32,
        pattern,
        pattern_inverse: (&pattern_inverse).into(),
        pattern_a: float4(pattern_a),
        pattern_b: float4(pattern_b),
    })
}

/// The `PATTERN_*` kind and two colors of one of the book's four patterns.
fn gpu_pattern(pattern: &dyn Pattern) -> Result<(u32, Tuple, Tuple), String> {
    // As in gpu_shape: upcast the `&dyn Pattern`, not the `Box`.
    let any: &dyn Any = pattern;
    if let Some(p) = any.downcast_ref::<StripePattern>() {
        Ok((PATTERN_STRIPE, p.a, p.b))
    } else if let Some(p) = any.downcast_ref::<GradientPattern>() {
        Ok((PATTERN_GRADIENT, p.a, p.b))
    } else if let Some(p) = any.downcast_ref::<RingPattern>() {
        Ok((PATTERN_RING, p.a, p.b))
    } else if let Some(p) = any.downcast_ref::<CheckersPattern>() {
        Ok((PATTERN_CHECKERS, p.a, p.b))
    } else {
        Err("only stripe, gradient, ring and checkers patterns are supported so far".into())
    }
}

/// Book chapter 14 without the tree: append the leaves under `shape` to
/// `out`, each carrying its whole world-to-object transform and world box,
/// so the shaders never see a group. `parent_transform` maps the space
/// `shape` sits in to world space, `parent_inverse` the other way. Composing
/// in f64 here, rather than per level in f32 on the GPU, also keeps deep
/// trees from accumulating rounding. Group bounds (and the BVH that
/// `World::divide` builds) are dropped: `bvh::build` makes its own over the
/// leaves (book chapter 15).
fn flatten(
    world: &World,
    shape: &dyn Shape,
    parent_transform: &Matrix4,
    parent_inverse: &Matrix4,
    out: &mut Vec<Leaf>,
) -> Result<(), String> {
    // world_to_object applies the parent chain first, then the shape's own;
    // object_to_world is the reverse.
    let full = *parent_transform * *shape.transform();
    let inverse = *shape.transform_inv() * *parent_inverse;
    match shape.children() {
        None => out.push(Leaf {
            shape: gpu_leaf(shape, &inverse)?,
            id: shape.id(),
            bounds: world_bounds(shape, full),
        }),
        Some(_) if shape.csg_operation().is_some() => {
            return Err(format!("shape {}: CSG is not supported yet", shape.id()));
        }
        Some(children) => {
            if shape.data().motion.is_some() {
                return Err(format!("shape {}: motion is not supported yet", shape.id()));
            }
            for &id in children {
                // World ids are 1-based arena indices.
                flatten(world, world.shapes[id - 1].as_ref(), &full, &inverse, out)?;
            }
        }
    }
    Ok(())
}

/// The f32 world box of a leaf whose object-to-world transform is `full`,
/// or None when it is unbounded. An infinite object box is unbounded as is:
/// transforming its corners would give inf * 0 = NaN under a rotation.
fn world_bounds(shape: &dyn Shape, full: Matrix4) -> Option<([f32; 3], [f32; 3])> {
    let b = shape.bounds();
    let xyz = |t: Tuple| [t.x, t.y, t.z];
    if !xyz(b.min).iter().chain(&xyz(b.max)).all(|v| v.is_finite()) {
        return None;
    }
    let b = b.transform(full);
    // Round outward, so the f32 box still holds the f64 one and a ray
    // grazing the shape is never culled.
    let lo = xyz(b.min).map(|v| {
        let f = v as f32;
        if f as f64 > v { f.next_down() } else { f }
    });
    let hi = xyz(b.max).map(|v| {
        let f = v as f32;
        if (f as f64) < v { f.next_up() } else { f }
    });
    // A finite f64 box can still overflow f32.
    lo.iter()
        .chain(&hi)
        .all(|v| v.is_finite())
        .then_some((lo, hi))
}

/// A shape outside any group, as the probes upload it.
pub fn gpu_shape(shape: &dyn Shape) -> Result<GpuShape, String> {
    if shape.data().parent.is_some() {
        return Err(format!(
            "shape {}: a group child needs its parents' transforms; use GpuScene",
            shape.id()
        ));
    }
    gpu_leaf(shape, shape.transform_inv())
}

/// A primitive seen through `inverse`, its world-to-object transform.
fn gpu_leaf(shape: &dyn Shape, inverse: &Matrix4) -> Result<GpuShape, String> {
    // Upcast to `&dyn Any` to learn the concrete type. This must be the
    // `&dyn Shape`, not the `Box`: a `&Box<dyn Shape>` would coerce to an
    // `Any` whose type is the Box itself, and `is::<Sphere>` would be false.
    let any: &dyn Any = shape;
    // (kind, minimum, maximum, closed); the extent only matters for
    // cylinders and cones. f64 infinities stay infinite in f32.
    let mut triangle = [[0.0; 4]; 6]; // p1, e1, e2, n1, n2, n3
    let (kind, minimum, maximum, closed) = if any.is::<Sphere>() {
        (SHAPE_SPHERE, 0.0, 0.0, false)
    } else if any.is::<Plane>() {
        (SHAPE_PLANE, 0.0, 0.0, false)
    } else if any.is::<Cube>() {
        (SHAPE_CUBE, 0.0, 0.0, false)
    } else if let Some(c) = any.downcast_ref::<Cylinder>() {
        (SHAPE_CYLINDER, c.minimum, c.maximum, c.closed)
    } else if let Some(c) = any.downcast_ref::<Cone>() {
        (SHAPE_CONE, c.minimum, c.maximum, c.closed)
    } else if let Some(t) = any.downcast_ref::<Triangle>() {
        // flat: every vertex normal is the face normal
        triangle = [t.p1, t.e1, t.e2, t.normal, t.normal, t.normal].map(float4);
        (SHAPE_TRIANGLE, 0.0, 0.0, false)
    } else if let Some(t) = any.downcast_ref::<TriangleUV>() {
        triangle = [t.p1, t.e1, t.e2, t.n1, t.n2, t.n3].map(float4);
        (SHAPE_TRIANGLE, 0.0, 0.0, false)
    } else {
        return Err(format!(
            "shape {}: only spheres, planes, cubes, cylinders, cones and triangles are supported so far",
            shape.id()
        ));
    };
    let data = shape.data();
    if data.motion.is_some() {
        return Err(format!("shape {}: motion is not supported yet", shape.id()));
    }
    Ok(GpuShape {
        inverse: inverse.into(),
        material: gpu_material(&data.material)?,
        kind,
        minimum: minimum as f32,
        maximum: maximum as f32,
        closed: closed as u32,
        p1: triangle[0],
        e1: triangle[1],
        e2: triangle[2],
        n1: triangle[3],
        n2: triangle[4],
        n3: triangle[5],
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{Float4, GpuMat4};
    use rtc_rust::shapes::csg::{Csg, CsgOperation};
    use rtc_rust::shapes::group::Group;
    use rtc_rust::shapes::torus::Torus;
    use std::f64::consts::PI;

    /// `m * v` with the f32 matrix the shader gets.
    fn apply(m: &GpuMat4, v: Float4) -> Float4 {
        [0, 1, 2, 3].map(|r| (0..4).map(|c| m.cols[c][r] * v[c]).sum())
    }

    #[test]
    fn unsupported_shapes_are_rejected() {
        let mut world = World::new();
        world.light = Some(Light::point_light(
            Tuple::point(0.0, 0.0, 0.0),
            Tuple::color(1.0, 1.0, 1.0),
        ));
        world.add_shape(Box::new(Sphere::new()));
        world.add_shape(Box::new(Plane::new()));
        world.add_shape(Box::new(Cube::new()));
        world.add_shape(Box::new(Cylinder::new()));
        world.add_shape(Box::new(Cone::new()));
        let g = world.add_shape(Box::new(Group::new()));
        world.add_child(g, Box::new(Sphere::new()));
        assert!(GpuScene::from_world(&world, &Camera::default()).is_ok());
        let child = world.add_child(g, Box::new(Torus::new()));
        let err = GpuScene::from_world(&world, &Camera::default()).err();
        assert!(err.is_some_and(|e| e.contains("only spheres, planes")));
        // Standalone upload would drop the group's transform.
        let err = gpu_shape(world.shapes[child - 1].as_ref()).err();
        assert!(err.is_some_and(|e| e.contains("use GpuScene")));

        let mut world = World::new();
        world.light = Some(Light::point_light(
            Tuple::point(0.0, 0.0, 0.0),
            Tuple::color(1.0, 1.0, 1.0),
        ));
        let c = world.add_shape(Box::new(Csg::new(CsgOperation::Union)));
        world.add_child(c, Box::new(Sphere::new()));
        world.add_child(c, Box::new(Cube::new()));
        let err = GpuScene::from_world(&world, &Camera::default()).err();
        assert!(err.is_some_and(|e| e.contains("CSG")));
    }

    #[test]
    fn test_chap_14_flattened_transforms_match_book() {
        // The book's nested groups: g1 rotates, g2 scales, s translates.
        let mut world = World::new();
        world.light = Some(Light::point_light(
            Tuple::point(0.0, 0.0, 0.0),
            Tuple::color(1.0, 1.0, 1.0),
        ));
        let mut g1 = Group::new();
        g1.set_transform(Matrix4::rotation_y(PI / 2.0));
        let g1 = world.add_shape(Box::new(g1));
        let mut g2 = Group::new();
        g2.set_transform(Matrix4::scaling(1.0, 2.0, 3.0));
        let g2 = world.add_child(g1, Box::new(g2));
        let mut s = Sphere::new();
        s.set_transform(Matrix4::translation(5.0, 0.0, 0.0));
        world.add_child(g2, Box::new(s));
        let scene = GpuScene::from_world(&world, &Camera::default()).unwrap();
        assert_eq!(scene.shapes.len(), 1, "only the leaf is uploaded");
        let inverse = &scene.shapes[0].inverse;

        // World to object space, as the shader transforms the ray: the
        // book's point, with g2 scaling (1, 2, 3) from its normal test
        // instead of 2: rotate (-2, 0, -10) to (10, 0, -2), scale to
        // (10, 0, -2/3), translate to (5, 0, -2/3).
        let p = apply(inverse, [-2.0, 0.0, -10.0, 1.0]);
        let want = [5.0, 0.0, -2.0 / 3.0, 1.0];
        assert!(
            p.iter().zip(want).all(|(a, b)| (a - b).abs() < 1e-5),
            "{p:?}"
        );

        // Object to world normal, as the shader does: transpose(inverse) * n.
        // The non-uniform scaling makes a wrong composition order show.
        let r = 3f32.sqrt() / 3.0;
        let n = [0, 1, 2].map(|c| (0..3).map(|k| inverse.cols[c][k] * r).sum::<f32>());
        let len = n.iter().map(|x| x * x).sum::<f32>().sqrt();
        let want = [0.2857, 0.4286, -0.8571];
        assert!(
            n.iter().zip(want).all(|(a, b)| (a / len - b).abs() < 1e-4),
            "{n:?}"
        );
    }

    #[test]
    fn test_chap_15_scene_bounds_and_order() {
        let mut world = World::new();
        world.light = Some(Light::point_light(
            Tuple::point(0.0, 0.0, 0.0),
            Tuple::color(1.0, 1.0, 1.0),
        ));
        let p1 = world.add_shape(Box::new(Plane::new()));
        // A naive transform of the plane's infinite box gives inf * 0 = NaN.
        let mut tilted = Plane::new();
        tilted.set_transform(Matrix4::rotation_x(0.3));
        let p2 = world.add_shape(Box::new(tilted));
        let mut cyl = Cylinder::new();
        cyl.set_transform(Matrix4::rotation_z(0.5));
        let c = world.add_shape(Box::new(cyl));
        let mut s = Sphere::new();
        s.set_transform(Matrix4::translation(5.0, 0.0, 0.0));
        world.add_shape(Box::new(s));
        let mut g = Group::new();
        g.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));
        let g = world.add_shape(Box::new(g));
        world.add_child(g, Box::new(Sphere::new()));

        let scene = GpuScene::from_world(&world, &Camera::default()).unwrap();
        assert_eq!(scene.info.counts.unbounded_count, 3);
        assert_eq!(scene.info.counts.shape_count, 5);
        assert_eq!(scene.shapes.len(), 5);
        assert_eq!(scene.ids[..3], [p1, p2, c]);
        assert_eq!(scene.nodes.len(), 1);
        assert_eq!(scene.info.counts.node_count, 1);
        let n = &scene.nodes[0];
        assert_eq!((n.start, n.count), (3, 2));
        let contains = |lo: [f64; 3], hi: [f64; 3]| {
            (0..3).all(|a| n.min[a] as f64 <= lo[a] && n.max[a] as f64 >= hi[a])
        };
        assert!(contains([4.0, -1.0, -1.0], [6.0, 1.0, 1.0]), "{n:?}");
        assert!(contains([-2.0; 3], [2.0; 3]), "{n:?}");
        // The outward cast keeps the f64 box inside the f32 one.
        for &id in &scene.ids[3..] {
            let shape = world.shapes[id - 1].as_ref();
            let parent = match shape.data().parent {
                Some(p) => *world.shapes[p - 1].transform(),
                None => Matrix4::identity(),
            };
            let b = shape.bounds().transform(parent * *shape.transform());
            assert!(
                contains([b.min.x, b.min.y, b.min.z], [b.max.x, b.max.y, b.max.z]),
                "shape {id}: {n:?}"
            );
        }
    }

    #[test]
    fn test_chap_15_scene_without_bounded_shapes_renders() {
        // Metal cannot bind an empty buffer: an empty node array must still
        // render.
        let mut camera = Camera::new(4, 2, PI / 2.0);
        camera.set_transform(&Matrix4::translation(0.0, 1.0, 0.0));
        let gpu = crate::gpu::Gpu::new().unwrap();
        for with_plane in [true, false] {
            let mut world = World::new();
            world.light = Some(Light::point_light(
                Tuple::point(0.0, 5.0, 0.0),
                Tuple::color(1.0, 1.0, 1.0),
            ));
            if with_plane {
                world.add_shape(Box::new(Plane::new()));
            }
            let scene = GpuScene::from_world(&world, &camera).unwrap();
            assert!(scene.nodes.is_empty());
            assert_eq!(scene.info.counts.node_count, 0);
            assert_eq!(gpu.render(&scene).unwrap().len(), 8);
        }
    }

    #[test]
    fn unsupported_patterns_are_rejected() {
        use rtc_rust::patterns::nestedpattern::NestedPattern;
        let (a, b) = (Tuple::color(1.0, 1.0, 1.0), Tuple::color(0.0, 0.0, 0.0));
        let mut m = Material::new();
        m.pattern = Some(Box::new(CheckersPattern::new(a, b)));
        assert!(gpu_material(&m).is_ok());
        m.pattern = Some(Box::new(NestedPattern::new(
            Box::new(StripePattern::new(a, b)),
            Box::new(RingPattern::new(a, b)),
        )));
        let err = gpu_material(&m).err();
        assert!(err.is_some_and(|e| e.contains("only stripe, gradient, ring and checkers")));
    }

    #[test]
    fn set_camera_replaces_only_the_camera() {
        let world = World::default_world();
        let mut scene = GpuScene::from_world(&world, &Camera::default()).unwrap();
        let (shapes, nodes, info) = (scene.shapes.clone(), scene.nodes.clone(), scene.info);
        let mut camera = Camera::default();
        camera.set_transform(&Matrix4::translation(1.0, 2.0, 3.0));
        scene.set_camera(&camera).unwrap();
        assert_ne!(scene.info.camera, info.camera);
        assert_eq!(scene.info.camera, gpu_camera(&camera).unwrap());
        assert_eq!(scene.shapes, shapes);
        assert_eq!(scene.nodes, nodes);
        assert_eq!(scene.info.light, info.light);
        assert_eq!(scene.info.counts, info.counts);
    }
}
