//! World
//!
//! This module defines what rays can hit.
//!

use rayon::prelude::*;
use std::path::Path;

use crate::bounds::BoundingBox;
use crate::camera::Camera;
use crate::canvas::Canvas;
use crate::intersection::{Intersection, Intersections};
use crate::light::Light;
use crate::material::Material;
use crate::math::approx_eq;
use crate::matrix::Matrix4;
use crate::obj::Parser;
use crate::ray::Ray;
use crate::shape::Shape;
use crate::shapes::csg::intersection_allowed;
use crate::shapes::cylinder::Cylinder;
use crate::shapes::group::Group;
use crate::shapes::sphere::Sphere;
use crate::tuple::Tuple;
use crate::{log::*, tuple};
use std::fmt;
use std::sync::atomic::{AtomicUsize, Ordering};

#[cfg(feature = "stats")]
pub static NODE_VISITS: AtomicUsize = AtomicUsize::new(0);
#[cfg(feature = "stats")]
pub static PRIM_TESTS: AtomicUsize = AtomicUsize::new(0);

// Zero-cost without `--features stats`: empty body, `#[inline(always)]` elides it.
#[inline(always)]
fn record_node_visit() {
    #[cfg(feature = "stats")]
    NODE_VISITS.fetch_add(1, Ordering::Relaxed);
}
#[inline(always)]
fn record_prim_test() {
    #[cfg(feature = "stats")]
    PRIM_TESTS.fetch_add(1, Ordering::Relaxed);
}

/// Reset traversal counters (no-op without `--features stats`).
pub fn reset_stats() {
    #[cfg(feature = "stats")]
    {
        NODE_VISITS.store(0, Ordering::Relaxed);
        PRIM_TESTS.store(0, Ordering::Relaxed);
    }
}

/// (node_visits, prim_tests) since last reset; (0, 0) without `--features stats`.
pub fn read_stats() -> (usize, usize) {
    #[cfg(feature = "stats")]
    {
        (
            NODE_VISITS.load(Ordering::Relaxed),
            PRIM_TESTS.load(Ordering::Relaxed),
        )
    }
    #[cfg(not(feature = "stats"))]
    {
        (0, 0)
    }
}

pub struct Computations<'a> {
    pub t: f64,
    pub object: &'a dyn Shape,
    pub point: Tuple,
    pub over_point: Tuple,
    pub under_point: Tuple,
    pub object_point: Tuple,
    pub eyev: Tuple,
    pub normalv: Tuple,
    pub reflectv: Tuple,
    pub inside: bool,
    pub n1: f64,
    pub n2: f64,
}

/// Encapsulating some precomputed information relating to the intersection.
///
/// Precomputes the following:
/// * the point (in world space) where the intersection occurred
/// * the eye vector (pointing back toward the eye or camera)
/// * the normal vector
///
pub fn prepare_computations_upto_chap10<'a>(
    intersection: Intersection,
    ray: &Ray,
    shape: &'a dyn Shape,
) -> Computations<'a> {
    let t = intersection.t;
    let point = ray.position(t);
    let eyev = -ray.direction;
    let mut normalv = shape.normal_at_no_intersect(point);
    let inside = if normalv.dot(eyev) < 0.0 {
        normalv = -normalv;
        true
    } else {
        false
    };
    let reflectv = ray.direction.reflect(normalv);
    let over_point = point + normalv * crate::math::EPSILON;
    let under_point = point - normalv * crate::math::EPSILON;
    let object_point = Tuple::point(0.0, 0.0, 0.0);

    Computations {
        t,
        object: shape,
        point,
        over_point,
        under_point,
        object_point,
        eyev,
        normalv,
        reflectv,
        inside,
        n1: 1.0,
        n2: 1.0,
    }
}

/// Schlick computation to approximate Fresnel's equation.
///
/// returns a number between 0 and 1, inclusive. This number is called the
/// reflectance and represents what fraction of the light is reflected, given the
/// surface information at the hit
///
pub fn schlick(comps: &Computations) -> f64 {
    // find the cosine of the angle between the eye and normal vectors
    let mut cos = comps.eyev.dot(comps.normalv);

    // total internal reflection can only occur if n1 > n2
    if comps.n1 > comps.n2 {
        let n = comps.n1 / comps.n2;
        let sin2_t = n * n * (1.0 - cos * cos);
        if sin2_t > 1.0 {
            return 1.0;
        }

        // compute cosine of theta_t using trig identity and
        // when n1 > n2, use cos(theta_t) instead
        cos = (1.0 - sin2_t).sqrt();
    }

    let r0 = ((comps.n1 - comps.n2) / (comps.n1 + comps.n2)).powf(2.0);

    r0 + (1.0 - r0) * (1.0 - cos).powf(5.0)
}

// O(1): every caller passes the World arena (self.shapes), where ids are
// assigned 1..=N in insertion order, so id == index + 1 (see add_shape).
fn shape_by_id(shapes: &[Box<dyn Shape>], id: usize) -> &dyn Shape {
    shapes[id - 1].as_ref()
}

pub fn normal_at(shapes: &[Box<dyn Shape>], shape_id: usize, world_point: Tuple) -> Tuple {
    let object_point = world_to_object(shapes, shape_id, world_point);
    let object_normal = shape_by_id(shapes, shape_id).local_normal_at_no_hit(object_point);
    normal_to_world(shapes, shape_id, object_normal)
}

pub fn normal_to_world(shapes: &[Box<dyn Shape>], shape_id: usize, normal: Tuple) -> Tuple {
    let shape = shape_by_id(shapes, shape_id);

    let mut normal = shape.transform_inv().transpose() * normal;
    normal.w = 0.0;

    match shape.data().parent {
        Some(parent_id) => normal_to_world(shapes, parent_id, normal), // recurse LAST
        // Normalize ONCE, at the root — not per level. Intermediate normalizes
        // are mathematically redundant (they only rescale, which the final one
        // undoes) and each adds ~1 ULP of drift, so a BVH's extra identity
        // levels would otherwise perturb the color. This keeps it drift-free
        // and cheaper (one sqrt), and still returns a unit vector.
        None => normal.normalize(),
    }
}

fn world_to_object(shapes: &[Box<dyn Shape>], id: usize, point: Tuple) -> Tuple {
    let shape = shape_by_id(shapes, id);

    let point = match shape.data().parent {
        Some(parent_id) => world_to_object(shapes, parent_id, point),
        None => point,
    };
    *shape.transform_inv() * point
}

pub fn prepare_computations<'a>(
    intersection: Intersection,
    ray: &Ray,
    shapes: &'a [Box<dyn Shape>],
    xs: &Intersections,
) -> Computations<'a> {
    let shape = shape_by_id(shapes, intersection.object_id);
    let refractive_index_of =
        |id: usize| -> f64 { shape_by_id(shapes, id).data().material.refractive_index };

    let t = intersection.t;
    let point = ray.position(t);
    let eyev = -ray.direction;
    // group-aware: walks the parent chain (world_to_object -> local_normal_at -> normal_to_world)
    let object_point = world_to_object(shapes, intersection.object_id, point);
    let object_normal = shape.local_normal_at(object_point, intersection);
    let mut normalv = normal_to_world(shapes, intersection.object_id, object_normal);

    let inside = if normalv.dot(eyev) < 0.0 {
        normalv = -normalv;
        true
    } else {
        false
    };

    let reflectv = ray.direction.reflect(normalv);
    let over_point = point + normalv * crate::math::EPSILON;
    let under_point = point - normalv * crate::math::EPSILON;
    // The pattern samples object_point. Derive it from over_point rather than
    // the raw hit: on a plane the raw y is ±1e-16, so floor()-based patterns
    // (checkers) flip between 0 and -1 and speckle. over_point sits a
    // consistent EPSILON above the surface.
    let object_point = world_to_object(shapes, intersection.object_id, over_point);

    let hit = &intersection;
    let mut n1 = 1.0;
    let mut n2 = 1.0;

    let mut containers = [0usize; 32]; // stack-allocated, no heap
    let last = |c: &[usize], len: usize| -> Option<usize> {
        if len == 0 { None } else { Some(c[len - 1]) }
    };
    let mut len = 0;

    for i in xs.iter() {
        let is_hit = approx_eq(i.t, hit.t) && i.object_id == hit.object_id;

        if is_hit {
            n1 = last(&containers, len).map_or(1.0, refractive_index_of);
        };

        // toggle membership: already inside => we're EXITING; otherwise ENTERING
        if let Some(pos) = containers[..len].iter().position(|&id| id == i.object_id) {
            containers.copy_within(pos + 1..len, pos); // shift left, preserve order
            len -= 1;
        } else {
            debug_assert!(len < containers.len(), "container overflow");
            containers[len] = i.object_id;
            len += 1;
        }

        // n2 = material the ray is ENTERING (last container, after the toggle)
        if is_hit {
            n2 = last(&containers, len).map_or(1.0, refractive_index_of);
            break; // (4) stop at the hit
        }
    }

    Computations {
        t,
        object: shape,
        point,
        over_point,
        under_point,
        object_point,
        eyev,
        normalv,
        reflectv,
        inside,
        n1,
        n2,
    }
}

pub fn hexagon_corner(material: Material) -> Box<dyn Shape> {
    let mut corner = Sphere::new();
    corner.set_material(material);
    corner.set_transform(
        Matrix4::translation(0.0, 0.0, -1.0) * Matrix4::scaling(1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0),
    );
    Box::new(corner)
}

pub fn hexagon_edge(material: Material) -> Box<dyn Shape> {
    let mut edge = Cylinder::new();
    edge.set_material(material);
    edge.minimum = 0.0;
    edge.maximum = 1.0;
    edge.set_transform(
        Matrix4::translation(0.0, 0.0, -1.0)
            * Matrix4::rotation_y(-std::f64::consts::PI / 6.0)
            * Matrix4::rotation_z(-std::f64::consts::PI / 2.0)
            * Matrix4::scaling(1.0 / 4.0, 1.0, 1.0 / 4.0),
    );
    Box::new(edge)
}

pub fn hexagon(
    world: &mut World,
    group_id: usize,
    transform: Matrix4,
    material: Material,
) -> usize {
    let mut g_hexagon = Group::new();
    g_hexagon.set_transform(transform);

    let g_id_hexagon = if group_id == 0 {
        let g = world.add_shape(Box::new(g_hexagon)); // group in arena first → real id
        g
    } else {
        let g = world.add_child(group_id, Box::new(g_hexagon)); // group in arena first → real id
        g
    };

    for side_n in 0..6 {
        let mut g_side = Group::new();
        g_side.set_transform(Matrix4::rotation_y(
            side_n as f64 * std::f64::consts::PI / 3.0,
        ));
        let g_id_side = world.add_child(g_id_hexagon, Box::new(g_side));
        world.add_child(g_id_side, hexagon_corner(material.clone()));
        world.add_child(g_id_side, hexagon_edge(material.clone()));
    }

    g_id_hexagon
}

pub struct World {
    pub shapes: Vec<Box<dyn Shape>>,
    pub light: Option<Light>,
    next_id: AtomicUsize,
}

impl World {
    pub fn add_child(&mut self, group_id: usize, child: Box<dyn Shape>) -> usize {
        let child_id = self.add_shape(child); // final world id, in arena
        self.shape_by_id_mut(child_id).unwrap().data_mut().parent = Some(group_id);
        self.shape_by_id_mut(group_id)
            .unwrap()
            .add_child_id(child_id);
        child_id
    }

    /// Increment the shape id and add the shape to world.
    /// # Examples
    /// ```
    /// use rtc_rust::world::World;
    /// use rtc_rust::shapes::sphere::Sphere;
    ///
    /// let mut w = World::default_world();
    /// let ball = Sphere::new();
    /// let ball_id = w.add_shape(Box::new(ball));
    ///
    /// assert!(ball_id > 0);
    /// ```
    pub fn add_shape(&mut self, mut shape: Box<dyn Shape>) -> usize {
        let id = self.next_id.fetch_add(1, Ordering::Relaxed);
        shape.set_id(id);
        self.shapes.push(shape);
        // Upholds the id == index + 1 invariant that shape_by_id relies on.
        debug_assert_eq!(
            self.shapes.len(),
            id,
            "id must equal 1-based insertion index"
        );
        id
    }

    /// Compute and cache every group's bounding box, ready for the cull.
    ///
    /// Walks the tree post-order from the roots (following parent/child links,
    /// NOT id order), so it is correct even when a parent's id exceeds its
    /// children's — as happens once divide() appends sub-groups.
    pub fn build_bounds(&mut self) {
        let roots: Vec<usize> = self
            .shapes
            .iter()
            .filter(|s| s.data().parent.is_none())
            .map(|s| s.id())
            .collect();
        for r in roots {
            self.store_subtree_bounds(r);
        }
    }

    /// Post-order: compute this subtree's own-space box, store it on the group,
    /// and return it so the parent can fold it in. Leaves return their own box
    /// without storing.
    fn store_subtree_bounds(&mut self, id: usize) -> BoundingBox {
        let children = match shape_by_id(&self.shapes, id).children() {
            None => return shape_by_id(&self.shapes, id).bounds(), // leaf
            Some(c) => c.to_vec(),
        };
        let mut bb = BoundingBox::empty();
        for cid in children {
            let child_tf = *shape_by_id(&self.shapes, cid).transform();
            let child_bb = self.store_subtree_bounds(cid); // recurse FIRST (post-order)
            bb.add_box(&child_bb.transform(child_tf));
        }
        self.shape_by_id_mut(id).unwrap().set_bounds(bb);
        bb
    }

    /// Every group's bounding box as `(depth, 8 world-space corners)`, for debug
    /// wireframe rendering. `depth` is 0 for a root group and increments per
    /// nesting level (useful for coloring). Call after build_bounds() so the
    /// boxes are current. Infinite boxes (e.g. a group holding a plane) are
    /// skipped.
    pub fn group_world_boxes(&self) -> Vec<(usize, [Tuple; 8])> {
        let mut out = Vec::new();
        let roots: Vec<usize> = self
            .shapes
            .iter()
            .filter(|s| s.data().parent.is_none())
            .map(|s| s.id())
            .collect();
        for r in roots {
            self.collect_boxes(r, Matrix4::identity(), 0, &mut out);
        }
        out
    }

    fn collect_boxes(
        &self,
        id: usize,
        parent_tf: Matrix4,
        depth: usize,
        out: &mut Vec<(usize, [Tuple; 8])>,
    ) {
        let shape = shape_by_id(&self.shapes, id);
        let world_tf = parent_tf * *shape.transform(); // group-local -> world
        if let Some(children) = shape.children() {
            let bb = shape.bounds();
            let finite = [bb.min.x, bb.min.y, bb.min.z, bb.max.x, bb.max.y, bb.max.z]
                .iter()
                .all(|v| v.is_finite());
            if finite {
                let (lo, hi) = (bb.min, bb.max);
                // corner index bits: 0=x, 1=y, 2=z (0 => lo, 1 => hi)
                out.push((
                    depth,
                    [
                        world_tf * Tuple::point(lo.x, lo.y, lo.z),
                        world_tf * Tuple::point(hi.x, lo.y, lo.z),
                        world_tf * Tuple::point(lo.x, hi.y, lo.z),
                        world_tf * Tuple::point(hi.x, hi.y, lo.z),
                        world_tf * Tuple::point(lo.x, lo.y, hi.z),
                        world_tf * Tuple::point(hi.x, lo.y, hi.z),
                        world_tf * Tuple::point(lo.x, hi.y, hi.z),
                        world_tf * Tuple::point(hi.x, hi.y, hi.z),
                    ],
                ));
            }
            let kids = children.to_vec();
            for k in kids {
                self.collect_boxes(k, world_tf, depth + 1, out);
            }
        }
    }

    pub fn set_light(&mut self, light: Light) {
        self.light = Some(light)
    }

    pub fn default_world() -> Self {
        let mut w = Self::new();

        let light = Light::point_light(
            Tuple::point(-10.0, 10.0, -10.0),
            Tuple::color(1.0, 1.0, 1.0),
        );

        w.light = Some(light);

        let mut m = Material::new();
        m.color = Tuple::color(0.8, 1.0, 0.6);
        m.diffuse = 0.7;
        m.specular = 0.2;
        let mut s1 = Sphere::new();
        s1.set_material(m);
        w.add_shape(Box::new(s1));

        let mut s2 = Sphere::new();
        s2.set_transform(Matrix4::scaling(0.5, 0.5, 0.5));
        w.add_shape(Box::new(s2));

        w
    }

    pub fn color_at(&self, ray: &Ray, remaining: i32) -> Tuple {
        let xs = self.intersect(ray);
        match xs.hit() {
            None => Tuple::color(0.0, 0.0, 0.0),
            Some(hit) => match self.shapes.iter().find(|s| s.id() == hit.object_id) {
                None => Tuple::color(0.0, 0.0, 0.0),
                Some(_shape) => {
                    let comps = prepare_computations(hit, ray, &self.shapes, &xs);
                    self.shade_hit(&comps, remaining)
                }
            },
        }
    }

    fn shape_by_id(&self, id: usize) -> Option<&dyn Shape> {
        // O(1): id == index + 1 (see add_shape). id 0 is the "none" sentinel:
        // wrapping_sub(1) -> usize::MAX -> get() -> None.
        self.shapes.get(id.wrapping_sub(1)).map(|b| b.as_ref())
    }

    fn shape_by_id_mut(&mut self, id: usize) -> Option<&mut Box<dyn Shape>> {
        self.shapes.iter_mut().find(|s| s.id() == id)
    }

    pub fn load_obj_file<P: AsRef<Path>>(
        &mut self,
        path: P,
        transform: Matrix4,
    ) -> std::io::Result<usize> {
        let path: &Path = path.as_ref(); // one conversion, type pinned to Path
        println!("loading {}", path.display()); // .display() → printable
        let p = Parser::parse_obj_file(path)?; // reuse: Path is AsRef<Path>

        println!("REPORT: Input file {}. {}", path.display(), p.report());

        let mut top = Group::new();
        top.set_transform(transform);
        let top_id = self.add_shape(Box::new(top));

        // default group (ungrouped faces) → children of the top group
        for t in &p.default_group {
            self.add_child(top_id, Box::new(t.clone()));
        }
        // named groups → a sub-group each
        for tris in p.named_groups.values() {
            let sub = self.add_child(top_id, Box::new(Group::new()));
            for t in tris {
                self.add_child(sub, Box::new(t.clone()));
            }
        }
        Ok(top_id)
    }

    /// Set `material` on the shape `id` and, when it is a group, on every
    /// shape below it. Typical use: give a whole OBJ model its material after
    /// `load_obj_file`.
    /// # Examples
    /// ```
    /// use rtc_rust::material::Material;
    /// use rtc_rust::shapes::group::Group;
    /// use rtc_rust::shapes::sphere::Sphere;
    /// use rtc_rust::world::World;
    ///
    /// let mut w = World::new();
    /// let g_id = w.add_shape(Box::new(Group::new()));
    /// let ball_id = w.add_child(g_id, Box::new(Sphere::new()));
    ///
    /// let mut m = Material::new();
    /// m.color = rtc_rust::tuple::Tuple::color(1.0, 0.0, 0.0);
    /// w.set_material_recursive(g_id, &m);
    ///
    /// assert!(w.shapes[ball_id - 1].material().color.approx_eq(m.color));
    /// ```
    pub fn set_material_recursive(&mut self, id: usize, material: &Material) {
        let kids = self
            .shape_by_id(id)
            .and_then(|s| s.children().map(|c| c.to_vec()));
        if let Some(shape) = self.shape_by_id_mut(id) {
            shape.set_material(material.clone());
        }
        if let Some(kids) = kids {
            for cid in kids {
                self.set_material_recursive(cid, material);
            }
        }
    }

    /// Is `target_id` the shape `root_id` itself, or anywhere below it?
    /// Used by CSG to decide which side of the tree an intersection hit.
    fn subtree_includes(&self, root_id: usize, target_id: usize) -> bool {
        if root_id == target_id {
            return true;
        }
        match self.shape_by_id(root_id).and_then(|s| s.children()) {
            Some(kids) => kids
                .iter()
                .any(|&cid| self.subtree_includes(cid, target_id)),
            None => false,
        }
    }

    /// Keep only the intersections that lie on the boundary of the CSG
    /// solid, per its operation. Walks the t-sorted list tracking whether
    /// the ray is currently inside the left/right child.
    fn filter_intersections(&self, csg: &dyn Shape, xs: &Intersections) -> Intersections {
        let op = csg
            .csg_operation()
            .expect("filter_intersections needs a CSG shape");
        let left = csg.children().expect("a CSG node has children")[0];

        let mut inl = false; // inside the left child?
        let mut inr = false; // inside the right child?
        let mut result = Intersections::new();
        for i in xs.iter() {
            let lhit = self.subtree_includes(left, i.object_id);
            if intersection_allowed(op, lhit, inl, inr) {
                result.push(i);
            }
            if lhit {
                inl = !inl;
            } else {
                inr = !inr;
            }
        }
        result
    }

    fn intersect_node(&self, shape: &dyn Shape, ray: &Ray, xs: &mut Intersections) {
        record_node_visit();
        // transform the ray into THIS shape's object space
        let local_ray = shape.transform_inv().transform_ray(ray);

        match shape.children() {
            Some(children) => {
                if !shape.bounds().intersects(&local_ray) {
                    return; // ray can't hit anything in this group → skip subtree
                }
                if shape.csg_operation().is_some() {
                    // CSG: collect BOTH children's hits, then keep only the
                    // ones on the boundary of the combined solid
                    let mut sub = Intersections::new();
                    for &cid in children {
                        if let Some(child) = self.shape_by_id(cid) {
                            self.intersect_node(child, &local_ray, &mut sub);
                        }
                    }
                    for i in self.filter_intersections(shape, &sub).iter() {
                        xs.push(i);
                    }
                    return;
                }
                for &cid in children {
                    if let Some(child) = self.shape_by_id(cid) {
                        self.intersect_node(child, &local_ray, xs); // recurse in group space
                    }
                }
            }
            None => {
                record_prim_test();
                for i in shape.local_intersect(&local_ray).iter() {
                    xs.push(i); // leaf hit
                }
            }
        }
    }

    pub fn intersect(&self, ray: &Ray) -> Intersections {
        let mut xs = Intersections::new();
        for shape in &self.shapes {
            // roots only
            if shape.data().parent.is_none() {
                self.intersect_node(shape.as_ref(), ray, &mut xs);
            }
        }
        xs
    }

    /// Is anything between `point` and `light_position`? Takes the position
    /// explicitly so an area light can test each of its sample points.
    ///
    /// Any-hit query: stops at the first occluder instead of collecting and
    /// sorting every intersection like `intersect` does.
    pub fn is_shadowed(&self, light_position: Tuple, point: Tuple) -> bool {
        let v = light_position - point;
        let distance = Tuple::magnitude(v);
        let direction = Tuple::normalize(v);
        let r = Ray::new(point, direction);
        self.shapes
            .iter()
            .filter(|shape| shape.data().parent.is_none())
            .any(|shape| self.occluded_node(shape.as_ref(), &r, distance))
    }

    /// True if `shape` (or anything under it) blocks `ray` before `distance`.
    /// `t` is preserved by the object-space transform because the ray
    /// direction is not re-normalised, so the comparison is valid in any space.
    fn occluded_node(&self, shape: &dyn Shape, ray: &Ray, distance: f64) -> bool {
        let blocks = |i: &Intersection| i.t >= 0.0 && i.t < distance;
        match shape.children() {
            Some(_) if shape.csg_operation().is_some() => {
                // CSG needs the full, ordered intersection list to decide
                // which hits lie on the boundary of the combined solid.
                let mut xs = Intersections::new();
                self.intersect_node(shape, ray, &mut xs);
                xs.iter().any(|i| blocks(&i))
            }
            Some(children) => {
                record_node_visit();
                let local_ray = shape.transform_inv().transform_ray(ray);
                if !shape.bounds().intersects(&local_ray) {
                    return false;
                }
                children.iter().any(|&cid| {
                    self.shape_by_id(cid)
                        .is_some_and(|child| self.occluded_node(child, &local_ray, distance))
                })
            }
            None => {
                record_node_visit();
                record_prim_test();
                let local_ray = shape.transform_inv().transform_ray(ray);
                shape.local_occludes(&local_ray, distance)
            }
        }
    }

    pub fn new() -> Self {
        Self {
            shapes: Vec::new(),
            light: None,
            next_id: AtomicUsize::new(1),
        }
    }

    pub fn shade_hit(&self, comps: &Computations, remaining: i32) -> Tuple {
        match &self.light {
            Some(light) => {
                let intensity = light.intensity_at(comps.over_point, self);
                let surface = light.lighting(
                    comps.object,
                    comps.over_point,
                    comps.object_point,
                    comps.eyev,
                    comps.normalv,
                    intensity,
                );
                let reflected = self.reflected_color(comps, remaining);
                let refracted = self.refracted_color(comps, remaining);
                if comps.object.material().reflective > 0.0
                    && comps.object.material().transparency > 0.0
                {
                    let reflectance = schlick(comps);
                    return surface + reflected * reflectance + refracted * (1.0 - reflectance);
                }
                surface + reflected + refracted
            }
            None => Tuple::color(0.0, 0.0, 0.0),
        }
    }

    pub fn render_single(&self, camera: Camera) -> Canvas {
        let mut image = Canvas::new(camera.hsize, camera.vsize);

        for y in 0..camera.vsize {
            for x in 0..camera.hsize {
                let ray = camera.ray_for_pixel(x, y);
                let color = self.color_at(&ray, 10);
                image.write_pixel(x, y, color);
            }
        }

        image
    }

    pub fn reflected_color(&self, comps: &Computations, remaining: i32) -> Tuple {
        if remaining <= 0 {
            return tuple::colors::BLACK;
        }

        if approx_eq(comps.object.material().reflective, 0.0) {
            return tuple::colors::BLACK;
        }

        let reflect_ray = Ray::new(comps.over_point, comps.reflectv);
        let color = self.color_at(&reflect_ray, remaining - 1);
        color * comps.object.material().reflective
    }

    pub fn refracted_color(&self, comps: &Computations, remaining: i32) -> Tuple {
        if remaining <= 0 {
            return tuple::colors::BLACK;
        }

        if approx_eq(comps.object.material().transparency, 0.0) {
            return Tuple::color(0.0, 0.0, 0.0);
        }

        // Handle total internal reflection.
        // Find the ratio of first index of refraction to the second.
        // (Yup, this is inverted from the definition of Snell's Law.)
        let n_ratio = comps.n1 / comps.n2;

        // cos(theta_i) is the same as the dot product of the two vectors
        let cos_i = comps.eyev.dot(comps.normalv);

        // Find sin(theta_t)^2 via trigonometric identity
        let sin2_t = n_ratio * n_ratio * (1.0 - cos_i * cos_i);

        // Return black when there is total internal reflection.
        if sin2_t > 1.0 {
            return tuple::colors::BLACK;
        }

        // Find cos(theta_t) via trigonometric identity
        let cos_t = (1.0 - sin2_t).sqrt();

        // Compute the direction of the refracted ray
        let direction = comps.normalv * (n_ratio * cos_i - cos_t) - comps.eyev * n_ratio;

        // Create the refracted ray
        let refracted_ray = Ray::new(comps.under_point, direction);

        // Find the color of the refracted ray, making sure to multiply
        // by the transparency value to account for any opacity
        self.color_at(&refracted_ray, remaining - 1) * comps.object.material().transparency
    }

    pub fn render_parallel(&self, camera: Camera) -> Canvas {
        let width = camera.hsize;
        let height = camera.vsize;

        let mut pixels = vec![Tuple::color(0.0, 0.0, 0.0); width * height];
        pixels.par_iter_mut().enumerate().for_each(|(i, pixel)| {
            let x = i % width;
            let y = i / width;
            *pixel = self.color_at(&camera.ray_for_pixel(x, y), 10);
        });

        if camera.aa_samples > 1 {
            // Edge-detected supersampling: only pixels that differ from a
            // neighbour get the full sub-pixel grid, so flat areas cost one
            // ray and edges cost n*n.
            let edges: Vec<usize> = (0..pixels.len())
                .into_par_iter()
                .filter(|&i| Self::is_edge(&pixels, width, height, i, camera.aa_threshold))
                .collect();
            let resampled: Vec<(usize, Tuple)> = edges
                .into_par_iter()
                .map(|i| (i, self.supersample(&camera, i % width, i / width)))
                .collect();
            for (i, color) in resampled {
                pixels[i] = color;
            }
        }

        let mut image = Canvas::new(width, height);
        for (i, color) in pixels.into_iter().enumerate() {
            image.write_pixel(i % width, i / width, color);
        }
        image
    }

    /// Does pixel `i` differ from any 4-neighbour by more than `threshold`
    /// in some channel?
    fn is_edge(pixels: &[Tuple], width: usize, height: usize, i: usize, threshold: f64) -> bool {
        let (x, y) = (i % width, i / width);
        let differs = |j: usize| {
            let d = pixels[i] - pixels[j];
            d.x.abs() > threshold || d.y.abs() > threshold || d.z.abs() > threshold
        };
        (x > 0 && differs(i - 1))
            || (x + 1 < width && differs(i + 1))
            || (y > 0 && differs(i - width))
            || (y + 1 < height && differs(i + width))
    }

    /// Average of an `n` x `n` grid of sub-pixel samples for pixel `(x, y)`.
    pub fn supersample(&self, camera: &Camera, x: usize, y: usize) -> Tuple {
        let n = camera.aa_samples;
        let mut sum = Tuple::color(0.0, 0.0, 0.0);
        for sy in 0..n {
            for sx in 0..n {
                let dx = (sx as f64 + 0.5) / n as f64;
                let dy = (sy as f64 + 0.5) / n as f64;
                sum = sum + self.color_at(&camera.ray_for_subpixel(x, y, dx, dy), 10);
            }
        }
        sum / (n * n) as f64
    }

    pub fn render(&self, camera: Camera) -> Canvas {
        // self.render_single(camera) // change this one line to switch
        self.render_parallel(camera) // change this one line to switch
    }

    fn set_children(&mut self, id: usize, ids: Vec<usize>) {
        self.shape_by_id_mut(id).unwrap().set_children(ids)
    }

    fn partition_children(&mut self, group_id: usize) -> (Vec<usize>, Vec<usize>) {
        let (left_box, right_box) = self.subtree_bounds(group_id).split();

        // read current children (clone ids so we can mutate the group below)
        let children: Vec<usize> = shape_by_id(&self.shapes, group_id)
            .children()
            .unwrap()
            .to_vec();

        let (mut left, mut right, mut stay) = (Vec::new(), Vec::new(), Vec::new());
        for cid in children {
            let child = shape_by_id(&self.shapes, cid);
            // child's box IN THE GROUP'S space = child.transform * child-subtree bounds
            let cbox = self.subtree_bounds(cid).transform(*child.transform());
            if left_box.contains_box(&cbox) {
                left.push(cid);
            } else if right_box.contains_box(&cbox) {
                right.push(cid);
            } else {
                stay.push(cid); // straddles the split plane → keep in the group
            }
        }

        // group keeps only the straddlers; left/right get pulled out
        self.set_children(group_id, stay);
        (left, right)
    }

    fn subtree_bounds(&self, id: usize) -> BoundingBox {
        let shape = shape_by_id(&self.shapes, id);
        match shape.children() {
            None => shape.bounds(), // leaf: object-space box
            Some(children) => {
                let mut bb = BoundingBox::empty();
                for &cid in children {
                    let child = shape_by_id(&self.shapes, cid);
                    bb.add_box(&self.subtree_bounds(cid).transform(*child.transform()));
                }
                bb
            }
        }
    }

    fn make_subgroup(&mut self, parent_id: usize, child_ids: Vec<usize>) -> usize {
        // 1. append a new, identity-transform group to the arena
        let sub_id = self.add_shape(Box::new(Group::new()));

        // 2. move each orphaned child under the sub-group
        for &cid in &child_ids {
            self.shape_by_id_mut(cid).unwrap().data_mut().parent = Some(sub_id);
            self.shape_by_id_mut(sub_id).unwrap().add_child_id(cid);
        }

        // 3. the sub-group becomes a child of the parent
        self.shape_by_id_mut(sub_id).unwrap().data_mut().parent = Some(parent_id);
        self.shape_by_id_mut(parent_id)
            .unwrap()
            .add_child_id(sub_id);

        sub_id
    }

    pub fn divide(&mut self, group_id: usize, threshold: usize) {
        // only groups divide; only when they have enough children to be worth splitting
        if let Some(children) = shape_by_id(&self.shapes, group_id).children() {
            // a CSG node's two children ARE its left/right operands — never
            // repartition them; still recurse below (they may be big groups)
            let is_csg = shape_by_id(&self.shapes, group_id)
                .csg_operation()
                .is_some();
            if !is_csg && children.len() >= threshold {
                let (left, right) = self.partition_children(group_id);
                if !left.is_empty() {
                    self.make_subgroup(group_id, left);
                }
                if !right.is_empty() {
                    self.make_subgroup(group_id, right);
                }
            }

            // recurse into ALL current children (snapshot the ids first — the list
            // just changed, and we're about to mutate deeper)
            let kids: Vec<usize> = shape_by_id(&self.shapes, group_id)
                .children()
                .unwrap()
                .to_vec();
            for kid in kids {
                self.divide(kid, threshold);
            }
        }
        // primitives: no children → the `if let` is None → no-op
    }
}

impl Default for World {
    fn default() -> Self {
        Self::new()
    }
}

impl fmt::Display for World {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.light {
            Some(l) => writeln!(
                f,
                "{}World{} light: {}{}{}",
                Color::Yellow,
                Color::Reset,
                Color::Green,
                l,
                Color::Reset
            )?,
            None => writeln!(
                f,
                "{}World{} light: {}none{}",
                Color::Yellow,
                Color::Reset,
                Color::Red,
                Color::Reset
            )?,
        }
        writeln!(
            f,
            "  shapes: {}{}{}",
            Color::Green,
            self.shapes.len(),
            Color::Reset
        )?;
        Ok(())
    }
}

/// A transformation matrix—like scaling, rotation, and translation—that orients the world relative
/// to your eye, thus allowing you to line everything up and get exactly the shot that you need
pub fn view_transform(from: Tuple, to: Tuple, up: Tuple) -> Matrix4 {
    let forward = (to - from).normalize();
    let upn = up.normalize();
    let left = forward.cross(upn);
    let true_up = left.cross(forward);

    let orientation = Matrix4::new([
        [left.x, left.y, left.z, 0.0],
        [true_up.x, true_up.y, true_up.z, 0.0],
        [-forward.x, -forward.y, -forward.z, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]);

    orientation * Matrix4::translation(-from.x, -from.y, -from.z)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bounds::BoundingBox;
    use crate::intersection::{Intersection, Intersections};
    use crate::light::Sequence;
    use crate::math::{EPSILON, approx_eq};
    // use crate::obj::Parser;
    use crate::pattern::Pattern;
    use crate::patterns::blendedpattern::BlendedPattern;
    use crate::patterns::checkerspattern::CheckersPattern;
    use crate::patterns::gradientpattern::GradientPattern;
    use crate::patterns::nestedpattern::NestedPattern;
    use crate::patterns::ringpattern::RingPattern;
    use crate::patterns::stripepattern::*;
    use crate::patterns::testpattern::TestPattern;
    use crate::shape::Shape;
    use crate::shapes::cone::Cone;
    use crate::shapes::csg::{Csg, CsgOperation};
    use crate::shapes::cube::Cube;
    use crate::shapes::cylinder::Cylinder;
    use crate::shapes::group::Group;
    use crate::shapes::plane::Plane;
    use crate::shapes::triangle::Triangle;
    use crate::shapes::triangleuv::TriangleUV;
    use crate::tuple::colors::*;
    use crate::{loge, logi, tuple::Tuple};

    /// Chap 7 - Creating a world
    #[test]
    fn test_chap_7_1() -> Result<(), String> {
        let w = World::new();
        let chk = w.shapes.is_empty() && w.light.is_none();
        if chk {
            Ok(())
        } else {
            Err("Creating a world".into())
        }
    }

    /// Chap 7 - The default world
    #[test]
    fn test_chap_7_2() -> Result<(), String> {
        let w = World::default_world();

        let chk = !w.shapes.is_empty();
        let chk = chk
            && w.light.as_ref().unwrap().approx_eq(&Light::point_light(
                Tuple::point(-10.0, 10.0, -10.0),
                Tuple::color(1.0, 1.0, 1.0),
            ));

        let chk = chk && w.shapes.len() == 2;
        let chk = chk
            && w.shapes[0]
                .material()
                .color
                .approx_eq(Tuple::color(0.8, 1.0, 0.6));

        // The transform should be a scaling of 0.5 in all directions.
        let xform = w.shapes[1].transform();
        let chk = chk && approx_eq(xform[(0, 0)], 0.5);
        let chk = chk && approx_eq(xform[(1, 1)], 0.5);
        let chk = chk && approx_eq(xform[(2, 2)], 0.5);

        // The id's should be sphere 1 and sphere 2.
        let chk = chk && w.shapes[0].id() == 1;
        let chk = chk && w.shapes[1].id() == 2;

        if chk {
            Ok(())
        } else {
            loge!("test_chap_7_2", "world:{}", w);
            loge!(
                "test_chap_7_2",
                "w.shapes[0].material().color:{}",
                w.shapes[0].material().color
            );
            loge!("test_chap_7_2", "shape len : {}", w.shapes.len());
            loge!("test_chap_7_2", "shape id 1: {}", w.shapes[0].id());
            loge!("test_chap_7_2", "shape id 2: {}", w.shapes[1].id());
            loge!("test_chap_7_2", "xform[(0,0)]: {}", xform[(0, 0)]);
            loge!("test_chap_7_2", "xform[(1,1)]: {}", xform[(1, 1)]);
            loge!("test_chap_7_2", "xform[(2,2)]: {}", xform[(2, 2)]);
            Err("The default world".into())
        }
    }

    /// Chap 7 - Intersect a world with a ray
    #[test]
    fn test_chap_7_3() -> Result<(), String> {
        let w = World::default_world();

        logi!("test_chap_7_3", "world:{}", w);

        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs = w.intersect(&r);

        let chk = xs.count() == 4;
        let chk = chk && approx_eq(xs[0].t, 4.0);
        let chk = chk && approx_eq(xs[1].t, 4.5);
        let chk = chk && approx_eq(xs[2].t, 5.5);
        let chk = chk && approx_eq(xs[3].t, 6.0);
        if chk {
            Ok(())
        } else {
            Err("Intersect a world with a ray".into())
        }
    }

    /// Chap 7 - Precomputing the state of an intersection
    #[test]
    fn test_chap_7_4() -> Result<(), String> {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Sphere::new();
        let i = Intersection::new(4.0, shape.id());
        let comps = prepare_computations_upto_chap10(i, &r, &shape as &dyn Shape);

        let chk = approx_eq(comps.t, i.t);
        let chk = chk && comps.object.id() == shape.id();
        let chk = chk && comps.point.approx_eq(Tuple::point(0.0, 0.0, -1.0));
        let chk = chk && comps.eyev.approx_eq(Tuple::vector(0.0, 0.0, -1.0));
        let chk = chk && comps.normalv.approx_eq(Tuple::vector(0.0, 0.0, -1.0));
        let chk = chk && !comps.inside;
        if chk {
            Ok(())
        } else {
            Err("Precomputing the state of an intersection".into())
        }
    }

    /// Chap x - The hit, when an intersection occurs on the outside
    #[test]
    fn test_chap_7_5() -> Result<(), String> {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Sphere::new();
        let i = Intersection::new(4.0, shape.id());
        let comps = prepare_computations_upto_chap10(i, &r, &shape as &dyn Shape);

        let chk = !comps.inside;
        if chk {
            Ok(())
        } else {
            Err("The hit, when an intersection occurs on the outside".into())
        }
    }

    /// Chap x - The hit, when an intersection occurs on the inside
    #[test]
    fn test_chap_7_6() -> Result<(), String> {
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Sphere::new();
        let i = Intersection::new(1.0, shape.id());
        let comps = prepare_computations_upto_chap10(i, &r, &shape as &dyn Shape);

        let chk = Tuple::point(0.0, 0.0, 1.0).approx_eq(comps.point);
        let chk = chk && Tuple::vector(0.0, 0.0, -1.0).approx_eq(comps.eyev);
        let chk = chk && comps.inside;
        // normal would have been (0, 0, 1), but is inverted!
        let chk = chk && Tuple::vector(0.0, 0.0, -1.0).approx_eq(comps.normalv);
        if chk {
            Ok(())
        } else {
            Err("The hit, when an intersection occurs on the inside".into())
        }
    }

    /// Chap 7 - Shading an intersection
    #[test]
    fn test_chap_7_7() -> Result<(), String> {
        let w = World::default_world();
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = w.shapes[0].as_ref();
        let i = Intersection::new(4.0, shape.id());
        let comps = prepare_computations_upto_chap10(i, &r, shape as &dyn Shape);
        let c = w.shade_hit(&comps, 10);

        let chk = c.approx_eq(Tuple::color(0.380661193081, 0.475826491351, 0.285495894811));
        if chk {
            Ok(())
        } else {
            loge!("test_chap_7_7", "c: {}", c);
            Err("Shading an intersection".into())
        }
    }

    /// Chap 7 - Shading an intersection from the inside
    #[test]
    fn test_chap_7_8() -> Result<(), String> {
        let mut w = World::default_world();

        w.light = Some(Light::point_light(
            Tuple::point(0.0, 0.25, 0.0),
            Tuple::color(1.0, 1.0, 1.0),
        ));

        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = w.shapes[1].as_ref();
        let i = Intersection::new(0.5, shape.id());
        let comps = prepare_computations_upto_chap10(i, &r, shape as &dyn Shape);
        let c = w.shade_hit(&comps, 10);

        let chk = c.approx_eq(Tuple::color(0.904984472083, 0.904984472083, 0.904984472083));
        if chk {
            Ok(())
        } else {
            loge!("test_chap_7_8", "c: {}", c);
            Err("Shading an intersection from the inside".into())
        }
    }

    /// Chap 7 - The color when a ray misses
    #[test]
    fn test_chap_7_9() -> Result<(), String> {
        let w = World::default_world();
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 1.0, 0.0));
        let c = w.color_at(&r, 10);
        let chk = c.approx_eq(Tuple::color(0.0, 0.0, 0.0));
        if chk {
            Ok(())
        } else {
            Err("The color when a ray misses".into())
        }
    }
    /// Chap 7 - The color when a ray hits
    #[test]
    fn test_chap_7_10() -> Result<(), String> {
        let w = World::default_world();
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let c = w.color_at(&r, 10);
        let chk = c.approx_eq(Tuple::color(0.380661193081, 0.475826491351, 0.285495894811));
        if chk {
            Ok(())
        } else {
            loge!("test_chap_7_10", "c: {}", c);
            Err("The color when a ray hits".into())
        }
    }

    /// Chap 7 - The color with an intersection behind the ray
    #[test]
    fn test_chap_7_11() -> Result<(), String> {
        let mut w = World::default_world();

        {
            let outer = w.shapes[0].as_mut();
            let mut mat = outer.material().clone();
            mat.ambient = 1.0;
            outer.set_material(mat);
        }

        {
            let inner = w.shapes[1].as_mut();
            let mut mat = inner.material().clone();
            mat.ambient = 1.0;
            inner.set_material(mat);
        }

        let r = Ray::new(Tuple::point(0.0, 0.0, 0.75), Tuple::vector(0.0, 0.0, -1.0));
        let c = w.color_at(&r, 10);

        let inner = w.shapes[1].as_ref();
        let chk = c.approx_eq(inner.material().color);
        if chk {
            Ok(())
        } else {
            Err("The color with an intersection behind the ray".into())
        }
    }

    /// Chap 7 - The transformation matrix for the default orientation
    #[test]
    fn test_chap_7_12() -> Result<(), String> {
        let from = Tuple::point(0.0, 0.0, 0.0);
        let to = Tuple::point(0.0, 0.0, -1.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let t = view_transform(from, to, up);
        let chk = t.approx_eq(Matrix4::identity());
        if chk {
            Ok(())
        } else {
            Err("The transformation matrix for the default orientation".into())
        }
    }

    /// Chap 7 - A view transformation matrix looking in positive z direction
    #[test]
    fn test_chap_7_13() -> Result<(), String> {
        let from = Tuple::point(0.0, 0.0, 0.0);
        let to = Tuple::point(0.0, 0.0, 1.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let t = view_transform(from, to, up);
        let m = Matrix4::scaling(-1.0, 1.0, -1.0);
        let chk = t.approx_eq(m);
        if chk {
            Ok(())
        } else {
            loge!("test_chap_7_13", "\nt:{}\nm:{}", t, m);
            Err("A view transformation matrix looking in positive z direction".into())
        }
    }

    /// Chap 7 - The view transformation moves the world
    #[test]
    fn test_chap_7_14() -> Result<(), String> {
        let from = Tuple::point(0.0, 0.0, 8.0);
        let to = Tuple::point(0.0, 0.0, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let t = view_transform(from, to, up);
        let m = Matrix4::translation(0.0, 0.0, -8.0);
        let chk = t.approx_eq(m);
        if chk {
            Ok(())
        } else {
            Err("The view transformation moves the world".into())
        }
    }

    /// Chap 7 - An arbitrary view transformation
    #[test]
    fn test_chap_7_15() -> Result<(), String> {
        let from = Tuple::point(1.0, 3.0, 2.0);
        let to = Tuple::point(4.0, -2.0, 8.0);
        let up = Tuple::vector(1.0, 1.0, 0.0);
        let t = view_transform(from, to, up);
        let m = Matrix4::new([
            [
                -0.507092552837110,
                0.507092552837110,
                0.676123403782813,
                -2.366431913239846,
            ],
            [
                0.767715933859680,
                0.606091526731326,
                0.121218305346265,
                -2.828427124746189,
            ],
            [
                -0.358568582800318,
                0.597614304667197,
                -0.717137165600636,
                0.000000000000000,
            ],
            [
                0.000000000000000,
                0.000000000000000,
                0.000000000000000,
                1.000000000000000,
            ],
        ]);
        let chk = t.approx_eq(m);
        if chk {
            Ok(())
        } else {
            loge!("test_chap_7_15", "t:{}", t);
            Err("An arbitrary view transformation".into())
        }
    }

    /// Chap 7 - Rendering a world with a camera
    #[test]
    fn test_chap_7_22() -> Result<(), String> {
        // A flat-shaded sphere: every interior pixel is exactly the same
        // colour, so only its silhouette can be an edge.
        let mut w = World::new();
        w.set_light(Light::point_light(Tuple::point(-10.0, 10.0, -10.0), WHITE));
        let mut s = Sphere::new();
        let mut m = Material::new();
        m.color = Tuple::color(0.8, 0.2, 0.2);
        m.ambient = 1.0;
        m.diffuse = 0.0;
        m.specular = 0.0;
        s.set_material(m);
        w.add_shape(Box::new(s));
        let from = Tuple::point(0.0, 0.0, -5.0);
        let to = Tuple::point(0.0, 0.0, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let transform = view_transform(from, to, up);

        // Option 1 - unmutable camera
        let c = Camera::new(11, 11, std::f64::consts::PI / 2.0).with_transform(transform);

        // Option 2 - mutable camera
        // let mut c = Camera::new(11, 11, std::f64::consts::PI / 2.0);
        // c.set_transform(&view_transform(from, to, up));

        let image = w.render(c);
        let rc = image.write_ppm("test_chap_7_22.ppm");
        let chk = rc.is_ok();

        if chk {
            Ok(())
        } else {
            Err("Rendering a world with a camera".into())
        }
    }

    /// Chap 7 - Chapter 7 Putting It  Together
    #[test]
    fn test_chap_7_23_putting_it_together() -> Result<(), String> {
        let mut world = World::new();
        let light = Light::point_light(
            Tuple::point(-10.0, 10.0, -10.0),
            Tuple::color(1.0, 1.0, 1.0),
        );
        world.light = Some(light);

        let mut material = Material::new();
        material.color = Tuple::color(1.0, 0.9, 0.9);

        let mut floor = Sphere::new();
        floor.set_transform(Matrix4::scaling(10.0, 0.01, 10.0));
        material.diffuse = 0.7;
        material.specular = 0.3;
        floor.set_material(material.clone());

        let mut left_wall = Sphere::new();
        left_wall.set_transform(
            Matrix4::translation(0.0, 0.0, 5.0)
                * Matrix4::rotation_y(-std::f64::consts::PI / 4.0)
                * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(10.0, 0.01, 10.0),
        );
        left_wall.set_material(floor.material().clone());

        let mut right_wall = Sphere::new();
        right_wall.set_transform(
            Matrix4::translation(0.0, 0.0, 5.0)
                * Matrix4::rotation_y(std::f64::consts::PI / 4.0)
                * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(10.0, 0.01, 10.0),
        );
        right_wall.set_material(floor.material().clone());

        let mut middle = Sphere::new();
        middle.set_transform(Matrix4::translation(-0.5, 1.0, 0.5));
        material.color = Tuple::color(0.1, 1.0, 0.5);
        material.diffuse = 0.7;
        material.specular = 0.3;
        middle.set_material(material.clone());

        let mut right = Sphere::new();
        right.set_transform(Matrix4::translation(1.5, 0.5, -0.5) * Matrix4::scaling(0.5, 0.5, 0.5));
        material.color = Tuple::color(0.5, 1.0, 0.1);
        material.diffuse = 0.7;
        material.specular = 0.3;
        right.set_material(material.clone());

        let mut left = Sphere::new();
        left.set_transform(
            Matrix4::translation(-1.5, 0.33, -0.75) * Matrix4::scaling(0.33, 0.33, 0.33),
        );
        material.color = Tuple::color(1.0, 0.8, 0.1);
        material.diffuse = 0.7;
        material.specular = 0.3;
        left.set_material(material);

        let from = Tuple::point(0.0, 1.5, -5.0);
        let to = Tuple::point(0.0, 1.0, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let transform = view_transform(from, to, up);

        let camera = Camera::new(100, 50, std::f64::consts::PI / 3.0).with_transform(transform);

        world.add_shape(Box::new(floor));
        world.add_shape(Box::new(left_wall));
        world.add_shape(Box::new(right_wall));
        world.add_shape(Box::new(middle));
        world.add_shape(Box::new(left));
        world.add_shape(Box::new(right));

        let image = world.render(camera);
        let rc = image.write_ppm("test_chap_7_23_putting_it_together.ppm");
        let chk = rc.is_ok();

        if chk {
            Ok(())
        } else {
            Err("Chapter 7 Putting It  Together".into())
        }
    }

    /// Chap 8 - There is no shadow when nothing is collinear with point and light
    #[test]
    fn test_chap_8_2() -> Result<(), String> {
        let w = World::default_world();
        let p = Tuple::point(0.0, 10.0, 0.0);
        let is_shadowed = w.is_shadowed(w.light.as_ref().unwrap().position, p);

        let chk = !is_shadowed;
        if chk {
            Ok(())
        } else {
            Err("There is no shadow when nothing is collinear with point and light".into())
        }
    }

    /// Chap 8 - The shadow when an object is between the point and the light
    #[test]
    fn test_chap_8_3() -> Result<(), String> {
        let w = World::default_world();
        let p = Tuple::point(10.0, -10.0, 10.0);
        let is_shadowed = w.is_shadowed(w.light.as_ref().unwrap().position, p);

        let chk = is_shadowed;
        if chk {
            Ok(())
        } else {
            Err("The shadow when an object is between the point and the light".into())
        }
    }

    /// Chap 8 - There is no shadow when an object is behind the light
    #[test]
    fn test_chap_8_5() -> Result<(), String> {
        let w = World::default_world();
        let p = Tuple::point(0.0, 10.0, 0.0);
        let is_shadowed = w.is_shadowed(w.light.as_ref().unwrap().position, p);

        let chk = !is_shadowed;
        if chk {
            Ok(())
        } else {
            Err("There is no shadow when an object is behind the light".into())
        }
    }

    /// Chap 8 - There is no shadow when an object is behind the point
    #[test]
    fn test_chap_8_6() -> Result<(), String> {
        let w = World::default_world();
        let p = Tuple::point(0.0, 10.0, 0.0);
        let is_shadowed = w.is_shadowed(w.light.as_ref().unwrap().position, p);

        let chk = !is_shadowed;
        if chk {
            Ok(())
        } else {
            Err("There is no shadow when an object is behind the point".into())
        }
    }

    /// Chap 8 - shade_hit() is given an intersection in shadow
    #[test]
    fn test_chap_8_7() -> Result<(), String> {
        let mut w = World::new();
        let light = Light::point_light(Tuple::point(0.0, 0.0, -10.0), Tuple::color(1.0, 1.0, 1.0));
        w.set_light(light);
        let s1 = Sphere::new();
        w.add_shape(Box::new(s1));
        let mut s2 = Sphere::new();
        s2.set_transform(Matrix4::translation(0.0, 0.0, 10.0));
        w.add_shape(Box::new(s2.clone()));
        let r = Ray::new(Tuple::point(0.0, 0.0, 5.0), Tuple::vector(0.0, 0.0, 1.0));
        let i = Intersection::new(4.0, s2.id());
        let comps = prepare_computations_upto_chap10(i, &r, &s2 as &dyn Shape);
        let c = w.shade_hit(&comps, 10);

        let chk = c.approx_eq(Tuple::color(0.1, 0.1, 0.1));
        if chk {
            Ok(())
        } else {
            Err("shade_hit() is given an intersection in shadow".into())
        }
    }

    /// Chap 8 - The hit should offset the point
    #[test]
    fn test_chap_8_8() -> Result<(), String> {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let mut s1 = Sphere::new();
        s1.set_transform(Matrix4::translation(0.0, 0.0, 1.0));
        let i = Intersection::new(5.0, s1.id());
        let comps = prepare_computations_upto_chap10(i, &r, &s1 as &dyn Shape);
        let chk = comps.over_point.z < -crate::math::EPSILON / 2.0;
        let chk = chk && comps.point.z > comps.over_point.z;

        if chk {
            Ok(())
        } else {
            Err("The hit should offset the point".into())
        }
    }

    /// Chap 9 - Chapter 9 Putting It  Together
    #[test]
    fn test_chap_9_6_putting_it_together() -> Result<(), String> {
        let mut world = World::new();
        let light = Light::point_light(
            Tuple::point(-10.0, 10.0, -10.0),
            Tuple::color(1.0, 1.0, 1.0),
        );
        world.light = Some(light);

        let mut material = Material::new();
        material.color = Tuple::color(1.0, 0.9, 0.9);

        let mut floor = Plane::new();
        floor.set_transform(Matrix4::scaling(10.0, 0.01, 10.0));
        material.diffuse = 0.7;
        material.specular = 0.3;
        floor.set_material(material.clone());

        let mut left_wall = Plane::new();
        left_wall.set_transform(
            Matrix4::translation(0.0, 0.0, 5.0)
                * Matrix4::rotation_y(-std::f64::consts::PI / 4.0)
                * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(10.0, 0.01, 10.0),
        );
        left_wall.set_material(floor.material().clone());

        let mut right_wall = Plane::new();
        right_wall.set_transform(
            Matrix4::translation(0.0, 0.0, 5.0)
                * Matrix4::rotation_y(std::f64::consts::PI / 4.0)
                * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(10.0, 0.01, 10.0),
        );
        right_wall.set_material(floor.material().clone());

        let mut middle = Sphere::new();
        middle.set_transform(Matrix4::translation(-0.5, 1.0, 0.5));
        material.color = Tuple::color(0.1, 1.0, 0.5);
        material.diffuse = 0.7;
        material.specular = 0.3;
        middle.set_material(material.clone());

        let mut right = Sphere::new();
        right.set_transform(Matrix4::translation(1.5, 0.5, -0.5) * Matrix4::scaling(0.5, 0.5, 0.5));
        material.color = Tuple::color(0.5, 1.0, 0.1);
        material.diffuse = 0.7;
        material.specular = 0.3;
        right.set_material(material.clone());

        let mut left = Sphere::new();
        left.set_transform(
            Matrix4::translation(-1.5, 0.33, -0.75) * Matrix4::scaling(0.33, 0.33, 0.33),
        );
        material.color = Tuple::color(1.0, 0.8, 0.1);
        material.diffuse = 0.7;
        material.specular = 0.3;
        left.set_material(material);

        let from = Tuple::point(0.0, 1.5, -50.0);
        let to = Tuple::point(0.0, 1.0, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let transform = view_transform(from, to, up);

        let camera = Camera::new(100, 50, std::f64::consts::PI / 1.1).with_transform(transform);

        world.add_shape(Box::new(floor));
        world.add_shape(Box::new(left_wall));
        world.add_shape(Box::new(right_wall));
        world.add_shape(Box::new(middle));
        world.add_shape(Box::new(left));
        world.add_shape(Box::new(right));

        let image = world.render(camera);
        let rc = image.write_ppm("test_chap_9_6_putting_it_together.ppm");
        let chk = rc.is_ok();

        if chk {
            Ok(())
        } else {
            Err("Chapter 7 Putting It  Together".into())
        }
    }

    /// Chap 10 - Stripes with an object transformation
    #[test]
    fn test_chap_10_6() -> Result<(), String> {
        let mut s1 = Sphere::new();
        s1.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));

        let pattern = StripePattern::new(WHITE, BLACK);
        let world_point = Tuple::point(1.5, 0.0, 0.0);

        let c = pattern.color_at_shape(&s1, world_point);

        let chk = c.approx_eq(WHITE);
        if chk {
            Ok(())
        } else {
            Err("Stripes with an object transformation".into())
        }
    }

    /// Chap 10 - Stripes with a pattern transformation
    #[test]
    fn test_chap_10_7() -> Result<(), String> {
        let s1 = Sphere::new();

        let mut pattern = StripePattern::new(WHITE, BLACK);
        pattern.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));

        let world_point = Tuple::point(1.5, 0.0, 0.0);

        let c = pattern.color_at_shape(&s1, world_point);

        let chk = c.approx_eq(WHITE);

        if chk {
            Ok(())
        } else {
            Err("Stripes with a pattern transformation".into())
        }
    }

    /// Chap 10 - Stripes with both an object and a pattern transformation
    #[test]
    fn test_chap_10_8() -> Result<(), String> {
        let mut s1 = Sphere::new();
        s1.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));

        let mut pattern = StripePattern::new(WHITE, BLACK);
        pattern.set_transform(Matrix4::translation(0.5, 0.0, 0.0));

        let world_point = Tuple::point(2.5, 0.0, 0.0);

        let c = pattern.color_at_shape(&s1, world_point);

        let chk = c.approx_eq(WHITE);
        if chk {
            Ok(())
        } else {
            Err("Stripes with both an object and a pattern transformation".into())
        }
    }

    /// Chap 10 - A pattern with an object transformation
    #[test]
    fn test_chap_10_11() -> Result<(), String> {
        let mut shape = Sphere::new();
        shape.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));
        let pattern = TestPattern::new();

        let world_point = Tuple::point(2.0, 3.0, 4.0);

        let c = pattern.color_at_shape(&shape, world_point);

        let chk = c.approx_eq(Tuple::color(1.0, 1.5, 2.0));
        if chk {
            Ok(())
        } else {
            Err("A pattern with an object transformation".into())
        }
    }

    /// Chap x - A pattern with a pattern transformation
    #[test]
    fn test_chap_10_12() -> Result<(), String> {
        let shape = Sphere::new();
        let mut pattern = TestPattern::new();
        pattern.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));

        let world_point = Tuple::point(2.0, 3.0, 4.0);

        let c = pattern.color_at_shape(&shape, world_point);

        let chk = c.approx_eq(Tuple::color(1.0, 1.5, 2.0));

        if chk {
            Ok(())
        } else {
            Err("A pattern with a pattern transformation".into())
        }
    }

    /// Chap 10 - A pattern with both an object and a pattern transformation
    #[test]
    fn test_chap_10_13() -> Result<(), String> {
        let mut shape = Sphere::new();
        shape.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));
        let mut pattern = TestPattern::new();
        pattern.set_transform(Matrix4::translation(0.5, 1.0, 1.5));

        let world_point = Tuple::point(2.5, 3.0, 3.5);

        let c = pattern.color_at_shape(&shape, world_point);

        let chk = c.approx_eq(Tuple::color(0.75, 0.5, 0.25));
        if chk {
            Ok(())
        } else {
            Err("A pattern with both an object and a pattern transformation".into())
        }
    }

    /// Chap 10 - A gradient linearly interpolates between colors
    #[test]
    fn test_chap_10_14() -> Result<(), String> {
        let pattern = GradientPattern::new(WHITE, BLACK);
        let chk_0 = pattern.color_at(Tuple::point(0.0, 0.0, 0.0));
        let chk_25 = pattern.color_at(Tuple::point(0.25, 0.0, 0.0));
        let chk_50 = pattern.color_at(Tuple::point(0.50, 0.0, 0.0));
        let chk_75 = pattern.color_at(Tuple::point(0.75, 0.0, 0.0));

        let chk = chk_0.approx_eq(WHITE);
        let chk = chk && chk_25.approx_eq(Tuple::color(0.75, 0.75, 0.75));
        let chk = chk && chk_50.approx_eq(Tuple::color(0.5, 0.5, 0.5));
        let chk = chk && chk_75.approx_eq(Tuple::color(0.25, 0.25, 0.25));
        if chk {
            Ok(())
        } else {
            loge!("test_chap_10_14", "chk_0: {}", chk_0);
            loge!("test_chap_10_14", "chk_25: {}", chk_25);
            loge!("test_chap_10_14", "chk_50: {}", chk_50);
            loge!("test_chap_10_14", "chk_75: {}", chk_75);
            Err("A gradient linearly interpolates between colors".into())
        }
    }

    /// Chap 10 - A ring should extend in both x and z
    #[test]
    fn test_chap_10_15() -> Result<(), String> {
        let pattern = RingPattern::new(WHITE, BLACK);
        let chk_0 = pattern.color_at(Tuple::point(0.0, 0.0, 0.0));
        let chk_1 = pattern.color_at(Tuple::point(1.0, 0.0, 0.0));
        let chk_2 = pattern.color_at(Tuple::point(0.0, 0.0, 1.0));
        let chk_3 = pattern.color_at(Tuple::point(0.708, 0.0, 0.708));

        let chk = chk_0.approx_eq(WHITE);
        let chk = chk && chk_1.approx_eq(BLACK);
        let chk = chk && chk_2.approx_eq(BLACK);
        let chk = chk && chk_3.approx_eq(BLACK);
        if chk {
            Ok(())
        } else {
            loge!("test_chap_10_15", "chk_0: {}", chk_0);
            Err("A ring should extend in both x and z".into())
        }
    }

    /// Chap 10 - Checkers should repeat in x
    #[test]
    fn test_chap_10_16() -> Result<(), String> {
        let pattern = CheckersPattern::new(WHITE, BLACK);
        let chk_0 = pattern.color_at(Tuple::point(0.0, 0.0, 0.0));
        let chk_1 = pattern.color_at(Tuple::point(0.99, 0.0, 0.0));
        let chk_2 = pattern.color_at(Tuple::point(1.01, 0.0, 0.0));

        let chk = chk_0.approx_eq(WHITE);
        let chk = chk && chk_1.approx_eq(WHITE);
        let chk = chk && chk_2.approx_eq(BLACK);
        if chk {
            Ok(())
        } else {
            Err("Checkers should repeat in x".into())
        }
    }

    /// Chap 10 - Checkers should repeat in y
    #[test]
    fn test_chap_10_17() -> Result<(), String> {
        let pattern = CheckersPattern::new(WHITE, BLACK);
        let chk_0 = pattern.color_at(Tuple::point(0.0, 0.0, 0.0));
        let chk_1 = pattern.color_at(Tuple::point(0.0, 0.99, 0.0));
        let chk_2 = pattern.color_at(Tuple::point(0.0, 1.01, 0.0));

        let chk = chk_0.approx_eq(WHITE);
        let chk = chk && chk_1.approx_eq(WHITE);
        let chk = chk && chk_2.approx_eq(BLACK);
        if chk {
            Ok(())
        } else {
            Err("Checkers should repeat in y".into())
        }
    }

    /// Chap 10 - Checkers should repeat in z
    #[test]
    fn test_chap_10_18() -> Result<(), String> {
        let pattern = CheckersPattern::new(WHITE, BLACK);
        let chk_0 = pattern.color_at(Tuple::point(0.0, 0.0, 0.0));
        let chk_1 = pattern.color_at(Tuple::point(0.0, 0.0, 0.99));
        let chk_2 = pattern.color_at(Tuple::point(0.0, 0.0, 1.01));

        let chk = chk_0.approx_eq(WHITE);
        let chk = chk && chk_1.approx_eq(WHITE);
        let chk = chk && chk_2.approx_eq(BLACK);
        if chk {
            Ok(())
        } else {
            Err("Checkers should repeat in z".into())
        }
    }

    /// Chap 10 - Chapter 10 Putting It  Together
    #[test]
    fn test_chap_10_19() -> Result<(), String> {
        let mut world = World::new();
        let light = Light::point_light(
            Tuple::point(-10.0, 10.0, -10.0),
            Tuple::color(1.0, 1.0, 1.0),
        );
        world.light = Some(light);

        let mut pattern = GradientPattern::new(WHITE, BLACK);
        pattern.set_transform(Matrix4::scaling(0.25, 0.25, 0.25));

        let mut material = Material::new();
        material.color = Tuple::color(1.0, 0.9, 0.9);
        material.pattern = Some(Box::new(pattern));
        material.reflective = 0.5;

        let mut floor = Plane::new();
        floor.set_transform(Matrix4::scaling(10.0, 0.01, 10.0));
        material.diffuse = 0.7;
        material.specular = 0.3;
        floor.set_material(material.clone());

        let mut left_wall = Plane::new();
        left_wall.set_transform(
            Matrix4::translation(0.0, 0.0, 5.0)
                * Matrix4::rotation_y(-std::f64::consts::PI / 4.0)
                * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(10.0, 0.01, 10.0),
        );
        left_wall.set_material(floor.material().clone());

        let mut right_wall = Plane::new();
        right_wall.set_transform(
            Matrix4::translation(0.0, 0.0, 5.0)
                * Matrix4::rotation_y(std::f64::consts::PI / 4.0)
                * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(10.0, 0.01, 10.0),
        );
        right_wall.set_material(floor.material().clone());

        let mut middle = Sphere::new();
        middle.set_transform(Matrix4::translation(-0.5, 1.0, 0.5));
        material.color = Tuple::color(0.1, 1.0, 0.5);
        material.diffuse = 0.7;
        material.specular = 0.3;
        middle.set_material(material.clone());

        let mut right = Sphere::new();
        right.set_transform(Matrix4::translation(1.5, 0.5, -0.5) * Matrix4::scaling(0.5, 0.5, 0.5));
        material.color = Tuple::color(0.5, 1.0, 0.1);
        material.diffuse = 0.7;
        material.specular = 0.3;
        right.set_material(material.clone());

        let mut left = Sphere::new();
        left.set_transform(
            Matrix4::translation(-1.5, 0.33, -0.75) * Matrix4::scaling(0.33, 0.33, 0.33),
        );
        material.color = Tuple::color(1.0, 0.8, 0.1);
        material.diffuse = 0.7;
        material.specular = 0.3;
        left.set_material(material);

        // VIEW TRANSFORM SETTINGS
        let from = Tuple::point(0.0, 1.5, -12.0);
        let to = Tuple::point(0.0, 1.0, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let transform = view_transform(from, to, up);

        let camera = Camera::new(100, 50, std::f64::consts::PI / 3.0).with_transform(transform);

        world.add_shape(Box::new(floor));
        world.add_shape(Box::new(left_wall));
        world.add_shape(Box::new(right_wall));
        world.add_shape(Box::new(middle));
        world.add_shape(Box::new(left));
        world.add_shape(Box::new(right));

        let image = world.render(camera);
        let rc = image.write_ppm("test_chap_10_19_putting_it_together.ppm");
        let chk = rc.is_ok();

        if chk {
            Ok(())
        } else {
            Err("Chapter 10 Putting It  Together".into())
        }
    }

    /// Chap 10 - Nested patterns showcase (4K, full rabbit hole)
    ///
    /// Renders a scene that exercises every pattern type AND nesting:
    /// - floor: checkers nested with a gradient (alternating tiles)
    /// - walls: stripes nested with rings
    /// - spheres: gradients, rings, stripes, and a doubly-nested pattern
    ///
    /// Rendered at 4K (3840x2160) via the parallel renderer.
    #[test]
    fn test_chap_10_20_nested_showcase() -> Result<(), String> {
        // --- colors -------------------------------------------------------
        let red = Tuple::color(0.9, 0.1, 0.1);
        let green = Tuple::color(0.1, 0.9, 0.2);
        let blue = Tuple::color(0.1, 0.2, 0.9);
        let cyan = Tuple::color(0.1, 0.9, 0.9);
        let magenta = Tuple::color(0.9, 0.1, 0.9);
        let yellow = Tuple::color(0.95, 0.85, 0.1);
        let orange = Tuple::color(1.0, 0.55, 0.0);

        let mut world = World::new();
        world.light = Some(Light::point_light(
            Tuple::point(-10.0, 10.0, -10.0),
            Tuple::color(1.0, 1.0, 1.0),
        ));

        // --- floor: checkers tiles, each tile filled by a gradient --------
        let mut floor_checkers = CheckersPattern::new(WHITE, BLACK);
        floor_checkers.set_transform(Matrix4::scaling(0.5, 0.5, 0.5));
        let mut floor_gradient = GradientPattern::new(blue, cyan);
        floor_gradient.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));
        let mut floor_pattern =
            NestedPattern::new(Box::new(floor_checkers), Box::new(floor_gradient));
        floor_pattern.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));

        let mut floor_mat = Material::new();
        floor_mat.pattern = Some(Box::new(floor_pattern));
        floor_mat.diffuse = 0.7;
        floor_mat.specular = 0.1;

        let mut floor = Plane::new();
        floor.set_material(floor_mat);

        // --- back wall: stripes alternating with rings --------------------
        let mut wall_stripes = StripePattern::new(magenta, WHITE);
        wall_stripes.set_transform(Matrix4::scaling(0.25, 0.25, 0.25));
        let mut wall_rings = RingPattern::new(yellow, orange);
        wall_rings.set_transform(Matrix4::scaling(0.5, 0.5, 0.5));
        let mut wall_pattern = NestedPattern::new(Box::new(wall_stripes), Box::new(wall_rings));
        wall_pattern.set_transform(Matrix4::rotation_y(std::f64::consts::PI / 6.0));

        let mut wall_mat = Material::new();
        wall_mat.pattern = Some(Box::new(wall_pattern));
        wall_mat.specular = 0.0;

        let mut back_wall = Plane::new();
        back_wall.set_transform(
            Matrix4::translation(0.0, 0.0, 8.0) * Matrix4::rotation_x(std::f64::consts::PI / 2.0),
        );
        back_wall.set_material(wall_mat);

        // --- middle sphere: doubly-nested (stripe-of-gradients vs ring) ---
        let mut inner_stripe = StripePattern::new(red, green);
        inner_stripe.set_transform(Matrix4::scaling(0.25, 0.25, 0.25));
        let inner_gradient = GradientPattern::new(yellow, magenta);
        // first nest: stripes alternating with a gradient
        let mut nest_a = NestedPattern::new(Box::new(inner_stripe), Box::new(inner_gradient));
        nest_a.set_transform(Matrix4::scaling(0.5, 0.5, 0.5));
        let mut ring_child = RingPattern::new(cyan, blue);
        ring_child.set_transform(Matrix4::scaling(0.3, 0.3, 0.3));
        // second nest: the whole thing above alternating with a ring
        let mut middle_pattern = NestedPattern::new(Box::new(nest_a), Box::new(ring_child));
        middle_pattern.set_transform(
            Matrix4::scaling(0.6, 0.6, 0.6) * Matrix4::rotation_z(std::f64::consts::PI / 4.0),
        );

        let mut middle_mat = Material::new();
        middle_mat.pattern = Some(Box::new(middle_pattern));
        middle_mat.diffuse = 0.7;
        middle_mat.specular = 0.3;

        let mut middle = Sphere::new();
        middle.set_transform(Matrix4::translation(-0.5, 1.0, 0.5));
        middle.set_material(middle_mat);

        // --- right sphere: rotated gradient -------------------------------
        let mut right_grad = GradientPattern::new(green, magenta);
        right_grad.set_transform(
            Matrix4::scaling(0.5, 0.5, 0.5) * Matrix4::rotation_y(std::f64::consts::PI / 4.0),
        );
        let mut right_mat = Material::new();
        right_mat.pattern = Some(Box::new(right_grad));
        right_mat.diffuse = 0.7;
        right_mat.specular = 0.3;
        let mut right = Sphere::new();
        right.set_transform(Matrix4::translation(1.5, 0.5, -0.5) * Matrix4::scaling(0.5, 0.5, 0.5));
        right.set_material(right_mat);

        // --- left sphere: fine rings --------------------------------------
        let mut left_rings = RingPattern::new(red, yellow);
        left_rings.set_transform(Matrix4::scaling(0.15, 0.15, 0.15));
        let mut left_mat = Material::new();
        left_mat.pattern = Some(Box::new(left_rings));
        left_mat.diffuse = 0.7;
        left_mat.specular = 0.3;
        let mut left = Sphere::new();
        left.set_transform(
            Matrix4::translation(-1.5, 0.33, -0.75) * Matrix4::scaling(0.33, 0.33, 0.33),
        );
        left.set_material(left_mat);

        // --- extra ball #1: nested stripe/ring, blue family ---------------
        let mut b1_stripe = StripePattern::new(blue, cyan);
        b1_stripe.set_transform(Matrix4::scaling(0.2, 0.2, 0.2));
        let mut b1_ring = RingPattern::new(WHITE, blue);
        b1_ring.set_transform(Matrix4::scaling(0.2, 0.2, 0.2));
        let mut b1_pattern = BlendedPattern::new(Box::new(b1_stripe), Box::new(b1_ring));
        b1_pattern.set_transform(Matrix4::rotation_z(std::f64::consts::PI / 3.0));
        let mut b1_mat = Material::new();
        b1_mat.pattern = Some(Box::new(b1_pattern));
        b1_mat.diffuse = 0.7;
        b1_mat.specular = 0.4;
        let mut ball1 = Sphere::new();
        ball1.set_transform(
            Matrix4::translation(2.6, 0.75, 1.2) * Matrix4::scaling(0.75, 0.75, 0.75),
        );
        ball1.set_material(b1_mat);

        // --- extra ball #2: warm gradient ---------------------------------
        let mut b2_grad = GradientPattern::new(orange, red);
        b2_grad.set_transform(Matrix4::scaling(0.5, 0.5, 0.5));
        let mut b2_mat = Material::new();
        b2_mat.pattern = Some(Box::new(b2_grad));
        b2_mat.diffuse = 0.8;
        b2_mat.specular = 0.5;
        b2_mat.shininess = 300.0;
        let mut ball2 = Sphere::new();
        ball2.set_transform(Matrix4::translation(-2.7, 0.5, 0.3) * Matrix4::scaling(0.5, 0.5, 0.5));
        ball2.set_material(b2_mat);

        // --- camera (4K) --------------------------------------------------
        let from = Tuple::point(0.0, 1.5, -12.0);
        let to = Tuple::point(0.0, 1.0, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let transform = view_transform(from, to, up);
        let camera = Camera::new(96, 54, std::f64::consts::PI / 3.0).with_transform(transform);

        world.add_shape(Box::new(floor));
        world.add_shape(Box::new(back_wall));
        world.add_shape(Box::new(middle));
        world.add_shape(Box::new(left));
        world.add_shape(Box::new(right));
        world.add_shape(Box::new(ball1));
        world.add_shape(Box::new(ball2));

        let image = world.render_parallel(camera);
        let rc = image.write_ppm("test_chap_10_20_nested_showcase.ppm");

        if rc.is_ok() {
            Ok(())
        } else {
            Err("Chapter 10 Nested patterns showcase".into())
        }
    }

    /// Chap 11 - Precomputing the reflection vector
    #[test]
    fn test_chap_11_2() -> Result<(), String> {
        let sqrt2_over_2 = (1.0 / std::f64::consts::FRAC_1_SQRT_2) / 2.0;

        let shape = Plane::new();
        let r = Ray::new(
            Tuple::point(0.0, 1.0, 1.0),
            Tuple::vector(0.0, -sqrt2_over_2, sqrt2_over_2),
        );
        let i = Intersection::new((2.0_f64).sqrt(), shape.id());
        let comps = prepare_computations_upto_chap10(i, &r, &shape);
        let chk = comps
            .reflectv
            .approx_eq(Tuple::vector(0.0, sqrt2_over_2, sqrt2_over_2));

        if chk {
            Ok(())
        } else {
            loge!("test_chap_11_2", "comps.reflectv:{}", comps.reflectv);
            Err("Precomputing the reflection vector ".into())
        }
    }

    /// Chap 11 - The reflected color for a nonreflective material
    #[test]
    fn test_chap_11_3() -> Result<(), String> {
        let mut w = World::default_world();

        let point = Tuple::point(0.0, 0.0, 0.0);
        let direction = Tuple::vector(0.0, 0.0, 1.0);
        let r = Ray::new(point, direction);
        w.shapes[1].data_mut().material.ambient = 1.0;
        let i = Intersection::new(1.0, w.shapes[1].as_ref().id());
        let comps = prepare_computations_upto_chap10(i, &r, w.shapes[1].as_ref());
        let color = w.reflected_color(&comps, 10);

        let chk = color.approx_eq(BLACK);
        if chk {
            Ok(())
        } else {
            Err("The reflected color for a nonreflective material".into())
        }
    }

    /// Chap 11 - The reflected color for a reflective material
    #[test]
    fn test_chap_11_4() -> Result<(), String> {
        let mut w = World::default_world();

        let mut material = Material::new();
        material.reflective = 0.5;

        let mut shape = Plane::new();
        shape.set_transform(Matrix4::translation(0.0, -1.0, 0.0));
        shape.set_material(material.clone());

        w.add_shape(Box::new(shape));

        let sqrt2_over_2 = (1.0 / std::f64::consts::FRAC_1_SQRT_2) / 2.0;
        let point = Tuple::point(0.0, 0.0, -3.0);
        let direction = Tuple::vector(0.0, -sqrt2_over_2, sqrt2_over_2);
        let r = Ray::new(point, direction);

        let shape = w.shapes.last().expect("world must have at least one shape");
        let i = Intersection::new(2.0_f64.sqrt(), shape.id());
        let comps = prepare_computations_upto_chap10(i, &r, shape.as_ref());
        let color = w.reflected_color(&comps, 10);

        let chk = color.approx_eq(Tuple::color(0.190330596701, 0.237913245876, 0.142747947526));
        if chk {
            Ok(())
        } else {
            loge!("test_chap_11_4", "color:{}", color);
            Err("The reflected color for a reflective material".into())
        }
    }

    /// Chap 11 - shade_hit() with a reflective material
    #[test]
    fn test_chap_11_5() -> Result<(), String> {
        let mut w = World::default_world();

        let mut material = Material::new();
        material.reflective = 0.5;

        let mut shape = Plane::new();
        shape.set_transform(Matrix4::translation(0.0, -1.0, 0.0));
        shape.set_material(material.clone());

        w.add_shape(Box::new(shape));

        let sqrt2_over_2 = (1.0 / std::f64::consts::FRAC_1_SQRT_2) / 2.0;
        let point = Tuple::point(0.0, 0.0, -3.0);
        let direction = Tuple::vector(0.0, -sqrt2_over_2, sqrt2_over_2);
        let r = Ray::new(point, direction);

        let shape = w.shapes.last().expect("world must have at least one shape");
        let i = Intersection::new(2.0_f64.sqrt(), shape.id());
        let comps = prepare_computations_upto_chap10(i, &r, shape.as_ref());
        let color = w.shade_hit(&comps, 10);

        let chk = color.approx_eq(Tuple::color(0.876755985652, 0.924338634827, 0.829173336477));

        if chk {
            Ok(())
        } else {
            loge!("test_chap_11_5", "color:{}", color);
            Err("shade_hit() with a reflective material".into())
        }
    }

    /// Chap 11 - color_at() with mutually reflective surfaces Test #6: Avoid Infinite Recursion
    /// Show that your code safely handles infinite recursion caused by two objects that mutually
    /// reflect rays between themselves. Create two parallel mirrors by positioning one plane above
    /// another and making them both reflective. Orient a ray so that it strikes one plane and
    /// bounces to the other. What will happen?
    #[test]
    fn test_chap_11_6() -> Result<(), String> {
        let mut w = World::new(); //World::default_world();

        let light = Light::point_light(Tuple::point(0.0, 0.0, 0.0), Tuple::color(1.0, 1.0, 1.0));
        w.set_light(light);

        let mut material = Material::new();
        material.reflective = 1.0;

        let mut lower = Plane::new();
        lower.set_transform(Matrix4::translation(0.0, -1.0, 0.0));
        lower.set_material(material.clone());

        w.add_shape(Box::new(lower));

        let mut upper = Plane::new();
        upper.set_transform(Matrix4::translation(0.0, 1.0, 0.0));
        upper.set_material(material.clone());

        w.add_shape(Box::new(upper));

        let point = Tuple::point(0.0, 0.0, 0.0);
        let direction = Tuple::vector(0.0, 1.0, 0.0);
        let r = Ray::new(point, direction);

        let color = w.color_at(&r, 10);
        let chk = color.x >= 1.0;

        if chk {
            Ok(())
        } else {
            loge!("test_chap_11_6", "color:{}", color);
            Err("shade_hit() with a reflective material".into())
        }
    }

    /// Chap 11 - Transparency and Refractive Index for the default material
    #[test]
    fn test_chap_11_7() -> Result<(), String> {
        let material = Material::new();
        let chk =
            approx_eq(material.transparency, 0.0) && approx_eq(material.refractive_index, 1.0);
        if chk {
            Ok(())
        } else {
            loge!("test_chap_11_7", "material:{}", material);
            Err("Transparency and Refractive Index for the default material".into())
        }
    }

    /// Chap 11 - A helper for producing a sphere with a glassy material
    #[test]
    fn test_chap_11_8() -> Result<(), String> {
        let s = Sphere::glass();
        let chk = s.transform().approx_eq(Matrix4::identity())
            && approx_eq(s.data.material.transparency, 1.0)
            && approx_eq(s.data.material.refractive_index, 1.5);
        if chk {
            Ok(())
        } else {
            Err("A helper for producing a sphere with a glassy material".into())
        }
    }

    /// Chap 11 - Finding n1 and n2 at various intersections
    #[test]
    fn test_chap_11_9() -> Result<(), String> {
        let mut w = World::new();

        let mut a = Sphere::glass();
        a.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));
        a.data.material.refractive_index = 1.5;

        let mut b = Sphere::glass();
        b.set_transform(Matrix4::translation(0.0, 0.0, -0.25));
        b.data.material.refractive_index = 2.0;

        let mut c = Sphere::glass();
        c.set_transform(Matrix4::translation(0.0, 0.0, 0.25));
        c.data.material.refractive_index = 2.5;

        // add to the arena → world ids 1,2,3 (== index + 1)
        let a_id = w.add_shape(Box::new(a));
        let b_id = w.add_shape(Box::new(b));
        let c_id = w.add_shape(Box::new(c));

        let r = Ray::new(Tuple::point(0.0, 0.0, -4.0), Tuple::vector(0.0, 0.0, 1.0));
        let mut xs = Intersections::new();
        xs.push(Intersection::new(2.0, a_id));
        xs.push(Intersection::new(2.75, b_id));
        xs.push(Intersection::new(3.25, c_id));
        xs.push(Intersection::new(4.75, b_id));
        xs.push(Intersection::new(5.25, c_id));
        xs.push(Intersection::new(6.0, a_id));

        let expected = [
            (1.0, 1.5),
            (1.5, 2.0),
            (2.0, 2.5),
            (2.5, 2.5),
            (2.5, 1.5),
            (1.5, 1.0),
        ];

        for (idx, (n1, n2)) in expected.iter().enumerate() {
            let comps = prepare_computations(xs[idx], &r, &w.shapes, &xs);
            if !approx_eq(comps.n1, *n1) || !approx_eq(comps.n2, *n2) {
                return Err(format!(
                    "idx {idx}: got ({}, {}), expected ({n1}, {n2})",
                    comps.n1, comps.n2
                ));
            }
        }
        Ok(())
    }

    /// Chap x - The under point is offset below the surface
    #[test]
    fn test_chap_11_10() -> Result<(), String> {
        let mut w = World::new();
        let ray = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let mut shape = Sphere::glass();
        shape.set_transform(Matrix4::translation(0.0, 0.0, 1.0));
        let shape_id = w.add_shape(Box::new(shape));

        let i = Intersection::new(5.0, shape_id);
        let mut xs = Intersections::new();
        xs.push(i);

        let comps = prepare_computations(i, &ray, &w.shapes, &xs);

        let chk = comps.under_point.z > EPSILON / 2_f64 && comps.point.z < comps.under_point.z;

        if chk {
            Ok(())
        } else {
            Err(" The under point is offset below the surface".into())
        }
    }

    /// Chap 11 - The refracted color with an opaque surface
    #[test]
    fn test_chap_11_11() -> Result<(), String> {
        let w = World::default_world();
        let shape = &w.shapes[0];
        let ray = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));

        let mut xs = Intersections::new();
        xs.push(Intersection::new(4.0, shape.id()));
        xs.push(Intersection::new(6.0, shape.id()));

        let comps = prepare_computations(xs[0], &ray, &w.shapes, &xs);
        let c = w.refracted_color(&comps, 0);

        let chk = c.approx_eq(Tuple::color(0.0, 0.0, 0.0));
        if chk {
            Ok(())
        } else {
            Err("The refracted color with an opaque surface".into())
        }
    }

    /// Chap 11 - The refracted color at the maximum recursive depth
    #[test]
    fn test_chap_11_12() -> Result<(), String> {
        let mut w = World::default_world();
        let shape = w.shapes[0].as_mut();

        let mut m = Material::new();
        m.transparency = 1.0;
        m.refractive_index = 1.5;
        shape.set_material(m.clone());

        let ray = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));

        let mut xs = Intersections::new();
        xs.push(Intersection::new(4.0, shape.id()));
        xs.push(Intersection::new(6.0, shape.id()));

        let comps = prepare_computations(xs[0], &ray, &w.shapes, &xs);
        let c = w.refracted_color(&comps, 0);

        let chk = c.approx_eq(Tuple::color(0.0, 0.0, 0.0));
        if chk {
            Ok(())
        } else {
            Err("The refracted color at the maximum recursive depth".into())
        }
    }

    /// Chap x - The refracted color under total internal reflection
    #[test]
    fn test_chap_11_13() -> Result<(), String> {
        let mut w = World::default_world();
        let shape = w.shapes[0].as_mut();

        let mut m = Material::new();
        m.transparency = 1.0;
        m.refractive_index = 1.5;
        shape.set_material(m.clone());

        let sqrt2_over_2 = (1.0 / std::f64::consts::FRAC_1_SQRT_2) / 2.0;
        let ray = Ray::new(
            Tuple::point(0.0, 0.0, sqrt2_over_2),
            Tuple::vector(0.0, 1.0, 0.0),
        );

        let mut xs = Intersections::new();
        xs.push(Intersection::new(-sqrt2_over_2, shape.id()));
        xs.push(Intersection::new(sqrt2_over_2, shape.id()));

        // NOTE: this time you're inside the sphere, so you need
        // to look at the second intersection, xs[1], not xs[0]
        let comps = prepare_computations(xs[1], &ray, &w.shapes, &xs);
        let c = w.refracted_color(&comps, 5);

        let chk = c.approx_eq(Tuple::color(0.0, 0.0, 0.0));
        if chk {
            Ok(())
        } else {
            loge!("test_chap_11_13", "c:{}", c);
            Err("The refracted color under total internal reflection".into())
        }
    }

    /// Chap x - The refracted color with a refracted ray
    #[test]
    fn test_chap_11_14() -> Result<(), String> {
        let mut w = World::default_world();

        // Avoid mutating twice on the vector by
        // split once; `left` holds shapes[0], `right` holds shapes[1..]
        let (left, right) = w.shapes.split_at_mut(1);

        let a = left[0].as_mut();
        let b = right[0].as_mut();

        let mut m = Material::new();
        m.ambient = 1.0;
        m.pattern = Some(Box::new(TestPattern::new()));
        a.set_material(m.clone());

        let mut m = Material::new();
        m.transparency = 1.0;
        m.refractive_index = 1.5;
        b.set_material(m.clone());

        let ray = Ray::new(Tuple::point(0.0, 0.0, 0.1), Tuple::vector(0.0, 1.0, 0.0));

        let mut xs = Intersections::new();
        xs.push(Intersection::new(-0.9899, a.id()));
        xs.push(Intersection::new(-0.4899, b.id()));
        xs.push(Intersection::new(0.4899, b.id()));
        xs.push(Intersection::new(0.9899, a.id()));

        let comps = prepare_computations(xs[2], &ray, &w.shapes, &xs);
        let c = w.refracted_color(&comps, 5);

        let chk = c.approx_eq(Tuple::color(0.000000000000, 0.998884681786, 0.047216421860));
        if chk {
            Ok(())
        } else {
            loge!("test_chap_11_14", "refracted_color c: {}", c);
            Err("The refracted color with a refracted ray".into())
        }
    }

    /// Chap 11 - shade_hit() with a transparent material
    #[test]
    fn test_chap_11_15() -> Result<(), String> {
        let mut w = World::default_world();
        let mut floor = Plane::new();
        floor.set_transform(Matrix4::translation(0.0, -1.0, 0.0));

        let mut m = Material::new();
        m.transparency = 0.5;
        m.refractive_index = 1.5;
        floor.set_material(m.clone());

        let mut ball = Sphere::new();
        let mut m = Material::new();
        m.color = Tuple::color(1.0, 0.0, 0.0);
        m.ambient = 0.5;
        ball.set_transform(Matrix4::translation(0.0, -3.5, -0.5));
        ball.set_material(m.clone());

        let floor_id = w.add_shape(Box::new(floor));
        w.add_shape(Box::new(ball));

        let sqrt2_over_2 = (1.0 / std::f64::consts::FRAC_1_SQRT_2) / 2.0;
        let ray = Ray::new(
            Tuple::point(0.0, 0.0, -3.0),
            Tuple::vector(0.0, -sqrt2_over_2, sqrt2_over_2),
        );

        let mut xs = Intersections::new();
        xs.push(Intersection::new(2_f64.sqrt(), floor_id));

        let comps = prepare_computations(xs[0], &ray, &w.shapes, &xs);
        let color = w.shade_hit(&comps, 5);

        let chk = color.approx_eq(Tuple::color(0.936425388951, 0.686425388951, 0.686425388951));
        if chk {
            Ok(())
        } else {
            loge!("test_chap_11_15", "shade hit produced color: {}", color);
            Err("shade_hit() with a transparent material".into())
        }
    }

    /// Chap 11 - The Schlick approximation under total internal reflection
    #[test]
    fn test_chap_11_16() -> Result<(), String> {
        let mut w = World::default_world();
        let shape = Sphere::glass();
        let shape_id = w.add_shape(Box::new(shape));

        let sqrt2_over_2 = (1.0 / std::f64::consts::FRAC_1_SQRT_2) / 2.0;
        let ray = Ray::new(
            Tuple::point(0.0, 0.0, sqrt2_over_2),
            Tuple::vector(0.0, 1.0, 0.0),
        );

        let mut xs = Intersections::new();
        xs.push(Intersection::new(-sqrt2_over_2, shape_id));
        xs.push(Intersection::new(sqrt2_over_2, shape_id));

        let comps = prepare_computations(xs[1], &ray, &w.shapes, &xs);
        let reflectance = schlick(&comps);

        let chk = approx_eq(reflectance, 1.0);
        if chk {
            Ok(())
        } else {
            Err("The Schlick approximation under total internal reflection".into())
        }
    }

    /// Chap 11 - The Schlick approximation with a perpendicular viewing angle
    #[test]
    fn test_chap_11_17() -> Result<(), String> {
        let mut w = World::default_world();
        let shape = Sphere::glass();
        let shape_id = w.add_shape(Box::new(shape));

        let ray = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 1.0, 0.0));

        let mut xs = Intersections::new();
        xs.push(Intersection::new(-1.0, shape_id));
        xs.push(Intersection::new(1.0, shape_id));

        let comps = prepare_computations(xs[1], &ray, &w.shapes, &xs);
        let reflectance = schlick(&comps);

        let chk = approx_eq(reflectance, 0.04);
        if chk {
            Ok(())
        } else {
            loge!("test_chap_11_17", "reflectance: {}", reflectance);
            Err("The Schlick approximation with a perpendicular viewing angle".into())
        }
    }

    /// Chap 11 - The Schlick approximation with small angle and n2 > n1
    #[test]
    fn test_chap_11_18() -> Result<(), String> {
        let mut w = World::default_world();
        let shape = Sphere::glass();
        let shape_id = w.add_shape(Box::new(shape));

        let ray = Ray::new(Tuple::point(0.0, 0.99, -2.0), Tuple::vector(0.0, 0.0, 1.0));

        let mut xs = Intersections::new();
        xs.push(Intersection::new(1.8589, shape_id));

        let comps = prepare_computations(xs[0], &ray, &w.shapes, &xs);
        let reflectance = schlick(&comps);

        let chk = approx_eq(reflectance, 0.48873081012212183);
        if chk {
            Ok(())
        } else {
            loge!("test_chap_11_18", "reflectance: {}", reflectance);
            Err("The Schlick approximation with small angle and n2 > n1".into())
        }
    }

    /// Chap 11 - shade_hit() with a reflective, transparent material
    #[test]
    fn test_chap_11_19() -> Result<(), String> {
        let mut w = World::default_world();
        let mut floor = Plane::new();
        floor.set_transform(Matrix4::translation(0.0, -1.0, 0.0));

        let mut m = Material::new();
        m.reflective = 0.5;
        m.transparency = 0.5;
        m.refractive_index = 1.5;
        floor.set_material(m.clone());

        let mut ball = Sphere::new();
        let mut m = Material::new();
        m.color = Tuple::color(1.0, 0.0, 0.0);
        m.ambient = 0.5;
        ball.set_transform(Matrix4::translation(0.0, -3.5, -0.5));
        ball.set_material(m.clone());

        let floor_id = w.add_shape(Box::new(floor));
        w.add_shape(Box::new(ball));

        let sqrt2_over_2 = (1.0 / std::f64::consts::FRAC_1_SQRT_2) / 2.0;
        let ray = Ray::new(
            Tuple::point(0.0, 0.0, -3.0),
            Tuple::vector(0.0, -sqrt2_over_2, sqrt2_over_2),
        );

        let mut xs = Intersections::new();
        xs.push(Intersection::new(2_f64.sqrt(), floor_id));

        let comps = prepare_computations(xs[0], &ray, &w.shapes, &xs);
        let color = w.shade_hit(&comps, 5);

        let chk = color.approx_eq(Tuple::color(0.933915140526, 0.696434226271, 0.692430691343));

        if chk {
            Ok(())
        } else {
            loge!("test_chap_11_19", "color: {}", color);
            Err("shade_hit() with a reflective, transparent material".into())
        }
    }

    /// Chap 11 - Chapter 11 Putting It  Together
    #[test]
    fn test_chap_11_20() -> Result<(), String> {
        let mut world = World::new();
        let light = Light::point_light(
            Tuple::point(-10.0, 10.0, -10.0),
            Tuple::color(1.0, 1.0, 1.0),
        );
        world.light = Some(light);

        let mut pattern = GradientPattern::new(WHITE, BLACK);
        pattern.set_transform(Matrix4::scaling(0.25, 0.25, 0.25));

        let mut material = Material::new();
        material.color = Tuple::color(1.0, 0.9, 0.9);
        // material.pattern = Some(Box::new(pattern));
        material.reflective = 0.25;

        let mut floor = Plane::new();
        floor.set_transform(Matrix4::scaling(10.0, 0.01, 10.0));
        material.diffuse = 0.7;
        material.specular = 0.3;
        floor.set_material(material.clone());

        let mut left_wall = Plane::new();
        left_wall.set_transform(
            Matrix4::translation(0.0, 0.0, 5.0)
                * Matrix4::rotation_y(-std::f64::consts::PI / 4.0)
                * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(10.0, 0.01, 10.0),
        );
        left_wall.set_material(floor.material().clone());

        let mut right_wall = Plane::new();
        right_wall.set_transform(
            Matrix4::translation(0.0, 0.0, 5.0)
                * Matrix4::rotation_y(std::f64::consts::PI / 4.0)
                * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(10.0, 0.01, 10.0),
        );
        right_wall.set_material(floor.material().clone());

        let mut material = Material::new();
        material.color = Tuple::color(1.0, 0.9, 0.9);
        let mut pattern = CheckersPattern::new(WHITE, BLACK);
        pattern.set_transform(Matrix4::scaling(0.1, 0.1, 0.1));
        material.pattern = Some(Box::new(pattern));

        let mut middle = Sphere::new();
        middle.set_transform(Matrix4::translation(-0.5, 1.0, 0.5));
        material.color = Tuple::color(0.1, 1.0, 0.5);
        material.diffuse = 0.7;
        material.specular = 0.3;
        middle.set_material(material.clone());

        let mut right = Sphere::glass();
        right.set_transform(Matrix4::translation(1.5, 0.5, -0.5) * Matrix4::scaling(0.5, 0.5, 0.5));
        let mut material = Material::new();
        material.pattern = Some(Box::new(pattern));
        material.color = Tuple::color(0.82, 0.0, 0.0);
        material.transparency = 0.95;
        material.reflective = 1.0;
        material.shininess = 300.0;
        material.specular = 1.0;
        right.set_material(material.clone());

        let mut left = Sphere::new();
        left.set_transform(
            Matrix4::translation(-1.5, 0.33, -0.75) * Matrix4::scaling(0.33, 0.33, 0.33),
        );
        material.color = Tuple::color(1.0, 0.8, 0.1);
        material.diffuse = 0.7;
        material.specular = 0.3;
        left.set_material(material);

        // VIEW TRANSFORM SETTINGS
        let from = Tuple::point(0.0, 1.5, -12.0);
        let to = Tuple::point(0.0, 1.0, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let transform = view_transform(from, to, up);

        let camera = Camera::new(60, 40, std::f64::consts::PI / 3.0).with_transform(transform);

        world.add_shape(Box::new(floor));
        world.add_shape(Box::new(left_wall));
        world.add_shape(Box::new(right_wall));
        world.add_shape(Box::new(middle));
        world.add_shape(Box::new(left));
        world.add_shape(Box::new(right));

        let image = world.render(camera);
        let rc = image.write_ppm("test_chap_11_20_putting_it_together.ppm");
        let chk = rc.is_ok();

        if chk {
            Ok(())
        } else {
            Err("Chapter 11_20 Putting It  Together".into())
        }
    }

    /// Chap 12 - A ray intersects a cube
    #[test]
    fn test_chap_12_1() -> Result<(), String> {
        let c = Cube::new();

        // +x
        let r = Ray::new(Tuple::point(5.0, 0.5, 0.0), Tuple::vector(-1.0, 0.0, 0.0));
        let xs = c.local_intersect(&r);
        let chk = xs.count() == 2;
        let chk = chk && approx_eq(xs[0].t, 4.0) && approx_eq(xs[1].t, 6.0);

        // -x
        let r = Ray::new(Tuple::point(-5.0, 0.5, 0.0), Tuple::vector(1.0, 0.0, 0.0));
        let xs = c.local_intersect(&r);
        let chk = chk && approx_eq(xs[0].t, 4.0) && approx_eq(xs[1].t, 6.0);

        // +y
        let r = Ray::new(Tuple::point(0.5, 5.0, 0.0), Tuple::vector(0.0, -1.0, 0.0));
        let xs = c.local_intersect(&r);
        let chk = chk && approx_eq(xs[0].t, 4.0) && approx_eq(xs[1].t, 6.0);

        // -y
        let r = Ray::new(Tuple::point(0.5, -5.0, 0.0), Tuple::vector(0.0, 1.0, 0.0));
        let xs = c.local_intersect(&r);
        let chk = chk && approx_eq(xs[0].t, 4.0) && approx_eq(xs[1].t, 6.0);

        // +z
        let r = Ray::new(Tuple::point(0.5, 0.0, 5.0), Tuple::vector(0.0, 0.0, -1.0));
        let xs = c.local_intersect(&r);
        let chk = chk && approx_eq(xs[0].t, 4.0) && approx_eq(xs[1].t, 6.0);

        // -z
        let r = Ray::new(Tuple::point(0.5, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs = c.local_intersect(&r);
        let chk = chk && approx_eq(xs[0].t, 4.0) && approx_eq(xs[1].t, 6.0);

        // inside
        let r = Ray::new(Tuple::point(0.0, 0.5, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs = c.local_intersect(&r);
        let chk = chk && approx_eq(xs[0].t, -1.0) && approx_eq(xs[1].t, 1.0);

        if chk {
            Ok(())
        } else {
            Err("A ray intersects a cube".into())
        }
    }

    /// Chap 12 - A ray misses a cube
    #[test]
    fn test_chap_12_2() -> Result<(), String> {
        let c = Cube::new();

        let r = Ray::new(
            Tuple::point(-2.0, 0.0, 0.0),
            Tuple::vector(0.26730, 0.5345, 0.8018),
        );
        let xs = c.local_intersect(&r);
        let chk = xs.count() == 0;

        let r = Ray::new(
            Tuple::point(0.0, -2.0, 0.0),
            Tuple::vector(0.8018, 0.2673, 0.5345),
        );
        let xs = c.local_intersect(&r);
        let chk = chk && xs.count() == 0;

        let r = Ray::new(
            Tuple::point(0.0, 0.0, -2.0),
            Tuple::vector(0.5345, 0.8018, 0.2673),
        );
        let xs = c.local_intersect(&r);
        let chk = chk && xs.count() == 0;

        let r = Ray::new(Tuple::point(2.0, 0.0, 2.0), Tuple::vector(0.0, 0.0, -1.0));
        let xs = c.local_intersect(&r);
        let chk = chk && xs.count() == 0;

        let r = Ray::new(Tuple::point(0.0, 2.0, 2.0), Tuple::vector(0.0, -1.0, 0.0));
        let xs = c.local_intersect(&r);
        let chk = chk && xs.count() == 0;

        let r = Ray::new(Tuple::point(2.0, 2.0, 0.0), Tuple::vector(-1.0, 0.0, 0.0));
        let xs = c.local_intersect(&r);
        let chk = chk && xs.count() == 0;

        if chk {
            Ok(())
        } else {
            Err("A ray misses a cube".into())
        }
    }

    /// Chap 12 - The normal on the surface of a cube
    #[test]
    fn test_chap_12_3() -> Result<(), String> {
        let c = Cube::new();

        let p = Tuple::point(1.0, 0.5, -0.8);
        let normal = c.local_normal_at_no_hit(p);
        let chk = Tuple::vector(1.0, 0.0, 0.0).approx_eq(normal);

        let p = Tuple::point(-1.0, -0.5, 0.9);
        let normal = c.local_normal_at_no_hit(p);
        let chk = chk && Tuple::vector(-1.0, 0.0, 0.0).approx_eq(normal);

        let p = Tuple::point(-0.4, 1.0, -0.1);
        let normal = c.local_normal_at_no_hit(p);
        let chk = chk && Tuple::vector(0.0, 1.0, 0.0).approx_eq(normal);

        let p = Tuple::point(0.3, -1.0, -0.7);
        let normal = c.local_normal_at_no_hit(p);
        let chk = chk && Tuple::vector(0.0, -1.0, 0.0).approx_eq(normal);

        let p = Tuple::point(-0.6, 0.3, 1.0);
        let normal = c.local_normal_at_no_hit(p);
        let chk = chk && Tuple::vector(0.0, 0.0, 1.0).approx_eq(normal);

        let p = Tuple::point(0.4, 0.4, -1.0);
        let normal = c.local_normal_at_no_hit(p);
        let chk = chk && Tuple::vector(0.0, 0.0, -1.0).approx_eq(normal);

        let p = Tuple::point(1.0, 1.0, 1.0);
        let normal = c.local_normal_at_no_hit(p);
        let chk = chk && Tuple::vector(1.0, 0.0, 0.0).approx_eq(normal);

        let p = Tuple::point(-1.0, -1.0, -1.0);
        let normal = c.local_normal_at_no_hit(p);
        let chk = chk && Tuple::vector(-1.0, 0.0, 0.0).approx_eq(normal);

        if chk {
            Ok(())
        } else {
            Err("The normal on the surface of a cube".into())
        }
    }

    /// Chap 12 - Chapter 12 Putting It  Together
    #[test]
    fn test_chap_12_4() -> Result<(), String> {
        let mut world = World::new();
        let light = Light::point_light(
            Tuple::point(-10.0, 10.0, -10.0),
            Tuple::color(1.0, 1.0, 1.0),
        );
        world.light = Some(light);

        let mut pattern = GradientPattern::new(WHITE, BLACK);
        pattern.set_transform(Matrix4::scaling(0.25, 0.25, 0.25));

        let mut material = Material::new();
        material.color = Tuple::color(1.0, 0.9, 0.9);
        // material.pattern = Some(Box::new(pattern));
        material.reflective = 0.05;

        let mut floor = Plane::new();
        floor.set_transform(Matrix4::scaling(10.0, 0.01, 10.0));
        material.diffuse = 0.7;
        material.specular = 0.3;
        floor.set_material(material.clone());

        let mut left_wall = Plane::new();
        left_wall.set_transform(
            Matrix4::translation(0.0, 0.0, 5.0)
                * Matrix4::rotation_y(-std::f64::consts::PI / 4.0)
                * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(10.0, 0.01, 10.0),
        );
        left_wall.set_material(floor.material().clone());

        let mut right_wall = Plane::new();
        right_wall.set_transform(
            Matrix4::translation(0.0, 0.0, 5.0)
                * Matrix4::rotation_y(std::f64::consts::PI / 4.0)
                * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(10.0, 0.01, 10.0),
        );
        right_wall.set_material(floor.material().clone());

        let mut table_top = Cube::new();
        table_top.set_transform(
            Matrix4::translation(-1.0, 0.5, -6.0)
                // * Matrix4::rotation_y(std::f64::consts::PI / 4.0)
                * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(1.0, 1.0, 0.05),
        );
        let mut material = Material::new();
        material.color = Tuple::color(1.0, 0.0, 0.0);
        table_top.set_material(material.clone());

        let mut leg1 = Cube::new();
        leg1.set_transform(
            Matrix4::translation(-1.95, 0.0, -6.95)
                // * Matrix4::rotation_y(std::f64::consts::PI / 4.0)
                // * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(0.05, 0.5, 0.05),
        );
        let mut material = Material::new();
        material.color = Tuple::color(1.0, 1.0, 0.0);
        leg1.set_material(material.clone());

        let mut leg2 = leg1.clone();
        leg2.set_transform(
            Matrix4::translation(-0.05, 0.0, -6.95)
                // * Matrix4::rotation_y(std::f64::consts::PI / 4.0)
                // * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(0.05, 0.5, 0.05),
        );
        let mut material = Material::new();
        material.color = Tuple::color(0.0, 1.0, 1.0);
        leg2.set_material(material.clone());

        let mut x_pos: f64 = 1.5;
        let mut y_pos: f64 = 2.0;
        let mut z_pos: f64 = -1.95;
        let pos_incr: f64 = 0.2;
        let mut b = Cube::new();
        b.set_transform(
            Matrix4::translation(x_pos, 1.0, z_pos)
                // * Matrix4::rotation_y(std::f64::consts::PI / 4.0)
                // * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(0.05, 0.05, 0.05),
        );

        let mut vbb: Vec<Box<dyn Shape>> = vec![Box::new(b.clone())];
        loop {
            loop {
                x_pos += pos_incr;
                if x_pos > 3.0 {
                    x_pos = 1.5;
                    y_pos -= pos_incr / 2.0;
                    z_pos -= pos_incr;

                    let mut material = Material::new();
                    material.color = Tuple::color(x_pos / 10.0, y_pos / 10.0, y_pos / x_pos);
                    b.set_material(material.clone());

                    break;
                }

                b.set_transform(
                    Matrix4::translation(x_pos + pos_incr, y_pos, z_pos - pos_incr)
                // * Matrix4::rotation_y(std::f64::consts::PI / 4.0)
                // * Matrix4::rotation_x(-std::f64::consts::PI / 2.0)
                * Matrix4::scaling(0.05, 0.05, 0.05),
                );
                vbb.push(Box::new(b.clone()));
            }

            if z_pos < -4.0 {
                break;
            }
        }

        let from = Tuple::point(0.0, 1.5, -12.0);
        let to = Tuple::point(0.0, 1.0, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let transform = view_transform(from, to, up);

        let (display_x, display_y) = (60, 40);
        // let (display_x, display_y) = (3456, 2234);
        let camera =
            Camera::new(display_x, display_y, std::f64::consts::PI / 3.0).with_transform(transform);

        world.add_shape(Box::new(floor));
        world.add_shape(Box::new(left_wall));
        world.add_shape(Box::new(right_wall));
        world.add_shape(Box::new(table_top));
        world.add_shape(Box::new(leg1));
        world.add_shape(Box::new(leg2));
        for b_elem in vbb {
            world.add_shape(b_elem);
        }

        let image = world.render(camera);
        let rc = image.write_ppm("test_chap_12_4_putting_it_together.ppm");
        let chk = rc.is_ok();

        if chk {
            Ok(())
        } else {
            Err("Chapter 12_4 Putting It  Together".into())
        }
    }

    /// Chap 13 - A ray misses a cylinder
    #[test]
    fn test_chap_13_1() -> Result<(), String> {
        let cyl = Cylinder::new();

        let direction = Tuple::vector(0.0, 1.0, 0.0).normalize();
        let r = Ray::new(Tuple::point(1.0, 0.0, 0.0), direction);
        let xs = cyl.local_intersect(&r);
        let chk = xs.count() == 0;

        let direction = Tuple::vector(0.0, 1.0, 0.0).normalize();
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), direction);
        let xs = cyl.local_intersect(&r);
        let chk = chk && xs.count() == 0;

        let direction = Tuple::vector(1.0, 1.0, 1.0).normalize();
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), direction);
        let xs = cyl.local_intersect(&r);
        let chk = chk && xs.count() == 0;

        if chk {
            Ok(())
        } else {
            Err("A ray misses a cylinder".into())
        }
    }

    /// Chap 13 - A ray strikes a cylinder
    #[test]
    fn test_chap_13_2() -> Result<(), String> {
        let cyl = Cylinder::new();

        let direction = Tuple::vector(0.0, 0.0, 1.0).normalize();
        let r = Ray::new(Tuple::point(1.0, 0.0, -5.0), direction);
        let xs = cyl.local_intersect(&r);
        let chk = xs.count() == 2 && approx_eq(xs[0].t, 5.0) && approx_eq(xs[1].t, 5.0);

        let direction = Tuple::vector(0.0, 0.0, 1.0).normalize();
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), direction);
        let xs = cyl.local_intersect(&r);
        let chk = chk && xs.count() == 2 && approx_eq(xs[0].t, 4.0) && approx_eq(xs[1].t, 6.0);

        let direction = Tuple::vector(0.1, 1.0, 1.0).normalize();
        let r = Ray::new(Tuple::point(0.5, 0.0, -5.0), direction);
        let xs = cyl.local_intersect(&r);
        let chk = chk
            && xs.count() == 2
            && approx_eq(xs[0].t, 6.80798191702732)
            && approx_eq(xs[1].t, 7.088723439378861);

        if chk {
            Ok(())
        } else {
            loge!("test_chap_13_2", "hit t's:: t0:{} t1:{}", xs[0].t, xs[1].t);
            Err("A ray misses a cylinder".into())
        }
    }

    /// Chap 13 - Normal vector on a cylinder
    #[test]
    fn test_chap_13_3() -> Result<(), String> {
        let cyl = Cylinder::new();

        let n = cyl.local_normal_at_no_hit(Tuple::point(1.0, 0.0, 0.0));
        let chk = n.approx_eq(Tuple::vector(1.0, 0.0, 0.0));

        let n = cyl.local_normal_at_no_hit(Tuple::point(0.0, 5.0, -1.0));
        let chk = chk && n.approx_eq(Tuple::vector(0.0, 0.0, -1.0));

        let n = cyl.local_normal_at_no_hit(Tuple::point(0.0, -2.0, 1.0));
        let chk = chk && n.approx_eq(Tuple::vector(0.0, 0.0, 1.0));

        let n = cyl.local_normal_at_no_hit(Tuple::point(-1.0, 1.0, 0.0));
        let chk = chk && n.approx_eq(Tuple::vector(-1.0, 0.0, 0.0));

        if chk {
            Ok(())
        } else {
            Err("Normal vector on a cylinder".into())
        }
    }

    /// Chap 13 - The default minimum and maximum for a cylinder
    #[test]
    fn test_chap_13_4() -> Result<(), String> {
        let cyl = Cylinder::new();

        let chk = cyl.minimum == f64::NEG_INFINITY && cyl.maximum == f64::INFINITY;
        if chk {
            Ok(())
        } else {
            loge!("test_chap_13_4", "minimum: {}", cyl.minimum);
            Err("The default minimum and maximum for a cylinder".into())
        }
    }

    /// Chap 13 - Intersecting a constrained cylinder
    #[test]
    fn test_chap_13_5() -> Result<(), String> {
        let mut cyl = Cylinder::new();
        cyl.minimum = 1.0;
        cyl.maximum = 2.0;

        let origin = Tuple::point(0.0, 1.5, 0.0);
        let direction = Tuple::vector(0.1, 1.0, 0.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cyl.local_intersect(&r);
        let chk = xs.count() == 0;

        let origin = Tuple::point(0.0, 3.0, -5.0);
        let direction = Tuple::vector(0.0, 0.0, 1.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cyl.local_intersect(&r);
        let chk = chk && xs.count() == 0;

        let origin = Tuple::point(0.0, 0.0, -5.0);
        let direction = Tuple::vector(0.0, 0.0, 1.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cyl.local_intersect(&r);
        let chk = chk && xs.count() == 0;

        let origin = Tuple::point(0.0, 2.0, 0.0);
        let direction = Tuple::vector(0.0, 0.0, 1.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cyl.local_intersect(&r);
        let chk = chk && xs.count() == 0;

        let origin = Tuple::point(0.0, 1.0, 0.0);
        let direction = Tuple::vector(0.0, 0.0, 1.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cyl.local_intersect(&r);
        let chk = chk && xs.count() == 0;

        let origin = Tuple::point(0.0, 1.5, -2.0);
        let direction = Tuple::vector(0.0, 0.0, 1.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cyl.local_intersect(&r);
        let chk = chk && xs.count() == 2;

        if chk {
            Ok(())
        } else {
            Err("Intersecting a constrained cylinder".into())
        }
    }

    /// Chap 13 - The default closed value for a cylinder
    #[test]
    fn test_chap_13_6() -> Result<(), String> {
        let cyl = Cylinder::new();

        let chk = !cyl.closed;
        if chk {
            Ok(())
        } else {
            Err("The default closed value for a cylinder".into())
        }
    }

    /// Chap 13 - Intersecting the caps of a closed cylinder
    #[test]
    fn test_chap_13_7() -> Result<(), String> {
        let mut cyl = Cylinder::new();
        cyl.minimum = 1.0;
        cyl.maximum = 2.0;
        cyl.closed = true;

        // | point        | direction   | count |
        let origin = Tuple::point(0.0, 3.0, 0.0);
        let direction = Tuple::vector(0.0, -1.0, 0.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cyl.local_intersect(&r);
        let chk = xs.count() == 2;

        let origin = Tuple::point(0.0, 3.0, -2.0);
        let direction = Tuple::vector(0.0, -1.0, 2.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cyl.local_intersect(&r);
        let chk = chk && xs.count() == 2;

        let origin = Tuple::point(0.0, 4.0, -2.0); // corner case
        let direction = Tuple::vector(0.0, -1.0, 1.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cyl.local_intersect(&r);
        let chk = chk && xs.count() == 2;

        let origin = Tuple::point(0.0, 0.0, -2.0);
        let direction = Tuple::vector(0.0, 1.0, 2.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cyl.local_intersect(&r);
        let chk = chk && xs.count() == 2;

        let origin = Tuple::point(0.0, -1.0, -2.0); // corner case
        let direction = Tuple::vector(0.0, 1.0, 1.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cyl.local_intersect(&r);
        let chk = chk && xs.count() == 2;

        if chk {
            Ok(())
        } else {
            Err("Intersecting the caps of a closed cylinder".into())
        }
    }

    /// Chap 13 - The normal vector on a cylinder's end caps
    #[test]
    fn test_chap_13_8() -> Result<(), String> {
        let mut cyl = Cylinder::new();
        cyl.minimum = 1.0;
        cyl.maximum = 2.0;
        cyl.closed = true;

        let point = Tuple::point(0.0, 1.0, 0.0);
        let n = cyl.local_normal_at_no_hit(point);
        let chk = n.approx_eq(Tuple::vector(0.0, -1.0, 0.0));

        let point = Tuple::point(0.5, 1.0, 0.0);
        let n = cyl.local_normal_at_no_hit(point);
        let chk = chk && n.approx_eq(Tuple::vector(0.0, -1.0, 0.0));

        let point = Tuple::point(0.0, 1.0, 0.5);
        let n = cyl.local_normal_at_no_hit(point);
        let chk = chk && n.approx_eq(Tuple::vector(0.0, -1.0, 0.0));

        let point = Tuple::point(0.0, 2.0, 0.0);
        let n = cyl.local_normal_at_no_hit(point);
        let chk = chk && n.approx_eq(Tuple::vector(0.0, 1.0, 0.0));

        let point = Tuple::point(0.5, 2.0, 0.0);
        let n = cyl.local_normal_at_no_hit(point);
        let chk = chk && n.approx_eq(Tuple::vector(0.0, 1.0, 0.0));

        let point = Tuple::point(0.0, 2.0, 0.5);
        let n = cyl.local_normal_at_no_hit(point);
        let chk = chk && n.approx_eq(Tuple::vector(0.0, 1.0, 0.0));

        if chk {
            Ok(())
        } else {
            Err("The normal vector on a cylinder's end caps".into())
        }
    }

    /// Chap 13 - Intersecting a cone with a ray
    #[test]
    fn test_chap_13_9() -> Result<(), String> {
        let cone = Cone::new();

        let origin = Tuple::point(0.0, 0.0, -5.0);
        let direction = Tuple::vector(0.0, 0.0, 1.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cone.local_intersect(&r);
        let chk = xs.count() == 2 && approx_eq(xs[0].t, 5.0) && approx_eq(xs[1].t, 5.0);

        let origin = Tuple::point(0.0, 0.0, -5.0);
        let direction = Tuple::vector(1.0, 1.0, 1.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cone.local_intersect(&r);
        let chk = chk
            && xs.count() == 2
            && approx_eq(xs[0].t, 8.660254037844386)
            && approx_eq(xs[1].t, 8.660254037844386);

        let origin = Tuple::point(1.0, 1.0, -5.0);
        let direction = Tuple::vector(-0.5, -1.0, 1.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = cone.local_intersect(&r);
        let chk = chk
            && xs.count() == 2
            && approx_eq(xs[0].t, 4.550055679356349)
            && approx_eq(xs[1].t, 49.449944320643645);

        if chk {
            Ok(())
        } else {
            if xs.count() > 1 {
                loge!("test_chap_13_9", "xs[0].t:{} xs[1].t:{}", xs[0].t, xs[1].t);
            }
            Err("A ray misses a cone".into())
        }
    }

    /// Chap 13 - Intersecting a cone with a ray parallel to one of its halves
    #[test]
    fn test_chap_13_10() -> Result<(), String> {
        let shape = Cone::new();

        let origin = Tuple::point(0.0, 0.0, -1.0);
        let direction = Tuple::vector(0.0, 1.0, 1.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = shape.local_intersect(&r);

        let chk = xs.count() == 1 && approx_eq(xs[0].t, 0.3535533905932738);
        if chk {
            Ok(())
        } else {
            if xs.count() > 0 {
                loge!("test_chap_13_10", "xs[0].t:{} ", xs[0].t);
            }
            Err("Intersecting a cone with a ray parallel to one of its halves".into())
        }
    }

    /// Chap 13 - Intersecting a cone's end caps
    #[test]
    fn test_chap_13_11() -> Result<(), String> {
        let mut shape = Cone::new();
        shape.minimum = -0.5;
        shape.maximum = 0.5;
        shape.closed = true;

        let origin = Tuple::point(0.0, 0.0, -5.0);
        let direction = Tuple::vector(0.0, 1.0, 0.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = shape.local_intersect(&r);
        let chk = xs.count() == 0;

        let origin = Tuple::point(0.0, 0.0, -0.25);
        let direction = Tuple::vector(0.0, 1.0, 1.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = shape.local_intersect(&r);
        let chk = chk && xs.count() == 2;

        let origin = Tuple::point(0.0, 0.0, -0.25);
        let direction = Tuple::vector(0.0, 1.0, 0.0).normalize();
        let r = Ray::new(origin, direction);
        let xs = shape.local_intersect(&r);
        let chk = chk && xs.count() == 4;

        if chk {
            Ok(())
        } else {
            Err("Intersecting a cone's end caps".into())
        }
    }

    /// Chap 13 - Computing the normal vector on a cone
    #[test]
    fn test_chap_13_12() -> Result<(), String> {
        let shape = Cone::new();

        let point = Tuple::point(0.0, 0.0, 0.0);
        let n = shape.local_normal_at_no_hit(point);
        let chk = n.approx_eq(Tuple::vector(0.0, 0.0, 0.0));

        let point = Tuple::point(1.0, 1.0, 1.0);
        let n = shape.local_normal_at_no_hit(point);
        let chk = chk && n.approx_eq(Tuple::vector(1.0, -2.0_f64.sqrt(), 1.0));

        let point = Tuple::point(-1.0, -1.0, 0.0);
        let n = shape.local_normal_at_no_hit(point);
        let chk = chk && n.approx_eq(Tuple::vector(-1.0, 1.0, 0.0));

        if chk {
            Ok(())
        } else {
            loge!("test_chap_13_12", "n:{}", n);
            Err("Computing the normal vector on a cone".into())
        }
    }

    /// Chap 13 - Putting It Together
    ///
    /// The book leaves this scene open-ended ("render cylinders and cones"), so
    /// this composes a showcase from the cylinder features built in this
    /// chapter: capped solids, a truncated open tube, and a transformed
    /// cylinder, following the same layout as the earlier "putting it together"
    /// tests. Renders to a PPM in the working directory.
    #[test]
    fn test_chap_13_putting_it_all_together() -> Result<(), String> {
        let mut world = World::new();

        let light = Light::point_light(
            Tuple::point(-10.0, 10.0, -10.0),
            Tuple::color(1.0, 1.0, 1.0),
        );
        world.light = Some(light);

        let mut floor_checkers = CheckersPattern::new(WHITE, BLACK);
        floor_checkers.set_transform(Matrix4::scaling(0.5, 0.5, 0.5));

        // Floor - a slightly reflective plane
        let mut floor = Plane::new();
        let mut floor_material = Material::new();
        floor_material.color = Tuple::color(0.8, 0.8, 0.85);
        floor_material.specular = 0.0;
        floor_material.reflective = 0.2;
        floor_material.pattern = Some(Box::new(floor_checkers));
        floor.set_material(floor_material);
        world.add_shape(Box::new(floor));

        // Tall capped cylinder (green)
        let mut tall = Cylinder::new();
        tall.minimum = 0.0;
        tall.maximum = 3.0;
        tall.closed = true;
        tall.set_transform(Matrix4::translation(-1.5, 0.0, 0.5) * Matrix4::scaling(0.5, 1.0, 0.5));
        let mut tall_material = Material::new();
        tall_material.color = Tuple::color(0.1, 0.8, 0.3);
        tall_material.diffuse = 0.7;
        tall_material.specular = 0.3;
        tall.set_material(tall_material);
        world.add_shape(Box::new(tall));

        // Capped drum (red, slightly reflective)
        let mut drum = Cylinder::new();
        drum.minimum = 0.0;
        drum.maximum = 2.5;
        drum.closed = true;
        drum.set_transform(
            Matrix4::translation(2.3, 0.0, -0.5) * Matrix4::scaling(0.29, 1.0, 0.29),
        );
        let mut drum_material = Material::new();
        drum_material.color = Tuple::color(0.9, 0.2, 0.2);
        drum_material.diffuse = 0.7;
        drum_material.specular = 0.3;
        drum_material.reflective = 0.2;
        drum.set_material(drum_material);
        world.add_shape(Box::new(drum));

        let mut cyl_radius = 0.39;
        let mut cyl_maximum = 2.0;
        let mut col_factor = 1.0;
        loop {
            // Thin open tube () - not closed, so you can see through it
            let mut tube = Cylinder::new();
            tube.minimum = 0.0;
            tube.maximum = cyl_maximum;
            tube.closed = false;
            tube.set_transform(
                Matrix4::translation(2.3, 0.0, -0.5)
                    * Matrix4::scaling(cyl_radius, 1.0, cyl_radius),
            );

            let mut mirror = Material::new();
            mirror.color = Tuple::color(0.0, 0.0, 0.0); // near-black base; reflection provides the look
            mirror.ambient = 0.0;
            mirror.diffuse = 0.0;
            mirror.specular = 1.0; // bright highlight where the light hits
            mirror.shininess = 300.0; // tight, sharp highlight (mirror-like, not matte)
            mirror.reflective = 0.94 * col_factor; // perfect mirror; 0.9 for "very polished but not perfect"

            // let mut tube_material = Material::new();
            // tube_material.color =
            //     Tuple::color(0.2 / col_factor, 0.4 / col_factor, 0.19 / col_factor);
            // // tube_material.diffuse = 0.7 * col_factor;
            // tube_material.transparency = 0.95;
            // tube_material.reflective = 1.0;
            // tube_material.shininess = 300.0 * col_factor;
            // tube_material.specular = 1.0 * col_factor;
            // tube.set_material(tube_material);
            tube.set_material(mirror.clone());
            world.add_shape(Box::new(tube));

            if cyl_radius > 1.5 {
                break;
            }

            cyl_radius += 0.4;
            cyl_maximum -= 0.4;
            col_factor *= 0.9;
        }

        // Thin open tube (blue) - not closed, so you can see through it
        let mut tube = Cylinder::new();
        tube.minimum = 0.0;
        tube.maximum = 2.0;
        tube.closed = false;
        tube.set_transform(Matrix4::translation(0.4, 0.0, 1.6) * Matrix4::scaling(0.3, 1.0, 0.3));
        let mut tube_material = Material::new();
        tube_material.color = Tuple::color(0.2, 0.4, 0.9);
        tube_material.diffuse = 0.7;
        tube_material.specular = 0.3;
        tube.set_material(tube_material);
        world.add_shape(Box::new(tube));

        // Tilted capped cylinder lying on its side (yellow)
        let mut tilted = Cylinder::new();
        tilted.minimum = 0.0;
        tilted.maximum = 2.0;
        tilted.closed = true;
        tilted.set_transform(
            Matrix4::translation(0.0, 0.25, -1.5)
                * Matrix4::rotation_y(std::f64::consts::PI / 4.0)
                * Matrix4::rotation_z(std::f64::consts::PI / 4.0)
                * Matrix4::scaling(0.25, 1.5, 0.25),
        );
        let mut tilted_material = Material::new();
        tilted_material.color = Tuple::color(0.9, 0.8, 0.1);
        tilted_material.diffuse = 0.7;
        tilted_material.specular = 0.3;
        tilted.set_material(tilted_material);
        world.add_shape(Box::new(tilted));

        // Thin cone ()
        let mut cone = Cone::new();
        cone.minimum = -1.4;
        cone.maximum = 1.4;
        cone.closed = false;
        cone.set_transform(
            Matrix4::translation(-2.0, 1.0, -1.6)
                * Matrix4::rotation_y(std::f64::consts::PI / 4.0)
                * Matrix4::rotation_z(std::f64::consts::PI / 4.0)
                * Matrix4::scaling(0.3, 1.0, 0.3),
        );
        let mut cone_material = Material::new();
        cone_material.color = Tuple::color(0.7, 0.9, 0.3);
        cone_material.diffuse = 0.7;
        cone_material.specular = 0.3;
        cone.set_material(cone_material);
        world.add_shape(Box::new(cone));

        // View transform / camera
        let from = Tuple::point(0.0, 2.5, -7.0);
        let to = Tuple::point(0.0, 1.0, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let transform = view_transform(from, to, up);
        let (display_x, display_y) = (60, 40);
        // let (display_x, display_y) = (3456, 2234);
        let camera =
            Camera::new(display_x, display_y, std::f64::consts::PI / 3.0).with_transform(transform);

        let image = world.render(camera);
        let rc = image.write_ppm("test_chap_13_putting_it_all_together.ppm");

        if rc.is_ok() {
            Ok(())
        } else {
            Err("Chapter 13 Putting It Together".into())
        }
    }

    /// Chap 14 - A shape has a parent attribute
    #[test]
    fn test_chap_14_2() -> Result<(), String> {
        let s = Sphere::new();

        let chk = s.data().parent.is_none();
        if chk {
            Ok(())
        } else {
            Err("A shape has a parent attribute".into())
        }
    }

    /// Chap 14 - Adding a child to a group
    #[test]
    fn test_chap_14_3() -> Result<(), String> {
        let mut w = World::new();
        let g_id = w.add_shape(Box::new(Group::new()));
        let s_id = w.add_child(g_id, Box::new(Sphere::new()));

        // group now lists the child, and the child points back at the group
        let chk = w.shape_by_id(g_id).unwrap().children().unwrap()[0] == s_id;
        let chk = chk && w.shape_by_id(s_id).unwrap().data().parent == Some(g_id);

        if chk {
            Ok(())
        } else {
            Err("A shape has a parent attribute".into())
        }
    }

    /// Chap 14 - Intersecting a ray with an empty group
    #[test]
    fn test_chap_14_4() -> Result<(), String> {
        let g = Group::new();
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs = g.local_intersect(&r);

        let chk = xs.is_empty();
        if chk {
            Ok(())
        } else {
            loge!("test_chap_14_4", "xs.count:{}", xs.count());
            Err("Intersecting a ray with an empty group".into())
        }
    }

    /// Chap 14 - Intersecting a ray with a nonempty group
    #[test]
    fn test_chap_14_5() -> Result<(), String> {
        let mut w = World::new();
        let g_id = w.add_shape(Box::new(Group::new())); // group in arena first → real id

        let s1 = Sphere::new();
        let mut s2 = Sphere::new();
        s2.set_transform(Matrix4::translation(0.0, 0.0, -3.0));
        let mut s3 = Sphere::new();
        s3.set_transform(Matrix4::translation(5.0, 0.0, 0.0));

        let s1_id = w.add_child(g_id, Box::new(s1));
        let s2_id = w.add_child(g_id, Box::new(s2)); // capture the FINAL id
        let _s3_id = w.add_child(g_id, Box::new(s3));

        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs = w.intersect(&r);

        let chk = xs.count() == 4
            && xs[0].object_id == s2_id
            && xs[1].object_id == s2_id
            && xs[2].object_id == s1_id
            && xs[3].object_id == s1_id;

        if chk {
            Ok(())
        } else {
            loge!("test_chap_14_5", "xs.count:{}", xs.count(),);
            Err("Intersecting a ray with a nonempty group".into())
        }
    }

    /// Chap 14 - Intersecting a transformed group
    #[test]
    fn test_chap_14_6() -> Result<(), String> {
        let mut w = World::new();
        let mut g = Group::new();
        g.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));
        let g_id = w.add_shape(Box::new(g)); // group in arena first → real id

        let mut s = Sphere::new();
        s.set_transform(Matrix4::translation(5.0, 0.0, 0.0));

        let _s_id = w.add_child(g_id, Box::new(s));

        let r = Ray::new(Tuple::point(10.0, 0.0, -10.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs = w.intersect(&r);

        let chk = xs.count() == 2;

        if chk {
            Ok(())
        } else {
            loge!("test_chap_14_6", "xs.count:{}", xs.count(),);
            Err("Intersecting a ray with a nonempty group".into())
        }
    }

    /// Chap 14 - Converting a point from world to object space
    #[test]
    fn test_chap_14_7() -> Result<(), String> {
        let mut w = World::new();

        let mut g1 = Group::new();
        g1.set_transform(Matrix4::rotation_y(std::f64::consts::PI / 2.0));
        let g1_id = w.add_shape(Box::new(g1)); // group in arena first → real id

        let mut g2 = Group::new();
        g2.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));
        let g2_id = w.add_child(g1_id, Box::new(g2));

        let mut s = Sphere::new();
        s.set_transform(Matrix4::translation(5.0, 0.0, 0.0));
        let s_id = w.add_child(g2_id, Box::new(s));

        let p = world_to_object(&w.shapes, s_id, Tuple::point(-2.0, 0.0, -10.0));
        let chk = p.approx_eq(Tuple::point(0.0, 0.0, -1.0));

        if chk {
            Ok(())
        } else {
            loge!("test_chap_14_7", "p: {}", p);
            Err("Converting a point from world to object space".into())
        }
    }

    /// Chap x - Converting a normal from object to world space
    #[test]
    fn test_chap_14_8() -> Result<(), String> {
        let mut w = World::new();

        let mut g1 = Group::new();
        g1.set_transform(Matrix4::rotation_y(std::f64::consts::PI / 2.0));
        let g1_id = w.add_shape(Box::new(g1)); // group in arena first → real id

        let mut g2 = Group::new();
        g2.set_transform(Matrix4::scaling(1.0, 2.0, 3.0));
        let g2_id = w.add_child(g1_id, Box::new(g2));

        let mut s = Sphere::new();
        s.set_transform(Matrix4::translation(5.0, 0.0, 0.0));
        let s_id = w.add_child(g2_id, Box::new(s));

        let sqrt3_3 = 3_f64.sqrt() / 3.0;
        let normal = Tuple::vector(sqrt3_3, sqrt3_3, sqrt3_3);
        let n = normal_to_world(&w.shapes, s_id, normal);

        let chk = n.approx_eq(Tuple::vector(
            2.0 / 7.0,  // 0.285714285714,
            3.0 / 7.0,  // 0.428571428571,
            -6.0 / 7.0, // -0.857142857143,
        ));

        if chk {
            Ok(())
        } else {
            loge!("test_chap_14_8", "n:{}", n);
            Err("Converting a normal from object to world space".into())
        }
    }

    /// Chap 14 - Finding the normal on a child object
    #[test]
    fn test_chap_14_9() -> Result<(), String> {
        let mut w = World::new();

        let mut g1 = Group::new();
        g1.set_transform(Matrix4::rotation_y(std::f64::consts::PI / 2.0));
        let g1_id = w.add_shape(Box::new(g1)); // group in arena first → real id

        let mut g2 = Group::new();
        g2.set_transform(Matrix4::scaling(1.0, 2.0, 3.0));
        let g2_id = w.add_child(g1_id, Box::new(g2));

        let mut s = Sphere::new();
        s.set_transform(Matrix4::translation(5.0, 0.0, 0.0));
        let s_id = w.add_child(g2_id, Box::new(s));

        let n = normal_at(&w.shapes, s_id, Tuple::point(1.7321, 1.1547, -5.5774));

        let chk = n.approx_eq(Tuple::vector(
            0.285703681841,
            0.428543151781,
            -0.857160529448,
        ));

        if chk {
            Ok(())
        } else {
            loge!("test_chap_14_9", "n:{}", n);
            Err("Finding the normal on a child object".into())
        }
    }

    /// Chap 14 - A pattern on a shape inside a transformed group sees the group's transform
    #[test]
    fn test_chap_14_10() -> Result<(), String> {
        let mut w = World::new();

        let mut g = Group::new();
        g.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));
        let g_id = w.add_shape(Box::new(g));

        let mut s = Sphere::new();
        s.set_transform(Matrix4::translation(2.0, 0.0, 0.0));
        let mut m = Material::new();
        m.pattern = Some(Box::new(StripePattern::new(WHITE, BLACK)));
        s.set_material(m);
        let s_id = w.add_child(g_id, Box::new(s));

        // Sphere ends up at world (4,0,0) with radius 2; hit its +x pole at (6,0,0).
        let r = Ray::new(Tuple::point(10.0, 0.0, 0.0), Tuple::vector(-1.0, 0.0, 0.0));
        let xs = w.intersect(&r);
        let hit = xs.hit().ok_or("expected a hit")?;
        let comps = prepare_computations(hit, &r, &w.shapes, &xs);

        // world (6,0,0) → group⁻¹ (scale ½) → (3,0,0) → sphere⁻¹ (translate -2) → (1,0,0)
        let chk = comps.object_point.approx_eq(Tuple::point(1.0, 0.0, 0.0));
        let chk = chk && comps.object.id() == s_id; // the hit is the nested sphere

        if chk {
            Ok(())
        } else {
            loge!("test_chap_14_10", "object_point: {}", comps.object_point);
            Err("Pattern on nested shape must see the group transform".into())
        }
    }

    /// Chap 14 - Create a hexagon
    #[test]
    fn test_chap_14_create_hexagon() -> Result<(), String> {
        let mut w = World::new();
        let light = Light::point_light(
            Tuple::point(-10.0, 10.0, -10.0),
            Tuple::color(1.0, 1.0, 1.0),
        );
        w.light = Some(light);

        let mut g_top = Group::new();
        g_top.set_transform(
            Matrix4::scaling(0.8, 0.8, 0.8) * Matrix4::rotation_x(-std::f64::consts::PI / 4.0),
        );
        let g_id_top = w.add_shape(Box::new(g_top));

        let mut material = Material::new();
        material.color = Tuple::color(1.0, 0.0, 0.0);

        let _g_id = hexagon(
            &mut w,
            0,
            Matrix4::translation(-1.0, 1.0, -3.0)
                * Matrix4::scaling(0.5, 0.5, 0.5)
                * Matrix4::rotation_x(std::f64::consts::PI / 2.0)
                * Matrix4::rotation_z(std::f64::consts::PI / 4.0),
            material.clone(),
        );

        material.color = Tuple::color(0.75, 0.5, 0.0);
        let _g_id = hexagon(
            &mut w,
            g_id_top,
            Matrix4::translation(1.0, 1.0, 3.0)
                * Matrix4::scaling(0.75, 0.75, 0.75)
                * Matrix4::rotation_x(0.0 * std::f64::consts::PI / 2.0)
                * Matrix4::rotation_z(0.0 * std::f64::consts::PI / 4.0),
            material.clone(),
        );

        let _g_id = hexagon(
            &mut w,
            g_id_top,
            Matrix4::translation(1.0, 1.5, 3.0)
                * Matrix4::scaling(0.75, 0.75, 0.75)
                * Matrix4::rotation_x(0.0 * std::f64::consts::PI / 2.0)
                * Matrix4::rotation_z(0.0 * std::f64::consts::PI / 4.0),
            material.clone(),
        );

        let _g_id = hexagon(
            &mut w,
            g_id_top,
            Matrix4::translation(1.0, 2.0, 3.0)
                * Matrix4::scaling(0.75, 0.75, 0.75)
                * Matrix4::rotation_x(0.0 * std::f64::consts::PI / 2.0)
                * Matrix4::rotation_z(0.0 * std::f64::consts::PI / 4.0),
            material.clone(),
        );

        let n = 5;
        for ix in 0..n {
            material.color = Tuple::color(ix as f64 / n as f64, 0.5, 0.7);
            let mut g_loop = Group::new();
            g_loop.set_transform(
                Matrix4::translation(-5.0 + 2.0 * ix as f64, 0.0, 0.0)
                    * Matrix4::scaling(0.4, 0.4, 0.4)
                    * Matrix4::rotation_x(-(ix as f64 / n as f64) * std::f64::consts::PI / 4.0),
            );
            let g_id_loop = w.add_shape(Box::new(g_loop));

            for iy in 0..n {
                let scale = 0.5 + iy as f64 / n as f64;
                let _g_id = hexagon(
                    &mut w,
                    g_id_loop,
                    Matrix4::translation(1.0, 1.0 + 0.5 * iy as f64, 3.0)
                        * Matrix4::scaling(scale, scale, scale)
                        * Matrix4::rotation_x(0.0 * std::f64::consts::PI / 2.0)
                        * Matrix4::rotation_z(0.0 * std::f64::consts::PI / 4.0),
                    material.clone(),
                );
            }
        }

        // View transform / camera
        let from = Tuple::point(0.0, 2.5, -7.0);
        let to = Tuple::point(0.0, 1.0, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let transform = view_transform(from, to, up);
        let (display_x, display_y) = (60, 40);
        // let (display_x, display_y) = (3456, 2234);
        let camera =
            Camera::new(display_x, display_y, std::f64::consts::PI / 3.0).with_transform(transform);

        w.divide(g_id_top, 4);
        w.build_bounds();
        let image = w.render(camera);
        let rc = image.write_ppm("test_chap_14_putting_it_all_together.ppm");
        if rc.is_ok() {
            Ok(())
        } else {
            Err("Create a hexagon".into())
        }
    }

    /// Chap 14 - Bounding box test 1 - add_point.
    #[test]
    fn test_chap_14_11() -> Result<(), String> {
        let mut bb = BoundingBox::empty();
        bb.add_point(Tuple::point(5.0, -2.0, 0.0));
        bb.add_point(Tuple::point(7.0, 0.0, -3.0));

        let chk = bb.contains_point(Tuple::point(5.0, -2.0, -3.0));
        let chk = chk && bb.contains_point(Tuple::point(7.0, 0.0, 0.0));
        let chk =
            chk && bb.min == Tuple::point(5.0, -2.0, -3.0) && bb.max == Tuple::point(7.0, 0.0, 0.0);

        if chk {
            Ok(())
        } else {
            Err("Bounding box test 1".into())
        }
    }

    /// Chap 14 - Bounding box test 2 - check intersects.
    #[test]
    fn test_chap_14_12() -> Result<(), String> {
        let bb = BoundingBox::new(Tuple::point(-1.0, -1.0, -1.0), Tuple::point(1.0, 1.0, 1.0));
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));

        let xs = bb.intersects(&r);
        let chk = xs;

        let r = Ray::new(Tuple::point(5.0, 0.5, 0.0), Tuple::vector(-1.0, 0.0, 0.0));
        let xs = bb.intersects(&r);
        let chk = chk && xs;

        let r = Ray::new(Tuple::point(2.0, 0.0, 2.0), Tuple::vector(0.0, 0.0, -1.0));
        let xs = bb.intersects(&r);
        let chk = chk && !xs;

        if chk {
            Ok(())
        } else {
            Err("Bounding box test 2".into())
        }
    }

    /// Chap 14 - Bounding box test 3 - check transform.
    #[test]
    fn test_chap_14_13() -> Result<(), String> {
        let bb = BoundingBox::new(Tuple::point(-1.0, -1.0, -1.0), Tuple::point(1.0, 1.0, 1.0));
        let bb = bb.transform(
            Matrix4::rotation_x(std::f64::consts::PI / 4.0)
                * Matrix4::rotation_y(std::f64::consts::PI / 4.0),
        );

        let one_plus_sqrt2 = 1.0 + std::f64::consts::SQRT_2 / 2.0;
        let min = Tuple::point(-std::f64::consts::SQRT_2, -one_plus_sqrt2, -one_plus_sqrt2);
        let max = Tuple::point(std::f64::consts::SQRT_2, one_plus_sqrt2, one_plus_sqrt2);
        let chk = bb.min.approx_eq(min) && bb.max.approx_eq(max);
        if chk {
            Ok(())
        } else {
            loge!("test_chap_14_13", "bb.min:{}. bb.max:{}", bb.min, bb.max);
            Err("Bounding box test 3 - check transform.".into())
        }
    }

    /// Chap 14 - Bounding box test 4 - check if box contains another box.
    #[test]
    fn test_chap_14_14() -> Result<(), String> {
        let bb = BoundingBox::new(Tuple::point(-1.0, -1.0, -1.0), Tuple::point(1.0, 1.0, 1.0));
        let chk = bb.contains_box(&BoundingBox::new(
            Tuple::point(-1.0, -1.0, -1.0),
            Tuple::point(1.0, 1.0, 1.0),
        ));

        let bb_translated = bb.transform(Matrix4::translation(2.0, 0.0, 0.0));
        let chk = chk && !bb_translated.contains_box(&bb);

        if chk {
            Ok(())
        } else {
            Err("Bounding box test 4 - check if box contains another box.".into())
        }
    }

    /// Chap 14 - Bounding box for a bounded cylinder.
    #[test]
    fn test_chap_14_15() -> Result<(), String> {
        let mut cyl = Cylinder::new();
        cyl.minimum = -5.0;
        cyl.maximum = 3.0;

        let bb = cyl.bounds();
        let chk = bb.min.approx_eq(Tuple::point(-1.0, -5.0, -1.0))
            && bb.max.approx_eq(Tuple::point(1.0, 3.0, 1.0));

        if chk {
            Ok(())
        } else {
            loge!("test_chap_14_15", "bb.min:{} bb.max:{}", bb.min, bb.max);
            Err("Bounding box for a bounded cylinder".into())
        }
    }

    /// Chap 14 - Bounding box for a bounded cone.
    #[test]
    fn test_chap_14_16() -> Result<(), String> {
        // symmetric truncation: radius grows to 5 at both extremes
        let mut cone = Cone::new();
        cone.minimum = -5.0;
        cone.maximum = 3.0;

        let bb = cone.bounds();
        let chk = bb.min.approx_eq(Tuple::point(-5.0, -5.0, -5.0))
            && bb.max.approx_eq(Tuple::point(5.0, 3.0, 5.0));

        // asymmetric case: limit = max(|-1|, |0|) = 1, exercising abs().max()
        let mut cone2 = Cone::new();
        cone2.minimum = -1.0;
        cone2.maximum = 0.0;

        let bb2 = cone2.bounds();
        let chk = chk
            && bb2.min.approx_eq(Tuple::point(-1.0, -1.0, -1.0))
            && bb2.max.approx_eq(Tuple::point(1.0, 0.0, 1.0));

        if chk {
            Ok(())
        } else {
            loge!("test_chap_14_16", "bb.min:{} bb.max:{}", bb.min, bb.max);
            Err("Bounding box for a bounded cone".into())
        }
    }

    /// Chap 14 - Bounding box split.
    #[test]
    fn test_chap_14_17() -> Result<(), String> {
        // dx is chosen over dy since x takes priority
        let bb = BoundingBox::new(Tuple::point(-1.0, -4.0, -5.0), Tuple::point(9.0, 6.0, 5.0));
        let (left, right) = bb.split();

        let chk = left.min.approx_eq(Tuple::point(-1.0, -4.0, -5.0));
        let chk = chk && left.max.approx_eq(Tuple::point(4.0, 6.0, 5.0));
        let chk = chk && right.min.approx_eq(Tuple::point(4.0, -4.0, -5.0));
        let chk = chk && right.max.approx_eq(Tuple::point(9.0, 6.0, 5.0));

        // dx is biggest - split on x
        let bb = BoundingBox::new(Tuple::point(-1.0, -2.0, -3.0), Tuple::point(9.0, 5.5, 3.0));
        let (left, right) = bb.split();

        let chk = chk && left.min.approx_eq(Tuple::point(-1.0, -2.0, -3.0));
        let chk = chk && left.max.approx_eq(Tuple::point(4.0, 5.5, 3.0));
        let chk = chk && right.min.approx_eq(Tuple::point(4.0, -2.0, -3.0));
        let chk = chk && right.max.approx_eq(Tuple::point(9.0, 5.5, 3.0));

        // dy is biggest - split on y
        let bb = BoundingBox::new(Tuple::point(-1.0, -2.0, -3.0), Tuple::point(5.0, 8.0, 3.0));
        let (left, right) = bb.split();

        let chk = chk && left.min.approx_eq(Tuple::point(-1.0, -2.0, -3.0));
        let chk = chk && left.max.approx_eq(Tuple::point(5.0, 3.0, 3.0));
        let chk = chk && right.min.approx_eq(Tuple::point(-1.0, 3.0, -3.0));
        let chk = chk && right.max.approx_eq(Tuple::point(5.0, 8.0, 3.0));

        // dz is biggest - split on z
        let bb = BoundingBox::new(Tuple::point(-1.0, -2.0, -3.0), Tuple::point(5.0, 3.0, 7.0));
        let (left, right) = bb.split();

        let chk = chk && left.min.approx_eq(Tuple::point(-1.0, -2.0, -3.0));
        let chk = chk && left.max.approx_eq(Tuple::point(5.0, 3.0, 2.0));
        let chk = chk && right.min.approx_eq(Tuple::point(-1.0, -2.0, 2.0));
        let chk = chk && right.max.approx_eq(Tuple::point(5.0, 3.0, 7.0));

        if chk {
            Ok(())
        } else {
            loge!(
                "test_chap_14_17",
                "left min:{} max:{}. right min:{} max:{}.",
                left.min,
                left.max,
                right.min,
                right.max
            );
            Err("Bounding box split.".into())
        }
    }

    /// Chap 14 - Split into subtree's.
    #[test]
    fn test_chap_14_18() -> Result<(), String> {
        let mut w = World::new();
        let g = Group::new();
        let g_id = w.add_shape(Box::new(g));

        let s0 = Sphere::new();
        let mut s1 = Sphere::new();
        let mut s2 = Sphere::new();
        s1.set_transform(Matrix4::translation(-2.0, 0.0, 0.0));
        s2.set_transform(Matrix4::translation(2.0, 0.0, 0.0));

        let s0_id = w.add_child(g_id, Box::new(s0));
        let s1_id = w.add_child(g_id, Box::new(s1));
        let s2_id = w.add_child(g_id, Box::new(s2));

        w.build_bounds();
        let (left, right) = w.partition_children(g_id);

        // the two buckets
        let chk = left == vec![s1_id] && right == vec![s2_id];

        // the group kept ONLY the straddler s0
        let group_children = w.shapes[g_id - 1].children().unwrap(); // id == index+1
        let chk = chk && group_children == [s0_id];

        // arena is UNCHANGED — nothing was removed (this is the mental-model fix)
        let chk = chk && w.shapes.len() == 4;

        if chk {
            Ok(())
        } else {
            Err("Split into subtree's.".into())
        }
    }

    /// Chap 14 - Test the subgroups after splitting.
    #[test]
    fn test_chap_14_19() -> Result<(), String> {
        let mut w = World::new();
        let g = Group::new();
        let g_id = w.add_shape(Box::new(g));

        let s0 = Sphere::new();
        let mut s1 = Sphere::new();
        let mut s2 = Sphere::new();
        s1.set_transform(Matrix4::translation(-2.0, 0.0, 0.0));
        s2.set_transform(Matrix4::translation(2.0, 0.0, 0.0));

        let _s0_id = w.add_child(g_id, Box::new(s0));
        let s1_id = w.add_child(g_id, Box::new(s1));
        let _s2_id = w.add_child(g_id, Box::new(s2));

        let (left, right) = w.partition_children(g_id);
        let sub_left = w.make_subgroup(g_id, left); // new sub-group id, holds s1
        let sub_right = w.make_subgroup(g_id, right); // new sub-group id, holds s2

        w.build_bounds();

        // arena grew by 2 (the two sub-groups appended)
        let chk = w.shapes.len() == 6;

        // group now has: straddler s0, plus the two sub-groups
        let gc = w.shapes[g_id - 1].children().unwrap();
        let chk = chk && gc.len() == 3; // [s0_id, sub_left, sub_right]

        // the LEFT sub-group holds exactly s1 (look up the sub-group by its id-1)
        let sub_left_children = w.shapes[sub_left - 1].children().unwrap();
        let chk = chk && sub_left_children == [s1_id];

        // s1's parent was rewired from the group to the sub-group
        let chk = chk && w.shapes[s1_id - 1].data().parent == Some(sub_left);
        // ...and s1 itself is still a leaf (it did NOT become a group)
        let chk = chk && w.shapes[s1_id - 1].children().is_none();

        if chk {
            Ok(())
        } else {
            loge!(
                "test_chap_14_19",
                "gc:{:?} sub_left:{} sub_right:{}",
                gc,
                sub_left,
                sub_right
            );
            Err("Test the subgroups after splitting.".into())
        }
    }

    /// Chap 14 - divide() reorganizes the tree WITHOUT changing what a ray hits.
    #[test]
    fn test_chap_14_20() -> Result<(), String> {
        let mut w = World::new();
        let g_id = w.add_shape(Box::new(Group::new()));

        // three spheres in a row on x; a ray along x hits all three (6 hits)
        let _s0 = w.add_child(g_id, Box::new(Sphere::new()));
        let mut s1 = Sphere::new();
        s1.set_transform(Matrix4::translation(-2.0, 0.0, 0.0));
        let _s1 = w.add_child(g_id, Box::new(s1));
        let mut s2 = Sphere::new();
        s2.set_transform(Matrix4::translation(2.0, 0.0, 0.0));
        let _s2 = w.add_child(g_id, Box::new(s2));

        let ray = Ray::new(Tuple::point(-5.0, 0.0, 0.0), Tuple::vector(1.0, 0.0, 0.0));

        // hits BEFORE dividing
        w.build_bounds();
        let before = w.intersect(&ray);

        // reorganize into a BVH, rebuild bounds, hits AFTER
        w.divide(g_id, 1);
        w.build_bounds();
        let after = w.intersect(&ray);

        // meaningful: the ray really does hit all three spheres
        let mut chk = before.count() == 6;
        // contract: divide preserves the hit set exactly (same t's and object ids)
        chk = chk && before.count() == after.count();
        for i in 0..before.count() {
            chk = chk && approx_eq(before[i].t, after[i].t);
            chk = chk && before[i].object_id == after[i].object_id;
        }

        if chk {
            Ok(())
        } else {
            loge!(
                "test_chap_14_20",
                "before.count:{} after.count:{}",
                before.count(),
                after.count()
            );
            Err("divide() must not change what a ray hits".into())
        }
    }

    /// Chap 15 - Converting an OBJ file to a group
    #[test]
    fn test_chap_15_15() -> std::io::Result<()> {
        let mut w = World::new();

        let light = Light::point_light(
            Tuple::point(-10.0, 10.0, -10.0),
            Tuple::color(1.0, 1.0, 1.0),
        );
        w.light = Some(light);

        // let contents = "v -1 1 0\nv -1 0 0\nv 1 0 0\nv 1 1 0\n \ng FirstGroup\nf 1 2 3\ng SecondGroup\nf 1 3 4";
        // Tetrahedron
        let _contents = "v 0 2 0\n\
                v -1 0 -1\n\
                v 1 0 -1\n\
                v 0 0 1\n\
                g Tetrahedron\n\
                f 1 2 3\n\
                f 1 3 4\n\
                f 1 4 2\n\
                f 2 4 3";

        // let path = std::env::temp_dir().join("test_chap_15_15.obj");
        // std::fs::write(&path, contents)?; // ? → I/O errors become the Err
        // let top_id = w.load_obj_file(&path, Matrix4::identity())?;
        // assert!(top_id > 0);
        // std::fs::remove_file(&path)?; // cleaned up BEFORE the asserts

        // let path = concat!(env!("CARGO_MANIFEST_DIR"), "/models/pumpkin.obj");
        // let s = 0.05; // 79-unit span × 0.05 ≈ 4 units tall
        // let transform = Matrix4::scaling(s, s, s) * Matrix4::translation(2.62, -0.87, 110.02); // move center → origin

        // let path = concat!(env!("CARGO_MANIFEST_DIR"), "/models/cow.obj");
        // let path = concat!(env!("CARGO_MANIFEST_DIR"), "/models/teapot.obj");
        // let s = 0.55; // 79-unit span × 0.05 ≈ 4 units tall

        // icosahedron.obj: circumradius 1, centered on origin; generated by
        // models/generate.py, so this test works from a fresh clone
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/models/icosahedron.obj");
        let s = 1.8;
        let transform = Matrix4::translation(0.0, 1.0, 0.0) * Matrix4::scaling(s, s, s);
        let top_id = w.load_obj_file(path, transform)?;
        assert!(top_id > 0);

        // View transform / camera
        let from = Tuple::point(4.0, 2.5, -7.0);
        let to = Tuple::point(0.0, 1.0, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let transform = view_transform(from, to, up);
        let (display_x, display_y) = (600, 400);
        // let (display_x, display_y) = (3456, 2234);
        let camera =
            Camera::new(display_x, display_y, std::f64::consts::PI / 3.0).with_transform(transform);

        w.divide(top_id, 4);
        w.build_bounds();
        let image = w.render(camera);
        let _rc = image.write_ppm("test_chap_15_putting_it_all_together.ppm");

        Ok(())
    }

    /// Chap 15 - Smooth shading: a torus with vertex normals. The mesh in
    /// models/torus_smooth.obj (generated by models/generate.py) carries `vn`
    /// records with analytic normals, so the parser emits smooth triangles and
    /// the silhouette of each facet disappears from the shading.
    #[test]
    fn test_chap_15_smooth_torus() -> std::io::Result<()> {
        let mut w = World::new();

        let light = Light::point_light(
            Tuple::point(-10.0, 10.0, -10.0),
            Tuple::color(1.0, 1.0, 1.0),
        );
        w.light = Some(light);

        // Floor - a slightly reflective plane
        let mut floor = Plane::new();
        let mut floor_material = Material::new();
        floor_material.color = Tuple::color(0.8, 0.8, 0.85);
        floor_material.specular = 0.0;
        floor_material.reflective = 0.2;
        floor.set_material(floor_material);
        w.add_shape(Box::new(floor));

        // Torus: xz plane, outer radius 1.4. Stand it up like a wheel
        // (rotation_x) resting on the floor, angled toward the camera.
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/models/torus_smooth.obj");
        let transform = Matrix4::translation(0.0, 1.4, 0.0)
            * Matrix4::rotation_y(std::f64::consts::PI / 6.0)
            * Matrix4::rotation_x(std::f64::consts::PI / 2.0);
        let top_id = w.load_obj_file(path, transform)?;
        assert!(top_id > 0);

        let mut torus_material = Material::new();
        torus_material.color = Tuple::color(0.2, 0.55, 0.75);
        torus_material.diffuse = 0.8;
        torus_material.specular = 0.4;
        torus_material.shininess = 40.0;
        w.set_material_recursive(top_id, &torus_material);

        // View transform / camera
        let from = Tuple::point(0.0, 2.2, -4.5);
        let to = Tuple::point(0.0, 0.8, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let transform = view_transform(from, to, up);
        let (display_x, display_y) = (600, 400);
        let camera =
            Camera::new(display_x, display_y, std::f64::consts::PI / 3.0).with_transform(transform);

        w.divide(top_id, 4);
        w.build_bounds();
        let image = w.render_parallel(camera);
        let _rc = image.write_ppm("test_chap_15_smooth_torus.ppm");

        Ok(())
    }

    /// Chap 15 - Constructing a smooth triangle
    #[test]
    fn test_chap_15_16() -> Result<(), String> {
        let p1 = Tuple::point(0.0, 1.0, 0.0);
        let p2 = Tuple::point(-1.0, 0.0, 0.0);
        let p3 = Tuple::point(1.0, 0.0, 0.0);
        let n1 = Tuple::vector(0.0, 1.0, 0.0);
        let n2 = Tuple::vector(-1.0, 0.0, 0.0);
        let n3 = Tuple::vector(1.0, 0.0, 0.0);
        let tri = TriangleUV::new(p1, p2, p3, n1, n2, n3);

        let chk = tri.p1.approx_eq(p1);
        let chk = chk && tri.p2.approx_eq(p2);
        let chk = chk && tri.p3.approx_eq(p3);
        let chk = chk && tri.n1.approx_eq(n1);
        let chk = chk && tri.n2.approx_eq(n2);
        let chk = chk && tri.n3.approx_eq(n3);
        if chk {
            Ok(())
        } else {
            Err("Constructing a smooth triangle".into())
        }
    }

    /// Chap 15 - An intersection can encapsulate `u` and `v`
    #[test]
    fn test_chap_15_17() -> Result<(), String> {
        let p1 = Tuple::point(0.0, 1.0, 0.0);
        let p2 = Tuple::point(-1.0, 0.0, 0.0);
        let p3 = Tuple::point(1.0, 0.0, 0.0);
        let s = Triangle::new(p1, p2, p3);

        let i = Intersection::new_with_uv(3.5, s.id(), 0.2, 0.4);

        let chk = approx_eq(i.u, 0.2);
        let chk = chk && approx_eq(i.v, 0.4);
        if chk {
            Ok(())
        } else {
            Err("An intersection can encapsulate `u` and `v`".into())
        }
    }

    /// Chap 15 - An intersection with a smooth triangle stores u/v
    #[test]
    fn test_chap_15_18() -> Result<(), String> {
        let origin = Tuple::point(-0.2, 0.3, -2.0);
        let direction = Tuple::vector(0.0, 0.0, 1.0);
        let r = Ray::new(origin, direction);

        let p1 = Tuple::point(0.0, 1.0, 0.0);
        let p2 = Tuple::point(-1.0, 0.0, 0.0);
        let p3 = Tuple::point(1.0, 0.0, 0.0);
        let n1 = Tuple::vector(0.0, 1.0, 0.0);
        let n2 = Tuple::vector(-1.0, 0.0, 0.0);
        let n3 = Tuple::vector(1.0, 0.0, 0.0);

        let tri = TriangleUV::new(p1, p2, p3, n1, n2, n3);
        let xs = tri.local_intersect(&r);

        let chk = tri.n1.approx_eq(n1);
        let chk = chk && xs.count() > 0;
        let chk = if chk && xs.count() > 0 {
            xs[0].t > 0.0
        } else {
            false
        };

        let chk = if chk && xs.count() > 0 {
            xs[0].u > 0.0 && approx_eq(xs[0].u, 0.45) && approx_eq(xs[0].v, 0.25)
        } else {
            false
        };

        if chk {
            Ok(())
        } else {
            Err("An intersection with a smooth triangle stores u/v".into())
        }
    }

    /// Chap 15 - A smooth triangle uses u/v to interpolate the normal
    #[test]
    fn test_chap_15_19() -> Result<(), String> {
        let i = Intersection::new_with_uv(1.0, 1, 0.45, 0.25);

        let p1 = Tuple::point(0.0, 1.0, 0.0);
        let p2 = Tuple::point(-1.0, 0.0, 0.0);
        let p3 = Tuple::point(1.0, 0.0, 0.0);
        let n1 = Tuple::vector(0.0, 1.0, 0.0);
        let n2 = Tuple::vector(-1.0, 0.0, 0.0);
        let n3 = Tuple::vector(1.0, 0.0, 0.0);

        let tri = TriangleUV::new(p1, p2, p3, n1, n2, n3);
        let n = tri.normal_at(Tuple::point(0.0, 0.0, 0.0), i);

        let expect = Tuple::vector(-0.554700196225, 0.832050294338, 0.000000000000);
        let chk = n.approx_eq(expect);

        if chk {
            Ok(())
        } else {
            loge!("test_chap_15_19", "n: {}, i: {:?}", n, i);
            Err("A smooth triangle uses u/v to interpolate the normal".into())
        }
    }

    /// Chap 15 - Preparing the normal on a smooth triangle
    #[test]
    fn test_chap_15_20() -> Result<(), String> {
        let mut w = World::new();

        let p1 = Tuple::point(0.0, 1.0, 0.0);
        let p2 = Tuple::point(-1.0, 0.0, 0.0);
        let p3 = Tuple::point(1.0, 0.0, 0.0);
        let n1 = Tuple::vector(0.0, 1.0, 0.0);
        let n2 = Tuple::vector(-1.0, 0.0, 0.0);
        let n3 = Tuple::vector(1.0, 0.0, 0.0);

        let tri = TriangleUV::new(p1, p2, p3, n1, n2, n3);

        let tri_id = w.add_shape(Box::new(tri));

        let i = Intersection::new_with_uv(1.0, tri_id, 0.45, 0.25);

        let origin = Tuple::point(-0.2, 0.3, -2.0);
        let direction = Tuple::vector(0.0, 0.0, 1.0);
        let r = Ray::new(origin, direction);

        let mut xs = Intersections::new();
        xs.push(i);
        let comps = prepare_computations(i, &r, &w.shapes, &xs);

        let chk = comps.normalv.approx_eq(Tuple::vector(
            -0.554700196225,
            0.832050294338,
            0.000000000000,
        ));

        if chk {
            Ok(())
        } else {
            loge!("test_chap_15_20", "comps.normalv: {}", comps.normalv);
            Err("Preparing the normal on a smooth triangle".into())
        }
    }

    /// Chap 16 - CSG is created with an operation and two shapes
    #[test]
    fn test_chap_16_3() -> Result<(), String> {
        let mut w = World::new();

        let c_id = w.add_shape(Box::new(Csg::new(CsgOperation::Union)));
        let s1_id = w.add_child(c_id, Box::new(Sphere::new()));
        let s2_id = w.add_child(c_id, Box::new(Cube::new()));

        let c = w.shapes[c_id - 1].as_ref();
        let chk = c.csg_operation() == Some(CsgOperation::Union)
            && c.children() == Some(&[s1_id, s2_id][..])
            && w.shapes[s1_id - 1].data().parent == Some(c_id)
            && w.shapes[s2_id - 1].data().parent == Some(c_id);
        if chk {
            Ok(())
        } else {
            loge!("test_chap_16_3", "children:{:?}", c.children());
            Err("CSG is created with an operation and two shapes".into())
        }
    }

    /// Chap 16 - Filtering a list of intersections
    #[test]
    fn test_chap_16_4() -> Result<(), String> {
        // (operation, expected xs indices after filtering)
        let cases = [
            (CsgOperation::Union, [0_usize, 3_usize]),
            (CsgOperation::Intersection, [1, 2]),
            (CsgOperation::Difference, [0, 1]),
        ];
        for (op, expected) in cases {
            let mut w = World::new();
            let c_id = w.add_shape(Box::new(Csg::new(op)));
            let s1_id = w.add_child(c_id, Box::new(Sphere::new()));
            let s2_id = w.add_child(c_id, Box::new(Cube::new()));

            let mut xs = Intersections::new();
            xs.push(Intersection::new(1.0, s1_id));
            xs.push(Intersection::new(2.0, s2_id));
            xs.push(Intersection::new(3.0, s1_id));
            xs.push(Intersection::new(4.0, s2_id));

            let result = w.filter_intersections(w.shapes[c_id - 1].as_ref(), &xs);

            let chk = result.count() == 2
                && approx_eq(result[0].t, xs[expected[0]].t)
                && result[0].object_id == xs[expected[0]].object_id
                && approx_eq(result[1].t, xs[expected[1]].t)
                && result[1].object_id == xs[expected[1]].object_id;
            if !chk {
                loge!("test_chap_16_4", "op:{:?} result:{:?}", op, result);
                return Err("Filtering a list of intersections".into());
            }
        }
        Ok(())
    }

    /// Chap 16 - A ray misses a CSG object
    #[test]
    fn test_chap_16_5() -> Result<(), String> {
        let mut w = World::new();
        let c_id = w.add_shape(Box::new(Csg::new(CsgOperation::Union)));
        let _s1_id = w.add_child(c_id, Box::new(Sphere::new()));
        let _s2_id = w.add_child(c_id, Box::new(Cube::new()));
        w.build_bounds();

        let r = Ray::new(Tuple::point(0.0, 2.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs = w.intersect(&r);

        if xs.is_empty() {
            Ok(())
        } else {
            loge!("test_chap_16_5", "xs.count:{}", xs.count());
            Err("A ray misses a CSG object".into())
        }
    }

    /// Chap 16 - A ray hits a CSG object
    #[test]
    fn test_chap_16_6() -> Result<(), String> {
        let mut w = World::new();
        let c_id = w.add_shape(Box::new(Csg::new(CsgOperation::Union)));
        let s1_id = w.add_child(c_id, Box::new(Sphere::new()));
        let mut s2 = Sphere::new();
        s2.set_transform(Matrix4::translation(0.0, 0.0, 0.5));
        let s2_id = w.add_child(c_id, Box::new(s2));
        w.build_bounds();

        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs = w.intersect(&r);

        let chk = xs.count() == 2
            && approx_eq(xs[0].t, 4.0)
            && xs[0].object_id == s1_id
            && approx_eq(xs[1].t, 6.5)
            && xs[1].object_id == s2_id;
        if chk {
            Ok(())
        } else {
            loge!("test_chap_16_6", "xs:{:?}", xs);
            Err("A ray hits a CSG object".into())
        }
    }

    /// Chap 16 - divide() must never repartition a CSG node's children
    /// (children[0]/children[1] ARE the left/right operands)
    #[test]
    fn test_chap_16_divide() -> Result<(), String> {
        let mut w = World::new();
        let c_id = w.add_shape(Box::new(Csg::new(CsgOperation::Difference)));
        let s1_id = w.add_child(c_id, Box::new(Sphere::new()));
        let s2_id = w.add_child(c_id, Box::new(Cube::new()));

        w.divide(c_id, 1); // threshold 1 would split any plain group

        let c = w.shapes[c_id - 1].as_ref();
        if c.children() == Some(&[s1_id, s2_id][..]) {
            Ok(())
        } else {
            loge!("test_chap_16_divide", "children:{:?}", c.children());
            Err("divide() must leave CSG children untouched".into())
        }
    }

    /// Chap 16 - Putting it together: the classic CSG widget, a rounded cube
    /// (cube ∩ sphere) with three cylindrical holes drilled through it
    /// (− union of three cylinders). Exercises nested CSG nodes.
    #[test]
    fn test_chap_16_putting_it_all_together() -> Result<(), String> {
        let mut w = World::new();

        let light = Light::point_light(
            Tuple::point(-10.0, 10.0, -10.0),
            Tuple::color(1.0, 1.0, 1.0),
        );
        w.light = Some(light);

        // Floor - a slightly reflective plane
        let mut floor = Plane::new();
        let mut floor_material = Material::new();
        floor_material.color = Tuple::color(0.8, 0.8, 0.85);
        floor_material.specular = 0.0;
        floor_material.reflective = 0.2;
        floor.set_material(floor_material);
        w.add_shape(Box::new(floor));

        // top = core − holes
        let mut top = Csg::new(CsgOperation::Difference);
        top.set_transform(
            Matrix4::translation(0.0, 1.5, 0.0)
                * Matrix4::rotation_y(std::f64::consts::PI / 4.0)
                * Matrix4::rotation_x(-0.4),
        );
        let top_id = w.add_shape(Box::new(top));

        // core = cube ∩ sphere (a rounded cube), the LEFT child
        let core_id = w.add_child(top_id, Box::new(Csg::new(CsgOperation::Intersection)));
        // holes = cyl_x ∪ (cyl_y ∪ cyl_z), the RIGHT child
        let holes_id = w.add_child(top_id, Box::new(Csg::new(CsgOperation::Union)));

        let _cube_id = w.add_child(core_id, Box::new(Cube::new()));
        let mut ball = Sphere::new();
        ball.set_transform(Matrix4::scaling(1.35, 1.35, 1.35));
        let _ball_id = w.add_child(core_id, Box::new(ball));

        let drill = |axis_rot: Matrix4| {
            let mut c = Cylinder::new();
            c.minimum = -2.0;
            c.maximum = 2.0;
            c.closed = true;
            c.set_transform(axis_rot * Matrix4::scaling(0.72, 1.0, 0.72));
            c
        };
        let cyl_x = drill(Matrix4::rotation_z(std::f64::consts::PI / 2.0));
        let _ = w.add_child(holes_id, Box::new(cyl_x));
        let yz_id = w.add_child(holes_id, Box::new(Csg::new(CsgOperation::Union)));
        let cyl_y = drill(Matrix4::identity());
        let _ = w.add_child(yz_id, Box::new(cyl_y));
        let cyl_z = drill(Matrix4::rotation_x(std::f64::consts::PI / 2.0));
        let _ = w.add_child(yz_id, Box::new(cyl_z));

        // materials: blue body, red hole walls
        let mut body = Material::new();
        body.color = Tuple::color(0.2, 0.55, 0.75);
        body.diffuse = 0.8;
        body.specular = 0.4;
        body.shininess = 40.0;
        w.set_material_recursive(core_id, &body);
        let mut hole = Material::new();
        hole.color = Tuple::color(0.85, 0.25, 0.2);
        hole.diffuse = 0.8;
        hole.specular = 0.3;
        w.set_material_recursive(holes_id, &hole);

        // View transform / camera
        let from = Tuple::point(0.0, 3.0, -5.5);
        let to = Tuple::point(0.0, 1.3, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let transform = view_transform(from, to, up);
        let camera = Camera::new(600, 400, std::f64::consts::PI / 3.0).with_transform(transform);

        w.build_bounds();
        let image = w.render_parallel(camera);
        let rc = image.write_ppm("test_chap_16_putting_it_all_together.ppm");
        if rc.is_err() {
            return Err("failed to write PPM".into());
        }

        Ok(())
    }

    /// Chap 16 - The chapter-opener figure: the same cube/sphere pair
    /// combined three ways — union, intersection, difference — side by side.
    #[test]
    fn test_chap_16_union_intersect_difference() -> Result<(), String> {
        let mut w = World::new();

        let light = Light::point_light(
            Tuple::point(-10.0, 10.0, -10.0),
            Tuple::color(1.0, 1.0, 1.0),
        );
        w.light = Some(light);

        // Checkered, slightly reflective floor
        // large squares: the tracer has no anti-aliasing, so small squares
        // moiré badly toward the horizon
        let floor_checkers = CheckersPattern::new(WHITE, BLACK);
        let mut floor = Plane::new();
        let mut floor_material = Material::new();
        floor_material.color = Tuple::color(0.8, 0.8, 0.85);
        floor_material.specular = 0.0;
        floor_material.reflective = 0.2;
        floor_material.pattern = Some(Box::new(floor_checkers));
        floor.set_material(floor_material);
        w.add_shape(Box::new(floor));

        let mut cube_material = Material::new();
        cube_material.color = Tuple::color(0.2, 0.55, 0.75);
        cube_material.diffuse = 0.8;
        cube_material.specular = 0.4;
        cube_material.shininess = 40.0;
        let mut ball_material = Material::new();
        ball_material.color = Tuple::color(0.85, 0.3, 0.25);
        ball_material.diffuse = 0.8;
        ball_material.specular = 0.4;
        ball_material.shininess = 40.0;

        // the same cube/sphere pair, one operation per column
        let ops = [
            (CsgOperation::Union, -3.2),
            (CsgOperation::Intersection, 0.0),
            (CsgOperation::Difference, 3.2),
        ];
        for (op, x) in ops {
            let mut top = Csg::new(op);
            top.set_transform(
                Matrix4::translation(x, 1.0, 0.0)
                    * Matrix4::rotation_y(-std::f64::consts::PI / 6.0),
            );
            let top_id = w.add_shape(Box::new(top));

            let mut cube = Cube::new();
            cube.set_material(cube_material.clone());
            let _ = w.add_child(top_id, Box::new(cube)); // left

            let mut ball = Sphere::new();
            ball.set_material(ball_material.clone());
            // overlap the cube's upper front corner
            ball.set_transform(
                Matrix4::translation(0.5, 0.5, -0.5) * Matrix4::scaling(0.9, 0.9, 0.9),
            );
            let _ = w.add_child(top_id, Box::new(ball)); // right
        }

        // View transform / camera. Tilted down far enough that the horizon
        // (and its checkerboard moiré) stays out of frame.
        let from = Tuple::point(0.0, 4.2, -7.5);
        let to = Tuple::point(0.0, 0.7, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let transform = view_transform(from, to, up);
        let camera = Camera::new(900, 450, std::f64::consts::PI / 3.0).with_transform(transform);

        w.build_bounds();
        let image = w.render_parallel(camera);
        let rc = image.write_ppm("test_chap_16_union_intersect_difference.ppm");
        if rc.is_err() {
            return Err("failed to write PPM".into());
        }

        Ok(())
    }

    /// Cornell box (Goral et al. 1984) — the classic ray tracing test scene:
    /// a 2x2x2 room (red left wall, green right wall, white floor/ceiling/back)
    /// with a tall and a short block. Walls at x = ±1, floor y = 0, ceiling
    /// y = 2, back wall z = 1; the front (z = -1) is open and the camera looks
    /// in from outside.
    fn cornell_box(light: Light) -> World {
        use std::f64::consts::PI;

        let matte = |r: f64, g: f64, b: f64| {
            let mut m = Material::new();
            m.color = Tuple::color(r, g, b);
            m.specular = 0.0;
            m
        };
        let white = matte(0.73, 0.73, 0.73);
        let red = matte(0.65, 0.05, 0.05);
        let green = matte(0.12, 0.45, 0.15);

        let mut w = World::new();
        w.set_light(light);

        let wall = |material: &Material, transform: Matrix4| -> Box<dyn Shape> {
            let mut p = Plane::new();
            p.set_material(material.clone());
            p.set_transform(transform);
            Box::new(p)
        };
        w.add_shape(wall(&white, Matrix4::identity()));
        w.add_shape(wall(&white, Matrix4::translation(0.0, 2.0, 0.0)));
        w.add_shape(wall(
            &white,
            Matrix4::translation(0.0, 0.0, 1.0) * Matrix4::rotation_x(PI / 2.0),
        ));
        w.add_shape(wall(
            &red,
            Matrix4::translation(-1.0, 0.0, 0.0) * Matrix4::rotation_z(PI / 2.0),
        ));
        w.add_shape(wall(
            &green,
            Matrix4::translation(1.0, 0.0, 0.0) * Matrix4::rotation_z(PI / 2.0),
        ));

        // Blocks (unit cube spans -1..1, so scale is half the block size).
        let block = |material: &Material, transform: Matrix4| -> Box<dyn Shape> {
            let mut c = Cube::new();
            c.set_material(material.clone());
            c.set_transform(transform);
            Box::new(c)
        };
        w.add_shape(block(
            &white,
            Matrix4::translation(-0.35, 0.6, 0.35)
                * Matrix4::rotation_y(PI * 17.0 / 180.0)
                * Matrix4::scaling(0.3, 0.6, 0.3),
        ));
        w.add_shape(block(
            &white,
            Matrix4::translation(0.35, 0.3, -0.3)
                * Matrix4::rotation_y(-PI * 17.0 / 180.0)
                * Matrix4::scaling(0.3, 0.3, 0.3),
        ));
        w.build_bounds();
        w
    }

    /// Render a Cornell box at 1000x1000, log wall-clock time, pixels/sec and
    /// (with `--features stats`) BVH counters, and write `<name>.ppm`.
    fn render_cornell_box(name: &str, w: &World) -> Result<(), String> {
        render_cornell_box_with(name, w, 1)
    }

    /// As `render_cornell_box`, with an `aa` x `aa` anti-aliasing grid.
    fn render_cornell_box_with(name: &str, w: &World, aa: usize) -> Result<(), String> {
        use std::f64::consts::PI;
        use std::time::Instant;

        let from = Tuple::point(0.0, 1.0, -3.5);
        let to = Tuple::point(0.0, 1.0, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let (hsize, vsize) = (1000, 1000);
        let camera = Camera::new(hsize, vsize, PI * 39.0 / 180.0)
            .with_transform(view_transform(from, to, up))
            .with_antialias(aa);

        reset_stats();
        let start = Instant::now();
        let image = w.render_parallel(camera);
        let elapsed = start.elapsed();
        let (node_visits, prim_tests) = read_stats();

        let pixels = (hsize * vsize) as f64;
        logi!(
            name,
            "{}x{} rendered in {:.3} s = {:.0} pixels/s (node visits: {}, prim tests: {})",
            hsize,
            vsize,
            elapsed.as_secs_f64(),
            pixels / elapsed.as_secs_f64(),
            node_visits,
            prim_tests
        );

        image
            .write_ppm(format!("{name}.ppm"))
            .map_err(|e| format!("failed to write PPM: {e}"))
    }

    /// Chap 17 - Anti-aliasing only touches edge pixels: flat interior and
    /// background pixels are identical with and without it, and the edge of
    /// the sphere is a blend of sphere and background
    #[test]
    fn test_chap_17_4() -> Result<(), String> {
        // A flat-shaded sphere: every interior pixel is exactly the same
        // colour, so only its silhouette can be an edge.
        let mut w = World::new();
        w.set_light(Light::point_light(Tuple::point(-10.0, 10.0, -10.0), WHITE));
        let mut s = Sphere::new();
        let mut m = Material::new();
        m.color = Tuple::color(0.8, 0.2, 0.2);
        m.ambient = 1.0;
        m.diffuse = 0.0;
        m.specular = 0.0;
        s.set_material(m);
        w.add_shape(Box::new(s));
        let from = Tuple::point(0.0, 0.0, -5.0);
        let to = Tuple::point(0.0, 0.0, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);
        let transform = view_transform(from, to, up);
        let plain = Camera::new(41, 41, std::f64::consts::PI / 3.0).with_transform(transform);
        let aa = plain.with_antialias(4);

        let img = w.render_parallel(plain);
        let img_aa = w.render_parallel(aa);

        // Centre of the sphere and a background corner are flat: unchanged.
        let flat = [(20, 20), (0, 0)];
        for (x, y) in flat {
            if !img[(x, y)].approx_eq(img_aa[(x, y)]) {
                loge!(
                    "test_chap_17_4",
                    "({x},{y}) {} vs {}",
                    img[(x, y)],
                    img_aa[(x, y)]
                );
                return Err("Anti-aliasing changed a flat pixel".into());
            }
        }

        // Walk from the centre to the right until the background; the last
        // sphere pixel is an edge and must differ after anti-aliasing.
        let black = Tuple::color(0.0, 0.0, 0.0);
        let edge_x = (20..41)
            .find(|&x| img[(x, 20)].approx_eq(black))
            .ok_or("no background found on row 20")?
            - 1;
        let (before, after) = (img[(edge_x, 20)], img_aa[(edge_x, 20)]);
        let blended = !before.approx_eq(after) && after.x < before.x && after.x > 0.0;
        if blended {
            Ok(())
        } else {
            loge!("test_chap_17_4", "edge x={edge_x}: {} vs {}", before, after);
            Err("Anti-aliasing did not blend the sphere edge".into())
        }
    }

    /// Cornell box with a point light: the fixed benchmark baseline that
    /// chapter-17 extensions (soft shadows, anti-aliasing, ...) are measured
    /// against. Run deliberately, in release:
    /// `cargo test --release --lib cornell_box -- --ignored --nocapture`
    #[test]
    #[ignore]
    fn test_cornell_box_benchmark() -> Result<(), String> {
        let light = Light::point_light(
            // Below and in front of the ceiling so it and the block fronts get
            // some direct light; a point light flush with the ceiling leaves
            // them at a grazing angle.
            Tuple::point(0.0, 1.8, -0.6),
            Tuple::color(1.0, 1.0, 1.0),
        );
        render_cornell_box("test_cornell_box_benchmark", &cornell_box(light))
    }

    /// Bonus (soft shadows) - Cornell box lit the way the original was: a
    /// square area light in the ceiling, giving penumbrae around the blocks.
    #[test]
    #[ignore]
    fn test_cornell_box_area_light() -> Result<(), String> {
        let mut light = Light::area_light(
            Tuple::point(-0.25, 1.95, -0.25),
            Tuple::vector(0.5, 0.0, 0.0),
            8,
            Tuple::vector(0.0, 0.0, 0.5),
            8,
            Tuple::color(1.0, 1.0, 1.0),
        );
        light.jitter_by = Sequence::new(vec![0.7, 0.3, 0.9, 0.1, 0.5, 0.2, 0.8, 0.4, 0.6]);
        render_cornell_box("test_cornell_box_area_light", &cornell_box(light))
    }

    /// Chap 17 - Cornell box with the area light and 4x4 edge-detected
    /// anti-aliasing, to measure what anti-aliasing costs on top.
    #[test]
    #[ignore]
    fn test_cornell_box_antialias() -> Result<(), String> {
        let mut light = Light::area_light(
            Tuple::point(-0.25, 1.95, -0.25),
            Tuple::vector(0.5, 0.0, 0.0),
            8,
            Tuple::vector(0.0, 0.0, 0.5),
            8,
            Tuple::color(1.0, 1.0, 1.0),
        );
        light.jitter_by = Sequence::new(vec![0.7, 0.3, 0.9, 0.1, 0.5, 0.2, 0.8, 0.4, 0.6]);
        render_cornell_box_with("test_cornell_box_antialias", &cornell_box(light), 4)
    }
}
