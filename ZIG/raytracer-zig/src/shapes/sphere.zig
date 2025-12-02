const std = @import("std");
const print = @import("std").debug.print;
const types = @import("../types.zig");
const canvas = @import("../canvas.zig");
const utils = @import("../utils.zig");
const shapes = @import("../shapes.zig");
const matrix = @import("../matrix.zig");
const mat_module = @import("../material.zig");

const S = types.S;
const Ray = types.Ray;
const Scalar = types.Scalar;
const Shape = shapes.Shape;
const Tuple = types.Tuple;
const Matrix = matrix.Mat4;
const Intersection = types.Intersection;
const Intersections = types.Intersections;
const LocalIntersections = types.LocalIntersections;
const point = types.Point;
const vector = types.Vector;
const approxEq = types.approxEq;
const log = utils.log;
const material = mat_module.Material;

var NEXT_SPHERE_ID: std.atomic.Value(usize) = .{ .raw = 0 };

pub const Sphere = struct {
    xs: Intersections,
    radius: Scalar,
    object_id: usize,
    transformed_m: matrix.Mat4,
    transformed_m_inv: matrix.Mat4,
    transposed_m: matrix.Mat4,
    transposed_m_inv: matrix.Mat4,
    material: mat_module.Material,

    pub fn init() Sphere {
        const obj_id = NEXT_SPHERE_ID.fetchAdd(1, .seq_cst);
        return .{
            .xs = Intersections.init(),
            .radius = S(1),
            .object_id = obj_id,
            .transformed_m = matrix.Mat4.identity(),
            .transformed_m_inv = matrix.Mat4.identity(),
            .transposed_m = matrix.Mat4.identity(),
            .transposed_m_inv = matrix.Mat4.identity(),
            .material = material.init(),
        };
    }

    pub inline fn id(self: *const Sphere) usize {
        return self.object_id;
    }

    pub fn intersection(t: Scalar, sphere: Sphere) !Intersection {
        const i: Intersection = .{ .t = t, .object_id = sphere.object_id };
        return i;
    }

    /// ---
    /// Compute a local ray by applying the inverse of the sphere
    /// transform. Use the local rays origin and direction to
    /// compute the Intersections.
    /// ---
    pub fn intersect(self: *const Sphere, ray: Ray) LocalIntersections {

        // ---
        // NOTE: Two options
        // 1. Compute the inverse of the sphere transform.
        // 2. Use the already computed invers.
        // ---
        // const lr = Ray{ .origin = self.inverse().mulT(ray.origin), .direction = self.inverse().mulT(ray.direction) };
        const lr = Ray{ .origin = self.transformed_m_inv.mulT(ray.origin), .direction = self.transformed_m_inv.mulT(ray.direction) };
        const sphere2ray = lr.origin.sub(point(0, 0, 0));
        const a = lr.direction.dot(lr.direction);
        const b = S(2) * lr.direction.dot(sphere2ray);
        const c = sphere2ray.dot(sphere2ray) - S(1);
        const discriminant = b * b - S(4) * a * c;

        if (discriminant < S(0)) {
            return LocalIntersections.init();
        }

        const t1 = (-b - std.math.sqrt(discriminant)) / (S(2) * a);
        const t2 = (-b + std.math.sqrt(discriminant)) / (S(2) * a);

        var xs = LocalIntersections.init();
        xs.add(Shape.intersection(S(t1), @ptrCast(self)), self.object_id);
        xs.add(Shape.intersection(S(t2), @ptrCast(self)), self.object_id);

        return xs;
    }

    /// ---
    /// NOTE: Compute the normal at the given world_point.
    ///       This function uses cached matrix's for the
    ///       inverted transform and the inverted transpose.
    /// ---
    pub fn normal_at(self: *const Sphere, world_point: Tuple) Tuple {
        const object_point = self.transformed_m_inv.mulT(world_point);
        const object_normal = (object_point.sub(Tuple.point(0, 0, 0)));
        var world_normal = self.transposed_m_inv.mulT(object_normal);
        world_normal.w = S(0);
        return world_normal.normalize();
    }

    /// ---
    /// Setting the transform matrix.
    /// NOTE: Also computes the inverse as a side effect...
    /// ---
    pub fn set_transform(self: *Sphere, m: *const matrix.Mat4) void {
        self.transformed_m = m.*;
        self.transformed_m_inv = (m.*).inverse();
        self.transposed_m = (m.*).transpose();
        self.transposed_m_inv = (m.*).transpose().inverse();
    }

    pub fn transform(self: *const Sphere) matrix.Mat4 {
        return self.transformed_m;
    }

    pub fn inverse(self: *const Sphere) matrix.Mat4 {
        return self.transformed_m_inv;
    }
};

test "Chap5 -A ray intersects a sphere at two points" {
    const r = Ray.init(point(-5, 0, 0), vector(1, 0, 0));
    const s = Sphere.init();
    const xs = s.intersect(r);

    try std.testing.expect(xs.count == 2);
    try std.testing.expect((xs.local_intersections_items[0]).t == S(4));
    try std.testing.expect((xs.local_intersections_items[1]).t == S(6));
}

test "Chap5 -A ray intersects a sphere at a tangent" {
    const r = Ray.init(point(0, 1, -5), vector(0, 0, 1));
    const s = Sphere.init();
    const xs = s.intersect(r);

    try std.testing.expect(xs.count == 2);
    try std.testing.expect((xs.local_intersections_items[0]).t == S(5));
    try std.testing.expect((xs.local_intersections_items[1]).t == S(5));
}

test "Chap5 -A ray misses a sphere" {
    const r = Ray.init(point(0, 2, -5), vector(0, 0, 1));
    const s = Sphere.init();
    const xs = s.intersect(r);

    try std.testing.expect(xs.count == 0);
}

test "Chap5 -A ray originates inside a sphere" {
    const r = Ray.init(point(0, 0, 0), vector(0, 0, 1));
    const s = Sphere.init();
    const xs = s.intersect(r);

    try std.testing.expect(xs.count == 2);
    try std.testing.expect((xs.local_intersections_items[0]).t == S(-1));
    try std.testing.expect((xs.local_intersections_items[1]).t == S(1));
}

test "Chap5 -A sphere is behind a ray" {
    const r = Ray.init(point(0, 0, 5), vector(0, 0, 1));
    const s = Sphere.init();
    const xs = s.intersect(r);

    try std.testing.expect(xs.count == 2);
    try std.testing.expect((xs.local_intersections_items[0]).t == S(-6));
    try std.testing.expect((xs.local_intersections_items[1]).t == S(-4));
}

test "Chap6 -A sphere has a default material" {
    const m = material.init();
    const s = Sphere.init();
    try std.testing.expect(m.equals(s.material));
}

test "Chap6 -A sphere may be assigned a material" {
    var m = material.init();
    m.ambient = S(1);
    var s = Sphere.init();
    s.material = m;
    try std.testing.expect(m.equals(s.material));
    try std.testing.expect(s.material.equals(m));
}
