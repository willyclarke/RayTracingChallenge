const std = @import("std");
const print = @import("std").debug.print;
const types = @import("../types.zig");
// const canvas = @import("../canvas.zig");
const utils = @import("../utils.zig");

const shapes = @import("shapes.zig");
const matrix = @import("../matrix.zig");
const mat_module = @import("../material.zig");
const intersection_mod = @import("intersections.zig");

const ShapeHeader = @import("shape_header.zig").ShapeHeader;

const S = types.S;
const Ray = types.Ray;
const Scalar = types.Scalar;
const Shape = shapes.Shape;
const Tuple = types.Tuple;
const Matrix = matrix.Mat4;
const Intersection = intersection_mod.Intersection;
const Intersections = intersection_mod.Intersections;
const LocalIntersections = intersection_mod.LocalIntersections;
const point = types.Point;
const vector = types.Vector;
const approxEq = types.approxEq;
const log = utils.log;
const Material = mat_module.Material;

var NEXT_SPHERE_ID: std.atomic.Value(usize) = .{ .raw = 1 };

pub const Sphere = struct {
    h: ShapeHeader,
    radius: Scalar,

    pub fn init() Sphere {
        const obj_id = NEXT_SPHERE_ID.fetchAdd(1, .seq_cst);
        // log(@src(), "\nNext Sphere Id:{}\n", .{obj_id});
        return .{
            .h = .{
                .object_id = obj_id,
                .transformed_m = matrix.Mat4.identity(),
                .transformed_m_inv = matrix.Mat4.identity(),
                .transposed_m = matrix.Mat4.identity(),
                .transposed_m_inv = matrix.Mat4.identity(),
                .material = Material.init(),
            },
            .radius = S(1),
        };
    }

    pub inline fn id(self: *const Sphere) usize {
        // log(@src(), "\nNext Sphere Id:{}\n", .{self.h.object_id});
        return self.h.object_id;
    }

    pub fn intersection(t: Scalar, sphere: *const Sphere) !Intersection {
        const i: Intersection = .{ .t = t, .object_id = sphere.h.object_id };
        return i;
    }

    /// ---
    /// A ray hitting a sphere can at most have two intersections.
    /// ---
    pub const LocalHits = struct {
        count: usize = 0,
        t: [2]Scalar = .{ S(0), S(0) },
    };

    /// ---
    /// Compute a local ray by applying the inverse of the sphere
    /// transform. Use the local rays origin and direction to
    /// compute the Intersections.
    /// ---
    pub fn intersect(self: *const Sphere, ray: Ray) LocalHits {
        const localray = Ray{
            .origin = self.h.transformed_m_inv.mulT(ray.origin),
            .direction = self.h.transformed_m_inv.mulT(ray.direction),
        };

        const sphere2ray = localray.origin.sub(point(0, 0, 0));
        const a = localray.direction.dot(localray.direction);
        const b = S(2) * localray.direction.dot(sphere2ray);
        const c = sphere2ray.dot(sphere2ray) - S(1);
        const discriminant = b * b - S(4) * a * c;

        if (discriminant < S(0)) return .{ .count = 0 };

        const sqrt_disc = std.math.sqrt(discriminant);
        const t1 = (-b - sqrt_disc) / (S(2) * a);
        const t2 = (-b + sqrt_disc) / (S(2) * a);

        return .{ .count = 2, .t = .{ t1, t2 } };
    }

    /// ---
    /// NOTE: Compute the normal at the given world_point.
    ///       This function uses cached matrix's for the
    ///       inverted transform and the inverted transpose.
    /// ---
    pub fn normal_at(self: *const Sphere, world_point: Tuple) Tuple {
        const object_point = self.h.transformed_m_inv.mulT(world_point);
        const object_normal = (object_point.sub(Tuple.point(0, 0, 0)));
        var world_normal = self.h.transposed_m_inv.mulT(object_normal);
        world_normal.w = S(0);
        return world_normal.normalize();
    }

    /// ---
    /// Setting the transform matrix.
    /// NOTE: Also computes the inverse as a side effect...
    /// ---
    pub fn set_transform(self: *Sphere, m: matrix.Mat4) void {
        self.h.transformed_m = m;
        self.h.transformed_m_inv = m.inverse();
        self.h.transposed_m = m.transpose();
        self.h.transposed_m_inv = m.transpose().inverse();
    }

    pub fn transform(self: *const Sphere) matrix.Mat4 {
        return self.h.transformed_m;
    }

    pub fn inverse(self: *const Sphere) matrix.Mat4 {
        return self.h.transformed_m_inv;
    }

    pub fn reset_id() void {
        // const obj_id = NEXT_SPHERE_ID.load(.seq_cst);
        // log(@src(), "\nNext Sphere Id:{}\n", .{obj_id});
        NEXT_SPHERE_ID.store(1, .seq_cst);
    }
};

test "Chap5 -A ray intersects a sphere at two points" {
    const r = Ray.init(point(-5, 0, 0), vector(1, 0, 0));
    var s = Sphere.init();
    const xs = s.intersect(r);

    try std.testing.expect(xs.count == 2);
    try std.testing.expect((xs.t[0]) == S(4));
    try std.testing.expect((xs.t[1]) == S(6));
}

test "Chap5 -A ray intersects a sphere at a tangent" {
    const r = Ray.init(point(0, 1, -5), vector(0, 0, 1));
    var s = Sphere.init();
    const xs = s.intersect(r);

    try std.testing.expect(xs.count == 2);
    try std.testing.expect((xs.t[0]) == S(5));
    try std.testing.expect((xs.t[1]) == S(5));
}

test "Chap5 -A ray misses a sphere" {
    const r = Ray.init(point(0, 2, -5), vector(0, 0, 1));
    var s = Sphere.init();
    const xs = s.intersect(r);

    try std.testing.expect(xs.count == 0);
}

test "Chap5 -A ray originates inside a sphere" {
    const r = Ray.init(point(0, 0, 0), vector(0, 0, 1));
    var s = Sphere.init();
    const xs = s.intersect(r);

    try std.testing.expect(xs.count == 2);
    try std.testing.expect((xs.t[0]) == S(-1));
    try std.testing.expect((xs.t[1]) == S(1));
}

test "Chap5 -A sphere is behind a ray" {
    const r = Ray.init(point(0, 0, 5), vector(0, 0, 1));
    var s = Sphere.init();
    const xs = s.intersect(r);

    try std.testing.expect(xs.count == 2);
    try std.testing.expect((xs.t[0]) == S(-6));
    try std.testing.expect((xs.t[1]) == S(-4));
}

test "Chap6 -A sphere has a default material" {
    const m = Material.init();
    const s = Sphere.init();
    try std.testing.expect(m.equals(s.h.material));
}

test "Chap6 -A sphere may be assigned a material" {
    var m = Material.init();
    m.ambient = S(1);
    var s = Sphere.init();
    s.h.material = m;
    try std.testing.expect(m.equals(s.h.material));
    try std.testing.expect(s.h.material.equals(m));
}

test "A total failure" {
    try std.testing.expect(3 == 3);
}
