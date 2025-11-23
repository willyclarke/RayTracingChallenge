const std = @import("std");
const print = @import("std").debug.print;
const types = @import("../types.zig");
const canvas = @import("../canvas.zig");
const utils = @import("../utils.zig");

const S = types.S;
const Ray = types.Ray;
const Scalar = types.Scalar;
const Tuple = types.Tuple;
const Intersection = types.Intersection;
const Intersections = types.Intersections;
const point = types.Point;
const vector = types.Vector;
const approxEq = types.approxEq;
const log = utils.log;

var NEXT_SPHERE_ID: std.atomic.Value(usize) = .{ .raw = 0 };

pub const Sphere = struct {
    xs: Intersections,
    radius: Scalar,
    object_id: usize,

    pub fn init() Sphere {
        const id = NEXT_SPHERE_ID.fetchAdd(1, .seq_cst);
        return .{
            .xs = Intersections.init(),
            .radius = S(1),
            .object_id = id,
        };
    }
};

pub fn intersect(sphere: Sphere, ray: Ray) !Intersections {
    const sphere2ray = ray.origin.sub(point(0, 0, 0));
    const a = ray.direction.dot(ray.direction);
    const b = S(2) * ray.direction.dot(sphere2ray);
    const c = sphere2ray.dot(sphere2ray) - S(1);
    const discriminant = b * b - S(4) * a * c;

    if (discriminant < S(0)) {
        return Intersections.init();
    }

    const t1 = (-b - std.math.sqrt(discriminant)) / (S(2) * a);
    const t2 = (-b + std.math.sqrt(discriminant)) / (S(2) * a);

    const intersection0 = types.Intersection{ .t = t1, .object_id = sphere.object_id };
    const intersection1 = types.Intersection{ .t = t2, .object_id = sphere.object_id };

    var intersections = Intersections.init();
    try intersections.append(intersection0);
    try intersections.append(intersection1);

    return intersections;
}

pub fn intersection(t: Scalar, sphere: Sphere) !Intersection {
    const xx: Intersection = .{ .t = t, .object_id = sphere.object_id };
    return xx;
}

test "Chap5 -A ray intersects a sphere at two points" {
    const r = Ray.init(point(-5, 0, 0), vector(1, 0, 0));
    const s = Sphere.init();
    const xs = try intersect(s, r);

    try std.testing.expect(xs.count == 2);
    try std.testing.expect((try xs.get(0)).t == S(4));
    try std.testing.expect((try xs.get(1)).t == S(6));
}

test "Chap5 -A ray intersects a sphere at a tangent" {
    const r = Ray.init(point(0, 1, -5), vector(0, 0, 1));
    const s = Sphere.init();
    const xs = try intersect(s, r);

    try std.testing.expect(xs.count == 2);
    try std.testing.expect((try xs.get(0)).t == S(5));
    try std.testing.expect((try xs.get(1)).t == S(5));
}

test "Chap5 -A ray misses a sphere" {
    const r = Ray.init(point(0, 2, -5), vector(0, 0, 1));
    const s = Sphere.init();
    const xs = try intersect(s, r);

    try std.testing.expect(xs.count == 0);
}

test "Chap5 -A ray originates inside a sphere" {
    const r = Ray.init(point(0, 0, 0), vector(0, 0, 1));
    const s = Sphere.init();
    const xs = try intersect(s, r);

    try std.testing.expect(xs.count == 2);
    try std.testing.expect((try xs.get(0)).t == S(-1));
    try std.testing.expect((try xs.get(1)).t == S(1));
}

test "Chap5 -A sphere is behind a ray" {
    const r = Ray.init(point(0, 0, 5), vector(0, 0, 1));
    const s = Sphere.init();
    const xs = try intersect(s, r);

    try std.testing.expect(xs.count == 2);
    try std.testing.expect((try xs.get(0)).t == S(-6));
    try std.testing.expect((try xs.get(1)).t == S(-4));
}

test "Chap5 -An intersection encapsulates t and object" {
    const r = Ray.init(point(0, 0, 5), vector(0, 0, 1));
    const s = Sphere.init();
    const xs = try intersect(s, r);

    try std.testing.expect(xs.count == 2);
    try std.testing.expect((try xs.get(0)).t == S(-6));
    try std.testing.expect((try xs.get(1)).t == S(-4));
}
