// src/shapes/shape.zig
const std = @import("std");
const print = @import("std").debug.print;
const utils = @import("utils.zig");
const log = utils.log;

const types = @import("types.zig");
const sphere_mod = @import("shapes/sphere.zig");
const Sphere = sphere_mod.Sphere;

pub const Scalar = types.Scalar;
pub const S = types.S;
pub const Ray = types.Ray;
pub const point = types.Point;
pub const vector = types.Vector;
pub const Intersection = types.Intersection;
pub const Intersections = types.Intersections;
pub const LocalIntersections = types.LocalIntersections;

// BEST PATTERN: union(enum) — no manual enum needed!
pub const Shape = union(enum) {
    sphere: Sphere,
    // cube: Cube,
    // plane: Plane,

    // Real methods — these work because they're inside the struct scope
    pub fn id(self: *const Shape) usize {
        return switch (self.*) {
            inline else => |obj| obj.id(),
        };
    }

    pub fn intersect(self: *const Shape, ray: Ray) LocalIntersections {
        return switch (self.*) {
            .sphere => |s| s.intersect(ray),
            // .cube => |c| c.intersect(ray),
        };
    }

    /// ---
    /// Factory: create an Intersection from any shape
    /// ---
    pub fn intersection(t: Scalar, shape: *const Shape) Intersection {
        return .{ .t = t, .object_id = shape.id() }; // shape.id() works!
    }

    // Helper: wrap a Sphere into a Shape
    pub fn fromSphere(s: Sphere) Shape {
        return .{ .sphere = s };
    }
};

// Your tests — now work perfectly!
test "sphere works" {
    const s = Sphere.init();
    try std.testing.expect(s.radius == 1.0);
}

test "Chap5 - An intersection encapsulates t and object" {
    const s = Shape.fromSphere(Sphere.init()); // wrap it
    const i = Shape.intersection(3.5, &s); // create intersection
    try std.testing.expect(s.id() == i.object_id);
    try std.testing.expect(types.approxEq(i.t, 3.5));
}

test "Chap5 -Aggregating intersections" {
    const s = Shape.fromSphere(Sphere.init()); // wrap it
    const i_1 = Shape.intersection(1, &s); // create intersection
    const i_2 = Shape.intersection(2, &s); // create intersection

    var xs = Intersections.aggregate(std.testing.allocator, .{ i_1, i_2 });
    defer xs.deinit(std.testing.allocator);

    try std.testing.expect(xs.count() == 2);
    try std.testing.expect((try xs.get(0)).t == S(1));
    try std.testing.expect((try xs.get(1)).t == S(2));
    try std.testing.expect((try xs.get(0)).t == i_1.t);
    try std.testing.expect((try xs.get(1)).t == i_2.t);
}

test "Chap5 -Intersect sets the object on the intersection" {
    const r = Ray.init(point(0, 0, -5), vector(0, 0, 1));
    const s = Shape.fromSphere(Sphere.init());
    // const xs = Shape.intersect(&s, r);
    const xs = s.intersect(r);
    try std.testing.expect(xs.count == 2);
    try std.testing.expect(xs.local_intersections_items[0].object_id == s.id());
    try std.testing.expect(xs.local_intersections_items[1].object_id == s.id());
}

test "Chap5 -The hit, when all intersections have positive t" {
    const s = Shape.fromSphere(Sphere.init());
    const i_1 = Shape.intersection(1, &s); // create intersection
    const i_2 = Shape.intersection(2, &s); // create intersection
    var xs = Intersections.aggregate(std.testing.allocator, .{ i_1, i_2 });
    defer xs.deinit(std.testing.allocator);

    const i = xs.hit();

    try std.testing.expect(i.?.object_id == i_1.object_id);
    try std.testing.expect(i != null);
    try std.testing.expect(std.meta.eql(i, i_1));
    try std.testing.expect(!std.meta.eql(i, i_2));
}

test "Chap5 -The hit, when some intersections have negative t" {
    const s = Shape.fromSphere(Sphere.init());
    const i_1 = Shape.intersection(-1, &s); // create intersection
    const i_2 = Shape.intersection(1, &s); // create intersection
    var xs = Intersections.aggregate(std.testing.allocator, .{ i_1, i_2 });
    defer xs.deinit(std.testing.allocator);

    const i = xs.hit();

    try std.testing.expect(i.?.object_id == i_2.object_id);
    try std.testing.expect(std.meta.eql(i, i_2));
    try std.testing.expect(i.?.eql(i_2));
    try std.testing.expect(!i.?.eql(i_1));
}

test "Chap5 -The hit, when all intersections have negative t" {
    const s = Shape.fromSphere(Sphere.init());
    const i_1 = Shape.intersection(-2, &s); // create intersection
    const i_2 = Shape.intersection(-1, &s); // create intersection
    var xs = Intersections.aggregate(std.testing.allocator, .{ i_1, i_2 });
    defer xs.deinit(std.testing.allocator);

    const i = xs.hit();

    try std.testing.expect(i == null);
    try std.testing.expect(!std.meta.eql(i, i_1));
    try std.testing.expect(!std.meta.eql(i, i_2));
}

test "Chap5 -The hit is always the lowest nonnegative intersection" {
    const s = Shape.fromSphere(Sphere.init());
    const i_1 = Shape.intersection(5, &s); // create intersection
    const i_2 = Shape.intersection(7, &s); // create intersection
    const i_3 = Shape.intersection(-3, &s); // create intersection
    const i_4 = Shape.intersection(2, &s); // create intersection
    var xs = Intersections.aggregate(std.testing.allocator, .{ i_1, i_2, i_3, i_4 });
    defer xs.deinit(std.testing.allocator);

    const i = xs.hit();

    try std.testing.expect(i.?.object_id == i_4.object_id);
    try std.testing.expect(std.meta.eql(i, i_4));
    try std.testing.expect(!i.?.eql(i_1));
    try std.testing.expect(!i.?.eql(i_2));
    try std.testing.expect(!i.?.eql(i_3));
    try std.testing.expect(i.?.eql(i_4));
    try std.testing.expect(std.meta.eql(xs.hit(), i_4));
}
