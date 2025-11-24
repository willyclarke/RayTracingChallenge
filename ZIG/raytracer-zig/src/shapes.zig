// src/shapes/shape.zig
const std = @import("std");
const print = @import("std").debug.print;

const types = @import("types.zig");
const sphere_mod = @import("shapes/sphere.zig");
const Sphere = sphere_mod.Sphere;

pub const Scalar = types.Scalar;
pub const Ray = types.Ray;
pub const Intersections = types.Intersections;
pub const Intersection = types.Intersection;

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

    pub fn intersect(self: *const Shape, ray: Ray) Intersections {
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
}
