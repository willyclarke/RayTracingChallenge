// src/shapes/shape.zig
const std = @import("std");
const print = @import("std").debug.print;
const utils = @import("utils.zig");
const log = utils.log;

const types = @import("types.zig");
const matrix = @import("matrix.zig");
pub const sphere_mod = @import("shapes/sphere.zig");
pub const Sphere = sphere_mod.Sphere;

pub const approxEq = types.approxEq;
pub const Mat4 = matrix.Mat4;
pub const Ray = types.Ray;
pub const S = types.S;
pub const Scalar = types.Scalar;
pub const Tuple = types.Tuple;

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

    // Helper: wrap a Sphere into a Shape
    pub fn fromSphere(s: Sphere) Shape {
        return .{ .sphere = s };
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

    pub fn inverse(self: *const Shape) matrix.Mat4 {
        return switch (self.*) {
            .sphere => |s| s.inverse(),
            // .cube => |c| c.inverse(m),
        };
    }

    pub fn normal_at(self: *const Shape, position: Tuple) Tuple {
        return switch (self.*) {
            .sphere => |s| s.normal_at(position),
            // .cube => |c| c.normal_at(position),
        };
    }

    pub fn set_transform(self: *Shape, m: *const matrix.Mat4) void {
        return switch (self.*) {
            .sphere => |*s| s.set_transform(m),
            // .cube => |*c| c.set_transform(m),
        };
    }

    pub fn transform(self: *const Shape) matrix.Mat4 {
        return switch (self.*) {
            .sphere => |s| s.transform(),
            // .cube => |c| c.transform(m),
        };
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

test "Chap5 -A sphere's default transformation" {
    const s = Shape.fromSphere(Sphere.init());
    try std.testing.expect(matrix.Mat4.identity().equals(&s.transform()));
}

test "Chap5 -Changing a sphere's transformation" {
    var s = Shape.fromSphere(Sphere.init());
    const t = matrix.Mat4.translation(2, 3, 4);
    s.set_transform(&t);
    try std.testing.expect(s.transform().equals(&t));
}

test "Chap5 -Intersecting a scaled sphere with a ray" {
    const r = Ray.init(point(0, 0, -5), vector(0, 0, 1));
    var s = Shape.fromSphere(Sphere.init());
    s.set_transform(&matrix.Mat4.scaling(2, 2, 2));
    const xs = s.intersect(r);
    try std.testing.expect(xs.count == 2);
    try std.testing.expect(approxEq(xs.local_intersections_items[0].t, S(3)));
    try std.testing.expect(approxEq(xs.local_intersections_items[1].t, S(7)));
}

test "Chap5 -Intersecting a translated sphere with a ray" {
    const r = Ray.init(point(0, 0, -5), vector(0, 0, 1));
    var s = Shape.fromSphere(Sphere.init());
    s.set_transform(&matrix.Mat4.translation(5, 0, 0));
    const xs = s.intersect(r);
    try std.testing.expect(xs.count == 0);
}

test "Chap6 -The normal on a sphere at a point on the x axis" {
    const s = Shape.fromSphere(Sphere.init());
    const n = s.normal_at(point(1, 0, 0));
    try std.testing.expect(n.equals(vector(1, 0, 0)));
}

test "Chap6 -The normal on a sphere at a point on the y axis" {
    const s = Shape.fromSphere(Sphere.init());
    const n = s.normal_at(point(0, 1, 0));
    try std.testing.expect(n.equals(vector(0, 1, 0)));
}

test "Chap6 -The normal on a sphere at a point on the z axis" {
    const s = Shape.fromSphere(Sphere.init());
    const n = s.normal_at(point(0, 0, 1));
    try std.testing.expect(n.equals(vector(0, 0, 1)));
}

test "Chap6 -The normal on a sphere at a nonaxial point" {
    const s = Shape.fromSphere(Sphere.init());
    const sqrt3 = std.math.sqrt(S(3));
    const x = sqrt3 / S(3);
    const y = x;
    const z = x;
    const n = s.normal_at(point(x, y, z));
    try std.testing.expect(n.equals(vector(x, x, y)));
}

test "Chap6 -The normal is a normalized vector" {
    const s = Shape.fromSphere(Sphere.init());
    const sqrt3 = std.math.sqrt(S(3));
    const x = sqrt3 / S(3);
    const y = x;
    const z = x;
    const n = s.normal_at(point(x, y, z));
    try std.testing.expect(n.equals(n.normalize()));
}

test "Chap6 -Computing the normal on a translated sphere" {
    var s = Shape.fromSphere(Sphere.init());
    s.set_transform(&Mat4.translation(0, 1, 0));
    const x = S(0);
    const y = S(1.70711);
    const z = S(-0.70711);
    const n = s.normal_at(point(x, y, z));
    try std.testing.expect(n.equals(vector(0, -z, z)));
}

test "Chap6 -Computing the normal on a transformed sphere" {
    var s = Shape.fromSphere(Sphere.init());
    const m = Mat4.scaling(1, 0.5, 1).mulM(&Mat4.rotz(std.math.pi / S(5)));
    s.set_transform(&m);
    const x = S(0);
    const y = std.math.sqrt2 / S(2);
    const z = -std.math.sqrt2 / S(2);
    const n = s.normal_at(point(x, y, z));
    try std.testing.expect(n.equals(vector(0, S(0.97014), S(-0.24254))));
}

test "Chap6 -Reflecting a vector approaching at 45°" {
    const v = vector(1, -1, 0);
    const n = vector(0, 1, 0);
    const r = v.reflect(n);
    try std.testing.expect(r.equals(vector(1, 1, 0)));
}
