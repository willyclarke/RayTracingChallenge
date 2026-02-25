// src/shapes/shape.zig
const std = @import("std");
const print = @import("std").debug.print;
const utils = @import("../utils.zig");
const log = utils.log;

const types = @import("../types.zig");
const matrix = @import("../matrix.zig");
const material_mod = @import("../material.zig");
const intersection_mod = @import("intersections.zig");
pub const sphere_mod = @import("sphere.zig");
pub const Sphere = sphere_mod.Sphere;

pub const Mat4 = matrix.Mat4;

pub const approxEq = types.approxEq;
pub const Ray = types.Ray;
pub const S = types.S;
pub const Scalar = types.Scalar;
pub const Tuple = types.Tuple;
pub const point = types.Point;
pub const vector = types.Vector;

pub const Intersection = intersection_mod.Intersection;
pub const Intersections = intersection_mod.Intersections;
pub const LocalIntersections = intersection_mod.LocalIntersections;

const material = material_mod.Material;
const Material = material_mod.Material;

// BEST PATTERN: union(enum) — no manual enum needed!
pub const Shape = union(enum) {
    sphere: *Sphere,
    // cube: Cube,
    // plane: Plane,

    // Real methods — these work because they're inside the struct scope
    pub fn id(self: *const Shape) usize {
        return switch (self.*) {
            .sphere => |s| s.id(),
            // .plane => |p| p.id(),
        };
    }

    // Helper: wrap a Sphere into a Shape using a mutable pointer
    pub fn fromSphere(s: *Sphere) Shape {
        return .{ .sphere = s };
    }

    pub fn intersect(self: *const Shape, ray: Ray) LocalIntersections {
        var xs = LocalIntersections.init();

        switch (self.*) {
            .sphere => |sp| {
                const lh = sp.intersect(ray);
                var i: usize = 0;
                while (i < lh.count) : (i += 1) {
                    xs.add(Shape.intersection(lh.t[i], self), self.id());
                }
            },
            // .cube => |c| c.intersect(ray),
        }

        return xs;
    }

    /// ---
    /// Factory: create an Intersection from any shape
    /// ---
    pub fn intersection(t: Scalar, ptrShape: *const Shape) Intersection {
        return .{ .t = t, .ptrShape = ptrShape };
    }

    /// ---
    /// return a pointer to the cached inverse. Does not calculate.
    /// ---
    pub fn inverse(self: *const Shape) *matrix.Mat4 {
        return switch (self.*) {
            .sphere => |s| &s.h.transformed_m_inv,
            // .cube => |c| c.inverse(m),
        };
    }

    pub fn normal_at(self: *const Shape, position: Tuple) Tuple {
        return switch (self.*) {
            .sphere => |s| s.normal_at(position),
            // .cube => |c| c.normal_at(position),
        };
    }

    pub fn material(self: *const Shape) *Material {
        return switch (self.*) {
            .sphere => |s| &s.h.material,
            // .cube => |c| c.material,
        };
    }

    pub fn setMaterial(self: *Shape, mat: Material) void {
        self.material().* = mat;
        // switch (self) {
        //     .sphere => |s| s.h.material = mat,
        //     // .cube => |c| c.material,
        // }
    }

    pub fn setTransform(self: *Shape, m: *const matrix.Mat4) void {
        switch (self.*) {
            .sphere => |s| {
                s.h.transformed_m = m.*;
                s.h.transposed_m = m.transpose();
                s.h.transformed_m_inv = m.inverse();
                s.h.transposed_m_inv = s.h.transformed_m_inv.transpose();
            },
            // .cube => |*c| c.set_transform(m),
        }
    }

    pub fn transform(self: *const Shape) *matrix.Mat4 {
        return switch (self.*) {
            .sphere => |s| &s.h.transformed_m,
            // .cube => |c| c.transform(m),
        };
    }

    pub fn reset_id(self: *const Shape) void {
        return switch (self.*) {
            .sphere => |s| s.reset_id(),
            // .cube => |c| c.transform(m),
        };
    }
};

// Your tests — now work perfectly!
test "sphere works" {
    var sphere = Sphere.init();
    const s = Shape.fromSphere(&sphere);
    try std.testing.expect(s.sphere.radius == 1.0);
}

test "Chap5 - An intersection encapsulates t and object" {
    var sphere = Sphere.init();
    const s = Shape.fromSphere(&sphere);
    const i = Shape.intersection(3.5, &s); // create intersection
    try std.testing.expect(s.id() == i.ptrShape.id());
    try std.testing.expect(types.approxEq(i.t, 3.5));
}

test "Chap5 -Aggregating intersections" {
    var sphere = Sphere.init();
    const s = Shape.fromSphere(&sphere); // wrap it
    const i_1 = Shape.intersection(1, &s); // create intersection
    const i_2 = Shape.intersection(2, &s); // create intersection

    var xs = Intersections.aggregate(std.testing.allocator, .{ i_1, i_2 });
    defer xs.deinit(std.testing.allocator);

    try std.testing.expect(xs.items().len == 2);
    try std.testing.expect(xs.intersections_items.items[0].t == S(1));
    try std.testing.expect(xs.intersections_items.items[1].t == S(2));
    try std.testing.expect(xs.intersections_items.items[0].t == i_1.t);
    try std.testing.expect(xs.intersections_items.items[1].t == i_2.t);
}

test "Chap5 -Intersect sets the object on the intersection" {
    const r = Ray.init(point(0, 0, -5), vector(0, 0, 1));
    var sphere = Sphere.init();
    const s = Shape.fromSphere(&sphere);
    // const xs = Shape.intersect(&s, r);
    const xs = s.intersect(r);
    try std.testing.expect(xs.count == 2);
    try std.testing.expect(xs.items[0].ptrShape.id() == s.id());
    try std.testing.expect(xs.items[1].ptrShape.id() == s.id());
}

test "Chap5 -The hit, when all intersections have positive t" {
    var sphere = Sphere.init();
    const s = Shape.fromSphere(&sphere);
    const i_1 = Shape.intersection(1, &s); // create intersection
    const i_2 = Shape.intersection(2, &s); // create intersection
    var xs = Intersections.aggregate(std.testing.allocator, .{ i_1, i_2 });
    defer xs.deinit(std.testing.allocator);

    const i = xs.hit();

    try std.testing.expect(i.?.ptrShape.id() == i_1.ptrShape.id());
    try std.testing.expect(i != null);
    try std.testing.expect(std.meta.eql(i, i_1));
    try std.testing.expect(!std.meta.eql(i, i_2));
}

test "Chap5 -The hit, when some intersections have negative t" {
    var sphere = Sphere.init();
    const s = Shape.fromSphere(&sphere);
    const i_1 = Shape.intersection(-1, &s); // create intersection
    const i_2 = Shape.intersection(1, &s); // create intersection
    var xs = Intersections.aggregate(std.testing.allocator, .{ i_1, i_2 });
    defer xs.deinit(std.testing.allocator);

    const i = xs.hit();

    try std.testing.expect(i.?.ptrShape.id() == i_2.ptrShape.id());
    try std.testing.expect(std.meta.eql(i, i_2));
    try std.testing.expect(i.?.eql(i_2));
    try std.testing.expect(!i.?.eql(i_1));
}

test "Chap5 -The hit, when all intersections have negative t" {
    var sphere = Sphere.init();
    const s = Shape.fromSphere(&sphere);
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
    var sphere = Sphere.init();
    const s = Shape.fromSphere(&sphere);
    const i_1 = Shape.intersection(5, &s); // create intersection
    const i_2 = Shape.intersection(7, &s); // create intersection
    const i_3 = Shape.intersection(-3, &s); // create intersection
    const i_4 = Shape.intersection(2, &s); // create intersection
    var xs = Intersections.aggregate(std.testing.allocator, .{ i_1, i_2, i_3, i_4 });
    defer xs.deinit(std.testing.allocator);

    const i = xs.hit();

    try std.testing.expect(i.?.ptrShape.id() == i_4.ptrShape.id());
    try std.testing.expect(std.meta.eql(i, i_4));
    try std.testing.expect(!i.?.eql(i_1));
    try std.testing.expect(!i.?.eql(i_2));
    try std.testing.expect(!i.?.eql(i_3));
    try std.testing.expect(i.?.eql(i_4));
    try std.testing.expect(std.meta.eql(xs.hit(), i_4));
}

test "Chap5 -A sphere's default transformation" {
    var sphere = Sphere.init();
    const s = Shape.fromSphere(&sphere);
    try std.testing.expect(matrix.Mat4.identity().equals(s.transform()));
}

test "Chap5 -Changing a sphere's transformation" {
    var sphere = Sphere.init();
    var s = Shape.fromSphere(&sphere);
    const t = matrix.Mat4.translation(2, 3, 4);
    s.setTransform(&t);
    try std.testing.expect(s.transform().equals(&t));
}

test "Chap5 -Intersecting a scaled sphere with a ray" {
    const r = Ray.init(point(0, 0, -5), vector(0, 0, 1));
    var sphere = Sphere.init();
    var s = Shape.fromSphere(&sphere);
    s.setTransform(&matrix.Mat4.scaling(2, 2, 2));
    const xs = s.intersect(r);

    try std.testing.expect(xs.count == 2);
    try std.testing.expect(approxEq(xs.items[0].t, S(3)));
    try std.testing.expect(approxEq(xs.items[1].t, S(7)));
}

test "Chap5 -Intersecting a translated sphere with a ray" {
    const r = Ray.init(point(0, 0, -5), vector(0, 0, 1));
    var sphere = Sphere.init();
    var s = Shape.fromSphere(&sphere);
    s.setTransform(&matrix.Mat4.translation(5, 0, 0));
    const xs = s.intersect(r);
    try std.testing.expect(xs.count == 0);
}

test "Chap6 -The normal on a sphere at a point on the x axis" {
    var sphere = Sphere.init();
    const s = Shape.fromSphere(&sphere);
    const n = s.normal_at(point(1, 0, 0));
    try std.testing.expect(n.equals(vector(1, 0, 0)));
}

test "Chap6 -The normal on a sphere at a point on the y axis" {
    var sphere = Sphere.init();
    const s = Shape.fromSphere(&sphere);
    const n = s.normal_at(point(0, 1, 0));
    try std.testing.expect(n.equals(vector(0, 1, 0)));
}

test "Chap6 -The normal on a sphere at a point on the z axis" {
    var sphere = Sphere.init();
    const s = Shape.fromSphere(&sphere);
    const n = s.normal_at(point(0, 0, 1));
    try std.testing.expect(n.equals(vector(0, 0, 1)));
}

test "Chap6 -The normal on a sphere at a nonaxial point" {
    var sphere = Sphere.init();
    const s = Shape.fromSphere(&sphere);
    const sqrt3 = std.math.sqrt(S(3));
    const x = sqrt3 / S(3);
    const y = x;
    const z = x;
    const n = s.normal_at(point(x, y, z));
    try std.testing.expect(n.equals(vector(x, x, y)));
}

test "Chap6 -The normal is a normalized vector" {
    var sphere = Sphere.init();
    const s = Shape.fromSphere(&sphere);
    const sqrt3 = std.math.sqrt(S(3));
    const x = sqrt3 / S(3);
    const y = x;
    const z = x;
    const n = s.normal_at(point(x, y, z));
    try std.testing.expect(n.equals(n.normalize()));
}

test "Chap6 -Computing the normal on a translated sphere" {
    var sphere = Sphere.init();
    var s = Shape.fromSphere(&sphere);
    s.setTransform(&Mat4.translation(0, 1, 0));
    const x = S(0);
    const y = S(1.707106781186548);
    const z = S(-0.707106781186548);
    const n = s.normal_at(point(x, y, z));
    // log(@src(), "n: {f} z: {}\n", .{ n, z });
    // log(@src(), "e: {f} z: {}\n", .{ vector(0, -z, z), z });
    try std.testing.expect(n.equals(vector(0, -z, z)));
}

test "Chap6 -Computing the normal on a transformed sphere" {
    var sphere = Sphere.init();
    var s = Shape.fromSphere(&sphere);
    const m = Mat4.scaling(1, 0.5, 1).mulM(&Mat4.rotz(std.math.pi / S(5)));
    s.setTransform(&m);
    const x = S(0);
    const y = std.math.sqrt2 / S(2);
    const z = -std.math.sqrt2 / S(2);
    const n = s.normal_at(point(x, y, z));
    // log(@src(), "n: {f} z: {}\n", .{ n, z });
    try std.testing.expect(n.equals(vector(0, 0.970142500145332, -0.242535625036333)));
}

test "Chap6 -Reflecting a vector approaching at 45°" {
    const v = vector(1, -1, 0);
    const n = vector(0, 1, 0);
    const r = v.reflect(n);
    try std.testing.expect(r.equals(vector(1, 1, 0)));
}
