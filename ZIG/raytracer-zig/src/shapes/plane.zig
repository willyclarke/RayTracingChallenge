const std = @import("std");
const print = @import("std").debug.print;
const tMod = @import("../types.zig");
// const canvas = @import("../canvas.zig");
const utils = @import("../utils.zig");

const shapes = @import("shapes.zig");
const matrix = @import("../matrix.zig");
const mat_module = @import("../material.zig");
const intersection_mod = @import("intersections.zig");
const local_hits_mod = @import("local_hits.zig");

pub const Mat4 = matrix.Mat4;
const ShapeHeader = @import("shape_header.zig").ShapeHeader;

const S = tMod.S;
const Ray = tMod.Ray;
const Scalar = tMod.Scalar;
const Shape = shapes.Shape;
const Tuple = tMod.Tuple;
const Matrix = matrix.Mat4;
const Intersection = intersection_mod.Intersection;
const Intersections = intersection_mod.Intersections;
const LocalIntersections = intersection_mod.LocalIntersections;
const point = tMod.Point;
const vector = tMod.Vector;
const approxEq = tMod.approxEq;
const log = utils.log;
const Material = mat_module.Material;
const LocalHits = local_hits_mod.LocalHits;
const local_hits = local_hits_mod.LocalHits;

var NEXT_PLANE_ID: std.atomic.Value(usize) = .{ .raw = 1 };

pub const Plane = struct {
    h: ShapeHeader,

    pub const PlaneHits = local_hits(1);

    pub fn init() Plane {
        const obj_id = NEXT_PLANE_ID.fetchAdd(1, .seq_cst);
        // log(@src(), "\nNext Plane Id:{}\n", .{obj_id});
        return .{
            .h = ShapeHeader.init(obj_id),
        };
    }

    pub inline fn id(self: *const Plane) usize {
        // log(@src(), "\nNext Plane Id:{}\n", .{self.h.object_id});
        return self.h.object_id;
    }

    // The logic to intersect a ray with a plane is the only other bit that
    // needs implementing, and it has four cases to consider:

    // 1. The ray is
    // parallel to the plane, and will thus never intersect it.

    // 2. The ray is coplanar with the plane, which is to say that the ray’s
    // origin is on the plane, and the ray’s direction is parallel to the
    // plane. You’re viewing the plane edge-on. In this case, every point on
    // the ray intersects the plane, resulting in an infinite number of
    // intersections. That’s unwieldy! But since a plane is infinitely thin,
    // it’s invisible when viewed like this, so we’ll assume the ray misses in
    // this case.

    // 3. The ray origin is above the plane.

    // 4. The ray origin is below the plane.

    pub fn local_intersect(localray: Ray) PlaneHits {
        if (@abs(localray.direction.y) < tMod.EPSILON) {
            return .{ .count = 0, .t = .{S(0)} };
        }
        const t1 = -localray.origin.y / localray.direction.y;
        return .{ .count = 1, .t = .{t1} };
    }

    /// ---
    /// Recieve a local ray by that has been processed by applying the
    /// inverse of the sphere transform. Use the local rays origin and
    /// direction to compute the Intersections.
    /// ---
    pub fn intersect(self: *const Plane, ray: Ray) PlaneHits {
        const localray = Ray{
            .origin = self.h.transformed_m_inv.mulT(ray.origin),
            .direction = self.h.transformed_m_inv.mulT(ray.direction),
        };

        return local_intersect(localray);
    }

    /// ---
    /// NOTE: returns the plane-specific normal
    /// ---
    pub fn local_normal_at(self: *const Plane, local_point: Tuple) Tuple {
        _ = self;
        _ = local_point;
        return vector(0, 1, 0);
    }

    /// ---
    /// Setting the transform matrix.
    /// NOTE: Also computes the inverse as a side effect...
    /// ---
    pub fn set_transform(self: *Plane, m: matrix.Mat4) void {
        return self.h.set_transform(m);
    }

    pub fn transform(self: *const Plane) matrix.Mat4 {
        return self.h.transformed_m;
    }

    pub fn inverse(self: *const Plane) matrix.Mat4 {
        return self.h.transformed_m_inv;
    }

    pub fn reset_id() void {
        const obj_id = NEXT_PLANE_ID.load(.seq_cst);
        log(@src(), "\nNext Plane Id:{}\n", .{obj_id});
        NEXT_PLANE_ID.store(1, .seq_cst);
    }
};

test "Make sure it works" {
    try std.testing.expect(3 == 3);
}

test "Chap9 -The normal of a plane is constant everywhere" {
    const pl = Plane.init();
    const n1 = pl.local_normal_at(point(0, 0, 0));
    const n2 = pl.local_normal_at(point(10, 0, -10));
    const n3 = pl.local_normal_at(point(-5, 0, 150));
    try std.testing.expect(pl.h.transformed_m.equals(&Mat4.identity()));
    try std.testing.expect(n1.equals(vector(0, 1, 0)));
    try std.testing.expect(n2.equals(vector(0, 1, 0)));
    try std.testing.expect(n3.equals(vector(0, 1, 0)));
}

test "Chap9 -Intersect with a ray parallel to the plane" {
    const pl = Plane.init();
    _ = pl;
    const r = Ray.init(point(0, 10, 0), vector(0, 0, 1));
    const xs = Plane.local_intersect(r);
    try std.testing.expect(xs.count == 0);
}

test "Chap9 -Intersect with a coplanar ray" {
    const r = Ray.init(point(0, 0, 0), vector(0, 0, 1));
    const xs = Plane.local_intersect(r);
    try std.testing.expect(xs.count == 0);
}

test "Chap9 -A ray intersecting a plane from above" {
    const r = Ray.init(point(0, -1, 0), vector(0, 1, 0));
    const xs = Plane.local_intersect(r);
    try std.testing.expect(xs.count == 1);
    try std.testing.expect(approxEq(S(1), xs.items()[0]));
    // Can not really do the id check based on local_intersect
    // since the object id is stored in the header and provided
    // by the Intersect.
}
