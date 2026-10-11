// src/shapes/intersections.zig
//
const std = @import("std");
const print = @import("std").debug.print;
const utils = @import("../utils.zig");
const types = @import("../types.zig");
const shapes = @import("shapes.zig");

const approxEq = types.approxEq;
const Scalar = types.Scalar;
const S = types.S;
const Shape = shapes.Shape;

/// ---
/// Intersection used by the various types of objects like spheres, cubes, etc...
/// ---
pub const Intersection = struct {
    t: Scalar,
    // object_id: usize,
    ptrShape: *const Shape,

    pub fn eql(a: Intersection, b: Intersection) bool {
        return approxEq(a.t, b.t) and a.ptrShape.id() == b.ptrShape.id();
    }

    pub fn init() Intersection {
        return .{
            .t = S(0),
            .object_id = undefined,
        };
    }
};

/// ---
/// Fixed size used by LocalIntersections.
/// ---
pub const MaxIntersectionsPerShape = 2;

/// ---
/// return: array of up to two intersections.
/// ---
pub const LocalIntersections = struct {
    count: usize,
    items: [MaxIntersectionsPerShape]Intersection,
    min_t_intersection: Intersection,
    object_id: usize,

    pub fn init() LocalIntersections {
        return .{
            .count = 0, //
            .items = undefined, //
            .min_t_intersection = undefined, //
            .object_id = undefined, //
        };
    }

    pub fn add(self: *LocalIntersections, i: Intersection, object_id: usize) void {
        if (self.count < MaxIntersectionsPerShape) {
            self.items[self.count] = i;
            self.count += 1;
            self.object_id = object_id;

            if (self.count == 1) {
                self.min_t_intersection = i;
            } else if (i.t < self.min_t_intersection.t) {
                self.min_t_intersection = i;
            }
        }
    }

    pub fn tmin(self: *const LocalIntersections) Scalar {
        return self.min_t_intersection.t;
    }

    pub fn hit(self: *const LocalIntersections) bool {
        return self.count > 0;
    }
};

/// ---
/// Intersections used by the various types of objects like spheres, cubes, etc...
/// ---
pub const Intersections = struct {
    intersections_items: std.ArrayListUnmanaged(Intersection) = .{},

    pub fn init() Intersections {
        return .{};
    }

    pub fn deinit(self: *Intersections, allocator: std.mem.Allocator) void {
        self.intersections_items.deinit(allocator);
    }

    // pub fn count(self: *const Intersections) usize {
    //     return self.intersections_items.items.len;
    // }

    pub fn add(self: *Intersections, allocator: std.mem.Allocator, hit_to_add: Intersection) !void {
        self.intersections_items.append(allocator, hit_to_add) catch @panic("OOM");
    }

    pub fn items(self: *const Intersections) []const Intersection {
        return self.intersections_items.items;
    }

    /// ---
    /// Convenience helper for tests: aggregate a fixed list of intersections.
    /// ---
    pub fn aggregate(allocator: std.mem.Allocator, ints: anytype) Intersections {
        var xs = Intersections.init();
        inline for (ints) |i| {
            try xs.add(allocator, i);
        }
        return xs;
    }

    /// Custom formatter so `{}` prints nicely.
    /// `fmt` and `options` let you add variants later; for now we ignore them.
    pub fn format(self: Intersections, writer: anytype) !void {
        try writer.print("Count: {}. ", .{self.items().len});
        for (self.items()) |i| {
            try writer.print("Shape {} at t [ {} ] :: ", .{ i.object_id, i.t });
        }
        try writer.print("\n", .{});
    }

    /// ---
    /// Return the Intersection item with the t > 0 when there
    /// has been a hit and the initialized version otherwise.
    /// i.e. no hit: t=0 and object_id=0.
    /// ---
    pub fn hit(self: *const Intersections) ?Intersection {
        var best: ?Intersection = null;

        for (self.intersections_items.items) |i| {
            // if (i.t >= 0) { // includes t=0!
            if (i.t > 0) {
                if (best) |b| {
                    if (i.t < b.t) best = i;
                } else {
                    best = i;
                }
            }
        }
        return best;
    }
};

test "basic add" {
    try std.testing.expect(3 + 7 == 10);
}
