const std = @import("std");
const print = @import("std").debug.print;
const utils = @import("utils.zig");
const log = utils.log;

pub const Scalar = f32; // switch to f64 later if you need it
//
/// Convert ints/floats (incl. comptime literals) to Scalar (f32).
/// Convert ints or floats (incl. comptime literals) to Scalar (f32).
pub inline fn scalar(x: anytype) Scalar {
    return switch (@TypeOf(x)) {
        comptime_int, comptime_float => @as(Scalar, x),
        // handle all builtin float widths
        f16, f32, f64, f128 => @floatCast(x),
        // anything else we care about here should be an integer type
        else => @as(Scalar, @floatFromInt(x)),
    };
}

/// ---
/// Alias to convert anytype to the Scalar type
/// ---
pub inline fn S(x: anytype) Scalar {
    return scalar(x);
}

pub inline fn Deg2Rad(x: anytype) Scalar {
    return scalar(std.math.pi / S(180) * x);
}

pub fn toByteSaturated(x: Scalar) u8 {
    const clamped = std.math.clamp(std.math.round(x), 0.0, 255.0);
    return @intFromFloat(clamped); // truncates toward 0
}

pub fn toUsizeSaturated(x: Scalar, min: Scalar, max: Scalar) usize {
    const clamped = std.math.clamp(std.math.round(x), min, max);
    return @intFromFloat(clamped); // truncates toward 0
}

pub const EPSILON: Scalar = 1e-5;

// pub inline fn approxEq(a: Scalar, b: Scalar) bool {
//     return @abs(a - b) < EPSILON;
// }

/// Returns true when |a - b| < EPSILON
/// Compare by using absolute tolerance for small values
/// and relative tolerance for big values.
pub inline fn approxEq(a: Scalar, b: Scalar) bool {
    const diff = @abs(a - b);
    if (diff <= EPSILON) return true;
    const norm = @max(@abs(a), @abs(b));
    return diff <= norm * EPSILON;
}

test "Chap1 -almostEqual works for close values" {
    try std.testing.expect(approxEq(1.000001, 1.000002));
    try std.testing.expect(!approxEq(1.0, 1.1));
}

pub const Tuple = struct {
    x: Scalar,
    y: Scalar,
    z: Scalar,
    w: Scalar,

    /// Construct from components.
    pub fn init(x: Scalar, y: Scalar, z: Scalar, w: Scalar) Tuple {
        return .{ .x = x, .y = y, .z = z, .w = w };
    }

    /// Convenience constructors used throughout the book.
    pub fn point(x: Scalar, y: Scalar, z: Scalar) Tuple {
        return .{ .x = x, .y = y, .z = z, .w = 1.0 };
    }

    pub fn vector(x: Scalar, y: Scalar, z: Scalar) Tuple {
        return .{ .x = x, .y = y, .z = z, .w = 0.0 };
    }

    pub fn equals(self: Tuple, other: Tuple) bool {
        return approxEq(self.x, other.x) and
            approxEq(self.y, other.y) and
            approxEq(self.z, other.z) and
            approxEq(self.w, other.w);
    }

    pub fn add(self: Tuple, other: Tuple) Tuple {
        return .{ .x = self.x + other.x, .y = self.y + other.y, .z = self.z + other.z, .w = self.w + other.w };
    }

    pub fn sub(self: Tuple, other: Tuple) Tuple {
        return .{ .x = self.x - other.x, .y = self.y - other.y, .z = self.z - other.z, .w = self.w - other.w };
    }

    pub fn neg(self: Tuple) Tuple {
        return .{ .x = -self.x, .y = -self.y, .z = -self.z, .w = -self.w };
    }

    /// Scalar multiplication
    pub fn muls(self: Tuple, s: Scalar) Tuple {
        return .{ .x = self.x * s, .y = self.y * s, .z = self.z * s, .w = self.w * s };
    }

    pub fn div(self: Tuple, s: Scalar) Tuple {
        return .{ .x = self.x / s, .y = self.y / s, .z = self.z / s, .w = self.w / s };
    }

    pub fn mag(self: Tuple) Scalar {
        return std.math.sqrt(mags(self));
    }

    /// Magnitude squared
    pub fn mags(self: Tuple) Scalar {
        return self.x * self.x + self.y * self.y + self.z * self.z + self.w * self.w;
    }

    /// Normalize to length 1
    pub fn normalize(self: Tuple) Tuple {
        return self.div(self.mag());
    }

    /// Inner product of a tuple
    pub fn dot(self: Tuple, other: Tuple) Scalar {
        return self.x * other.x + self.y * other.y + self.z * other.z + self.w * other.w;
    }

    /// Tuple multiplication
    /// Also called Hadamard product or Schur product
    pub fn mult(self: Tuple, other: Tuple) Tuple {
        return .{ .x = self.x * other.x, .y = self.y * other.y, .z = self.z * other.z, .w = self.w * other.w };
    }

    pub fn cross(self: Tuple, other: Tuple) Tuple {
        return vector(self.y * other.z - self.z * other.y, self.z * other.x - self.x * other.z, self.x * other.y - self.y * other.x);
    }

    pub fn reflect(in: Tuple, normal: Tuple) Tuple {
        return in.sub(normal.muls(S(2) * in.dot(normal)));
    }

    /// Custom formatter so `{}` prints nicely.
    /// `fmt` and `options` let you add variants later; for now we ignore them.
    pub fn format(self: Tuple, writer: anytype) !void {
        if (self.w != S(0)) {
            try writer.print("P({d:8.15}, {d:8.15}, {d:8.15}, {d:8.15})", .{ self.x, self.y, self.z, self.w });
        } else {
            try writer.print("V({d:8.15}, {d:8.15}, {d:8.15}, {d:8.15})", .{ self.x, self.y, self.z, self.w });
        }
    }

    /// Handle printing via a *const Tuple pointer
    pub fn formatPtr(self: *const Tuple, writer: anytype) !void {
        // Forward to the value formatter
        return Tuple.format(self, writer);
    }

    /// Color - red channel
    pub inline fn r(self: Tuple) Scalar {
        return self.x;
    }

    /// Color - green channel
    pub inline fn g(self: Tuple) Scalar {
        return self.y;
    }

    /// Color - blue channel
    pub inline fn b(self: Tuple) Scalar {
        return self.z;
    }

    /// Color - alpha channel
    pub inline fn alpha(self: Tuple) Scalar {
        return self.w;
    }
};

pub const Ray = struct {
    const Self = @This();

    origin: Tuple,
    direction: Tuple,

    pub fn init(origin: Tuple, direction: Tuple) Ray {
        return .{
            .origin = origin,
            .direction = direction,
        };
    }

    pub fn position(self: Self, t: Scalar) Tuple {
        const p = self.origin.add(self.direction.muls(t));
        return p;
    }
};

/// Alias: Color *is* Tuple (same type)
pub const Color = Tuple;
pub const Point = Tuple.point;
pub const Vector = Tuple.vector;

/// Helper constructors & accessors for color semantics
pub inline fn color(red: Scalar, green: Scalar, blue: Scalar) Color {
    // store in x,y,z; keep w = 0 since it's a “vector-like” quantity
    return .{ .x = red, .y = green, .z = blue, .w = S(0) };
}

pub const Projectile = struct {
    position: Tuple,
    velocity: Tuple,

    pub fn init(position: Tuple, velocity: Tuple) Projectile {
        return .{ .position = position, .velocity = velocity };
    }
};

pub const Environment = struct {
    /// gravitation
    gravition: Tuple,
    /// wind
    wind: Tuple,

    pub fn init(gravitation: Tuple, wind: Tuple) Environment {
        return .{ .gravition = gravitation, .wind = wind };
    }
};

/// Utility function for printing a tuple.
pub fn format(t: Tuple, comptime fmt: []const u8, options: std.fmt.FormatOptions, writer: anytype) !void {
    _ = fmt;
    _ = options;
    try writer.print("Tuple({d:.6}, {d:.6}, {d:.6}, {d:.6})", .{ t.x, t.y, t.z, t.w });
}

/// ---
/// Intersection used by the various types of objects like spheres, cubes, etc...
/// ---
pub const Intersection = struct {
    t: Scalar,
    object_id: usize,

    pub fn eql(a: Intersection, b: Intersection) bool {
        return approxEq(a.t, b.t) and a.object_id == b.object_id;
    }

    pub fn init() Intersection {
        return .{
            .t = S(0),
            .object_id = 0,
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
    local_intersections_items: [MaxIntersectionsPerShape]Intersection,

    pub fn init() LocalIntersections {
        return .{ .count = 0, .local_intersections_items = undefined };
    }

    pub fn add(self: *LocalIntersections, i: Intersection) void {
        if (self.count < MaxIntersectionsPerShape) {
            self.local_intersections_items[self.count] = i;
            self.count += 1;
        }
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

    pub fn count(self: *const Intersections) usize {
        return self.intersections_items.items.len;
    }

    pub fn add(self: *Intersections, allocator: std.mem.Allocator, hit_to_add: Intersection) void {
        self.intersections_items.append(allocator, hit_to_add) catch @panic("OOM");
    }

    pub fn get(self: *const Intersections, index: usize) !Intersection {
        if (index >= self.intersections_items.items.len)
            return error.OutOfBounds;
        return self.intersections_items.items[index];
    }

    pub fn aggregate(allocator: std.mem.Allocator, ints: anytype) Intersections {
        var xs = Intersections.init();
        inline for (ints) |i| {
            xs.add(allocator, i);
        }
        return xs;
    }

    /// ---
    /// Return the Intersection item with the t > 0 when there
    /// has been a hit and the initialized version otherwise.
    /// i.e. no hit: t=0 and object_id=0.
    /// ---
    pub fn hit(self: *const Intersections) ?Intersection {
        var best: ?Intersection = null;

        for (self.intersections_items.items) |i| {
            if (i.t >= 0) { // includes t=0!
                if (best) |b| {
                    if (i.t < b.t) best = i;
                } else {
                    best = i;
                }
            }
        }
        return best;
        // var result = Intersection.init();
        //
        // for (self.intersections_items.items) |item| {
        //     if (item.t >= S(0) and result.object_id == 0) {
        //         result = item;
        //     } else if (item.t < result.t and item.t > 0) {
        //         result = item;
        //     }
        // }
        // return result;
    }
};

test "Chap1 -tuple initialization (Scalar)" {
    const eps: Scalar = 1e-6;

    const t = Tuple.init(4.3, -4.2, 3.1, 1.0);

    try std.testing.expectApproxEqAbs(@as(Scalar, 4.3), t.x, eps);
    try std.testing.expectApproxEqAbs(@as(Scalar, -4.2), t.y, eps);
    try std.testing.expectApproxEqAbs(@as(Scalar, 3.1), t.z, eps);
    try std.testing.expectApproxEqAbs(@as(Scalar, 1.0), t.w, eps);
}

test "Chap1 -point and vector constructors set w correctly" {
    const p = Tuple.point(1, 2, 3);
    const v = Tuple.vector(1, 2, 3);

    try std.testing.expectEqual(@as(Scalar, 1.0), p.w);
    try std.testing.expectEqual(@as(Scalar, 0.0), v.w);
    try std.testing.expect(!Tuple.equals(p, v));
}

test "Chap1 -point equality uses EPSILON" {
    const p1 = Tuple.point(1.0, 2.0, 3.0);
    const p2 = Tuple.point(1.0 + 1e-6, 2.0 - 5e-6, 3.0 + 2e-6);
    try std.testing.expect(Tuple.equals(p1, p2));
}

test "Chap1 -vector equality uses EPSILON" {
    const v1 = Tuple.vector(0.0, -1.0, 4.5);
    const v2 = Tuple.vector(0.0 + 9e-6, -1.0 - 2e-6, 4.5 + 1e-6);
    try std.testing.expect(Tuple.equals(v1, v2));
}

test "Chap1 -inequality when outside EPSILON" {
    const p1 = Tuple.point(1, 2, 3);
    const p2 = Tuple.point(1 + 1e-3, 2, 3); // too far on x for EPSILON=1e-5
    try std.testing.expect(!Tuple.equals(p1, p2));
}

test "Chap1 -adding two tuples" {
    const a1 = Tuple.init(3, -2, 5, 1);
    const a2 = Tuple.init(-2, 3, 1, 0);
    const a = Tuple.add(a1, a2);
    const e = Tuple.init(1, 1, 6, 1);
    try std.testing.expect(Tuple.equals(a, e));
}

test "Chap1 -subtracting two points" {
    const p1 = Tuple.point(3, 2, 1);
    const p2 = Tuple.point(5, 6, 7);
    const v = Tuple.sub(p1, p2);
    const e = Tuple.vector(-2, -4, -6);
    try std.testing.expect(Tuple.equals(v, e));
}

test "Chap1 -subtracting a vector from a point" {
    const p = Tuple.point(3, 2, 1);
    const v = Tuple.vector(5, 6, 7);
    const p_minus_v = Tuple.sub(p, v);
    const e = Tuple.point(-2, -4, -6);
    try std.testing.expect(Tuple.equals(p_minus_v, e));
}

test "Chap1 -subtractin two vectors" {
    const v1 = Tuple.vector(3, 2, 1);
    const v2 = Tuple.vector(5, 6, 7);
    const v = Tuple.sub(v1, v2);
    const e = Tuple.vector(-2, -4, -6);
    try std.testing.expect(Tuple.equals(v, e));
}

test "Chap1 -negate a tuple" {
    const a = Tuple.init(1, -2, 3, -4);
    const aneg = Tuple.neg(a);
    const e = Tuple.init(-1, 2, -3, 4);
    try std.testing.expect(Tuple.equals(aneg, e));
}

test "Chap1 -multiplying a tuple by scalar" {
    const a = Tuple.init(1, -2, 3, -4);
    const amult = Tuple.muls(a, 3.5);
    const e = Tuple.init(3.5, -7, 10.5, -14);
    try std.testing.expect(Tuple.equals(amult, e));
}

test "Chap1 -multiplying a tuple by a fraction" {
    const a = Tuple.init(1, -2, 3, -4);
    const amult = Tuple.muls(a, 0.5);
    const e = Tuple.init(0.5, -1, 1.5, -2);
    try std.testing.expect(Tuple.equals(amult, e));
}

test "Chap1 -dividing a tuple by a scalar" {
    const a = Tuple.init(1, -2, 3, -4);
    const amult = Tuple.div(a, 2);
    const e = Tuple.init(0.5, -1, 1.5, -2);
    try std.testing.expect(Tuple.equals(amult, e));
}

test "Chap1 -computing the magnitude of vector(1, 0, 0)" {
    const v = Tuple.vector(1, 0, 0);
    const mag = Tuple.mag(v);
    const e = 1;
    try std.testing.expect(approxEq(mag, e));
}

test "Chap1 -computing the magnitude of vector(0, 1, 0)" {
    const v = Tuple.vector(0, 1, 0);
    const mag = Tuple.mag(v);
    const e = 1;
    try std.testing.expect(approxEq(mag, e));
}

test "Chap1 -computing the magnitude of vector(0, 0, 1)" {
    const v = Tuple.vector(0, 0, 1);
    const mag = Tuple.mag(v);
    const e = 1;
    try std.testing.expect(approxEq(mag, e));
}

test "Chap1 -computing the magnitude of vector(1, 2, 3)" {
    const v = Tuple.vector(1, 2, 3);
    const mag = Tuple.mag(v);
    const e = std.math.sqrt(S(14));
    try std.testing.expect(approxEq(mag, e));
}

test "Chap1 -computing the magnitude of vector(-1, -2, -3)" {
    const v = Tuple.vector(-1, -2, -3);
    const mag = Tuple.mag(v);
    const e = std.math.sqrt(S(14));
    try std.testing.expect(approxEq(mag, e));
}

test "Chap1 -Normalizing vector(4, 0, 0)" {
    const v = Tuple.vector(4, 0, 0);
    const norm = Tuple.normalize(v);
    const e = Tuple.vector(1, 0, 0);
    try std.testing.expect(Tuple.equals(norm, e));
}

test "Chap1 -Normalizing vector(1, 2, 3)" {
    const v = Tuple.vector(1, 2, 3);
    const norm = Tuple.normalize(v);
    const e = Tuple.vector(S(1) / std.math.sqrt(S(14)), S(2) / std.math.sqrt(S(14)), S(3) / std.math.sqrt(S(14)));
    try std.testing.expect(Tuple.equals(norm, e));
}

test "Chap1 -The magnitude of a normalized vector" {
    const v = Tuple.vector(1, 2, 3);
    const norm = Tuple.normalize(v);
    const magnitude = Tuple.mag(norm);
    try std.testing.expect(approxEq(S(1), magnitude));
}

test "Chap1 -The dot product of two tuples" {
    const a = Tuple.vector(1, 2, 3);
    const b = Tuple.vector(2, 3, 4);
    const dot = Tuple.dot(a, b);
    try std.testing.expect(approxEq(dot, S(20)));
}

test "Chap1 -The cross product of two vectors" {
    const a = Tuple.vector(1, 2, 3);
    const b = Tuple.vector(2, 3, 4);
    const crossab = Tuple.cross(a, b);
    const crossba = Tuple.cross(b, a);
    try std.testing.expect(Tuple.equals(crossab, Tuple.vector(-1, 2, -1)));
    try std.testing.expect(Tuple.equals(crossba, Tuple.vector(1, -2, 1)));
}

pub fn tick(env: Environment, proj: Projectile) Projectile {
    const position = proj.position.add(proj.velocity);
    const velocity = proj.velocity.add(env.gravition.add(env.wind));

    // const position = Tuple.add(proj.p, proj.v);
    // const velocity = Tuple.add(proj.v, Tuple.add(env.g, env.w));

    const projectile = Projectile.init(position, velocity);
    return projectile;
}

test "Chap1 -Putting it together" {
    const r = error.SkipZigTest;
    if (r == error.SkipZigTest) return; // make not equal to for running this one

    // Projectile starts one unit above the origin.
    // Velocity is normalized to 1 unit/tick.
    var projectile = Projectile.init(Tuple.point(0, 1, 0), Tuple.normalize(Tuple.vector(1, 1, 0)));
    const e = Environment.init(Tuple.vector(0, -0.1, 0), Tuple.vector(-0.01, 0, 0));

    while (projectile.position.y > S(0)) {
        print("projectile.position: {f} projectile.velocity: {f}\n", .{ &projectile.position, &projectile.velocity });
        projectile = tick(e, projectile);
    }
}

test "Chap2 -Colors are (red, green, blue) Tuples" {
    const c = color(-0.5, 0.4, 1.7);
    try std.testing.expect(approxEq(c.r(), -0.5));
    try std.testing.expect(approxEq(c.g(), 0.4));
    try std.testing.expect(approxEq(c.b(), 1.7));
}

test "Chap2 -Adding colors" {
    const c1 = color(0.9, 0.6, 0.75);
    const c2 = color(0.7, 0.1, 0.25);
    const csum = c1.add(c2);
    const expect = color(c1.x + c2.x, c1.y + c2.y, c1.z + c2.z);
    try std.testing.expect(Tuple.equals(csum, expect));
}

test "Chap2 -Subtracting colors" {
    const c1 = color(0.9, 0.6, 0.75);
    const c2 = color(0.7, 0.1, 0.25);
    const csum = c1.sub(c2);
    const expect = color(c1.x - c2.x, c1.y - c2.y, c1.z - c2.z);
    try std.testing.expect(Tuple.equals(csum, expect));
}

test "Chap2 -Multiplying color by a scalar" {
    const c1 = color(0.2, 0.3, 0.4);
    const val: Scalar = S(2);
    const result = c1.muls(val);
    const expect = color(0.4, 0.6, 0.8);
    try std.testing.expect(Tuple.equals(result, expect));
}

test "Chap2 -Multiplying colors" {
    const c1 = color(1, 0.2, 0.4);
    const c2 = color(0.9, 1, 0.1);
    const result = c1.mult(c2);
    const expect = color(0.9, 0.2, 0.04);
    try std.testing.expect(Tuple.equals(result, expect));
}

test "Chap5 -Creating and querying a ray" {
    const origin = Tuple.point(1, 2, 3);
    const direction = Tuple.vector(4, 5, 6);
    const r = Ray.init(origin, direction);

    // log(@src(), "\norigin:{f}\ndirection:{f}\n", .{ origin, direction });
    // log(@src(), "\nray.origin:{f}\nray.direction:{f}\n", .{ ray.origin, ray.direction });

    try std.testing.expect(r.origin.equals(origin));
    try std.testing.expect(r.direction.equals(direction));
}

test "Chap5 -Computing a point from a distance" {
    const r = Ray.init(Tuple.point(2, 3, 4), Tuple.vector(1, 0, 0));
    const positionr0 = r.position(0);
    const positionrplus1 = r.position(1);
    const positionrminus1 = r.position(-1);
    const positionr25 = r.position(2.5);

    try std.testing.expect(positionr0.equals(Tuple.point(2, 3, 4)));
    try std.testing.expect(positionrplus1.equals(Tuple.point(3, 3, 4)));
    try std.testing.expect(positionrminus1.equals(Tuple.point(1, 3, 4)));
    try std.testing.expect(positionr25.equals(Tuple.point(4.5, 3, 4)));
}
