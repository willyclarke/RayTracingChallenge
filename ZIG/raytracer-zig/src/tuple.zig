const std = @import("std");
const print = @import("std").debug.print;

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

/// Alias to convert anytype to the Scalar type
pub inline fn S(x: anytype) Scalar {
    return scalar(x);
}

pub const EPSILON: Scalar = 1e-5;

/// Returns true when |a - b| < EPSILON
pub inline fn almostEqual(a: Scalar, b: Scalar) bool {
    return @abs(a - b) < EPSILON;
}

test "Chap1 -almostEqual works for close values" {
    try std.testing.expect(almostEqual(1.000001, 1.000002));
    try std.testing.expect(!almostEqual(1.0, 1.1));
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

    pub fn equals(a: Tuple, b: Tuple) bool {
        return almostEqual(a.x, b.x) and
            almostEqual(a.y, b.y) and
            almostEqual(a.z, b.z) and
            almostEqual(a.w, b.w);
    }

    pub fn add(a: Tuple, b: Tuple) Tuple {
        return .{ .x = a.x + b.x, .y = a.y + b.y, .z = a.z + b.z, .w = a.w + b.w };
    }

    pub fn sub(a: Tuple, b: Tuple) Tuple {
        return .{ .x = a.x - b.x, .y = a.y - b.y, .z = a.z - b.z, .w = a.w - b.w };
    }

    pub fn neg(a: Tuple) Tuple {
        return .{ .x = -a.x, .y = -a.y, .z = -a.z, .w = -a.w };
    }

    pub fn mul(a: Tuple, s: Scalar) Tuple {
        return .{ .x = a.x * s, .y = a.y * s, .z = a.z * s, .w = a.w * s };
    }

    pub fn div(a: Tuple, s: Scalar) Tuple {
        return .{ .x = a.x / s, .y = a.y / s, .z = a.z / s, .w = a.w / s };
    }

    pub fn mag(a: Tuple) Scalar {
        return std.math.sqrt(a.x * a.x + a.y * a.y + a.z * a.z + a.w * a.w);
    }

    pub fn normalize(a: Tuple) Tuple {
        return a.div(a.mag());
    }

    pub fn dot(a: Tuple, b: Tuple) Scalar {
        return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
    }

    pub fn cross(a: Tuple, b: Tuple) Tuple {
        return vector(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
    }

    /// Custom formatter so `{}` prints nicely.
    /// `fmt` and `options` let you add variants later; for now we ignore them.
    pub fn format(self: Tuple, writer: anytype) !void {
        if (self.w != S(0)) {
            try writer.print("P({d:8.5}, {d:8.5}, {d:8.5}, {d:1.0})", .{ self.x, self.y, self.z, self.w });
        } else {
            try writer.print("V({d:8.5}, {d:8.5}, {d:8.5}, {d:1.0})", .{ self.x, self.y, self.z, self.w });
        }
    }

    /// Handle printing via a *const Tuple pointer
    pub fn formatPtr(self: *const Tuple, writer: anytype) !void {
        // Forward to the value formatter
        return Tuple.format(self.*, writer);
    }
};

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
    const amult = Tuple.mul(a, 3.5);
    const e = Tuple.init(3.5, -7, 10.5, -14);
    try std.testing.expect(Tuple.equals(amult, e));
}

test "Chap1 -multiplying a tuple by a fraction" {
    const a = Tuple.init(1, -2, 3, -4);
    const amult = Tuple.mul(a, 0.5);
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
    try std.testing.expect(almostEqual(mag, e));
}

test "Chap1 -computing the magnitude of vector(0, 1, 0)" {
    const v = Tuple.vector(0, 1, 0);
    const mag = Tuple.mag(v);
    const e = 1;
    try std.testing.expect(almostEqual(mag, e));
}

test "Chap1 -computing the magnitude of vector(0, 0, 1)" {
    const v = Tuple.vector(0, 0, 1);
    const mag = Tuple.mag(v);
    const e = 1;
    try std.testing.expect(almostEqual(mag, e));
}

test "Chap1 -computing the magnitude of vector(1, 2, 3)" {
    const v = Tuple.vector(1, 2, 3);
    const mag = Tuple.mag(v);
    const e = std.math.sqrt(S(14));
    try std.testing.expect(almostEqual(mag, e));
}

test "Chap1 -computing the magnitude of vector(-1, -2, -3)" {
    const v = Tuple.vector(-1, -2, -3);
    const mag = Tuple.mag(v);
    const e = std.math.sqrt(S(14));
    try std.testing.expect(almostEqual(mag, e));
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
    try std.testing.expect(almostEqual(S(1), magnitude));
}

test "Chap1 -The dot product of two tuples" {
    const a = Tuple.vector(1, 2, 3);
    const b = Tuple.vector(2, 3, 4);
    const dot = Tuple.dot(a, b);
    try std.testing.expect(almostEqual(dot, S(20)));
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
        print("projectile.position: {f} projectile.velocity: {f}\n", .{ &projectile.position, projectile.velocity });
        projectile = tick(e, projectile);
    }
}
