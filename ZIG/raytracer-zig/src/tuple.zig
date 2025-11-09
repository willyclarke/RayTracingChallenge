const std = @import("std");

pub const Scalar = f32; // switch to f64 later if you need it
pub const EPSILON: Scalar = 1e-5;

/// Returns true when |a - b| < EPSILON
pub inline fn almostEqual(a: Scalar, b: Scalar) bool {
    return @abs(a - b) < EPSILON;
}

test "almostEqual works for close values" {
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
};

test "tuple initialization (Scalar)" {
    const eps: Scalar = 1e-6;

    const t = Tuple.init(4.3, -4.2, 3.1, 1.0);

    try std.testing.expectApproxEqAbs(@as(Scalar, 4.3), t.x, eps);
    try std.testing.expectApproxEqAbs(@as(Scalar, -4.2), t.y, eps);
    try std.testing.expectApproxEqAbs(@as(Scalar, 3.1), t.z, eps);
    try std.testing.expectApproxEqAbs(@as(Scalar, 1.0), t.w, eps);
}

test "point and vector constructors set w correctly" {
    const p = Tuple.point(1, 2, 3);
    const v = Tuple.vector(1, 2, 3);

    try std.testing.expectEqual(@as(Scalar, 1.0), p.w);
    try std.testing.expectEqual(@as(Scalar, 0.0), v.w);
    try std.testing.expect(!Tuple.equals(p, v));
}

test "point equality uses EPSILON" {
    const p1 = Tuple.point(1.0, 2.0, 3.0);
    const p2 = Tuple.point(1.0 + 1e-6, 2.0 - 5e-6, 3.0 + 2e-6);
    try std.testing.expect(Tuple.equals(p1, p2));
}

test "vector equality uses EPSILON" {
    const v1 = Tuple.vector(0.0, -1.0, 4.5);
    const v2 = Tuple.vector(0.0 + 9e-6, -1.0 - 2e-6, 4.5 + 1e-6);
    try std.testing.expect(Tuple.equals(v1, v2));
}

test "inequality when outside EPSILON" {
    const p1 = Tuple.point(1, 2, 3);
    const p2 = Tuple.point(1 + 1e-3, 2, 3); // too far on x for EPSILON=1e-5
    try std.testing.expect(!Tuple.equals(p1, p2));
}
