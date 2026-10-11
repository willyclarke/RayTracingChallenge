const std = @import("std");
const print = @import("std").debug.print;
const utils = @import("../utils.zig");
const tMod = @import("../types.zig");
const matrix = @import("../matrix.zig");
pub const Mat4 = matrix.Mat4;
const S = tMod.S;

const Color = tMod.Color;
const Tuple = tMod.Tuple;
const color = tMod.Vector;
const point = tMod.Point;
const vector = tMod.Vector;
const approxEq = tMod.approxEq;
const log = utils.log;

const StripePattern = @import("stripe_pattern.zig").StripePattern;
const PatternHeader = @import("pattern_header.zig").PatternHeader;

pub const Pattern = union(enum) {
    stripe: StripePattern,
    // gradient: GradientPattern,
    // ring: RingPattern,
    // checker: CheckerPattern,

    pub fn fromStripe(ptrStripePattern: *const StripePattern) Pattern {
        return .{ .stripe = ptrStripePattern.* };
    }

    pub fn init_stripe(a: Color, b: Color) Pattern {
        return .{ .stripe = StripePattern.init(a, b) };
    }

    pub fn pattern_at(self: *const Pattern, p: Tuple) Color {
        return switch (self.*) {
            .stripe => |*s| s.local_color_at(p),
        };
    }

    pub fn setTransform(self: *Pattern, m: Mat4) void {
        switch (self.*) {
            inline else => |*pat| pat.h.setTransform(m),
        }
    }

    pub fn header(self: *const Pattern) *const PatternHeader {
        return switch (self.*) {
            // .stripe => |*s| s.header(),
            inline else => |pat| &pat.h,
        };
    }

    pub fn local_color_at(self: *const Pattern, p: Tuple) Color {
        return switch (self.*) {
            inline else => |pat| pat.local_color_at(p),
        };
    }
};

test "Chap10 -Ensure it works" {
    try std.testing.expect(7 == 7);
}

test "Chap10 -Creating a stripe pattern" {
    const pattern = StripePattern.init(tMod.white(), tMod.black());
    try std.testing.expect(pattern.a.equals(tMod.white()));
    try std.testing.expect(pattern.b.equals(tMod.black()));
}

test "Chap10 -A stripe pattern is constant in y" {
    const pattern = StripePattern.init(tMod.white(), tMod.black());
    try std.testing.expect(pattern.local_color_at(point(0, 0, 0)).equals(tMod.white()));
    try std.testing.expect(pattern.local_color_at(point(0, 1, 0)).equals(tMod.white()));
    try std.testing.expect(pattern.local_color_at(point(0, 2, 0)).equals(tMod.white()));
    try std.testing.expect(pattern.b.equals(tMod.black()));
}

test "Chap10 -A stripe pattern is constant in z" {
    const pattern = StripePattern.init(tMod.white(), tMod.black());
    try std.testing.expect(pattern.local_color_at(point(0, 0, 0)).equals(tMod.white()));
    try std.testing.expect(pattern.local_color_at(point(0, 0, 1)).equals(tMod.white()));
    try std.testing.expect(pattern.local_color_at(point(0, 0, 2)).equals(tMod.white()));
    try std.testing.expect(pattern.b.equals(tMod.black()));
}

test "Chap10 -A stripe pattern alternates in x" {
    const pattern = StripePattern.init(tMod.white(), tMod.black());
    try std.testing.expect(pattern.local_color_at(point(0, 0, 0)).equals(tMod.white()));
    try std.testing.expect(pattern.local_color_at(point(0.9, 0, 0)).equals(tMod.white()));
    try std.testing.expect(pattern.local_color_at(point(1, 0, 0)).equals(tMod.black()));
    try std.testing.expect(pattern.local_color_at(point(-0.1, 0, 0)).equals(tMod.black()));
    try std.testing.expect(pattern.local_color_at(point(-1.0, 0, 0)).equals(tMod.black()));
    try std.testing.expect(pattern.local_color_at(point(-1.1, 0, 0)).equals(tMod.white()));
    try std.testing.expect(pattern.b.equals(tMod.black()));
}
