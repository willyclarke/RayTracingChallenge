const std = @import("std");
const print = @import("std").debug.print;
const utils = @import("../utils.zig");
const tMod = @import("../types.zig");
const matrix = @import("../matrix.zig");
pub const Mat4 = matrix.Mat4;
const approxEq = tMod.approxEq;
const S = tMod.S;

const Color = tMod.Color;
const Tuple = tMod.Tuple;
const point = tMod.Point;
const vector = tMod.Vector;
const log = utils.log;

const PatternHeader = @import("pattern_header.zig").PatternHeader;

pub const StripePattern = struct {
    h: PatternHeader = .{},
    a: Color,
    b: Color,

    pub fn init(a: Color, b: Color) StripePattern {
        return .{ .a = a, .b = b };
    }

    pub fn local_color_at(self: *const StripePattern, p: Tuple) Color {
        return if (@mod(@as(i64, @intFromFloat(std.math.floor(p.x))), 2) == 0)
            self.a
        else
            self.b;
    }

    pub fn header(self: *const StripePattern) *const PatternHeader {
        return &self.h;
    }
};
