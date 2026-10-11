const std = @import("std");
const print = @import("std").debug.print;
const utils = @import("../utils.zig");
const tMod = @import("../types.zig");
const matrix = @import("../matrix.zig");
pub const Mat4 = matrix.Mat4;

const S = tMod.S;
const point = tMod.Point;
const vector = tMod.Vector;
const approxEq = tMod.approxEq;
const log = utils.log;

pub const PatternHeader = struct {
    m: Mat4 = Mat4.identity(),
    inv: Mat4 = Mat4.identity(),

    pub fn setTransform(self: *PatternHeader, m: Mat4) void {
        self.m = m;
        self.inv = m.inverse();
    }
};
