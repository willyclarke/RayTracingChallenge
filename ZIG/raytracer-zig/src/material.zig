const std = @import("std");
const print = @import("std").debug.print;
const types = @import("types.zig");
const utils = @import("utils.zig");

const S = types.S;
const Ray = types.Ray;
const Scalar = types.Scalar;
const Tuple = types.Tuple;
const color = types.color;
const Color = types.Color;
const point = types.Point;
const vector = types.Vector;
const approxEq = types.approxEq;
const log = utils.log;

pub const Material = struct {
    col: Tuple,
    ambient: Scalar,
    diffuse: Scalar,
    specular: Scalar,
    shininess: Scalar,

    pub fn color(self: *const Material) Tuple {
        return self.col;
    }

    pub fn init() Material {
        return .{
            .col = types.color(1, 1, 1),
            .ambient = S(0.1),
            .diffuse = S(0.9),
            .specular = S(0.9),
            .shininess = S(200),
        };
    }

    pub fn equals(self: *const Material, other: Material) bool {
        const IsEqual = other.col.equals(self.col) //
            and approxEq(self.ambient, other.ambient) //
            and approxEq(self.diffuse, other.diffuse) //
            and approxEq(self.specular, other.specular) //
            and approxEq(self.shininess, other.shininess);
        return IsEqual;
    }
};

test "Chap6 -The default material" {
    const m = Material.init();

    try std.testing.expect(m.color().equals(color(1, 1, 1)));
    try std.testing.expect(approxEq(m.ambient, S(0.1)));
    try std.testing.expect(approxEq(m.diffuse, S(0.9)));
    try std.testing.expect(approxEq(m.specular, S(0.9)));
    try std.testing.expect(approxEq(m.shininess, S(200)));
}
