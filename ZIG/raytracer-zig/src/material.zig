const std = @import("std");
const print = @import("std").debug.print;
const types = @import("types.zig");
const utils = @import("utils.zig");
const Pattern = @import("patterns/pattern.zig").Pattern;

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
    /// Color
    col: Tuple,
    ambient: Scalar,
    diffuse: Scalar,
    specular: Scalar,
    shininess: Scalar,
    pattern: ?Pattern = null, // patterns are optional

    pub fn color(self: *const Material) Tuple {
        return self.col;
    }

    pub fn equals(self: *const Material, other: Material) bool {
        const IsEqual = other.col.equals(self.col) //
            and approxEq(self.ambient, other.ambient) //
            and approxEq(self.diffuse, other.diffuse) //
            and approxEq(self.specular, other.specular) //
            and approxEq(self.shininess, other.shininess);
        return IsEqual;
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

    /// Computes the "surface color" at a point, using pattern if present.
    // pub fn surfaceColor(self: *const Material, shape: *const Shape, world_point: Tuple) Color {
    //     if (self.pattern) |pat| {
    //         return pat.color_at_shape(shape, world_point);
    //     }
    //     return self.col;
    // }

    /// Fluent-ish helper: returns a modified copy
    pub fn withPattern(self: Material, p: Pattern) Material {
        var m = self;
        m.pattern = p;
        return m;
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
