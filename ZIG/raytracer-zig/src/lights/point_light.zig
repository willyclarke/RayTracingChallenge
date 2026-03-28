const std = @import("std");
const print = @import("std").debug.print;
const types = @import("../types.zig");
const utils = @import("../utils.zig");
const matrix = @import("../matrix.zig");

const S = types.S;
const Ray = types.Ray;
const Scalar = types.Scalar;
const Tuple = types.Tuple;
const Matrix = matrix.Mat4;
const Intersection = types.Intersection;
const Intersections = types.Intersections;
const LocalIntersections = types.LocalIntersections;
const color = types.color;
const Color = types.Color;
const point = types.Point;
const vector = types.Vector;
const approxEq = types.approxEq;
const log = utils.log;

pub const PointLight = struct {
    pos: Tuple,
    intsty: Tuple,

    pub fn equals(self: *const PointLight, other: PointLight) bool {
        return self.pos.approxEq(other.pos) and self.intsty.approxEq(other.intsty);
    }

    pub fn init() PointLight {
        return .{ .pos = Tuple.point(0, 0, 0), .intsty = color(1, 1, 1) };
    }

    pub fn init_at(in_position: Tuple, in_intensity: Tuple) PointLight {
        return PointLight{ .pos = in_position, .intsty = in_intensity };
    }

    pub fn get(self: *const PointLight) PointLight {
        return .{ self.pos, self.intsty };
    }

    pub fn intensity(self: *const PointLight) Tuple {
        return self.intsty;
    }

    pub fn position(self: *const PointLight) Tuple {
        return self.pos;
    }
};
