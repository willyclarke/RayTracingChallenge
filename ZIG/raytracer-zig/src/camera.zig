const std = @import("std");
const print = @import("std").debug.print;
const utils = @import("utils.zig");
const log = utils.log;

const types_mod = @import("types.zig");
const matrix_mod = @import("matrix.zig");
const Scalar = types_mod.Scalar;
const S = types_mod.S;
const Tuple = types_mod.Tuple;
const approxEq = types_mod.approxEq;
const point = types_mod.Point;
const vector = types_mod.Vector;
const color = types_mod.Vector;
const matrix = matrix_mod.Mat4;
const Matrix = matrix_mod.Mat4;
const Ray = types_mod.Ray;

pub const Camera = struct {
    hsize: usize,
    vsize: usize,
    field_of_view: Scalar,
    half_width: Scalar,
    half_height: Scalar,
    pixel_size: Scalar,
    transform: Matrix,

    pub fn init(hsize: usize, vsize: usize, field_of_view: Scalar) Camera {
        var half_width: Scalar = undefined;
        var half_height: Scalar = undefined;

        const half_view = std.math.tan(field_of_view / S(2));
        const aspect = @as(Scalar, @floatFromInt(hsize)) / @as(Scalar, @floatFromInt(vsize));

        if (aspect >= S(1)) {
            half_width = half_view;
            half_height = half_view / aspect;
        } else {
            half_width = half_view * aspect;
            half_height = half_view;
        }

        const pixel_size = (half_width * S(2)) / @as(Scalar, @floatFromInt(hsize));

        return .{
            .hsize = hsize,
            .vsize = vsize,
            .field_of_view = field_of_view,
            .half_width = half_width,
            .half_height = half_height,
            .pixel_size = pixel_size,
            .transform = Matrix.identity(),
        };
    }

    pub fn equals(self: *const Camera, other: Camera) bool {
        const IsEqual =
            approxEq(self.hsize, other.hsize) and
            approxEq(self.vsize, other.vsize) and
            approxEq(self.field_of_view, other.field_of_view) and
            other.transform.equals(&self.transform);
        return IsEqual;
    }
};

pub fn ray_for_pixel(camera: *const Camera, px: usize, py: usize) Ray {

    // the offset from the edge of the canvas to the pixel's center
    const xoffset = (@as(Scalar, @floatFromInt(px)) + S(0.5)) * camera.pixel_size;
    const yoffset = (@as(Scalar, @floatFromInt(py)) + S(0.5)) * camera.pixel_size;

    // the untransformed coordinates of the pixel in world space.
    // (remember that the camera looks toward -z, so +x is to the *left*.)
    const world_x = camera.half_width - xoffset;
    const world_y = camera.half_height - yoffset;

    // using the camera matrix, transform the canvas point and the origin,
    // and then compute the ray's direction vector.
    // (remember that the canvas is at z=-1)
    const pixel = camera.transform.inverse().mulT(point(world_x, world_y, S(-1)));
    const origin = camera.transform.inverse().mulT(point(0, 0, 0));
    const direction = pixel.sub(origin).normalize();
    return Ray.init(origin, direction);
}

test "Chap7 -Make sure it works" {
    try std.testing.expect(7 == 7);
}

test "Chap7 -Constructing a camera" {
    const hsize: usize = 160;
    const vsize: usize = 120;
    const field_of_view = std.math.pi / S(2);
    const c = Camera.init(hsize, vsize, field_of_view);
    try std.testing.expect(hsize == c.hsize);
    try std.testing.expect(vsize == c.vsize);
    try std.testing.expect(approxEq(field_of_view, c.field_of_view));
    const i = matrix.identity();
    try std.testing.expect(c.transform.equals(&i));
}

test "Chap7 -The pixel size for a horizontal canvas" {
    const c = Camera.init(200, 125, std.math.pi / S(2));
    try std.testing.expect(approxEq(c.pixel_size, 0.01));
}

test "Chap7 -The pixel size for a vertical canvas" {
    const c = Camera.init(125, 200, std.math.pi / S(2));
    try std.testing.expect(approxEq(c.pixel_size, 0.01));
}

test "Chap7 -Constructing a ray through the center of the canvas" {
    const c = Camera.init(201, 101, std.math.pi / S(2));
    const r = ray_for_pixel(&c, 100, 50);
    try std.testing.expect(r.origin.equals(point(0, 0, 0)));
    try std.testing.expect(r.direction.equals(vector(0, 0, -1)));
}

test "Chap7 -Constructing a ray through a corner of the canvas" {
    const c = Camera.init(201, 101, std.math.pi / S(2));
    const r = ray_for_pixel(&c, 0, 0);
    try std.testing.expect(r.origin.equals(point(0, 0, 0)));
    try std.testing.expect(r.direction.equals(vector(0.66519, 0.33259, -0.66851)));
}

test "Chap7 -Constructing a ray when the camera is transformed" {
    var c = Camera.init(201, 101, std.math.pi / S(2));
    const translation = Matrix.translation(0, -2, 5);
    c.transform = Matrix.roty(std.math.pi / S(4)).mulM(&translation);
    const r = ray_for_pixel(&c, 100, 50);
    try std.testing.expect(r.origin.equals(point(0, 2, -5)));
    try std.testing.expect(r.direction.equals(vector(std.math.sqrt2 / S(2), 0, -std.math.sqrt2 / S(2))));
}
