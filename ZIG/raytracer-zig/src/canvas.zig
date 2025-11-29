const std = @import("std");
const print = @import("std").debug.print;
const utils = @import("utils.zig");
const log = utils.log;

const types = @import("types.zig");
const shapes = @import("shapes.zig");
const matrix = @import("matrix.zig");
const Color = types.Color;
const Projectile = types.Projectile;
const rgb = types.Color;
const Scalar = types.Scalar;
const S = types.S;
const Tuple = types.Tuple;
const Ray = types.Ray;
const approxEq = types.approxEq;

pub const Canvas = struct {
    width: usize,
    height: usize,
    pixels: []Color,
    count: usize,

    /// Construct from components.
    pub fn init(alloc: std.mem.Allocator, width: usize, height: usize) !Canvas {
        const count = width * height;
        const pixels = try alloc.alloc(Color, count);
        const black = rgb.init(0, 0, 0, 0);
        @memset(pixels, black); // mutate the slice contents
        return .{ .width = width, .height = height, .pixels = pixels, .count = count };
    }

    pub inline fn index(self: Canvas, x: usize, y: usize) usize {
        return y * self.width + x;
    }

    pub fn writePixel(self: *Canvas, x: usize, y: usize, c: Color) void {
        const idx = self.index(x, y);
        if (idx < self.count)
            self.pixels[idx] = c;
    }

    pub fn pixelAt(self: Canvas, x: usize, y: usize) Color {
        const idx = self.index(x, y);
        if (idx < self.count)
            return self.pixels[idx];
        return Color.init(0.0, 0.0, 0.0, 0.0);
    }

    pub fn deinit(self: *Canvas, alloc: std.mem.Allocator) void {
        alloc.free(self.pixels);
        self.* = undefined;
    }
};

// Common/safer: clamp to [0,255] then cast to u8
pub fn createCanvasFile(canvas: *const Canvas, filename: []const u8) !void {
    const a = std.heap.page_allocator;

    var buffer: std.ArrayList(u8) = .empty;
    try buffer.writer(a).print("P3\n{} {}\n255\n", .{ canvas.width, canvas.height });

    // ---
    // NOTE: Get three chars per byte and there are three
    //       of them plus the three spaces.
    // ---
    const Increment = 3 * 3 + 3;
    const MaxCharsPerLine = 50;
    var NumCharsOnLine: u32 = 0;
    for (0..canvas.height) |y| {
        for (0..canvas.width) |x| {
            const R = types.toByteSaturated(types.S(255) * canvas.pixelAt(x, y).r());
            const G = types.toByteSaturated(types.S(255) * canvas.pixelAt(x, y).g());
            const B = types.toByteSaturated(types.S(255) * canvas.pixelAt(x, y).b());
            try buffer.writer(a).print("{:03} {:03} {:03} ", .{ R, G, B });

            NumCharsOnLine = NumCharsOnLine + Increment;
            if (MaxCharsPerLine < NumCharsOnLine) {
                try buffer.writer(a).print("\n", .{});
                NumCharsOnLine = 0;
            }
        }
    }

    // ---
    // NOTE: The file need to end with a new line.
    // ---
    try buffer.writer(a).print("\n", .{});

    const file = try std.fs.cwd().createFile(filename, .{ .truncate = true });
    defer file.close();
    try file.writeAll(buffer.items);
}

test "Chap2 -Creating a canvas" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const alloc = gpa.allocator();

    var c = try Canvas.init(alloc, 10, 20);
    defer c.deinit(alloc);

    try std.testing.expectEqual(@as(usize, 10), c.width);
    try std.testing.expectEqual(@as(usize, 20), c.height);

    // default to black
    for (0..c.height) |y| {
        for (0..c.width) |x| {
            try std.testing.expectEqual(c.pixels[y * c.width + x].r(), types.S(0));
            try std.testing.expectEqual(c.pixels[y * c.width + x].g(), types.S(0));
            try std.testing.expectEqual(c.pixels[y * c.width + x].b(), types.S(0));
        }
    }

    for (0..c.height) |y| {
        for (0..c.width) |x| {
            try std.testing.expectEqual(c.pixelAt(x, y).r(), types.S(0));
            try std.testing.expectEqual(c.pixelAt(x, y).g(), types.S(0));
            try std.testing.expectEqual(c.pixelAt(x, y).b(), types.S(0));
        }
    }
}

test "Chap2 -Writing pixels to a canvas" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const alloc = gpa.allocator();

    var c = try Canvas.init(alloc, 10, 20);
    defer c.deinit(alloc);

    // write and read
    const red = rgb.init(1, 0, 0, 0);
    c.writePixel(2, 3, red);
    try std.testing.expect(types.Tuple.equals(c.pixelAt(2, 3), red));
}

test "Chap2 -Contructing the PPM header" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const alloc = gpa.allocator();

    var c = try Canvas.init(alloc, 500, 300);
    defer c.deinit(alloc);

    for (0..c.height) |y| {
        for (0..c.width) |x| {
            const color = rgb.init(types.S(y) / types.S(c.height), types.S(x) / types.S(c.width), types.S(y) / types.S(c.height), 0);
            c.writePixel(x, y, color);
        }
    }

    try createCanvasFile(&c, "image.ppm");
}

test "Chap2 -Putting it together" {
    const r = error.SkipZigTest;
    if (r == error.SkipZigTest) return; // make not equal to for running this one

    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const alloc = gpa.allocator();

    var c = try Canvas.init(alloc, 900, 550);
    defer c.deinit(alloc);

    // Projectile starts one unit above the origin.
    // Velocity is normalized to 1 unit/tick.
    const start = Tuple.point(0, 1, 0);
    const velocity = Tuple.muls(Tuple.normalize(Tuple.vector(1, 1.8, 0)), 11.25);
    var projectile = Projectile.init(start, velocity);

    const gravity = Tuple.vector(0, -0.1, 0);
    const wind = Tuple.vector(-0.01, 0, 0);
    const e = types.Environment.init(gravity, wind);

    while (projectile.position.y > types.S(0)) {
        projectile = types.tick(e, projectile);

        const x = types.toUsizeSaturated(projectile.position.x, types.S(0), types.S(c.width));
        const y = c.height - types.toUsizeSaturated(projectile.position.y, types.S(0), types.S(c.height));
        // print("projectile.position: {f} projectile.velocity: {f} canvaspos x:{} y:{}\n", .{ &projectile.position, projectile.velocity, x, y });
        const color = rgb(types.S(y) / types.S(c.height), types.S(x) / types.S(c.width), types.S(y) / types.S(c.height));
        c.writePixel(x, y, color);
    }

    try createCanvasFile(&c, "projectile.ppm");
}

test "Chap5 -Putting it together" {
    try std.testing.expect(3 == 3);

    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const alloc = gpa.allocator();

    const canvas_pixels: usize = 100;
    var canvas = try Canvas.init(alloc, canvas_pixels, canvas_pixels);
    defer canvas.deinit(alloc);

    const ray_origin = types.Point(0, 0, -5);
    const wall_z = S(10);
    const wall_size = S(7);
    const pixel_size = S(wall_size) / S(canvas_pixels);
    const half = wall_size / S(2);
    const color = types.color(1, 0, 0);
    var s = shapes.Shape.fromSphere(shapes.Sphere.init());

    // const mxform = matrix.Mat4.rotz(std.math.pi / S(4)).mulM(&matrix.Mat4.scaling(0.5, 1, 1));
    // s.set_transform(&mxform);

    s.set_transform(&matrix.Mat4.shearing(1, 0, 0, 0, 0, 0).mulM(&matrix.Mat4.scaling(0.5, 1, 1)));

    for (0..canvas_pixels) |y| {
        const world_y = half - pixel_size * S(y);
        for (0..canvas_pixels) |x| {
            const world_x = -half + pixel_size * S(x);
            const position = types.Point(world_x, world_y, wall_z);
            const r = types.Ray.init(ray_origin, position.sub(ray_origin).normalize());
            const xs = s.intersect(r);
            if (xs.hit()) {
                canvas.writePixel(x, y, color);
            }
        }
    }

    try std.testing.expect(approxEq(ray_origin.z, S(-5)));
    try std.testing.expect(pixel_size > S(0));

    try createCanvasFile(&canvas, "chap5puttingtogether.ppm");
}
