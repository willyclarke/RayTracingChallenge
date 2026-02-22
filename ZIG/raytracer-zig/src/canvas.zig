const std = @import("std");
const print = @import("std").debug.print;
const utils = @import("utils.zig");
const log = utils.log;

const types = @import("types.zig");
const shapes = @import("shapes/shapes.zig");
const matrix = @import("matrix.zig");
const Mat4 = matrix.Mat4;

const material_mod = @import("material.zig");
const material = material_mod.Material;
const Material = material_mod.Material;

const light_mod = @import("lights.zig");
const point_light_mod = @import("lights/point_light.zig");
const light = light_mod.Light;
const Light = light_mod.Light;
const PointLight = point_light_mod.PointLight;
const lighting = light_mod.lighting;

const Color = types.Color;
const Projectile = types.Projectile;
const rgb = types.Color;
const Scalar = types.Scalar;
const S = types.S;
const Tuple = types.Tuple;
const Ray = types.Ray;
const approxEq = types.approxEq;
const point = types.Point;
const vector = types.Vector;

pub const Canvas = struct {
    width: usize,
    height: usize,
    pixels: []Color,

    pub fn init(alloc: std.mem.Allocator, width: usize, height: usize) !Canvas {
        const count = width * height;
        const pixels = try alloc.alloc(Color, count);
        const black = Color.init(0, 0, 0, 0);
        @memset(pixels, black);
        return .{ .width = width, .height = height, .pixels = pixels };
    }

    pub inline fn index(self: *const Canvas, x: usize, y: usize) usize {
        return y * self.width + x;
    }

    pub fn writePixel(self: *Canvas, x: usize, y: usize, c: Color) void {
        if (x >= self.width or y >= self.height) return;
        self.pixels[self.index(x, y)] = c;
    }

    pub fn pixelAt(self: *const Canvas, x: usize, y: usize) Color {
        if (x >= self.width or y >= self.height) return Color.init(0, 0, 0, 0);
        return self.pixels[self.index(x, y)];
    }

    pub fn pixelsSlice(self: *const Canvas) []const Color {
        return self.pixels;
    }

    pub fn deinit(self: *Canvas, alloc: std.mem.Allocator) void {
        alloc.free(self.pixels);
        self.* = undefined;
    }
};

// Common/safer: clamp to [0,255] then cast to u8
pub fn createCanvasFileP3(canvas: *const Canvas, filename: []const u8) !void {
    utils.log(@src(), "Starting\n", .{});
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
    utils.log(@src(), "Ending\n", .{});
}

pub fn createCanvasFile(canvas: *const Canvas, filename: []const u8) !void {
    utils.log(@src(), "Starting\n", .{});
    var file = try std.fs.cwd().createFile(filename, .{ .truncate = true });
    defer file.close();

    // Zig 0.15.x: the buffer is provided to the writer
    var buf: [64 * 1024]u8 = undefined; // tweak size if you like
    var fw = file.writer(&buf);
    const out = &fw.interface;

    try out.print("P6\n{} {}\n255\n", .{ canvas.width, canvas.height });

    for (0..canvas.height) |y| {
        for (0..canvas.width) |x| {
            const px = canvas.pixelAt(x, y);

            const r: u8 = types.toByteSaturated(types.S(255) * px.r());
            const g: u8 = types.toByteSaturated(types.S(255) * px.g());
            const b: u8 = types.toByteSaturated(types.S(255) * px.b());

            // Write 3 bytes (binary PPM)
            try out.writeAll(&.{ r, g, b });
        }
    }

    try out.flush(); // important
    utils.log(@src(), "Ending\n", .{});
}

pub fn createCanvasFile2(
    alloc: std.mem.Allocator,
    canvas: *const Canvas,
    filename: []const u8,
) !void {
    utils.log(@src(), "Starting\n", .{});

    var file = try std.fs.cwd().createFile(filename, .{ .truncate = true });
    defer file.close();

    // Zig 0.15.x: buffered output is done by providing a buffer to the writer.
    var out_buf: [1024 * 1024]u8 = undefined;
    var fw = file.writer(&out_buf);
    const out = &fw.interface;

    // P6 header (binary PPM)
    try out.print("P6\n{} {}\n255\n", .{ canvas.width, canvas.height });

    // Row buffer: width * 3 bytes (RGB)
    var row = try alloc.alloc(u8, canvas.width * 3);
    defer alloc.free(row);

    const pixels = canvas.pixelsSlice();

    // Safety: make sure the slice matches width*height
    // (optional in ReleaseFast; great for Debug)
    std.debug.assert(pixels.len == canvas.width * canvas.height);

    for (0..canvas.height) |y| {
        const row_start = y * canvas.width;
        const row_px = pixels[row_start .. row_start + canvas.width];

        var i: usize = 0;
        for (row_px) |px| {
            row[i + 0] = types.toByteSaturated(types.S(255) * px.r());
            row[i + 1] = types.toByteSaturated(types.S(255) * px.g());
            row[i + 2] = types.toByteSaturated(types.S(255) * px.b());
            i += 3;
        }

        try out.writeAll(row);
    }

    try out.flush();

    utils.log(@src(), "Ending\n", .{});
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

fn drawSquare(
    ptrCanvas: *Canvas,
    center_x: usize,
    center_y: usize,
    half_size: usize,
    ptrcolor: *const types.Color,
) void {
    var y: usize = center_y - half_size;
    while (y <= center_y + half_size) : (y += 1) {
        var x: usize = center_x - half_size;
        while (x <= center_x + half_size) : (x += 1) {
            ptrCanvas.writePixel(x, y, ptrcolor.*);
        }
    }
}

test "Chap4 -Putting It Together" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const alloc = gpa.allocator();

    var c = try Canvas.init(alloc, 900, 900);
    defer c.deinit(alloc);

    const origin = point(0, 0, 0);
    var p = point(0, 0, 0);
    try std.testing.expect(p.equals(origin));

    var color = types.color(1, 0, 0);

    const canvastranslate = Mat4.translation(S(c.width) / S(2), S(c.height) / S(2), 0);

    // Use the translation number x a percentage t set scaling.
    const scale = Mat4.scaling(canvastranslate.get(0, 3) * S(0.75), canvastranslate.get(1, 3) * S(0.75), 0);

    var canvpoint = canvastranslate.mulM(&scale).mulT(p);
    var centerx = types.toUsizeSaturated(canvpoint.x, S(0), S(c.width));
    var centery = types.toUsizeSaturated(canvpoint.y, S(0), S(c.height));
    const squaresize = types.toUsizeSaturated(S(c.width) / S(20), S(3), S(10));
    drawSquare(&c, centerx, centery, squaresize, &color);

    color.z = 1;
    p.x = 1;

    for (0..12) |idx| {
        const alfa = S(idx) / S(12) * std.math.tau;
        color.y = S(1) - color.z;
        color.z = S(idx) / S(12);

        canvpoint = canvastranslate.mulM(&scale.mulM(&Mat4.rotz(alfa))).mulT(p);
        centerx = types.toUsizeSaturated(canvpoint.x, S(0), S(c.width));
        centery = types.toUsizeSaturated(canvpoint.y, S(0), S(c.height));
        drawSquare(&c, centerx, centery, squaresize, &color);
    }

    try createCanvasFile(&c, "chap4clock.ppm");
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
    var sph = shapes.Sphere.init();
    var s = shapes.Shape.fromSphere(&sph);

    const sxform = matrix.Mat4.shearing(1, 0, 0, 0, 0, 0).mulM(&matrix.Mat4.scaling(0.5, 1, 1));
    s.setTransform(&sxform);

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

test "Chap6 -Putting it together" {
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
    // const color = types.color(1, 0, 0);
    var sph = shapes.Sphere.init();
    var s = shapes.Shape.fromSphere(&sph);
    s.sphere.h.material.col = types.color(1, 0.2, 1);

    const light_position = point(-10, 10, -10);
    const light_color = types.color(1, 1, 1);
    const pl = Light.fromPointLight(PointLight.init_at(light_position, light_color));
    // s.set_transform(&matrix.Mat4.shearing(1, 0, 0, 0, 0, 0).mulM(&matrix.Mat4.scaling(0.5, 1, 1)));

    for (0..canvas_pixels) |y| {
        const world_y = half - pixel_size * S(y);
        for (0..canvas_pixels) |x| {
            const world_x = -half + pixel_size * S(x);
            const position = types.Point(world_x, world_y, wall_z);
            const r = Ray.init(ray_origin, position.sub(ray_origin).normalize());
            const xs = s.intersect(r);

            if (xs.hit()) {
                const hitpoint = r.position(xs.tmin());
                const normal = s.normal_at(hitpoint);
                const eye = r.direction.muls(S(-1));
                const color = lighting(s.material().*, pl, hitpoint, eye, normal);
                canvas.writePixel(x, y, color);
            }
        }
    }

    try createCanvasFile(&canvas, "chap6puttingtogether.ppm");
}
