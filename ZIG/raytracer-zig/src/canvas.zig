const std = @import("std");
const print = @import("std").debug.print;

const tuple = @import("tuple.zig");
const Color = tuple.Color;
const Projectile = tuple.Projectile;
const rgb = tuple.color;
const Scalar = tuple.Scalar;
const Tuple = tuple.Tuple;

pub const Canvas = struct {
    width: usize,
    height: usize,
    pixels: []Color,
    count: usize,

    /// Construct from components.
    pub fn init(alloc: std.mem.Allocator, width: usize, height: usize) !Canvas {
        const count = width * height;
        const pixels = try alloc.alloc(Color, count);
        @memset(pixels, rgb(0, 0, 0)); // mutate the slice contents
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
            const R = tuple.toByteSaturated(tuple.S(255) * canvas.pixelAt(x, y).r());
            const G = tuple.toByteSaturated(tuple.S(255) * canvas.pixelAt(x, y).g());
            const B = tuple.toByteSaturated(tuple.S(255) * canvas.pixelAt(x, y).b());
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
            try std.testing.expectEqual(c.pixels[y * c.width + x].r(), tuple.S(0));
            try std.testing.expectEqual(c.pixels[y * c.width + x].g(), tuple.S(0));
            try std.testing.expectEqual(c.pixels[y * c.width + x].b(), tuple.S(0));
        }
    }

    for (0..c.height) |y| {
        for (0..c.width) |x| {
            try std.testing.expectEqual(c.pixelAt(x, y).r(), tuple.S(0));
            try std.testing.expectEqual(c.pixelAt(x, y).g(), tuple.S(0));
            try std.testing.expectEqual(c.pixelAt(x, y).b(), tuple.S(0));
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
    const red = rgb(1, 0, 0);
    c.writePixel(2, 3, red);
    try std.testing.expect(tuple.Tuple.equals(&c.pixelAt(2, 3), &red));
}

test "Chap2 -Contructing the PPM header" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const alloc = gpa.allocator();

    var c = try Canvas.init(alloc, 500, 300);
    defer c.deinit(alloc);

    for (0..c.height) |y| {
        for (0..c.width) |x| {
            const color = rgb(tuple.S(y) / tuple.S(c.height), tuple.S(x) / tuple.S(c.width), tuple.S(y) / tuple.S(c.height));
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
    const e = tuple.Environment.init(gravity, wind);

    while (projectile.position.y > tuple.S(0)) {
        projectile = tuple.tick(e, projectile);

        const x = tuple.toUsizeSaturated(projectile.position.x, tuple.S(0), tuple.S(c.width));
        const y = c.height - tuple.toUsizeSaturated(projectile.position.y, tuple.S(0), tuple.S(c.height));
        // print("projectile.position: {f} projectile.velocity: {f} canvaspos x:{} y:{}\n", .{ &projectile.position, projectile.velocity, x, y });
        const color = rgb(tuple.S(y) / tuple.S(c.height), tuple.S(x) / tuple.S(c.width), tuple.S(y) / tuple.S(c.height));
        c.writePixel(x, y, color);
    }

    try createCanvasFile(&c, "projectile.ppm");
}
