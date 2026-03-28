const std = @import("std");
const print = @import("std").debug.print;

const utils = @import("utils.zig");
const log = utils.log;

const camera_mod = @import("camera.zig");
const Camera = camera_mod.Camera;
const ray_for_pixel = camera_mod.ray_for_pixel;

const light_mod = @import("lights.zig");
const Light = light_mod.Light;
const PointLight = light_mod.PointLight;

const world_mod = @import("world.zig");
const World = world_mod.World;
const default_world = world_mod.default_world;
const view_transform = world_mod.view_transform;
const color_at = world_mod.color_at;

const canvas_mod = @import("canvas.zig");
const Canvas = canvas_mod.Canvas;
const createCanvasFile = canvas_mod.createCanvasFile;
const createCanvasFile2 = canvas_mod.createCanvasFile2;

const shapesMod = @import("shapes/shapes.zig");

const material = @import("material.zig").Material;
const Material = @import("material.zig").Material;
const Pattern = @import("patterns/pattern.zig").Pattern;
const StripePattern = @import("patterns/stripe_pattern.zig").StripePattern;

const tMod = @import("types.zig");
const mMod = @import("matrix.zig");

const Scalar = tMod.Scalar;
const S = tMod.S;
const Tuple = tMod.Tuple;
const approxEq = tMod.approxEq;
const point = tMod.Point;
const vector = tMod.Vector;
const color = tMod.Vector;
const Ray = tMod.Ray;
const matrix = mMod.Mat4;
const Matrix = mMod.Mat4;

pub fn renderSingleThread(alloc: std.mem.Allocator, camera: *const Camera, world: *const World) !Canvas {
    var image = try Canvas.init(alloc, camera.hsize, camera.vsize);
    errdefer image.deinit(alloc);

    // scratch arena for per-ray allocations (intersections, etc.)
    var scratch = std.heap.ArenaAllocator.init(alloc);
    defer scratch.deinit();
    const temp_alloc = scratch.allocator();

    var y: usize = 0;
    while (y < camera.vsize) : (y += 1) {
        _ = scratch.reset(.retain_capacity); // key: reuse memory without hitting OS
        var x: usize = 0;
        while (x < camera.hsize) : (x += 1) {
            const r = ray_for_pixel(camera, x, y);
            const c = color_at(world, r, temp_alloc);
            image.writePixel(x, y, c);
        }
    }

    return image;
}

pub fn render(alloc: std.mem.Allocator, camera: *const Camera, world: *const World) !Canvas {
    var image = try Canvas.init(alloc, camera.hsize, camera.vsize);
    errdefer image.deinit(alloc);

    // Decide number of workers.
    const cpu_count = (std.Thread.getCpuCount() catch 1);
    const worker_count = @min(cpu_count, camera.vsize);
    // utils.log(@src(), "Worker count: {}...\n", .{worker_count});

    if (worker_count <= 1) {
        utils.log(@src(), "Rendering single threaded...\n", .{});

        // scratch arena for per-ray allocations (intersections, etc.)
        var scratch = std.heap.ArenaAllocator.init(alloc);
        defer scratch.deinit();
        const temp_alloc = scratch.allocator();

        var y: usize = 0;
        while (y < camera.vsize) : (y += 1) {
            _ = scratch.reset(.retain_capacity); // key: reuse memory without hitting OS
            var x: usize = 0;
            while (x < camera.hsize) : (x += 1) {
                const r = ray_for_pixel(camera, x, y);
                const c = color_at(world, r, temp_alloc);
                image.writePixel(x, y, c);
            }
        }

        return image;
    }

    // utils.log(@src(), "Rendering multi threaded. Starting...\n", .{});
    const Worker = struct {
        camera: *const Camera,
        world: *const World,
        image: *Canvas,
        y0: usize,
        y1: usize,
        parent_alloc: std.mem.Allocator,

        fn run(self: @This()) void {
            var arena = std.heap.ArenaAllocator.init(self.parent_alloc);
            defer arena.deinit();
            const temp_alloc = arena.allocator();

            var y: usize = self.y0;
            while (y < self.y1) : (y += 1) {

                // optional: clear between rows
                _ = arena.reset(.retain_capacity);
                var x: usize = 0;
                while (x < self.camera.hsize) : (x += 1) {
                    const r = ray_for_pixel(self.camera, x, y);
                    const c = color_at(@constCast(self.world), r, temp_alloc);
                    self.image.writePixel(x, y, c);
                }
            }
        }
    };

    // Spawn workers
    var threads = try alloc.alloc(std.Thread, worker_count);
    defer alloc.free(threads);

    // Divide scanlines into chunks
    const rows_per = camera.vsize / worker_count;
    const rem = camera.vsize % worker_count;

    var y_start: usize = 0;
    var i: usize = 0;
    while (i < worker_count) : (i += 1) {
        const extra: usize = if (i < rem) 1 else 0;
        const y_end = y_start + rows_per + extra;

        const w = Worker{
            .camera = camera,
            .world = world,
            .image = &image,
            .y0 = y_start,
            .y1 = y_end,
            .parent_alloc = alloc,
        };

        threads[i] = try std.Thread.spawn(.{}, Worker.run, .{w});
        y_start = y_end;
    }

    // Join all
    for (threads) |t| t.join();

    // utils.log(@src(), "Rendering multi threaded. Ended...\n", .{});
    return image;
}

test "Chap7 -Make sure it works" {
    try std.testing.expect(7 == 7);
}

test "Chap7 -Rendering a world with a camera" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var w = try default_world(gpa.allocator());
    defer w.deinit(gpa.allocator());

    var c = Camera.init(11, 11, std.math.pi / S(2));
    const from = point(0, 0, -5);
    const to = point(0, 0, 0);
    const up = vector(0, 1, 0);
    c.transform = view_transform(from, to, up);

    var image = try renderSingleThread(gpa.allocator(), &c, &w);
    defer image.deinit(gpa.allocator());

    const px5y5 = image.pixelAt(5, 5);

    // print("px5y5: {f}\n", .{px5y5});
    try std.testing.expect(color(0.380661190703326, 0.475826488379158, 0.285495893027495).equals(px5y5));
}

test "Chap7 -Putting it together" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var w = try default_world(gpa.allocator());
    defer w.deinit(gpa.allocator());

    try w.setSingleLight(Light.fromPointLight(PointLight.init_at(point(-10, 10, -10), color(1, 1, 1))));

    var floor = try w.allocator().create(shapesMod.Sphere);
    floor.* = shapesMod.Sphere.init();
    floor.h.material.col = tMod.color(1.0, 0.9, 0.9);
    floor.h.material.diffuse = S(1.0);
    floor.h.material.specular = S(0.0);

    {
        floor.set_transform(mMod.Mat4.scaling(10, 0.01, 10));
        try w.setSingleShape(shapesMod.Shape.fromSphere(floor));
    }

    {
        var left_wall = try w.allocator().create(shapesMod.Sphere);
        left_wall.* = shapesMod.Sphere.init();
        left_wall.h.material = floor.h.material;

        const mTr = mMod.Mat4.translation(0, 0, 5);
        const mRoty = mMod.Mat4.roty(-std.math.pi / S(4));
        const mRotx = mMod.Mat4.rotx(std.math.pi / S(2));
        const mScale = mMod.Mat4.scaling(10, 0.01, 10);
        const mXform = mTr.mulM(&mRoty).mulM(&mRotx).mulM(&mScale);
        left_wall.set_transform(mXform);
        try w.addShape(shapesMod.Shape.fromSphere(left_wall));
    }

    {
        var right_wall = try w.allocator().create(shapesMod.Sphere);
        right_wall.* = shapesMod.Sphere.init();
        right_wall.h.material = floor.h.material;

        const m2a = mMod.Mat4.translation(0, 0, 5);
        const m2b = mMod.Mat4.roty(std.math.pi / S(4));
        const m2c = mMod.Mat4.rotx(std.math.pi / S(2));
        const m2d = mMod.Mat4.scaling(10, 0.01, 10);
        const m3 = m2a.mulM(&m2b).mulM(&m2c).mulM(&m2d);
        right_wall.set_transform(m3);
        const sh1 = shapesMod.Shape.fromSphere(right_wall);
        try w.addShape(sh1);
    }

    {
        var middle = try w.allocator().create(shapesMod.Sphere);
        middle.* = shapesMod.Sphere.init();
        middle.h.material.col = tMod.color(0.1, 1, 0.5);
        middle.h.material.diffuse = S(0.7);
        middle.h.material.specular = S(0.3);

        const xform = mMod.Mat4.translation(-0.5, 1, 0.5);
        middle.set_transform(xform);
        try w.addShape(shapesMod.Shape.fromSphere(middle));
    }

    {
        var right = try w.allocator().create(shapesMod.Sphere);
        right.* = shapesMod.Sphere.init();
        right.h.material.col = tMod.color(0.5, 1, 0.1);
        right.h.material.diffuse = S(0.7);
        right.h.material.specular = S(0.3);

        const xform = mMod.Mat4.translation(1.5, 0.5, -0.5).mulM(&mMod.Mat4.scaling(0.5, 0.5, 0.5));
        right.set_transform(xform);
        try w.addShape(shapesMod.Shape.fromSphere(right));
    }

    {
        var left = try w.allocator().create(shapesMod.Sphere);
        left.* = shapesMod.Sphere.init();
        left.h.material.col = tMod.color(1.0, 0.8, 0.1);
        left.h.material.diffuse = S(0.7);
        left.h.material.specular = S(0.3);

        const xform = mMod.Mat4.translation(-1.5, 0.33, -0.75).mulM(&mMod.Mat4.scaling(0.33, 0.33, 0.33));
        left.set_transform(xform);
        try w.addShape(shapesMod.Shape.fromSphere(left));
    }

    var c = Camera.init(600, 400, std.math.pi / S(3));
    // var c = Camera.init(3456, 2234, std.math.pi / S(3));
    const from = point(0, 1.5, -5);
    const to = point(0, 1, 0);
    const up = vector(0, 1, 0);
    c.transform = view_transform(from, to, up);

    // var image = try render(gpa.allocator(), &c, &w);
    var image = try renderSingleThread(gpa.allocator(), &c, &w);
    defer image.deinit(gpa.allocator());

    try createCanvasFile2(gpa.allocator(), &image, "chap7puttingtogether.ppm");
}

test "Chap9 -Putting it together" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var w = try default_world(gpa.allocator());
    defer w.deinit(gpa.allocator());

    try w.setSingleLight(Light.fromPointLight(PointLight.init_at(point(-5, 5, -10), color(1, 1, 1))));

    var floor = try w.allocator().create(shapesMod.Plane);
    floor.* = shapesMod.Plane.init();
    floor.h.material.col = tMod.color(1.0, 0.9, 0.9);
    floor.h.material.diffuse = S(1.0);
    floor.h.material.specular = S(0.0);

    const deg = std.math.pi / S(180);

    {
        try w.setSingleShape(shapesMod.Shape.fromPlane(floor));
    }

    {
        var leftWall = try w.allocator().create(shapesMod.Plane);
        leftWall.* = shapesMod.Plane.init();
        leftWall.h.material.col = tMod.color(0.5, 0.5, 0.5);
        leftWall.h.material.ambient = S(0.9);
        leftWall.h.material.diffuse = S(0.8);
        leftWall.h.material.specular = S(0.9);

        const rz = mMod.Mat4.rotz(S(90) * deg);
        const ry = mMod.Mat4.roty(S(45) * deg);
        const tr = mMod.Mat4.translation(S(-10), S(0), S(10));
        const m = tr.mulM(&ry).mulM(&rz);
        leftWall.set_transform(m);

        try w.addShape(shapesMod.Shape.fromPlane(leftWall));
    }

    {
        var rightWall = try w.allocator().create(shapesMod.Plane);
        rightWall.* = shapesMod.Plane.init();
        rightWall.h.material.col = tMod.color(0, 0.2, 0);
        rightWall.h.material.ambient = S(0.9);
        rightWall.h.material.diffuse = S(0.9);
        rightWall.h.material.specular = S(0.9);
        rightWall.h.material.shininess = S(900);

        const rz = mMod.Mat4.rotz(S(90) * deg);
        const ry = mMod.Mat4.roty(S(-45) * deg);
        const tr = mMod.Mat4.translation(S(10), S(0), S(-2.1));
        const m = tr.mulM(&ry).mulM(&rz);
        rightWall.set_transform(m);
        try w.addShape(shapesMod.Shape.fromPlane(rightWall));
    }

    {
        var middle = try w.allocator().create(shapesMod.Sphere);
        middle.* = shapesMod.Sphere.init();
        middle.h.material.col = tMod.color(0.1, 1, 0.5);
        middle.h.material.diffuse = S(0.7);
        middle.h.material.specular = S(0.3);

        const xform = mMod.Mat4.translation(-0.5, 1, 0.5);
        middle.set_transform(xform);
        try w.addShape(shapesMod.Shape.fromSphere(middle));
    }

    {
        var right = try w.allocator().create(shapesMod.Sphere);
        right.* = shapesMod.Sphere.init();
        right.h.material.col = tMod.color(0.5, 1, 0.1);
        right.h.material.diffuse = S(0.7);
        right.h.material.specular = S(0.3);

        const xform = mMod.Mat4.translation(1.5, 0.5, -0.5).mulM(&mMod.Mat4.scaling(0.5, 0.5, 0.5));
        right.set_transform(xform);
        try w.addShape(shapesMod.Shape.fromSphere(right));
    }

    {
        var left = try w.allocator().create(shapesMod.Sphere);
        left.* = shapesMod.Sphere.init();
        left.h.material.col = tMod.color(1.0, 0.8, 0.1);
        left.h.material.diffuse = S(0.7);
        left.h.material.specular = S(0.3);

        const xform = mMod.Mat4.translation(-1.5, 0.33, -0.75).mulM(&mMod.Mat4.scaling(0.33, 0.33, 0.33));
        left.set_transform(xform);
        try w.addShape(shapesMod.Shape.fromSphere(left));
    }

    var c = Camera.init(600, 400, std.math.pi / S(3));
    // var c = Camera.init(3456, 2234, std.math.pi / S(3));
    const from = point(0, 1.5, -5);
    const to = point(0, 1, 0);
    const up = vector(0, 1, 0);
    c.transform = view_transform(from, to, up);

    // var image = try render(gpa.allocator(), &c, &w);
    var image = try renderSingleThread(gpa.allocator(), &c, &w);
    defer image.deinit(gpa.allocator());

    try createCanvasFile2(gpa.allocator(), &image, "chap9puttingtogether.ppm");
}

test "Chap10 -Putting it together" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var w = try default_world(gpa.allocator());
    defer w.deinit(gpa.allocator());

    try w.setSingleLight(Light.fromPointLight(PointLight.init_at(point(-5, 5, -10), color(1, 1, 1))));

    var floor = try w.allocator().create(shapesMod.Plane);
    floor.* = shapesMod.Plane.init();
    floor.h.material.col = tMod.color(1.0, 0.9, 0.9);
    floor.h.material.diffuse = S(1.0);
    floor.h.material.specular = S(0.0);
    const stripe = StripePattern.init(color(0.5, 1, 1), color(1, 0, 0));
    floor.h.material.pattern = Pattern.fromStripe(&stripe);

    const deg = std.math.pi / S(180);

    {
        try w.setSingleShape(shapesMod.Shape.fromPlane(floor));
    }

    {
        var leftWall = try w.allocator().create(shapesMod.Plane);
        leftWall.* = shapesMod.Plane.init();
        leftWall.h.material.col = tMod.color(0.5, 0.5, 0.5);
        leftWall.h.material.ambient = S(0.9);
        leftWall.h.material.diffuse = S(0.8);
        leftWall.h.material.specular = S(0.9);

        const rz = mMod.Mat4.rotz(S(90) * deg);
        const ry = mMod.Mat4.roty(S(45) * deg);
        const tr = mMod.Mat4.translation(S(-10), S(0), S(10));
        const m = tr.mulM(&ry).mulM(&rz);
        leftWall.set_transform(m);

        try w.addShape(shapesMod.Shape.fromPlane(leftWall));
    }

    {
        var rightWall = try w.allocator().create(shapesMod.Plane);
        rightWall.* = shapesMod.Plane.init();
        rightWall.h.material.col = tMod.color(0, 0.2, 0);
        rightWall.h.material.ambient = S(0.9);
        rightWall.h.material.diffuse = S(0.9);
        rightWall.h.material.specular = S(0.9);
        rightWall.h.material.shininess = S(900);

        const rz = mMod.Mat4.rotz(S(90) * deg);
        const ry = mMod.Mat4.roty(S(-45) * deg);
        const tr = mMod.Mat4.translation(S(10), S(0), S(-2.1));
        const m = tr.mulM(&ry).mulM(&rz);
        rightWall.set_transform(m);
        try w.addShape(shapesMod.Shape.fromPlane(rightWall));
    }

    {
        var middle = try w.allocator().create(shapesMod.Sphere);
        middle.* = shapesMod.Sphere.init();
        middle.h.material.col = tMod.color(0.1, 1, 0.5);
        middle.h.material.diffuse = S(0.7);
        middle.h.material.specular = S(0.3);

        const xform = mMod.Mat4.translation(-0.5, 1, 0.5);
        middle.set_transform(xform);
        try w.addShape(shapesMod.Shape.fromSphere(middle));
    }

    {
        var right = try w.allocator().create(shapesMod.Sphere);
        right.* = shapesMod.Sphere.init();
        right.h.material.col = tMod.color(0.5, 1, 0.1);
        right.h.material.diffuse = S(0.7);
        right.h.material.specular = S(0.3);

        const xform = mMod.Mat4.translation(1.5, 0.5, -0.5).mulM(&mMod.Mat4.scaling(0.5, 0.5, 0.5));
        right.set_transform(xform);
        try w.addShape(shapesMod.Shape.fromSphere(right));
    }

    {
        var left = try w.allocator().create(shapesMod.Sphere);
        left.* = shapesMod.Sphere.init();
        left.h.material.col = tMod.color(1.0, 0.8, 0.1);
        left.h.material.diffuse = S(0.7);
        left.h.material.specular = S(0.3);

        const xform = mMod.Mat4.translation(-1.5, 0.33, -0.75).mulM(&mMod.Mat4.scaling(0.33, 0.33, 0.33));
        left.set_transform(xform);
        try w.addShape(shapesMod.Shape.fromSphere(left));
    }

    var c = Camera.init(600, 400, std.math.pi / S(3));
    // var c = Camera.init(3456, 2234, std.math.pi / S(3));
    const from = point(0, 1.5, -5);
    const to = point(0, 1, 0);
    const up = vector(0, 1, 0);
    c.transform = view_transform(from, to, up);

    // var image = try render(gpa.allocator(), &c, &w);
    var image = try renderSingleThread(gpa.allocator(), &c, &w);
    defer image.deinit(gpa.allocator());

    try createCanvasFile2(gpa.allocator(), &image, "chap10puttingtogether.ppm");
}
