const std = @import("std");
const print = @import("std").debug.print;

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

const sMod = @import("shapes/shapes.zig");

const tMod = @import("types.zig");
const mMod = @import("matrix.zig");

const Scalar = tMod.Scalar;
const S = tMod.S;
const Tuple = tMod.Tuple;
const approxEq = tMod.approxEq;
const point = tMod.Point;
const vector = tMod.Vector;
const color = tMod.Vector;
const matrix = mMod.Mat4;
const Matrix = mMod.Mat4;
const Ray = tMod.Ray;

pub fn render(alloc: std.mem.Allocator, camera: *const Camera, world: *World) !Canvas {
    var image = try Canvas.init(alloc, camera.hsize, camera.vsize);
    errdefer image.deinit(alloc);

    var y: usize = 0;
    while (y < camera.vsize) : (y += 1) {
        var x: usize = 0;
        while (x < camera.hsize) : (x += 1) {
            const r = ray_for_pixel(camera, x, y);
            const c = color_at(world, r);
            image.writePixel(x, y, c);
        }
    }

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

    var image = try render(gpa.allocator(), &c, &w);
    defer image.deinit(gpa.allocator());

    const px5y5 = image.pixelAt(5, 5);

    // print("px5y5: {f}\n", .{px5y5});
    try std.testing.expect(color(0.38066, 0.47583, 0.2855).equals(px5y5));
}

test "Chap7 -Putting it together" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var w = try default_world(gpa.allocator());
    defer w.deinit(gpa.allocator());

    try w.setSingleLight(Light.fromPointLight(PointLight.init_at(point(-10, 10, -10), color(1, 1, 1))));

    var floor = try w.allocator().create(sMod.Sphere);
    floor.* = sMod.Sphere.init();
    floor.h.material.col = tMod.color(1.0, 0.9, 0.9);
    floor.h.material.diffuse = S(1.0);
    floor.h.material.specular = S(0.0);

    {
        floor.set_transform(mMod.Mat4.scaling(10, 0.01, 10));
        try w.setSingleShape(sMod.Shape.fromSphere(floor));
    }

    {
        var left_wall = try w.allocator().create(sMod.Sphere);
        left_wall.* = sMod.Sphere.init();
        left_wall.h.material = floor.h.material;

        const mTr = mMod.Mat4.translation(0, 0, 5);
        const mRoty = mMod.Mat4.roty(-std.math.pi / S(4));
        const mRotx = mMod.Mat4.rotx(std.math.pi / S(2));
        const mScale = mMod.Mat4.scaling(10, 0.01, 10);
        const mXform = mTr.mulM(&mRoty).mulM(&mRotx).mulM(&mScale);
        left_wall.set_transform(mXform);
        try w.addShape(sMod.Shape.fromSphere(left_wall));
    }

    {
        var right_wall = try w.allocator().create(sMod.Sphere);
        right_wall.* = sMod.Sphere.init();
        right_wall.h.material = floor.h.material;

        const m2a = mMod.Mat4.translation(0, 0, 5);
        const m2b = mMod.Mat4.roty(std.math.pi / S(4));
        const m2c = mMod.Mat4.rotx(std.math.pi / S(2));
        const m2d = mMod.Mat4.scaling(10, 0.01, 10);
        const m3 = m2a.mulM(&m2b).mulM(&m2c).mulM(&m2d);
        right_wall.set_transform(m3);
        const sh1 = sMod.Shape.fromSphere(right_wall);
        try w.addShape(sh1);
    }

    {
        var middle = try w.allocator().create(sMod.Sphere);
        middle.* = sMod.Sphere.init();
        middle.h.material.col = tMod.color(0.1, 1, 0.5);
        middle.h.material.diffuse = S(0.7);
        middle.h.material.specular = S(0.3);

        const xform = mMod.Mat4.translation(-0.5, 1, 0.5);
        middle.set_transform(xform);
        try w.addShape(sMod.Shape.fromSphere(middle));
    }

    {
        var right = try w.allocator().create(sMod.Sphere);
        right.* = sMod.Sphere.init();
        right.h.material.col = tMod.color(0.5, 1, 0.1);
        right.h.material.diffuse = S(0.7);
        right.h.material.specular = S(0.3);

        const xform = mMod.Mat4.translation(1.5, 0.5, -0.5).mulM(&mMod.Mat4.scaling(0.5, 0.5, 0.5));
        right.set_transform(xform);
        try w.addShape(sMod.Shape.fromSphere(right));
    }

    {
        var left = try w.allocator().create(sMod.Sphere);
        left.* = sMod.Sphere.init();
        left.h.material.col = tMod.color(1.0, 0.8, 0.1);
        left.h.material.diffuse = S(0.7);
        left.h.material.specular = S(0.3);

        const xform = mMod.Mat4.translation(-1.5, 0.33, -0.75).mulM(&mMod.Mat4.scaling(0.33, 0.33, 0.33));
        left.set_transform(xform);
        try w.addShape(sMod.Shape.fromSphere(left));
    }

    var c = Camera.init(100, 40, std.math.pi / S(3));
    const from = point(0, 1.5, -5);
    const to = point(0, 1, 0);
    const up = vector(0, 1, 0);
    c.transform = view_transform(from, to, up);

    var image = try render(gpa.allocator(), &c, &w);
    defer image.deinit(gpa.allocator());

    try createCanvasFile(&image, "chap7puttingtogether.ppm");
}
