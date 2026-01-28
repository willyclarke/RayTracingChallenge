const std = @import("std");
const print = @import("std").debug.print;
const utils = @import("utils.zig");
const log = utils.log;

const types_mod = @import("types.zig");
const shapes_mod = @import("shapes.zig");
const matrix_mod = @import("matrix.zig");
const material_mod = @import("material.zig");
const light_mod = @import("lights.zig");
const point_light_mod = @import("lights/point_light.zig");
const canvas_mod = @import("canvas.zig");

const Canvas = canvas_mod.Canvas;

const material = material_mod.Material;
const Material = material_mod.Material;

const light = light_mod.Light;
const Light = light_mod.Light;
const PointLight = point_light_mod.PointLight;
const lighting = light_mod.lighting;

const Shape = shapes_mod.Shape;

const Projectile = types_mod.Projectile;
const Scalar = types_mod.Scalar;
const S = types_mod.S;
const Tuple = types_mod.Tuple;
const Ray = types_mod.Ray;
const approxEq = types_mod.approxEq;
const point = types_mod.Point;
const vector = types_mod.Vector;

pub const Intersection = types_mod.Intersection;
pub const Intersections = types_mod.Intersections;
pub const LocalIntersections = types_mod.LocalIntersections;

pub const World = struct {
    arena: std.heap.ArenaAllocator,

    lights: []const Light,
    shapes: []const Shape, // Shape is union holding pointers to Sphere, Cube etc.
    canvas: Canvas,

    pub fn allocator(self: *World) std.mem.Allocator {
        return self.arena.allocator();
    }

    pub fn init(
        parent_alloc: std.mem.Allocator,
        canvas_width: usize,
        canvas_height: usize,
    ) !World {
        var world: World = undefined;

        world.arena = std.heap.ArenaAllocator.init(parent_alloc);

        world.lights = &[_]Light{}; // will replace
        world.shapes = &[_]Shape{}; //

        //    // ^ I strongly recommend canvas uses parent_alloc, not arena,
        //   so you can deinit it explicitly.
        world.canvas = try Canvas.init(parent_alloc, canvas_width, canvas_height); //

        return world;
    }

    pub fn deinit(self: *World, parent_alloc: std.mem.Allocator) void {
        self.canvas.deinit(parent_alloc);
        self.arena.deinit(); // will free shapes, lights etc...
    }
};

pub fn default_world(parent_alloc: std.mem.Allocator) !World {
    var world = try World.init(parent_alloc, 800, 600);
    const alloc = world.allocator();

    shapes_mod.Sphere.reset_id();

    var s1 = try alloc.create(shapes_mod.Sphere); // shapes_mod.Shape.fromSphere(sphere0);
    s1.* = shapes_mod.Sphere.init();
    s1.h.material.col = types_mod.color(0.8, 1.0, 0.6);
    s1.h.material.diffuse = S(0.7);
    s1.h.material.specular = S(0.2);

    var s2 = try alloc.create(shapes_mod.Sphere); //  shapes_mod.Shape.fromSphere(sphere1);
    s2.* = shapes_mod.Sphere.init();
    s2.set_transform(matrix_mod.Mat4.scaling(0.5, 0.5, 0.5));

    // NOTE: Shapes[] is immutable, build complete array, then assign.
    const shapes = try alloc.alloc(Shape, 2);
    shapes[0] = Shape.fromSphere(s1);
    shapes[1] = Shape.fromSphere(s2);
    world.shapes = shapes;

    // NOTE: Light[] is immutable, build complet array, then assign.
    const lights = try alloc.alloc(Light, 1);
    lights[0] = Light.fromPointLight(PointLight.init_at(point(-10, 10, -10), types_mod.color(1, 1, 1)));
    world.lights = lights;

    return world;
}

pub fn intersect_world(world: *World, r: Ray) Intersections {
    var xs = Intersections.init();
    const alloc = world.allocator();

    for (world.shapes) |shape| {
        const locint = shape.intersect(r);

        var i: usize = 0;
        while (i < locint.count) : (i += 1) {
            try xs.add(alloc, locint.items[i]);
        }
    }

    std.sort.block(
        Intersection,
        xs.intersections_items.items,
        {},
        struct {
            fn lessThan(_: void, a: Intersection, b: Intersection) bool {
                return a.t < b.t;
            }
        }.lessThan,
    );

    return xs;
}

test "Chap7 -Creating a world" {
    try std.testing.expect(3 == 3);
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();

    var world = try World.init(gpa.allocator(), 0, 0);
    defer world.deinit(gpa.allocator());

    try std.testing.expect(world.lights.len == 0);
    try std.testing.expect(world.shapes.len == 0);
}

test "Chap7 -The default world" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var world = try default_world(gpa.allocator());
    defer world.deinit(gpa.allocator());

    try std.testing.expect(800 == world.canvas.width);
    try std.testing.expect(600 == world.canvas.height);
    try std.testing.expect(2 == world.shapes.len);
    try std.testing.expect(1 == world.lights.len);
    try std.testing.expect(1 == world.shapes[0].id());
    try std.testing.expect(2 == world.shapes[1].id());
    try std.testing.expect(world.lights[0].intensity().equals(vector(1, 1, 1)));
    try std.testing.expect(world.lights[0].position().equals(point(-10, 10, -10)));
}

test "Chap7 -Intersect a world with a ray" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var w = try default_world(gpa.allocator());
    defer w.deinit(gpa.allocator());

    const origin = Tuple.point(0, 0, -5);
    const direction = Tuple.vector(0, 0, 1);
    const r = Ray.init(origin, direction);

    // Find all intersections between the ray and objects in the world
    const xs = intersect_world(&w, r);

    // log(@src(), "\nIntersections xs: {f}\n", .{xs});
    try std.testing.expect(3 == 3);
    try std.testing.expect(4 == xs.items().len);
    try std.testing.expect(approxEq(S(4), xs.items()[0].t));
    try std.testing.expect(approxEq(S(4.5), xs.items()[1].t));
    try std.testing.expect(approxEq(S(5.5), xs.items()[2].t));
    try std.testing.expect(approxEq(S(6), xs.items()[3].t));
}
