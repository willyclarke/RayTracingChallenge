const std = @import("std");
const print = @import("std").debug.print;
const utils = @import("utils.zig");
const log = utils.log;

const camera_mod = @import("canvas.zig");
const tMod = @import("types.zig");
const mMod = @import("matrix.zig");
const material_mod = @import("material.zig");
const light_mod = @import("lights.zig");
const point_light_mod = @import("lights/point_light.zig");
const canvas_mod = @import("canvas.zig");
const shapes_mod = @import("shapes/shapes.zig");
const intersection_mod = @import("shapes/intersections.zig");
const sphereMod = @import("shapes/sphere.zig");
const shapesMod = @import("shapes/shapes.zig");

const Camera = camera_mod.Camera;
const Canvas = canvas_mod.Canvas;

const material = material_mod.Material;
const Material = material_mod.Material;

const light = light_mod.Light;
const Light = light_mod.Light;
const PointLight = point_light_mod.PointLight;
const lighting = light_mod.lighting;

const Shape = shapes_mod.Shape;
pub const Sphere = sphereMod.Sphere;

const Projectile = tMod.Projectile;
const Scalar = tMod.Scalar;
const S = tMod.S;
const Tuple = tMod.Tuple;
const Ray = tMod.Ray;
const approxEq = tMod.approxEq;
const point = tMod.Point;
const vector = tMod.Vector;
const color = tMod.Vector;
const matrix = mMod.Mat4;
const Matrix = mMod.Mat4;

pub const Intersection = intersection_mod.Intersection;
pub const Intersections = intersection_mod.Intersections;
pub const LocalIntersections = intersection_mod.LocalIntersections;

pub const World = struct {
    arena: std.heap.ArenaAllocator,

    lights: std.ArrayListUnmanaged(Light) = .{}, //[]Light,
    shapes: std.ArrayListUnmanaged(Shape) = .{}, // []const Shape, // Shape is union holding pointers to Sphere, Cube etc.
    //
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

        world.lights = .{}; //&[_]Light{}; // will replace
        world.shapes = .{}; //&[_]Shape{}; //

        //    // ^ I strongly recommend canvas uses parent_alloc, not arena,
        //   so you can deinit it explicitly.
        world.canvas = try Canvas.init(parent_alloc, canvas_width, canvas_height); //

        return world;
    }

    pub fn deinit(self: *World, parent_alloc: std.mem.Allocator) void {
        // lists were allocated from arena, so no need to deinit them individually
        self.canvas.deinit(parent_alloc);
        self.arena.deinit(); // will free shapes, lights etc...
    }

    pub fn addShape(self: *World, s: Shape) !void {
        try self.shapes.append(self.allocator(), s);
    }

    pub fn setSingleShape(self: *World, s: Shape) !void {
        self.shapes.clearRetainingCapacity();
        try self.shapes.append(self.allocator(), s);
    }

    pub fn setSingleLight(self: *World, l: Light) !void {
        self.lights.clearRetainingCapacity();
        try self.lights.append(self.allocator(), l);
    }

    pub fn addLight(self: *World, l: Light) !void {
        try self.lights.append(self.allocator(), l);
    }

    pub fn shapesSlice(self: *const World) []const Shape {
        return self.shapes.items;
    }

    pub fn lightsSlice(self: *const World) []const Light {
        return self.lights.items;
    }

    pub fn shapePtr(self: *const World, index: usize) *const Shape {
        return &self.shapes.items[index];
    }
};

pub const PrepareComputations = struct {
    const Self = @This();

    t: Scalar,
    ptrShape: *const Shape,
    point: Tuple,
    eyev: Tuple,
    normalv: Tuple,
    over_point: Tuple,
    inside: bool,

    pub fn init(t: Scalar, ptrShape: *const Shape) PrepareComputations {
        return .{
            .t = t,
            .ptrShape = ptrShape,
            .point = Tuple.point(0, 0, -1),
            .eyev = Tuple.vector(0, 0, -1),
            .normalv = Tuple.vector(0, 0, 1),
            .over_point = Tuple.point(0, 0, -1).add(Tuple.vector(0, 0, 1).muls(tMod.EPSILON)), // point+normalv*EPSILON
            .inside = false,
        };
    }

    pub fn normal_at(self: PrepareComputations, point_xs: Tuple) Tuple {
        return self.ptrShape.normal_at(point_xs);
    }
};

pub fn prepare_computations(i: Intersection, r: Ray) PrepareComputations {
    // instantiate a data structure for storing some precomputed values
    // copy the intersection's properties, for convenience
    var comps = PrepareComputations.init(i.t, i.ptrShape);

    // precompute some useful values
    comps.point = r.position(i.t);
    comps.eyev = r.direction.neg();
    comps.normalv = i.ptrShape.normal_at(comps.point);

    if (comps.normalv.dot(comps.eyev) < S(0)) {
        comps.inside = true;
        comps.normalv = comps.normalv.neg();
    }

    comps.over_point = comps.point.add(comps.normalv.muls(tMod.EPSILON));

    return comps;
}

pub fn shade_hit(world: *const World, comps: PrepareComputations, temp_alloc: std.mem.Allocator) Tuple {
    const shadowed = is_shadowed(world, comps.over_point, temp_alloc);
    return lighting(comps.ptrShape.material().*, world.lightsSlice()[0], comps.over_point, comps.eyev, comps.normalv, shadowed);
}

pub fn is_shadowed(world: *const World, p: Tuple, temp_alloc: std.mem.Allocator) bool {
    // 1. Measure the distance from point to the light source by subtracting point
    // from the light position, and taking the magnitude of the resulting vector.
    // Call this distance.
    const v = world.lightsSlice()[0].point_light.pos.sub(p);
    const distance = v.mag();

    // 2. Create a ray from point toward the light source by normalizing the vector
    // from step 1.
    const direction = v.normalize();
    const r = Ray.init(p, direction);

    // 3. Intersect the world with that ray.
    var xs = intersect_world(world, r, temp_alloc);
    defer xs.deinit(temp_alloc); // or world.allocator(), whichever you use

    // 4. Check to see if there was a hit, and if so, whether t is less than distance. If
    // so, the hit lies between the point and the light source, and the point is in
    // shadow.
    const hit_opt = xs.hit();
    if (hit_opt == null) return false;
    const hit = hit_opt.?; // safe now
    if (hit.t < distance) return true;

    return false;
}

pub fn default_world(parent_alloc: std.mem.Allocator) !World {
    var world = try World.init(parent_alloc, 800, 600);
    const alloc = world.allocator();

    shapes_mod.Sphere.reset_id();

    var s1 = try alloc.create(shapes_mod.Sphere); // shapes_mod.Shape.fromSphere(sphere0);
    s1.* = shapes_mod.Sphere.init();
    s1.h.material.col = tMod.color(0.8, 1.0, 0.6);
    s1.h.material.diffuse = S(0.7);
    s1.h.material.specular = S(0.2);

    var s2 = try alloc.create(shapes_mod.Sphere); //  shapes_mod.Shape.fromSphere(sphere1);
    s2.* = shapes_mod.Sphere.init();
    s2.set_transform(mMod.Mat4.scaling(0.5, 0.5, 0.5));

    try world.addShape(Shape.fromSphere(s1));
    try world.addShape(Shape.fromSphere(s2));

    // NOTE: Shapes[] is immutable, build complete array, then assign.
    // const shapes = try alloc.alloc(Shape, 2);
    // shapes[0] = Shape.fromSphere(s1);
    // shapes[1] = Shape.fromSphere(s2);
    // world.shapes = shapes;

    // NOTE: Light[] is immutable, build complet array, then assign.
    // const lights = try alloc.alloc(Light, 1);
    // lights[0] = Light.fromPointLight(PointLight.init_at(point(-10, 10, -10), types_mod.color(1, 1, 1)));
    // world.lights = lights;
    try world.addLight(Light.fromPointLight(PointLight.init_at(point(-10, 10, -10), tMod.color(1, 1, 1))));

    return world;
}

pub fn intersect_world(world: *const World, r: Ray, temp_alloc: std.mem.Allocator) Intersections {
    var xs = Intersections.init();
    // const alloc = world.allocator();

    // IMPORTANT: iterate by pointer so ptrShape points into the world's stable storage
    for (world.shapesSlice()) |*shape| {
        const locint = shape.intersect(r);

        var i: usize = 0;

        while (i < locint.count) : (i += 1) {
            try xs.add(temp_alloc, .{
                .t = locint.items[i].t,
                .ptrShape = shape,
            });
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

pub fn color_at(world: *const World, r: Ray, temp_alloc: std.mem.Allocator) Tuple {
    // 1. Call intersect_world to find the intersections of the given ray with the given world.
    var xs = intersect_world(world, r, temp_alloc);
    defer xs.deinit(temp_alloc); // or world.allocator(), whichever you use

    // 2. Find the hit from the resulting intersections.
    const hit_opt = xs.hit();

    // 3. Return the color black if there is no such intersection.
    if (hit_opt == null) return color(0, 0, 0);

    // 4. Otherwise, precompute the necessary values with prepare_computations.
    const hit = hit_opt.?; // safe now
    const comps = prepare_computations(hit, r);

    // 5. Finally, call shade_hit to find the color at the hit.
    return shade_hit(world, comps, temp_alloc);
}

pub fn view_transform(from: Tuple, to: Tuple, up: Tuple) Matrix {
    //
    // 1. Compute the forward vector by subtracting from from to. Normalize the result.
    // 2. Compute the left vector by taking the cross product of forward and the
    // normalized up vector.
    // 3. Compute the true_up vector by taking the cross product of left and forward.
    // This allows your original up vector to be only approximately up, which
    // makes framing scenes a lot easier, since you don’t need to personally
    // break out a calculator to figure out the precise upward direction.
    // 4. With these left, true_up, and forward vectors, you can now construct a matrix
    // that represents the orientation transformation:
    // 5. All that’s left is to append a translation to that transformation to move the
    // scene into place before orienting it. Multiply orientation by translation(-from.x,
    // -from.y, -from.z), and you’re golden!
    //
    const forward = to.sub(from).normalize();
    const forward_neg = forward.neg();
    const upn = up.normalize();
    const left = forward.cross(upn);
    const true_up = left.cross(forward);
    var orientation = mMod.mat4FromRows(vector(left.x, left.y, left.z), vector(true_up.x, true_up.y, true_up.z), vector(forward_neg.x, forward_neg.y, forward_neg.z), point(0, 0, 0));
    const translation = matrix.translation(-from.x, -from.y, -from.z);
    return orientation.mulM(&translation);
}

test "Chap7 -Creating a world" {
    try std.testing.expect(3 == 3);
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();

    var world = try World.init(gpa.allocator(), 0, 0);
    defer world.deinit(gpa.allocator());

    try std.testing.expect(world.lightsSlice().len == 0);
    try std.testing.expect(world.shapesSlice().len == 0);
}

test "Chap7 -The default world" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var world = try default_world(gpa.allocator());
    defer world.deinit(gpa.allocator());

    try std.testing.expect(800 == world.canvas.width);
    try std.testing.expect(600 == world.canvas.height);
    try std.testing.expect(2 == world.shapesSlice().len);
    try std.testing.expect(1 == world.lightsSlice().len);
    try std.testing.expect(1 == world.shapesSlice()[0].id());
    try std.testing.expect(2 == world.shapesSlice()[1].id());
    try std.testing.expect(world.lightsSlice()[0].intensity().equals(vector(1, 1, 1)));
    try std.testing.expect(world.lightsSlice()[0].position().equals(point(-10, 10, -10)));
}

test "Chap7 -Intersect a world with a ray" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var w = try default_world(gpa.allocator());
    defer w.deinit(gpa.allocator());

    var temp_arena = std.heap.ArenaAllocator.init(gpa.allocator());
    defer temp_arena.deinit();
    const temp_alloc = temp_arena.allocator();

    const origin = point(0, 0, -5);
    const direction = vector(0, 0, 1);
    const r = Ray.init(origin, direction);

    // Find all intersections between the ray and objects in the world
    const xs = intersect_world(&w, r, temp_alloc);

    // log(@src(), "\nIntersections xs: {f}\n", .{xs});
    try std.testing.expect(3 == 3);
    try std.testing.expect(4 == xs.items().len);
    try std.testing.expect(approxEq(S(4), xs.items()[0].t));
    try std.testing.expect(approxEq(S(4.5), xs.items()[1].t));
    try std.testing.expect(approxEq(S(5.5), xs.items()[2].t));
    try std.testing.expect(approxEq(S(6), xs.items()[3].t));
}

test "Chap7 -Precomputing the state of an intersection" {
    const origin = point(0, 0, -5);
    const direction = vector(0, 0, 1);
    const r = Ray.init(origin, direction);

    var sphere = Sphere.init();
    var s = Shape.fromSphere(&sphere);

    const i = Intersection{
        .t = S(4),
        .ptrShape = &s,
    };

    var comps = PrepareComputations.init(i.t, &s);

    comps = prepare_computations(i, r);
    try std.testing.expect(comps.t == i.t);
    try std.testing.expect(comps.ptrShape.id() == i.ptrShape.id());
    try std.testing.expect(comps.point.equals(point(0, 0, -1)) == true);
    try std.testing.expect(comps.eyev.equals(vector(0, 0, -1)) == true);
    try std.testing.expect(comps.normalv.equals(vector(0, 0, -1)) == true);
}

test "Chap7 -The hit, when an intersection occurs on the outside" {
    const r = Ray.init(point(0, 0, -5), vector(0, 0, 1));

    var sphere = Sphere.init();
    var s = Shape.fromSphere(&sphere);

    const i = Intersection{
        .t = S(4),
        .ptrShape = &s,
    };

    var comps = PrepareComputations.init(i.t, &s);

    comps = prepare_computations(i, r);
    try std.testing.expect(comps.t == i.t);
    try std.testing.expect(comps.inside == false);
}

test "Chap7 -The hit, when an intersection occurs on the inside" {
    const r = Ray.init(point(0, 0, 0), vector(0, 0, 1));

    var sphere = Sphere.init();
    var s = Shape.fromSphere(&sphere);

    const i = Intersection{
        .t = S(1),
        .ptrShape = &s,
    };

    var comps = PrepareComputations.init(i.t, &s);
    comps = prepare_computations(i, r);

    // log(@src(), "\ncomps.point : {f}\n", .{comps.point});
    try std.testing.expect(comps.point.equals(point(0, 0, 1)) == true);
    try std.testing.expect(comps.eyev.equals(vector(0, 0, -1)) == true);
    try std.testing.expect(comps.inside == true);
    // normal would have been (0, 0, 1), but is inverted!
    try std.testing.expect(comps.normalv.equals(vector(0, 0, -1)) == true);
}

test "Chap7 -Shading an intersection" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var w = try default_world(gpa.allocator());
    defer w.deinit(gpa.allocator());
    const r = Ray.init(point(0, 0, -5), vector(0, 0, 1));
    const shape = w.shapesSlice()[0];

    const i = Intersection{
        .t = S(4),
        .ptrShape = &shape,
    };

    var comps = PrepareComputations.init(i.t, &shape);
    comps = prepare_computations(i, r);

    const c = shade_hit(&w, comps, w.allocator());

    // log(@src(), "\nc : {f}\n", .{c});
    try std.testing.expect(c.equals(color(0.380661190703326, 0.475826488379158, 0.285495893027495)) == true);
}

test "Chap7 -Shading an intersection from the inside" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var w = try default_world(gpa.allocator());
    defer w.deinit(gpa.allocator());

    try w.setSingleLight(Light.fromPointLight(PointLight.init_at(point(0, 0.25, 0), color(1, 1, 1))));

    const r = Ray.init(point(0, 0, 0), vector(0, 0, 1));
    const ptrShape = w.shapePtr(1);

    const i = Intersection{
        .t = S(0.5),
        .ptrShape = ptrShape,
    };

    const comps = prepare_computations(i, r);
    const c = shade_hit(&w, comps, w.allocator());

    // log(@src(), "\nc : {f}\n", .{c});
    try std.testing.expect(c.equals(color(0.904984439883870, 0.904984439883870, 0.904984439883870)) == true);
}

test "Chap7 -The color when a ray misses" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var w = try default_world(gpa.allocator());
    defer w.deinit(gpa.allocator());
    const r = Ray.init(point(0, 0, -5), vector(0, 1, 0));
    const c = color_at(&w, r, gpa.allocator());
    try std.testing.expect(c.equals(color(0, 0, 0)) == true);
}

test "Chap7 -The color when a ray hits" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var w = try default_world(gpa.allocator());
    defer w.deinit(gpa.allocator());

    const r = Ray.init(point(0, 0, -5), vector(0, 0, 1));
    const c = color_at(&w, r, gpa.allocator());
    // log(@src(), "c: {f}\n", .{c});
    try std.testing.expect(c.equals(color(0.380661190703326, 0.475826488379158, 0.285495893027495)) == true);
}

test "Chap7 -The color with an intersection behind the ray" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var w = try default_world(gpa.allocator());
    defer w.deinit(gpa.allocator());

    var outer = &w.shapesSlice()[0]; // The outer sphere
    outer.material().ambient = S(1);

    var inner = &w.shapesSlice()[1]; // The inner sphere
    inner.material().ambient = S(1);

    const r = Ray.init(point(0, 0, 0.75), vector(0, 0, -1));

    const c = color_at(&w, r, gpa.allocator());
    try std.testing.expect(c.equals(inner.material().color()) == true);
}

test "Chap7 -The transformation matrix for the default orientation" {
    const from = point(0, 0, 0);
    const to = point(0, 0, 1);
    const up = vector(0, 1, 0);
    const t = view_transform(from, to, up);
    const scaling = vector(t.data[0][0], t.data[1][1], t.data[2][2]);
    // print("vt={f}\n", .{t});
    try std.testing.expect(true == (scaling.equals(vector(-1, 1, -1))));
}

test "Chap7 -The view transformation moves the world" {
    const from = point(0, 0, 8);
    const to = point(0, 0, 0);
    const up = vector(0, 1, 0);
    const t = view_transform(from, to, up);
    const translation = vector(t.data[0][3], t.data[1][3], t.data[2][3]);
    // print("vt={f}\n", .{t});
    try std.testing.expect(true == (translation.equals(vector(0, 0, -8))));
}

test "Chap7 -An arbitrary view transformation" {
    const from = point(1, 3, 2);
    const to = point(4, -2, 8);
    const up = vector(1, 1, 0);
    const t = view_transform(from, to, up);

    const e = mMod.Mat4{
        .data = .{
            .{ -0.507092552837109900000, 0.507092552837109900000, 0.676123403782813200000, -2.366431913239846000000 },
            .{ 0.767715933859680100000, 0.606091526731326300000, 0.121218305346265240000, -2.828427124746189400000 },
            .{ -0.358568582800318060000, 0.597614304667196800000, -0.717137165600636100000, 0.000000000000000000000 },
            .{ 0.000000000000000000000, 0.000000000000000000000, 0.000000000000000000000, 1.000000000000000000000 },
        },
    };

    // print("t={f}\n", .{t});
    // print("expect={f}\n", .{e});
    try std.testing.expect(true == (t.equals(&e)));
}

test "Chap8 -There is no shadow when nothing is collinear with point and light" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var w = try default_world(gpa.allocator());
    defer w.deinit(gpa.allocator());

    var temp_arena = std.heap.ArenaAllocator.init(gpa.allocator());
    defer temp_arena.deinit();
    const temp_alloc = temp_arena.allocator();

    const p = point(0, 10, 0);
    const result = is_shadowed(&w, p, temp_alloc);

    try std.testing.expect(result == false);
}

test "Chap8 -The shadow when an object is between the point and the light" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var w = try default_world(gpa.allocator());
    defer w.deinit(gpa.allocator());

    var temp_arena = std.heap.ArenaAllocator.init(gpa.allocator());
    defer temp_arena.deinit();
    const temp_alloc = temp_arena.allocator();

    const p = point(10, -10, 10);
    const result = is_shadowed(&w, p, temp_alloc);

    try std.testing.expect(result == true);
}

test "Chap8 -There is no shadow when an object is behind the light" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var w = try default_world(gpa.allocator());
    defer w.deinit(gpa.allocator());

    var temp_arena = std.heap.ArenaAllocator.init(gpa.allocator());
    defer temp_arena.deinit();
    const temp_alloc = temp_arena.allocator();

    const p = point(-20, 20, -20);
    const result = is_shadowed(&w, p, temp_alloc);

    try std.testing.expect(result == false);
}

test "Chap8 -There is no shadow when an object is behind the point" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var w = try default_world(gpa.allocator());
    defer w.deinit(gpa.allocator());

    var temp_arena = std.heap.ArenaAllocator.init(gpa.allocator());
    defer temp_arena.deinit();
    const temp_alloc = temp_arena.allocator();

    const p = point(-2, 2, -2);
    const result = is_shadowed(&w, p, temp_alloc);

    try std.testing.expect(result == false);
}

test "Chap8 -shade_hit() is given an intersection in shadow" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var w = try default_world(gpa.allocator());
    defer w.deinit(gpa.allocator());
    try w.setSingleLight(Light.fromPointLight(PointLight.init_at(point(0, 0, -10), color(1, 1, 1))));

    var temp_arena = std.heap.ArenaAllocator.init(gpa.allocator());
    defer temp_arena.deinit();
    // const temp_alloc = temp_arena.allocator();

    const s1 = try w.allocator().create(sphereMod.Sphere);
    s1.* = sphereMod.Sphere.init();
    try w.setSingleShape(shapesMod.Shape.fromSphere(s1));

    const s2 = try w.allocator().create(sphereMod.Sphere);
    s2.* = sphereMod.Sphere.init();
    s2.set_transform(mMod.Mat4.translation(0, 0, 10));
    try w.addShape(shapesMod.Shape.fromSphere(s2));

    const r = Ray.init(point(0, 0, 5), vector(0, 0, 1));

    const ptrShape = w.shapePtr(1);
    const i = Intersection{
        .t = S(4),
        .ptrShape = ptrShape,
    };
    const comps = prepare_computations(i, r);
    // const shadowed = is_shadowed(&w, r.origin, temp_alloc);
    const c = shade_hit(&w, comps, w.allocator());

    // print("shade_hit: {}\n", .{shadowed});
    // print("s1.transform: {f}\n", .{s1.transform()});
    // print("s2.transform: {f}\n", .{s2.transform()});
    // log(@src(), "\nc : {f}\n", .{c});
    try std.testing.expect(c.equals(color(0.1, 0.1, 0.1)) == true);
}

test "Chap8 -The hit should offset the point" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var w = try default_world(gpa.allocator());
    defer w.deinit(gpa.allocator());
    try w.setSingleLight(Light.fromPointLight(PointLight.init_at(point(0, 0, -10), color(1, 1, 1))));

    const s1 = try w.allocator().create(sphereMod.Sphere);
    s1.* = sphereMod.Sphere.init();
    s1.set_transform(mMod.Mat4.translation(0, 0, 1));
    try w.setSingleShape(shapesMod.Shape.fromSphere(s1));
    const ptrShape = w.shapePtr(0);

    const r = Ray.init(point(0, 0, -5), vector(0, 0, 1));

    const i = Intersection{
        .t = S(5),
        .ptrShape = ptrShape,
    };

    const comps = prepare_computations(i, r);

    // The point has been moved in the negative direction of the normalv.
    // So the z-value should be slightly bigger the compared to the over_point.z
    try std.testing.expect(comps.over_point.z < -tMod.EPSILON / S(2));
    try std.testing.expect(comps.point.z > comps.over_point.z);
    // utils.log(@src(), "comps.point: {f}\n", .{comps.point});
    // utils.log(@src(), "comps_over.point: {f}\n", .{comps.over_point});
}
