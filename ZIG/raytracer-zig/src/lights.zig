// src/lights/lights.zig
const std = @import("std");
const print = @import("std").debug.print;
const utils = @import("utils.zig");
const log = utils.log;

const types = @import("types.zig");
const matrix = @import("matrix.zig");
const material_mod = @import("material.zig");

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
const material = material_mod.Material;
const Material = material_mod.Material;

pub const point_light_mod = @import("lights/point_light.zig");
pub const PointLight = point_light_mod.PointLight;

pub const Light = union(enum) {
    point_light: PointLight,

    pub fn equals(self: *const Light, other: Light) bool {
        return switch (self.*) {
            .point_light => |pl| pl.equals(other.point_light),
        };
    }

    pub fn fromPointLight(pl: PointLight) Light {
        return .{ .point_light = pl };
    }

    pub fn init_at(self: *Light, in_position: Tuple, in_intensity: Tuple) Light {
        return switch (self.*) {
            .point_light => |pl| pl.init_at(in_position, in_intensity),
        };
    }

    /// ---
    /// intensity is the color with a magnitude i.e. a vector.
    /// ---
    pub fn intensity(self: *const Light) Tuple {
        return switch (self.*) {
            .point_light => |pl| pl.intensity(),
        };
    }

    pub fn position(self: *const Light) Tuple {
        return switch (self.*) {
            .point_light => |pl| pl.position(),
        };
    }
};

/// ---
/// This lighting() function is what will shade your objects so that they appear
/// three-dimensional.
/// ---
pub fn lighting(matrial: Material, lght: Light, pnt: Tuple, eyev: Tuple, normv: Tuple, in_shadow: bool) Tuple {
    // combine the surface color with the light's color/intensity
    const effective_color = matrial.color().mult(lght.intensity());

    // find the direction to the light source
    const lightv = lght.position().sub(pnt).normalize();

    // compute the ambient contribution
    const ambient = effective_color.muls(matrial.ambient);

    // ---
    // NOTE: Initialize colors to black to ease the if select below.
    // ---
    const black = color(0, 0, 0);
    var diffuse = black;
    var specular = black;

    const shadow_factor = S(!in_shadow);

    // light_dot_normal represents the cosine of the angle between the
    // light vector and the normal vector. A negative number means the
    // light is on the other side of the surface.
    const light_dot_normal = lightv.dot(normv);
    if (light_dot_normal >= S(0)) {
        // compute the diffuse contribution
        diffuse = effective_color.muls(matrial.diffuse).muls(light_dot_normal).muls(shadow_factor);

        // reflect_dot_eye represents the cosine of the angle between the
        // reflection vector and the eye vector. A negative number means the
        // light reflects away from the eye.
        const reflectv = lightv.reflect(normv).muls(S(-1));
        const reflect_dot_eye = reflectv.dot(eyev);
        if (reflect_dot_eye > S(0)) {
            // compute the specular contribution
            const factor = std.math.pow(Scalar, reflect_dot_eye, matrial.shininess);
            specular = lght.intensity().muls(matrial.specular).muls(factor).muls(shadow_factor);
        }
    }

    // Add the three contributions together to get the final shading
    return ambient.add(diffuse.add(specular));
}

test "Chap6 -A point light has a position and intensity" {
    const intensity = color(1, 1, 1);
    const position = point(0, 0, 0);
    const light = Light.fromPointLight(PointLight.init_at(position, intensity));
    try std.testing.expect(light.position().equals(position));
    try std.testing.expect(light.intensity().equals(intensity));
}

test "Chap6 -Lighting with the eye between the light and the surface" {
    const m = material.init();
    const position = point(0, 0, 0);
    const eyev = vector(0, 0, -1);
    const normalv = vector(0, 0, -1);
    const light = Light.fromPointLight(PointLight.init_at(point(0, 0, -10), color(1, 1, 1)));
    const result = lighting(m, light, position, eyev, normalv, false);
    try std.testing.expect(result.equals(color(1.9, 1.9, 1.9)));
}

test "Chap6 -Lighting with the eye between light and surface, eye offset 45°" {
    const m = material.init();
    const position = point(0, 0, 0);
    const eyev = vector(0, std.math.sqrt2 / S(2), -std.math.sqrt2 / S(2));
    const normalv = vector(0, 0, -1);
    const light = Light.fromPointLight(PointLight.init_at(point(0, 0, -10), color(1, 1, 1)));
    const result = lighting(m, light, position, eyev, normalv, false);
    try std.testing.expect(result.equals(color(1, 1, 1)));
}

test "Chap6 -Lighting with eye opposite surface, light offset 45°" {
    const m = material.init();
    const position = point(0, 0, 0);
    const eyev = vector(0, 0, -1);
    const normalv = vector(0, 0, -1);
    const light = Light.fromPointLight(PointLight.init_at(point(0, 10, -10), color(1, 1, 1)));
    const result = lighting(m, light, position, eyev, normalv, false);
    try std.testing.expect(result.equals(color(0.7364, 0.7364, 0.7364)));
}

test "Chap6 -Lighting with eye in the path of the reflection vector" {
    const m = material.init();
    const position = point(0, 0, 0);
    const eyev = vector(0, -std.math.sqrt2 / S(2), -std.math.sqrt2 / S(2));
    const normalv = vector(0, 0, -1);
    const light = Light.fromPointLight(PointLight.init_at(point(0, 10, -10), color(1, 1, 1)));
    const result = lighting(m, light, position, eyev, normalv, false);
    try std.testing.expect(result.equals(color(1.6364, 1.6364, 1.6364)));
}

test "Chap6 -Lighting with the light behind the surface" {
    const m = material.init();
    const position = point(0, 0, 0);
    const eyev = vector(0, 0, 1);
    const normalv = vector(0, 0, -1);
    const light = Light.fromPointLight(PointLight.init_at(point(0, 0, 10), color(1, 1, 1)));
    const result = lighting(m, light, position, eyev, normalv, false);
    try std.testing.expect(result.equals(color(0.1, 0.1, 0.1)));
}

test "Chap8 -Lighting with the surface in shadow" {
    const m = material.init();
    const position = point(0, 0, 0);
    const eyev = vector(0, 0, 1);
    const normalv = vector(0, 0, -1);
    const light = Light.fromPointLight(PointLight.init_at(point(0, 0, -10), color(1, 1, 1)));
    const in_shadow = true;
    const result = lighting(m, light, position, eyev, normalv, in_shadow);
    try std.testing.expect(result.equals(color(0.1, 0.1, 0.1)));
}
