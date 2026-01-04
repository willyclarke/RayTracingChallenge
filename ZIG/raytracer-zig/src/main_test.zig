//! This file runs ALL tests in the project
//! Run with: zig build test
//! Or directly: zig test src/main_test.zig

const std = @import("std");

// Import all your modules that contain tests
const types = @import("types.zig");
const matrix = @import("matrix.zig");
const material = @import("material.zig");
const shapes = @import("shapes.zig");
const lights = @import("lights.zig");
const canvas = @import("canvas.zig");
const world = @import("world.zig");

// These lines are CRUCIAL — they force Zig to include and run the tests
test {
    // This tells Zig: "run all tests in the current file AND all referenced files"
    std.testing.refAllDecls(@This());

    // Explicitly reference each module's tests
    std.testing.refAllDecls(types);
    std.testing.refAllDecls(matrix);
    std.testing.refAllDecls(material);
    std.testing.refAllDecls(shapes);
    std.testing.refAllDecls(lights);
    std.testing.refAllDecls(canvas);
    std.testing.refAllDecls(world);

    const sphere = @import("shapes/sphere.zig");
    std.testing.refAllDecls(sphere);

    const point_light = @import("lights/point_light.zig");
    std.testing.refAllDecls(point_light);
}

// Optional: keep your simple add test here too
test "basic add" {
    try std.testing.expect(3 + 7 == 10);
}
