const std = @import("std");
const print = @import("std").debug.print;
const types = @import("types.zig");
const utils = @import("utils.zig");
const Sphere = @import("shapes/sphere.zig").Sphere;
const log = utils.log;

pub const Shape = union(enum) {
    sphere: Sphere,
    // cube: Cube,

    pub fn intersect(self: *const Shape, ray: types.Ray) types.Intersections {
        return switch (self.*) {
            .sphere => |s| s.intersect(ray),
            // .cube => |c| c.intersect(ray),
        };
    }
};

test "sphere works" {
    const s = Sphere.init();
    try std.testing.expect(s.radius == 1);
}
