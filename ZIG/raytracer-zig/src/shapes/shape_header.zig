// src/shapes/shape_header.zig
const Mat4 = @import("../matrix.zig").Mat4;
const Material = @import("../material.zig").Material;
const Types = @import("../types.zig");
const Intersections = Types.Intersections;

pub const ShapeHeader = struct {
    object_id: usize,
    transformed_m: Mat4,
    transformed_m_inv: Mat4,
    transposed_m: Mat4,
    transposed_m_inv: Mat4,
    material: Material,
};
