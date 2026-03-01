// src/shapes/shape_header.zig
const Mat4 = @import("../matrix.zig").Mat4;
const Material = @import("../material.zig").Material;
const Types = @import("../types.zig");
const matrix = @import("../matrix.zig");
const Intersections = Types.Intersections;

pub const ShapeHeader = struct {
    object_id: usize,
    transformed_m: Mat4,
    transformed_m_inv: Mat4,
    transposed_m: Mat4,
    transposed_m_inv: Mat4,
    material: Material,

    pub fn init(obj_id: usize) ShapeHeader {
        return .{
            .object_id = obj_id,
            .transformed_m = matrix.Mat4.identity(),
            .transformed_m_inv = matrix.Mat4.identity(),
            .transposed_m = matrix.Mat4.identity(),
            .transposed_m_inv = matrix.Mat4.identity(),
            .material = Material.init(),
        };
    }

    pub fn set_material(self: *ShapeHeader, mat: Material) void {
        self.material = mat;
    }

    /// ---
    /// Setting the transform matrix.
    /// NOTE: Also computes the inverse as a side effect...
    /// ---
    pub fn set_transform(self: *ShapeHeader, m: matrix.Mat4) void {
        self.transformed_m = m;
        self.transformed_m_inv = m.inverse();
        self.transposed_m = m.transpose();
        self.transposed_m_inv = m.transpose().inverse();
    }
};
