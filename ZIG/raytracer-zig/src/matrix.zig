const std = @import("std");
const print = @import("std").debug.print;
const canvas = @import("canvas.zig");
const types = @import("types.zig");
const utils = @import("utils.zig");

const S = types.S;
const Scalar = types.Scalar;
const Tuple = types.Tuple;
const Ray = types.Ray;
const point = types.Point;
const vector = types.Vector;
const approxEq = types.approxEq;
const log = utils.log;

// ───── INVERSE STRATEGY — change this one line to switch! ─────
const UseFast4x4Inverse = true; // ← set to false → uses generic cofactor version
// const UseFast4x4Inverse = false;

/// Generic N×N matrix of Scalars
pub fn Matrix(comptime N: usize) type {
    return struct {
        const Self = @This();

        data: [N][N]Scalar,

        /// CTOR - all zero's
        pub fn zero() Self {
            return .{ .data = .{.{0} ** N} ** N };
        }

        /// CTOR - identity matrix
        pub fn identity() Self {
            var m = Self.zero();
            inline for (0..N) |i| {
                m.data[i][i] = 1.0;
            }
            return m;
        }

        /// Immutable read
        /// r - row
        /// c - col
        pub inline fn get(self: *const Self, r: usize, c: usize) Scalar {
            return self.data[r][c];
        }

        /// Mutable set
        /// r - row
        /// c - col
        pub inline fn set(self: *Self, r: usize, c: usize, v: Scalar) void {
            self.data[r][c] = v;
        }

        /// Pointer to element (so you can assign with `.*`)
        /// r - row
        /// c - col
        pub inline fn ptr(self: *Self, r: usize, c: usize) *Scalar {
            return &self.data[r][c];
        }

        /// Comparison - use approxEq
        pub inline fn equals(self: *const Self, other: *const Self) bool {
            comptime var i = 0;
            inline while (i < N * N) : (i += 1) {
                const r = i / N;
                const c = i % N;
                if (!approxEq(self.data[r][c], other.data[r][c])) {
                    return false;
                }
            }
            return true;
        }

        /// Row view (useful for bulk ops)
        /// r - row
        pub inline fn row(self: *Self, r: usize) *[N]Scalar {
            return &self.data[r];
        }

        /// Const row view
        /// r - row
        pub inline fn rowConst(self: *const Self, r: usize) *const [N]Scalar {
            return &self.data[r];
        }

        /// CTOR from Tuple rows
        pub inline fn rows(r: [N]Tuple) Self {
            comptime if (N != 4)
                @compileError("rows([N]Tuple) is only defined for N == 4");
            return .{
                .data = .{
                    .{ r[0].x, r[0].y, r[0].z, r[0].w },
                    .{ r[1].x, r[1].y, r[1].z, r[1].w },
                    .{ r[2].x, r[2].y, r[2].z, r[2].w },
                    .{ r[3].x, r[3].y, r[3].z, r[3].w },
                },
            };
        }

        // Generic matrix transpose
        pub inline fn transpose(self: *const Self) Self {
            var result = Self.zero();

            inline for (0..N) |r| {
                const a_row = self.rowConst(r);
                inline for (0..N) |c| {
                    result.data[c][r] = a_row[c];
                }
            }
            return result;
        }

        // Generic matrix multiplication
        pub inline fn mulM(self: *const Self, other: *const Self) Self {
            var result = Self.zero();

            inline for (0..N) |i| {
                const a_row = self.rowConst(i);
                inline for (0..N) |j| {
                    var sum: Scalar = 0.0;
                    inline for (0..N) |k| {
                        sum += a_row[k] * other.data[k][j];
                    }
                    result.data[i][j] = sum;
                }
            }
            return result;
        }

        /// CTOR for translation of x, y, z.
        pub inline fn translation(x: Scalar, y: Scalar, z: Scalar) Self {
            return .{
                .data = .{
                    .{ S(1), S(0), S(0), x },
                    .{ S(0), S(1), S(0), y },
                    .{ S(0), S(0), S(1), z },
                    .{ S(0), S(0), S(0), S(1) },
                },
            };
        }
        /// CTOR for translation of x, y, z.
        pub inline fn scaling(x: Scalar, y: Scalar, z: Scalar) Self {
            return .{
                .data = .{
                    .{ x, S(0), S(0), S(0) },
                    .{ S(0), y, S(0), S(0) },
                    .{ S(0), S(0), z, S(0) },
                    .{ S(0), S(0), S(0), S(1) },
                },
            };
        }

        /// Multiplication of 4×4 matrix × tuple
        pub inline fn mulT(self: *const Self, tpl: Tuple) Tuple {
            comptime if (N != 4)
                @compileError("multiplyTuple is only defined for N == 4");

            const x = self.data[0][0] * tpl.x + self.data[0][1] * tpl.y +
                self.data[0][2] * tpl.z + self.data[0][3] * tpl.w;
            const y = self.data[1][0] * tpl.x + self.data[1][1] * tpl.y +
                self.data[1][2] * tpl.z + self.data[1][3] * tpl.w;
            const z = self.data[2][0] * tpl.x + self.data[2][1] * tpl.y +
                self.data[2][2] * tpl.z + self.data[2][3] * tpl.w;
            const w = self.data[3][0] * tpl.x + self.data[3][1] * tpl.y +
                self.data[3][2] * tpl.z + self.data[3][3] * tpl.w;
            return Tuple.init(x, y, z, w);
        }

        /// Returns an (N-1)×(N-1) submatrix by removing one row and one column.
        /// Works for any N (2, 3, 4, ...) → returns Matrix(N-1)
        /// Runtime indices → fully safe, no runtime branching cost (inline loops unroll)
        pub inline fn subMatrix(self: *const Self, comptime remove_row: usize, comptime remove_col: usize) Matrix(N - 1) {
            var result: Matrix(N - 1) = undefined;

            var dst_row: usize = 0;
            comptime {
                @setEvalBranchQuota(100_000);
            }
            inline for (0..N) |src_row| {
                if (src_row == remove_row) continue;

                var dst_col: usize = 0;
                inline for (0..N) |src_col| {
                    if (src_col == remove_col) continue;

                    result.data[dst_row][dst_col] = self.data[src_row][src_col];
                    dst_col += 1;
                }
                dst_row += 1;
            }
            return result;
        }

        /// Compute determinant by cofactor expansion
        pub inline fn determinant(self: *const Self) Scalar {
            if (2 == N) {
                return self.data[0][0] * self.data[1][1] - self.data[0][1] * self.data[1][0];
            }

            var det: Scalar = 0;
            inline for (0..N) |c| {
                det += self.data[0][c] * self.cofactor(0, c);
            }
            return det;
        }

        /// Compute the minor of a matrix.
        /// The minor is the determinant of the
        /// matrix of which remove_row and remove_col
        /// are removed.
        pub inline fn minor(self: *const Self, comptime r: usize, comptime c: usize) Scalar {
            return self.subMatrix(r, c).determinant();
        }

        /// Compute the cofactor of a matrix
        /// The cofactor change sign of the minor when the sum of row + col is odd.
        pub inline fn cofactor(self: *const Self, comptime r: usize, comptime c: usize) Scalar {
            const sign = if (((r + c) % 2) == 0) S(1) else S(-1);
            return sign * self.minor(r, c);
        }

        /// Check if determinant is not 0 to verify that it is possible to invert
        pub inline fn isInvertible(self: *const Self) bool {
            const det = determinant(self);
            return (!(approxEq(det, S(0))));
        }

        /// Slow version with cofactors
        pub fn inverseGeneric(self: *const Self) Self {
            const det = self.determinant();

            // fail when not invertible.
            if (approxEq(det, S(0))) @panic("singular matrix");
            if (approxEq(det, S(0))) return Self.zero();

            var inv = Self.zero();

            inline for (0..N) |r| {
                inline for (0..N) |c| {
                    inv.data[c][r] = self.cofactor(r, c) / det; // NOTE: col row order accomplish the transpose.
                }
            }
            log(@src(), "\ninv:{f}\ndeterminant:{}\n", .{ inv, det });
            return inv;
        }

        pub fn inverse4x4Fast(self: *const Self) Self {
            if (N != 4) {
                // 2x2 and 3x3: generic is fine
                const det = self.determinant();
                if (approxEq(det, 0)) @panic("singular");
                var inv = Self.zero();
                inline for (0..N) |r| {
                    inline for (0..N) |c| {
                        inv.data[c][r] = self.cofactor(r, c) / det;
                    }
                }
                return inv;
            }

            // === FAST 4x4 INVERSE (analytical, ~50 FLOPs) ===
            // Source: https://github.com/g-truc/glm/blob/master/glm/detail/func_matrix.inl
            // Adapted to Zig — this is what everyone uses in production
            const m = self.data;

            var inv: [4][4]Scalar = undefined;

            // Row 0
            inv[0][0] = m[1][1] * m[2][2] * m[3][3] - m[1][1] * m[2][3] * m[3][2] - m[2][1] * m[1][2] * m[3][3] + m[2][1] * m[1][3] * m[3][2] + m[3][1] * m[1][2] * m[2][3] - m[3][1] * m[1][3] * m[2][2];
            inv[0][1] = -m[0][1] * m[2][2] * m[3][3] + m[0][1] * m[2][3] * m[3][2] + m[2][1] * m[0][2] * m[3][3] - m[2][1] * m[0][3] * m[3][2] - m[3][1] * m[0][2] * m[2][3] + m[3][1] * m[0][3] * m[2][2];
            inv[0][2] = m[0][1] * m[1][2] * m[3][3] - m[0][1] * m[1][3] * m[3][2] - m[1][1] * m[0][2] * m[3][3] + m[1][1] * m[0][3] * m[3][2] + m[3][1] * m[0][2] * m[1][3] - m[3][1] * m[0][3] * m[1][2];
            inv[0][3] = -m[0][1] * m[1][2] * m[2][3] + m[0][1] * m[1][3] * m[2][2] + m[1][1] * m[0][2] * m[2][3] - m[1][1] * m[0][3] * m[2][2] - m[2][1] * m[0][2] * m[1][3] + m[2][1] * m[0][3] * m[1][2];

            // Row 1
            inv[1][0] = -m[1][0] * m[2][2] * m[3][3] + m[1][0] * m[2][3] * m[3][2] + m[2][0] * m[1][2] * m[3][3] - m[2][0] * m[1][3] * m[3][2] - m[3][0] * m[1][2] * m[2][3] + m[3][0] * m[1][3] * m[2][2];
            inv[1][1] = m[0][0] * m[2][2] * m[3][3] - m[0][0] * m[2][3] * m[3][2] - m[2][0] * m[0][2] * m[3][3] + m[2][0] * m[0][3] * m[3][2] + m[3][0] * m[0][2] * m[2][3] - m[3][0] * m[0][3] * m[2][2];
            inv[1][2] = -m[0][0] * m[1][2] * m[3][3] + m[0][0] * m[1][3] * m[3][2] + m[1][0] * m[0][2] * m[3][3] - m[1][0] * m[0][3] * m[3][2] - m[3][0] * m[0][2] * m[1][3] + m[3][0] * m[0][3] * m[1][2];
            inv[1][3] = m[0][0] * m[1][2] * m[2][3] - m[0][0] * m[1][3] * m[2][2] - m[1][0] * m[0][2] * m[2][3] + m[1][0] * m[0][3] * m[2][2] + m[2][0] * m[0][2] * m[1][3] - m[2][0] * m[0][3] * m[1][2];

            // Row 2
            inv[2][0] = m[1][0] * m[2][1] * m[3][3] - m[1][0] * m[2][3] * m[3][1] - m[2][0] * m[1][1] * m[3][3] + m[2][0] * m[1][3] * m[3][1] + m[3][0] * m[1][1] * m[2][3] - m[3][0] * m[1][3] * m[2][1];
            inv[2][1] = -m[0][0] * m[2][1] * m[3][3] + m[0][0] * m[2][3] * m[3][1] + m[2][0] * m[0][1] * m[3][3] - m[2][0] * m[0][3] * m[3][1] - m[3][0] * m[0][1] * m[2][3] + m[3][0] * m[0][3] * m[2][1];
            inv[2][2] = m[0][0] * m[1][1] * m[3][3] - m[0][0] * m[1][3] * m[3][1] - m[1][0] * m[0][1] * m[3][3] + m[1][0] * m[0][3] * m[3][1] + m[3][0] * m[0][1] * m[1][3] - m[3][0] * m[0][3] * m[1][1];
            inv[2][3] = -m[0][0] * m[1][1] * m[2][3] + m[0][0] * m[1][3] * m[2][1] + m[1][0] * m[0][1] * m[2][3] - m[1][0] * m[0][3] * m[2][1] - m[2][0] * m[0][1] * m[1][3] + m[2][0] * m[0][3] * m[1][1];

            // Row 3
            inv[3][0] = -m[1][0] * m[2][1] * m[3][2] + m[1][0] * m[2][2] * m[3][1] + m[2][0] * m[1][1] * m[3][2] - m[2][0] * m[1][2] * m[3][1] - m[3][0] * m[1][1] * m[2][2] + m[3][0] * m[1][2] * m[2][1];
            inv[3][1] = m[0][0] * m[2][1] * m[3][2] - m[0][0] * m[2][2] * m[3][1] - m[2][0] * m[0][1] * m[3][2] + m[2][0] * m[0][2] * m[3][1] + m[3][0] * m[0][1] * m[2][2] - m[3][0] * m[0][2] * m[2][1];
            inv[3][2] = -m[0][0] * m[1][1] * m[3][2] + m[0][0] * m[1][2] * m[3][1] + m[1][0] * m[0][1] * m[3][2] - m[1][0] * m[0][2] * m[3][1] - m[3][0] * m[0][1] * m[1][2] + m[3][0] * m[0][2] * m[1][1];
            inv[3][3] = m[0][0] * m[1][1] * m[2][2] - m[0][0] * m[1][2] * m[2][1] - m[1][0] * m[0][1] * m[2][2] + m[1][0] * m[0][2] * m[2][1] + m[2][0] * m[0][1] * m[1][2] - m[2][0] * m[0][2] * m[1][1];

            const det = m[0][0] * inv[0][0] + m[0][1] * inv[1][0] + m[0][2] * inv[2][0] + m[0][3] * inv[3][0];
            if (approxEq(det, 0)) @panic("singular");

            // Then scale all entries by 1/det
            inline for (0..4) |r| {
                inline for (0..4) |c| {
                    inv[r][c] /= det;
                }
            }

            return .{ .data = inv };
        }

        pub fn inverse(self: *const Self) Self {
            if (N != 4) {
                // 2x2 and 3x3: always use generic (it's fast enough)
                return self.inverseGeneric();
            }

            // 4x4: choose strategy at compile time
            if (comptime UseFast4x4Inverse) {
                return self.inverse4x4Fast();
            } else {
                return self.inverseGeneric();
            }
        }

        pub fn transform_ray(self: *const Self, ray: Ray) Ray {
            return Ray{ .origin = self.mulT(ray.origin), .direction = self.mulT(ray.direction) };
        }

        pub fn transform_ray2objspace(self: *const Self, ray: Ray) Ray {
            const inv = self.inverse();
            return Ray{ .origin = inv.mulT(ray.origin), .direction = inv.mulT(ray.direction) };
        }

        /// Rotation matrix around the x-axis.
        pub fn rotx(r: Scalar) Self {
            return .{
                .data = .{
                    .{ S(1), S(0), S(0), S(0) },
                    .{ S(0), std.math.cos(r), -std.math.sin(r), S(0) },
                    .{ S(0), std.math.sin(r), std.math.cos(r), S(0) },
                    .{ S(0), S(0), S(0), S(1) },
                },
            };
        }

        /// Rotation matrix around the y-axis.
        pub fn roty(r: Scalar) Self {
            return .{
                .data = .{
                    .{ std.math.cos(r), S(0), std.math.sin(r), S(0) },
                    .{ S(0), S(1), S(0), S(0) },
                    .{ -std.math.sin(r), S(0), std.math.cos(r), S(0) },
                    .{ S(0), S(0), S(0), S(1) },
                },
            };
        }

        /// Rotation matrix around the z-axis.
        pub fn rotz(r: Scalar) Self {
            return .{
                .data = .{
                    .{ std.math.cos(r), -std.math.sin(r), S(0), S(0) },
                    .{ std.math.sin(r), std.math.cos(r), S(0), S(0) },
                    .{ S(0), S(0), S(1), S(0) },
                    .{ S(0), S(0), S(0), S(1) },
                },
            };
        }

        /// Rotation matrix around the z-axis.
        pub fn shearing(xy: Scalar, xz: Scalar, yx: Scalar, yz: Scalar, zx: Scalar, zy: Scalar) Self {
            return .{
                .data = .{
                    .{ S(1), xy, xz, S(0) },
                    .{ yx, S(1), yz, S(0) },
                    .{ zx, zy, S(1), S(0) },
                    .{ S(0), S(0), S(0), S(1) },
                },
            };
        }

        /// Pretty printing for `{f}` – works for all N.
        pub fn format(self: Self, w: anytype) !void {
            try w.print("Matrix{}(\n", .{N});
            inline for (0..N) |r| {
                const printrow = self.rowConst(r);
                try w.print("  |", .{});
                inline for (0..N) |c| {
                    // width 10, 4 decimals
                    try w.print(" {:>10.4}", .{printrow[c]});
                }
                try w.print(" |\n", .{});
            }
            try w.print(")", .{});
        }

        pub fn formatPtr(self: *const Self, w: anytype) !void {
            return self.*.format(w);
        }
    };
}

pub const Mat2 = Matrix(2);
pub const Mat3 = Matrix(3);
pub const Mat4 = Matrix(4);

test "matrix:Mat4 getters/setters" {
    var M = Mat4.zero();

    // setter/getter
    M.set(1, 2, 42);
    try std.testing.expect(M.get(1, 2) == 42);

    // pointer for in-place update
    M.ptr(1, 2).* = 7;
    try std.testing.expect(M.get(1, 2) == 7);

    // row view, then normal []
    M.row(0)[3] = 3.14;
    try std.testing.expect(M.get(0, 3) == 3.14);
}

test "matrix:Chap3 -Constructing and inspecting a 4x4 matrix" {
    const M = Mat4.rows(.{ Tuple.init(1, 2, 3, 4), Tuple.init(5.5, 6.5, 7.5, 8.5), Tuple.init(9, 10, 11, 12), Tuple.init(13.5, 14.5, 15.5, 16.5) });

    // print(" {f}\n", .{M});

    try std.testing.expect(M.get(0, 0) == S(1));
    try std.testing.expect(M.get(0, 3) == S(4));
    try std.testing.expect(M.get(1, 0) == S(5.5));
    try std.testing.expect(M.get(1, 2) == S(7.5));
    try std.testing.expect(M.get(2, 2) == S(11));
    try std.testing.expect(M.get(3, 0) == S(13.5));
    try std.testing.expect(M.get(3, 2) == S(15.5));
}

test "matrix: Chap3 -A 2x2 matrix ought to be representable" {
    const M = Mat2{ .data = .{ .{ S(-3), S(5) }, .{ S(1), S(-2) } } };
    // print("{f}\n", .{M});
    try std.testing.expect(M.get(0, 0) == S(-3));
    try std.testing.expect(M.get(0, 1) == S(5));
    try std.testing.expect(M.get(1, 0) == S(1));
    try std.testing.expect(M.get(1, 1) == S(-2));
}

test "matrix: Chap3 -A 3x3 matrix ought to be representable" {
    const M = Mat3{ .data = .{ .{ S(-3), S(5), S(0) }, .{ S(1), S(-2), S(-7) }, .{ S(0), S(1), S(1) } } };
    // print("{f}\n", .{M});
    try std.testing.expect(M.get(0, 0) == S(-3));
    try std.testing.expect(M.get(0, 1) == S(5));
    try std.testing.expect(M.get(0, 2) == S(0));
    try std.testing.expect(M.get(1, 0) == S(1));
    try std.testing.expect(M.get(1, 1) == S(-2));
    try std.testing.expect(M.get(1, 2) == S(-7));
    try std.testing.expect(M.get(2, 0) == S(0));
    try std.testing.expect(M.get(2, 1) == S(1));
    try std.testing.expect(M.get(2, 2) == S(1));
}

test "matrix: Chap3 -Matrix equality with identical matrices" {
    const M1 = Mat2{ .data = .{ .{ S(-3), S(5) }, .{ S(1), S(-2) } } };
    const M2 = Mat2{ .data = .{ .{ S(-3), S(5) }, .{ S(1), S(-2) } } };
    const ok = M1.equals(&M2);
    // log(@src(), "\n{f}\n{f}\nAre equal: {}\n", .{ M1, M2, ok });
    try std.testing.expect(ok == true);
}

test "log with timing and sleep" {
    const r = error.SkipZigTest;
    if (r == error.SkipZigTest) return; // make not equal to for running this one

    utils.log(@src(), "Before sleep\n", .{});
    std.Thread.sleep(100 * std.time.ns_per_ms);
    utils.log(@src(), "After 100 ms sleep\n", .{});
}

test "matrix: Chap3 -A 4x4 equality with identical matrices" {
    const M1 = Mat4{ .data = .{ .{ S(1), S(2), S(3), S(4) }, .{ S(5), S(6), S(7), S(8) }, .{ S(9), S(8), S(7), S(6) }, .{ S(5), S(4), S(3), S(2) } } };
    const M2 = M1;
    // log(@src(), "\n{f}\n{f}\n", .{ M1, M2 });
    try std.testing.expect(M2.equals(&M1));
}

test "matrix: Chap3 -A 4x4 equality with different matrices" {
    const M1 = Mat4{ .data = .{ .{ S(1), S(2), S(3), S(4) }, .{ S(5), S(6), S(7), S(8) }, .{ S(9), S(8), S(7), S(6) }, .{ S(5), S(4), S(3), S(2) } } };
    const M2 = Mat4{ .data = .{ .{ S(2), S(3), S(4), S(5) }, .{ S(6), S(7), S(8), S(9) }, .{ S(8), S(7), S(6), S(5) }, .{ S(4), S(3), S(2), S(1) } } };
    // log(@src(), "\n{f}\n{f}\n", .{ M1, M2 });
    try std.testing.expect(!M2.equals(&M1));
}

test "matrix: Chap3 -Multiplying two 4×4 matrices" {
    const A = Mat4{
        .data = .{
            .{ S(1), S(2), S(3), S(4) },
            .{ S(5), S(6), S(7), S(8) },
            .{ S(9), S(8), S(7), S(6) },
            .{ S(5), S(4), S(3), S(2) },
        },
    };
    const B = Mat4{
        .data = .{
            .{ S(-2), S(1), S(2), S(3) },
            .{ S(3), S(2), S(1), S(-1) },
            .{ S(4), S(3), S(6), S(5) },
            .{ S(1), S(2), S(7), S(8) },
        },
    };
    const expected = Mat4{
        .data = .{
            .{ S(20), S(22), S(50), S(48) },
            .{ S(44), S(54), S(114), S(108) },
            .{ S(40), S(58), S(110), S(102) },
            .{ S(16), S(26), S(46), S(42) },
        },
    };
    const C = A.mulM(&B);
    try std.testing.expect(C.equals(&expected));
}

test "matrix: Chap3 -A matrix multiplied by a tuple" {
    const M = Mat4{
        .data = .{
            .{ S(1), S(2), S(3), S(4) },
            .{ S(2), S(4), S(4), S(2) },
            .{ S(8), S(6), S(4), S(1) },
            .{ S(0), S(0), S(0), S(1) },
        },
    };
    const t = Tuple.init(1, 2, 3, 1);
    const result = M.mulT(t);
    const expected = Tuple.init(18, 24, 33, 1);
    try std.testing.expect(result.equals(expected));
}

test "matrix: generic multiplication works for 2×2" {
    const A = Mat2{ .data = .{ .{ S(1), S(2) }, .{ S(3), S(4) } } };
    const B = Mat2{ .data = .{ .{ S(5), S(6) }, .{ S(7), S(8) } } };
    const C = A.mulM(&B);
    const expected = Mat2{ .data = .{ .{ S(19), S(22) }, .{ S(43), S(50) } } };
    try std.testing.expect(C.equals(&expected));
}

test "matrix: Chap3 -Multiplying a matrix by the identity matrix" {
    const M = Mat4{
        .data = .{
            .{ S(0), S(1), S(2), S(4) },
            .{ S(1), S(2), S(4), S(8) },
            .{ S(2), S(4), S(8), S(16) },
            .{ S(4), S(8), S(16), S(32) },
        },
    };
    const result = M.mulM(&Mat4.identity());
    const expected = Mat4{
        .data = .{
            .{ S(0), S(1), S(2), S(4) },
            .{ S(1), S(2), S(4), S(8) },
            .{ S(2), S(4), S(8), S(16) },
            .{ S(4), S(8), S(16), S(32) },
        },
    };

    try std.testing.expect(result.equals(&expected));
}

test "matrix: Chap3 -Identity matrix multiplied by a tuple" {
    const t = Tuple.init(1, 2, 3, 4);
    const result = Mat4.identity().mulT(t);
    const expected = Tuple.init(1, 2, 3, 4);
    try std.testing.expect(result.equals(expected));
}

test "matrix: Chap3 -Transposing a matrix" {
    const M = Mat4{
        .data = .{
            .{ S(0), S(9), S(3), S(0) },
            .{ S(9), S(8), S(0), S(8) },
            .{ S(1), S(8), S(5), S(3) },
            .{ S(0), S(0), S(5), S(8) },
        },
    };
    const result = M.transpose();
    const expected = Mat4{
        .data = .{
            .{ S(0), S(9), S(1), S(0) },
            .{ S(9), S(8), S(8), S(0) },
            .{ S(3), S(0), S(5), S(5) },
            .{ S(0), S(8), S(3), S(8) },
        },
    };

    // log(@src(), "\nM:{f}\nexpected:{f}\n", .{ M, expected });
    try std.testing.expect(result.equals(&expected));
}

test "matrix: Chap3 -Transposing the identity matrix" {
    const M = Mat4.identity();
    const expected = Mat4.identity();
    const transposed = M.transpose();

    // log(@src(), "\nM:{f}\ntransposed:{f}\nexpected:{f}\n", .{ M, transposed, expected });
    try std.testing.expect(transposed.equals(&expected));
}

test "matrix: Chap3 -Calculate the determinant of a 2x2 matrix" {
    const M = Mat2{
        .data = .{ .{ S(1), S(5) }, .{ S(-3), S(2) } },
    };
    const determinant = M.determinant();
    // log(@src(), "\n{f}\nDeterminant:{}\n", .{ M, determinant });
    try std.testing.expect(approxEq(S(17), determinant));
}

test "matrix: subMatrix - 4x4 → 3x3" {
    const M = Mat4{
        .data = .{
            .{ 1, 2, 3, 4 },
            .{ 5, 6, 7, 8 },
            .{ 9, 10, 11, 12 },
            .{ 13, 14, 15, 16 },
        },
    };

    const sub = M.subMatrix(1, 2); // remove row 1, col 2

    const expected = Mat3{
        .data = .{
            .{ 1, 2, 4 },
            .{ 9, 10, 12 },
            .{ 13, 14, 16 },
        },
    };

    try std.testing.expect(sub.equals(&expected));
}

test "matrix: subMatrix - 3x3 → 2x2" {
    const M = Mat3{
        .data = .{
            .{ -3, 5, 0 },
            .{ 1, -2, -7 },
            .{ 0, 1, 1 },
        },
    };

    const sub = M.subMatrix(0, 1); // remove row 0, col 1

    const expected = Mat2{
        .data = .{
            .{ 1, -7 },
            .{ 0, 1 },
        },
    };

    try std.testing.expect(sub.equals(&expected));
}

test "Chap3 -Calculating a minor of a 3x3 matrix" {
    const A = Mat3{
        .data = .{
            .{ 3, 5, 0 },
            .{ 2, -1, -7 },
            .{ 6, -1, 5 },
        },
    };
    const B = A.subMatrix(1, 0);
    const determinantB = B.determinant();
    const minorA = A.minor(1, 0);
    // log(@src(), "\nA:{f}\nB:{f}\ndeterminantB:{}\nminorA:{}\n", .{ A, B, determinantB, minorA });
    try std.testing.expect(approxEq(S(25), determinantB));
    try std.testing.expect(approxEq(S(25), minorA));
}

test "Chap3 -Calculating a cofactor of a 3x3 matrix" {
    const A = Mat3{
        .data = .{
            .{ 3, 5, 0 },
            .{ 2, -1, -7 },
            .{ 6, -1, 5 },
        },
    };
    const minorA0 = A.minor(0, 0);
    const cofactorA0 = A.cofactor(0, 0);
    const minorA1 = A.minor(1, 0);
    const cofactorA1 = A.cofactor(1, 0);

    // log(@src(), "\nminorA:{}\ncofactorA:{}\nminorA1:{}\ncofactorA1:{}\n", .{ minorA0, cofactorA0, minorA1, cofactorA1 });
    try std.testing.expect(approxEq(S(-12), minorA0));
    try std.testing.expect(approxEq(S(-12), cofactorA0));
    try std.testing.expect(approxEq(S(25), minorA1));
    try std.testing.expect(approxEq(S(-25), cofactorA1));
}

test "Chap3 -Calculating the determinant of a 3x3 matrix" {
    const A = Mat3{
        .data = .{
            .{ 1, 2, 6 },
            .{ -5, 8, -4 },
            .{ 2, 6, 4 },
        },
    };
    const cofactorA00 = A.cofactor(0, 0);
    const cofactorA01 = A.cofactor(0, 1);
    const cofactorA02 = A.cofactor(0, 2);
    const determinantA = A.determinant();

    // log(@src(), "\ncofactorA00:{}\n cofactorA01:{}\n cofactorA02:{}\n determinantA:{}\n", .{ cofactorA00, cofactorA01, cofactorA02, determinantA });

    try std.testing.expect(approxEq(S(56), cofactorA00));
    try std.testing.expect(approxEq(S(12), cofactorA01));
    try std.testing.expect(approxEq(S(-46), cofactorA02));
    try std.testing.expect(approxEq(S(-196), determinantA));
}

test "Chap3 -Calculating the determinant of a 4x4 matrix" {
    const A = Mat4{
        .data = .{
            .{ -2, -8, 3, 5 },
            .{ -3, 1, 7, 3 },
            .{ 1, 2, -9, 6 },
            .{ -6, 7, 7, -9 },
        },
    };
    const cofactorA00 = A.cofactor(0, 0);
    const cofactorA01 = A.cofactor(0, 1);
    const cofactorA02 = A.cofactor(0, 2);
    const cofactorA03 = A.cofactor(0, 3);
    const determinantA = A.determinant();

    // log(@src(), "\nA:{f}\n", .{A});
    // log(@src(), "\ncofactorA00:{}\ncofactorA01:{}\ncofactorA02:{}\ncofactorA03:{}\ndeterminantA:{}\n", .{ cofactorA00, cofactorA01, cofactorA02, cofactorA03, determinantA });

    try std.testing.expect(approxEq(S(690), cofactorA00));
    try std.testing.expect(approxEq(S(447), cofactorA01));
    try std.testing.expect(approxEq(S(210), cofactorA02));
    try std.testing.expect(approxEq(S(51), cofactorA03));
    try std.testing.expect(approxEq(S(-4071), determinantA));
    try std.testing.expect(true == A.isInvertible());
}

test "Chap3 -Testing an invertible matrix for invertibility" {
    const A = Mat4{
        .data = .{
            .{ 6, 4, 4, 4 },
            .{ 5, 5, 7, 6 },
            .{ 4, -9, 3, -7 },
            .{ 9, 1, 7, -6 },
        },
    };
    const cofactorA00 = A.cofactor(0, 0);
    const cofactorA01 = A.cofactor(0, 1);
    const cofactorA02 = A.cofactor(0, 2);
    const cofactorA03 = A.cofactor(0, 3);
    const determinantA = A.determinant();

    // log(@src(), "\nA:{f}\n", .{A});
    // log(@src(), "\ncofactorA00:{}\ncofactorA01:{}\ncofactorA02:{}\ncofactorA03:{}\ndeterminantA:{}\n", .{ cofactorA00, cofactorA01, cofactorA02, cofactorA03, determinantA });

    try std.testing.expect(approxEq(S(-668), cofactorA00));
    try std.testing.expect(approxEq(S(112), cofactorA01));
    try std.testing.expect(approxEq(S(620), cofactorA02));
    try std.testing.expect(approxEq(S(-260), cofactorA03));
    try std.testing.expect(approxEq(S(-2120), determinantA));
    try std.testing.expect(true == A.isInvertible());
}

test "Chap3 -Testing a noninvertible matrix for invertibility" {
    const A = Mat4{
        .data = .{
            .{ -4, 2, -2, -3 },
            .{ 9, 6, 2, 6 },
            .{ 0, -5, 1, -5 },
            .{ 0, 0, 0, 0 },
        },
    };
    const cofactorA00 = A.cofactor(0, 0);
    const cofactorA01 = A.cofactor(0, 1);
    const cofactorA02 = A.cofactor(0, 2);
    const cofactorA03 = A.cofactor(0, 3);
    const determinantA = A.determinant();

    // log(@src(), "\nA:{f}\n", .{A});
    // log(@src(), "\ncofactorA00:{}\ncofactorA01:{}\ncofactorA02:{}\ncofactorA03:{}\ndeterminantA:{}\n", .{ cofactorA00, cofactorA01, cofactorA02, cofactorA03, determinantA });

    try std.testing.expect(approxEq(S(0), cofactorA00));
    try std.testing.expect(approxEq(S(0), cofactorA01));
    try std.testing.expect(approxEq(S(0), cofactorA02));
    try std.testing.expect(approxEq(S(0), cofactorA03));
    try std.testing.expect(approxEq(S(0), determinantA));
    try std.testing.expect(false == A.isInvertible());
}

test "Chap3 -Calculating the inverse of a matrix" {
    const A = Mat4{
        .data = .{
            .{ -5, 2, 6, -8 },
            .{ 1, -5, 1, 8 },
            .{ 7, 7, -6, -7 },
            .{ 1, -3, 7, 4 },
        },
    };
    const B = Mat4{
        .data = .{
            .{ 0.21805, 0.45113, 0.24060, -0.04511 },
            .{ -0.80827, -1.45677, -0.44361, 0.52068 },
            .{ -0.07895, -0.22368, -0.05263, 0.19737 },
            .{ -0.52256, -0.81391, -0.30075, 0.30639 },
        },
    };
    const cofactorA23 = A.cofactor(2, 3);
    const cofactorA32 = A.cofactor(3, 2);
    const determinantA = A.determinant();

    // const inverseA = A.inverse();
    // log(@src(), "\n       A:{f}\n", .{A});
    // log(@src(), "\n       B:{f}\n", .{B});
    // log(@src(), "\ninverseA:{f}\n", .{inverseA});
    // log(@src(), "\ncofactorA23:{}\ncofactorA32:{}\ndeterminantA:{}\n", .{ cofactorA23, cofactorA32, determinantA });

    try std.testing.expect(approxEq(S(-160), cofactorA23));
    try std.testing.expect(approxEq(S(105), cofactorA32));
    try std.testing.expect(approxEq(S(532), determinantA));
    try std.testing.expect(approxEq(S(-160) / S(532), cofactorA23 / determinantA));
    try std.testing.expect(approxEq(B.get(3, 2), cofactorA23 / determinantA));

    try std.testing.expect(true == A.isInvertible());
    try std.testing.expect(A.inverse().equals(&B));
}

test "Chap3 -Calculating the inverse of another matrix" {
    const A = Mat4{
        .data = .{
            .{ 8, -5, 9, 2 },
            .{ 7, 5, 6, 1 },
            .{ -6, 0, 9, 6 },
            .{ -3, 0, -9, -4 },
        },
    };
    const B = Mat4{
        .data = .{
            .{ -0.15385, -0.15385, -0.28205, -0.53846 },
            .{ -0.07692, 0.12308, 0.02564, 0.03077 },
            .{ 0.35897, 0.35897, 0.43590, 0.92308 },
            .{ -0.69231, -0.69231, -0.76923, -1.92308 },
        },
    };

    try std.testing.expect(true == A.isInvertible());
    try std.testing.expect(A.inverse().equals(&B));
}

test "Chap3 -Calculating the inverse of a third matrix" {
    const A = Mat4{
        .data = .{
            .{ 9, 3, 0, 9 },
            .{ -5, -2, -6, -3 },
            .{ -4, 9, 6, 4 },
            .{ -7, 6, 6, 2 },
        },
    };
    const B = Mat4{
        .data = .{
            .{ -0.04074, -0.07778, 0.14444, -0.22222 },
            .{ -0.07778, 0.03333, 0.36667, -0.33333 },
            .{ -0.02901, -0.14630, -0.10926, 0.12963 },
            .{ 0.17778, 0.06667, -0.26667, 0.33333 },
        },
    };

    try std.testing.expect(true == A.isInvertible());
    try std.testing.expect(A.inverse().equals(&B));
}

test "Chap3 -Multiplying a product by its invers" {
    const A = Mat4{
        .data = .{
            .{ 3, -9, 7, 3 },
            .{ 3, -8, 2, -9 },
            .{ -4, 4, 4, 1 },
            .{ -6, 5, -1, 1 },
        },
    };

    const B = Mat4{
        .data = .{
            .{ 8, 2, 2, 2 },
            .{ 3, -1, 7, 0 },
            .{ 7, 0, 5, 4 },
            .{ 6, -2, 0, 5 },
        },
    };

    const C = A.mulM(&B);

    // Explanation:
    //“One last thing to note about the inverse: at the beginning of this section,
    // you read that “if you multiply some matrix A by another matrix B,
    // producing C, you can multiply C by the inverse of B to get A again.”
    // Well, we can’t let such a statement slide by unproven!
    // Add one more test to show that the inverse does, in truth,
    // behave as described.”

    try std.testing.expect(C.mulM(&B.inverse()).equals(&A));
}

test "Chap4 -Multiplying by a translation matrix" {
    const transform = Mat4.translation(5, -3, 2);
    const p = point(-3, 4, 5);
    const translp = transform.mulT(p);
    const expected = point(2, 1, 7);
    // log(@src(), "\ntransform:{f}\np:{f}\ntranslp:{f}\n", .{ transform, p, translp });
    try std.testing.expect(expected.equals(translp));
}

// Further, if you take the inverse of a translation matrix, you get another
// translation matrix that moves points in reverse.
test "Chap4 -Multiplying by the inverse of a translation matrix" {
    const transform = Mat4.translation(5, -3, 2);
    const inv = transform.inverse();
    const p = point(-3, 4, 5);
    const translp = inv.mulT(p);
    const expected = point(-8, 7, 3);
    // log(@src(), "\ntransform:{f}\np:{f}\ntranslp:{f}\n", .{ transform, p, translp });
    try std.testing.expect(expected.equals(translp));
}

test "Chap4 -Translation does not affect vectors" {
    const transform = Mat4.translation(5, -3, 2);
    const v = Tuple.vector(-3, 4, 5);
    const translv = transform.mulT(v);
    // log(@src(), "\ntransform:{f}\nv      :{f}\ntranslv:{f}\n", .{ transform, v, translv });
    try std.testing.expect(v.equals(translv));
}

test "Chap4 -A scaling matrix applied to a point" {
    const transform = Mat4.scaling(2, 3, 4);
    const p = point(-4, 6, 8);
    const scaledp = transform.mulT(p);
    const expectedp = point(-8, 18, 32);
    // log(@src(), "\ntransform:{f}\np:{f}\nscaledp:{f}\n", .{ transform, p, scaledp });
    try std.testing.expect(scaledp.equals(expectedp));
}

test "Chap4 -A scaling matrix applied to a vector" {
    const transform = Mat4.scaling(2, 3, 4);
    const v = Tuple.vector(-4, 6, 8);
    const scaledv = transform.mulT(v);
    const expectedv = Tuple.vector(-8, 18, 32);
    // log(@src(), "\ntransform:{f}\nv      :{f}\nscaledv:{f}\n", .{ transform, v, scaledv });
    try std.testing.expect(scaledv.equals(expectedv));
}

test "Chap4 -Multiplying by the inverse of a scaling matrix" {
    const transform = Mat4.scaling(2, 3, 4);
    const inv = transform.inverse();
    const v = Tuple.vector(-4, 6, 8);
    const invv = inv.mulT(v);
    const expectedv = Tuple.vector(-2, 2, 2);
    // log(@src(), "\ntransform:{f}\nv          :{f}\ninv scaledv:{f}\n", .{ transform, v, invv });
    try std.testing.expect(invv.equals(expectedv));
}

test "Chap4 -Reflection is scaling by a negative value" {
    const transform = Mat4.scaling(-1, 1, 1);
    const p = point(2, 3, 4);
    const expectedp = point(-2, 3, 4);
    const reflectedp = transform.mulT(p);
    try std.testing.expect(reflectedp.equals(expectedp));
}

test "Chap4 -Convert from Degrees to Radians" {
    const angle = S(180);
    const anglerad = types.Deg2Rad(angle);
    try std.testing.expect(approxEq(anglerad, std.math.pi));
}

test "Chap4 -Rotating a point around the x axis" {
    const p = point(0, 1, 0);
    const half_quarter = Mat4.rotx(std.math.pi / S(4));
    const full_quarter = Mat4.rotx(std.math.pi / S(2));
    const point_hq = half_quarter.mulT(p);
    const point_fq = full_quarter.mulT(p);
    const expect_hq = point(0, std.math.sqrt2 / S(2), std.math.sqrt2 / S(2));
    const expect_fq = point(0, 0, 1);
    try std.testing.expect(point_hq.equals(expect_hq));
    try std.testing.expect(point_fq.equals(expect_fq));
}

test "Chap4 -The inverse of an x-rotation rotates in the opposite direction" {
    const p = point(0, 1, 0);
    const half_quarter = Mat4.rotx(std.math.pi / S(4));
    const inv = half_quarter.inverse();
    const point_opposite = inv.mulT(p);
    const expectedp = point(S(0), std.math.sqrt2 / S(2), -std.math.sqrt2 / S(2));
    try std.testing.expect(point_opposite.equals(expectedp));
}

test "Chap4 -Rotating a point around the y axis" {
    const p = point(0, 0, 1);
    const half_quarter = Mat4.roty(std.math.pi / S(4));
    const full_quarter = Mat4.roty(std.math.pi / S(2));
    const point_hq = half_quarter.mulT(p);
    const point_fq = full_quarter.mulT(p);
    const expect_hq = point(std.math.sqrt2 / S(2), 0, std.math.sqrt2 / S(2));
    const expect_fq = point(1, 0, 0);
    try std.testing.expect(point_hq.equals(expect_hq));
    try std.testing.expect(point_fq.equals(expect_fq));
}

test "Chap4 -Rotating a point around the z axis" {
    const p = point(0, 1, 0);
    const half_quarter = Mat4.rotz(std.math.pi / S(4));
    const full_quarter = Mat4.rotz(std.math.pi / S(2));
    const point_hq = half_quarter.mulT(p);
    const point_fq = full_quarter.mulT(p);
    const expect_hq = point(-std.math.sqrt2 / S(2), std.math.sqrt2 / S(2), 0);
    const expect_fq = point(-1, 0, 0);
    try std.testing.expect(point_hq.equals(expect_hq));
    try std.testing.expect(point_fq.equals(expect_fq));
}

test "Rotating a point around the z axis" {
    const p = point(0, 1, 0);
    const half_quarter = Mat4.rotz(std.math.pi / S(4));
    const full_quarter = Mat4.rotz(std.math.pi / S(2));
    const point_hq = half_quarter.mulT(p);
    const point_fq = full_quarter.mulT(p);
    const expect_hq = point(-std.math.sqrt2 / S(2), std.math.sqrt2 / S(2), 0);
    const expect_fq = point(-1, 0, 0);
    try std.testing.expect(point_hq.equals(expect_hq));
    try std.testing.expect(point_fq.equals(expect_fq));
}

test "Chap4 -A shearing transformation moves x in proportion to y" {
    const transform = Mat4.shearing(1, 0, 0, 0, 0, 0);
    const p = point(2, 3, 4);
    const transformedp = transform.mulT(p);
    const expectedp = Tuple.point(5, 3, 4);
    // log(@src(), "\ninput point:{f}\ntransform:{f}\ntransformedp:{f}\n", .{ p, transform, transformedp });
    try std.testing.expect(transformedp.equals(expectedp));
}

test "Chap4 -A shearing transformation moves x in proportion to z" {
    const transform = Mat4.shearing(0, 1, 0, 0, 0, 0);
    const p = point(2, 3, 4);
    const transformedp = transform.mulT(p);
    const expectedp = Tuple.point(6, 3, 4);
    // log(@src(), "\ninput point:{f}\ntransform:{f}\ntransformedp:{f}\n", .{ p, transform, transformedp });
    try std.testing.expect(transformedp.equals(expectedp));
}

test "Chap4 -A shearing transformation moves y in proportion to x" {
    const transform = Mat4.shearing(0, 0, 1, 0, 0, 0);
    const p = point(2, 3, 4);
    const transformedp = transform.mulT(p);
    const expectedp = Tuple.point(2, 5, 4);
    // log(@src(), "\ninput point:{f}\ntransform:{f}\ntransformedp:{f}\n", .{ p, transform, transformedp });
    try std.testing.expect(transformedp.equals(expectedp));
}

test "Chap4 -A shearing transformation moves y in proportion to z" {
    const transform = Mat4.shearing(0, 0, 0, 1, 0, 0);
    const p = point(2, 3, 4);
    const transformedp = transform.mulT(p);
    const expectedp = Tuple.point(2, 7, 4);
    // log(@src(), "\ninput point:{f}\ntransform:{f}\ntransformedp:{f}\n", .{ p, transform, transformedp });
    try std.testing.expect(transformedp.equals(expectedp));
}

test "Chap4 -A shearing transformation moves z in proportion to x" {
    const transform = Mat4.shearing(0, 0, 0, 0, 1, 0);
    const p = point(2, 3, 4);
    const transformedp = transform.mulT(p);
    const expectedp = Tuple.point(2, 3, 6);
    // log(@src(), "\ninput point:{f}\ntransform:{f}\ntransformedp:{f}\n", .{ p, transform, transformedp });
    try std.testing.expect(transformedp.equals(expectedp));
}

test "Chap4 -A shearing transformation moves z in proportion to y" {
    const transform = Mat4.shearing(0, 0, 0, 0, 0, 1);
    const p = point(2, 3, 4);
    const transformedp = transform.mulT(p);
    const expectedp = Tuple.point(2, 3, 7);
    // log(@src(), "\ninput point:{f}\ntransform:{f}\ntransformedp:{f}\n", .{ p, transform, transformedp });
    try std.testing.expect(transformedp.equals(expectedp));
}

test "Chap4 -Individual transformations are applied in sequence" {
    const p = point(1, 0, 1);
    const A = Mat4.rotx(std.math.pi / S(2));
    const B = Mat4.scaling(5, 5, 5);
    const C = Mat4.translation(10, 5, 7);

    // Apply rotation around x first.
    const p2 = A.mulT(p);
    try std.testing.expect(p2.equals(point(1, -1, 0)));

    // then apply scaling.
    const p3 = B.mulT(p2);
    try std.testing.expect(p3.equals(point(5, -5, 0)));

    // then apply translation.
    const p4 = C.mulT(p3);
    try std.testing.expect(p4.equals(point(15, 0, 7)));
}

test "Chap4 -Chained transformations must be applied in reverse order" {
    const p = point(1, 0, 1);
    const A = Mat4.rotx(std.math.pi / S(2));
    const B = Mat4.scaling(5, 5, 5);
    const C = Mat4.translation(10, 5, 7);
    const T = C.mulM(&B.mulM(&A));
    const transformedp = T.mulT(p);
    const exptectedp = point(15, 0, 7);
    try std.testing.expect(transformedp.equals(exptectedp));
}

test "Chap4 -Fluent Chained transformations must be applied in reverse order" {
    const p = point(1, 0, 1);
    const T = Mat4.translation(10, 5, 7).mulM(&Mat4.scaling(5, 5, 5).mulM(&Mat4.rotx(std.math.pi / S(2))));
    const transformedp = T.mulT(p);
    const exptectedp = point(15, 0, 7);
    try std.testing.expect(transformedp.equals(exptectedp));
}

fn drawSquare(
    ptrCanvas: *canvas.Canvas,
    center_x: usize,
    center_y: usize,
    half_size: usize,
    ptrcolor: *const types.Color,
) void {
    var y: usize = center_y - half_size;
    while (y <= center_y + half_size) : (y += 1) {
        var x: usize = center_x - half_size;
        while (x <= center_x + half_size) : (x += 1) {
            ptrCanvas.writePixel(x, y, ptrcolor.*);
        }
    }
}

test "Chap4 -Putting It Together" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const alloc = gpa.allocator();

    var c = try canvas.Canvas.init(alloc, 900, 900);
    defer c.deinit(alloc);

    const origin = point(0, 0, 0);
    var p = point(0, 0, 0);
    try std.testing.expect(p.equals(origin));

    var color = types.color(1, 0, 0);

    const canvastranslate = Mat4.translation(S(c.width) / S(2), S(c.height) / S(2), 0);

    // Use the translation number x a percentage t set scaling.
    const scale = Mat4.scaling(canvastranslate.get(0, 3) * S(0.75), canvastranslate.get(1, 3) * S(0.75), 0);

    var canvpoint = canvastranslate.mulM(&scale).mulT(p);
    var centerx = types.toUsizeSaturated(canvpoint.x, S(0), S(c.width));
    var centery = types.toUsizeSaturated(canvpoint.y, S(0), S(c.height));
    const squaresize = types.toUsizeSaturated(S(c.width) / S(20), S(3), S(10));
    drawSquare(&c, centerx, centery, squaresize, &color);

    color.z = 1;
    p.x = 1;

    for (0..12) |idx| {
        const alfa = S(idx) / S(12) * std.math.tau;
        color.y = S(1) - color.z;
        color.z = S(idx) / S(12);

        canvpoint = canvastranslate.mulM(&scale.mulM(&Mat4.rotz(alfa))).mulT(p);
        centerx = types.toUsizeSaturated(canvpoint.x, S(0), S(c.width));
        centery = types.toUsizeSaturated(canvpoint.y, S(0), S(c.height));
        drawSquare(&c, centerx, centery, squaresize, &color);
    }

    try canvas.createCanvasFile(&c, "chap4clock.ppm");
}

test "Chap5 -Translating a ray" {
    const r = Ray.init(Tuple.point(1, 2, 3), Tuple.vector(0, 1, 0));
    const m = Mat4.translation(3, 4, 5);
    const r2 = m.transform_ray(r);
    const ray_expected = Ray.init(Tuple.point(4, 6, 8), Tuple.vector(0, 1, 0));
    try std.testing.expect(r2.origin.equals(ray_expected.origin));
    try std.testing.expect(r2.direction.equals(ray_expected.direction));
}

test "Chap5 -Scaling a ray" {
    const r = Ray.init(Tuple.point(1, 2, 3), Tuple.vector(0, 1, 0));
    const m = Mat4.scaling(2, 3, 4);
    const r2 = m.transform_ray(r);
    const ray_expected = Ray.init(Tuple.point(2, 6, 12), Tuple.vector(0, 3, 0));
    try std.testing.expect(r2.origin.equals(ray_expected.origin));
    try std.testing.expect(r2.direction.equals(ray_expected.direction));
}
