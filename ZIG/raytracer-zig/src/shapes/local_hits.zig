const std = @import("std");
const print = @import("std").debug.print;

const Types = @import("../types.zig");
const Scalar = Types.Scalar;
const S = Types.S;

pub fn LocalHits(comptime N: usize) type {
    return struct {
        const Self = @This();
        count: usize = 0,
        t: [N]Scalar = [_]Scalar{S(0)} ** N,

        pub inline fn add(self: *Self, value: Scalar) void {
            // debug safety
            std.debug.assert(self.count < N);
            self.t[self.count] = value;
            self.count += 1;
        }

        pub inline fn clear(self: *Self) void {
            self.count = 0;
        }

        pub inline fn items(self: *const Self) []const Scalar {
            return self.t[0..self.count];
        }

        pub inline fn sort2(self: *Self) void {
            if (self.count == 2 and self.t[1] < self.t[0]) {
                const tmp = self.t[0];
                self.t[0] = self.t[1];
                self.t[1] = tmp;
            }
        }
    };
}
