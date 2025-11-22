const std = @import("std");
const print = @import("std").debug.print;
const raytracer_zig = @import("raytracer_zig");
const tuple = @import("tuple"); // comes from build.zig registration

pub fn main() !void {
    // Prints to stderr, ignoring potential errors.
    print("Xll your {s} are belong to us.\n", .{"codebase"});
    try raytracer_zig.bufferedPrint();

    const t = tuple.Tuple.init(4.3, -4.2, 3.1, 1.0);
    print("t = ({}, {}, {}, {})\n", .{ t.x, t.y, t.z, t.w });
    print("tt = {f}\n", .{t});
}

test "simple test" {
    const gpa = std.testing.allocator;
    var list: std.ArrayList(i32) = .empty;
    defer list.deinit(gpa); // Try commenting this out and see if zig detects the memory leak!
    try list.append(gpa, 42);
    try std.testing.expectEqual(@as(i32, 42), list.pop());
}

test "fuzz example" {
    const Context = struct {
        fn testOne(context: @This(), input: []const u8) anyerror!void {
            _ = context;
            // Try passing `--fuzz` to `zig build test` and see if it manages to fail this test case!
            try std.testing.expect(!std.mem.eql(u8, "canyoufindme", input));
        }
    };
    try std.testing.fuzz(Context{}, Context.testOne, .{});
}
