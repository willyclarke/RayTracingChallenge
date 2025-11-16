const std = @import("std");

var g_start_time_ns: i128 = 0;

/// Simple debug logger with relative time, source location and message.
/// Usage: log(@src(), "value = {}\n", .{v});
pub fn log(src: anytype, comptime fmt: []const u8, args: anytype) void {
    const ns_per_us = std.time.ns_per_us;

    // Init start time on first use
    if (g_start_time_ns == 0) {
        g_start_time_ns = std.time.nanoTimestamp();
    }

    const now_ns = std.time.nanoTimestamp();
    const delta_us: i128 = @divTrunc(now_ns - g_start_time_ns, ns_per_us);

    std.debug.print(
        "[+{d:>6} µs {s}:{d}:{d} {s}] ",
        .{ delta_us, src.file, src.line, src.column, src.fn_name },
    );
    std.debug.print(fmt, args);
}
