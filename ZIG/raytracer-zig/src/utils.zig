const std = @import("std");

var g_start_time_ns: i128 = 0;

/// ---
/// If you redirect output to a file, raw escape codes are ugly. You can check if stderr is a TTY (since std.debug.print writes to stderr):
/// ---
fn useColor() bool {
    const stderr = std.fs.File.stderr();
    return std.posix.isatty(stderr.handle);
}

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

    if (useColor()) {
        const green = "\x1b[32m";
        const acolor = "\x1b[33m";
        const cyan = "\x1b[36m";
        const dim = "\x1b[2m";
        const reset = "\x1b[0m";
        std.debug.print(
            "{s}[+{d:>6} µs{s} {s}- {s} - {s}{d}:{d}{s} ] {s}{s}{s}{s} : ",
            .{ dim, delta_us, reset, green, src.file, acolor, src.line, src.column, dim, reset, cyan, src.fn_name, reset },
        );
    } else {
        std.debug.print("[+{d:>6} µs - {s} - {d}:{d} ] {s} : ", .{
            delta_us,
            src.file,
            src.line,
            src.column,
            src.fn_name,
        });
    }

    std.debug.print(fmt, args);
}
