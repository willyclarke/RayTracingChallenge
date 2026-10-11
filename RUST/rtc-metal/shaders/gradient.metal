#include <metal_stdlib>
using namespace metal;

// One GPU thread per pixel. `gid` is this thread's (x, y) in the 2D grid
// the host dispatched; the host passes the image size so we can index the
// flat, row-major output buffer.
kernel void gradient(device float4 *pixels  [[buffer(0)]],
                     constant uint2 &size   [[buffer(1)]],
                     uint2 gid              [[thread_position_in_grid]])
{
    if (gid.x >= size.x || gid.y >= size.y) {
        return;
    }
    float2 uv = float2(gid) / float2(size - 1);
    pixels[gid.y * size.x + gid.x] = float4(uv.x, uv.y, 0.25, 1.0);
}
