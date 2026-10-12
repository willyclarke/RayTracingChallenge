
// One thread per pixel: camera ray -> color_at -> row-major RGBA output.
// scene.sampling == 0: rtc-rust's image, one ray through each pixel centre,
// written to `pixels`. scene.sampling == 1: a random point in the pixel and
// random light samples, added into `pixels` as an accumulation buffer
// (frame 0 overwrites it); resolve_* divide by the frame count.
kernel void render(device float4 *pixels          [[buffer(0)]],
                   constant SceneInfo &scene      [[buffer(1)]],
                   device const Shape *shapes     [[buffer(2)]],
                   device const BvhNode *nodes    [[buffer(3)]],
                   device const uint *tables      [[buffer(4)]],
                   uint2 gid                      [[thread_position_in_grid]])
{
    Camera c = scene.camera;
    if (gid.x >= c.width || gid.y >= c.height) {
        return;
    }
    uint pixel = gid.y * c.width + gid.x;
    Sampler s = Sampler{pixel, scene.frame, scene.sampling, 0};
    if (s.on == 0) {
        Ray r = ray_for_pixel(c, gid.x, gid.y, 0.5, 0.5);
        pixels[pixel] = float4(color_at(r, scene, shapes, nodes, tables, s), 1.0);
        return;
    }
    Ray r = ray_for_pixel(c, gid.x, gid.y, random01(pixel, s.frame, 0),
                          random01(pixel, s.frame, 1));
    float4 color = float4(color_at(r, scene, shapes, nodes, tables, s), 1.0);
    pixels[pixel] = s.frame == 0 ? color : pixels[pixel] + color;
}

// size_frames = (width, height, frames accumulated); frames > 0.
kernel void resolve_buffer(device const float4 *accum      [[buffer(0)]],
                           device float4 *out              [[buffer(1)]],
                           constant uint3 &size_frames     [[buffer(2)]],
                           uint2 gid                       [[thread_position_in_grid]])
{
    if (gid.x >= size_frames.x || gid.y >= size_frames.y) {
        return;
    }
    uint i = gid.y * size_frames.x + gid.x;
    out[i] = accum[i] / float(size_frames.z);
}

// The accumulated mean, clamped for display.
kernel void resolve_texture(device const float4 *accum               [[buffer(0)]],
                            texture2d<float, access::write> out      [[texture(0)]],
                            constant uint3 &size_frames              [[buffer(2)]],
                            uint2 gid                                [[thread_position_in_grid]])
{
    if (gid.x >= size_frames.x || gid.y >= size_frames.y ||
        gid.x >= out.get_width() || gid.y >= out.get_height()) {
        return;
    }
    uint i = gid.y * size_frames.x + gid.x;
    float3 mean = accum[i].rgb / float(size_frames.z);
    out.write(float4(clamp(mean, 0.0, 1.0), 1.0), gid);
}
