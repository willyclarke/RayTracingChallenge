
// One thread per pixel: camera ray -> color_at -> row-major RGBA output.
// scene.sampling == 0: rtc-rust's image, one ray through each pixel centre,
// written to `pixels`; with a lens (focal blur) each pixel averages
// rtc-rust's lens_samples rays instead. scene.sampling == 1: a random point in the pixel and
// random light samples, added into `pixels` as an accumulation buffer
// (frame 0 overwrites it); resolve_* divide by the frame count.
kernel void render(device float4 *pixels          [[buffer(0)]],
                   constant SceneInfo &scene      [[buffer(1)]],
                   device const Shape *shapes     [[buffer(2)]],
                   device const BvhNode *nodes    [[buffer(3)]],
                   device const uint *tables      [[buffer(4)]],
                   device const uint *rt_ids      [[buffer(5), function_constant(METAL_RT)]],
                   raytracing::primitive_acceleration_structure triangles
                                                  [[buffer(6), function_constant(METAL_RT)]],
                   raytracing::primitive_acceleration_structure boxes
                                                  [[buffer(7), function_constant(METAL_RT)]],
                   uint2 gid                      [[thread_position_in_grid]])
{
    RtScene rt;
    if (METAL_RT) {
        rt = RtScene{rt_ids, triangles, boxes};
    }
    Camera c = scene.camera;
    if (gid.x >= c.width || gid.y >= c.height) {
        return;
    }
    uint pixel = gid.y * c.width + gid.x;
    Sampler s = Sampler{pixel, scene.frame, scene.sampling, 0};
    // Deterministic with a lens: rtc-rust's lens_samples rays per pixel.
    // Otherwise one ray, through the centre or a random point of the pixel.
    // Keep color_at to this one call site: with three inlined copies the
    // METAL_RT pipeline returned garbage (values ~4e9) while the BVH one
    // was right (milestone 11).
    uint n = (s.on == 0 && c.aperture > 0.0) ? c.lens_samples : 1;
    float3 sum = 0.0;
    for (uint i = 0; i < n; i++) {
        Ray r;
        if (s.on == 0 && c.aperture > 0.0) {
            float4 l = lens_sample(gid.x, gid.y, i, n);
            r = ray_for_lens(c, gid.x, gid.y, l.z, l.w, l.x, l.y);
        } else if (s.on == 0) {
            r = ray_for_pixel(c, gid.x, gid.y, 0.5, 0.5);
        } else {
            r = ray_for_pixel(c, gid.x, gid.y, random01(pixel, s.frame, 0),
                              random01(pixel, s.frame, 1));
        }
        sum += color_at(r, scene, shapes, nodes, tables, rt, s);
    }
    float4 color = float4(sum / float(n), 1.0);
    pixels[pixel] = (s.on == 0 || s.frame == 0) ? color : pixels[pixel] + color;
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
