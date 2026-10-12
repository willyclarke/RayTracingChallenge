
// One thread per pixel: camera ray -> color_at -> row-major RGBA output.
kernel void render(device float4 *pixels          [[buffer(0)]],
                   constant SceneInfo &scene      [[buffer(1)]],
                   device const Shape *shapes     [[buffer(2)]],
                   device const BvhNode *nodes    [[buffer(3)]],
                   uint2 gid                      [[thread_position_in_grid]])
{
    Camera c = scene.camera;
    if (gid.x >= c.width || gid.y >= c.height) {
        return;
    }
    uint pixel = gid.y * c.width + gid.x;
    Ray r = ray_for_pixel(c, gid.x, gid.y);
    Sampler s = Sampler{pixel, 0, 0, 0};  // deterministic: rtc-rust's cell centres
    pixels[pixel] = float4(color_at(r, scene, shapes, nodes, s), 1.0);
}
