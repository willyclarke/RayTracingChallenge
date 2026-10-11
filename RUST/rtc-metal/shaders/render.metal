
// One thread per pixel: camera ray -> color_at -> row-major RGBA output.
kernel void render(device float4 *pixels          [[buffer(0)]],
                   constant SceneInfo &scene      [[buffer(1)]],
                   device const Shape *shapes     [[buffer(2)]],
                   uint2 gid                      [[thread_position_in_grid]])
{
    Camera c = scene.camera;
    if (gid.x >= c.width || gid.y >= c.height) {
        return;
    }
    Ray r = ray_for_pixel(c, gid.x, gid.y);
    pixels[gid.y * c.width + gid.x] = float4(color_at(r, scene, shapes), 1.0);
}
