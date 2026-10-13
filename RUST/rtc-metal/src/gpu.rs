//! Metal host side: device, command queue and the compiled shader library.

use crate::scene::GpuScene;
use crate::types::{GpuBvhNode, GpuSceneInfo, GpuShape};
use objc2::rc::Retained;
use objc2::runtime::ProtocolObject;
use objc2_foundation::NSString;
use objc2_metal::{
    MTLBuffer, MTLCommandBuffer, MTLCommandBufferStatus, MTLCommandEncoder, MTLCommandQueue,
    MTLComputeCommandEncoder, MTLComputePipelineState, MTLCreateSystemDefaultDevice, MTLDevice,
    MTLLibrary, MTLResourceOptions, MTLSize, MTLTexture,
};
use std::ptr::NonNull;

/// All MSL sources as one translation unit, compiled by the Metal driver when
/// `Gpu::new` runs. `rt.metal` comes first: the kernels below it call into it.
const SHADER_SRC: &str = concat!(
    include_str!("../shaders/rt.metal"),
    include_str!("../shaders/gradient.metal"),
    include_str!("../shaders/render.metal"),
    include_str!("../shaders/probes.metal"),
);

type Buffer = Retained<ProtocolObject<dyn MTLBuffer>>;
type Pipeline = Retained<ProtocolObject<dyn MTLComputePipelineState>>;
pub(crate) type Encoder = ProtocolObject<dyn MTLComputeCommandEncoder>;

/// Owns the Metal objects that live for the whole program.
pub struct Gpu {
    device: Retained<ProtocolObject<dyn MTLDevice>>,
    queue: Retained<ProtocolObject<dyn MTLCommandQueue>>,
    library: Retained<ProtocolObject<dyn MTLLibrary>>,
}

impl Gpu {
    pub fn new() -> Result<Self, String> {
        let device = MTLCreateSystemDefaultDevice().ok_or("no Metal device")?;
        let queue = device
            .newCommandQueue()
            .ok_or("could not create command queue")?;
        let library = device
            .newLibraryWithSource_options_error(&NSString::from_str(SHADER_SRC), None)
            .map_err(|e| format!("shader compile failed: {}", e.localizedDescription()))?;
        Ok(Self {
            device,
            queue,
            library,
        })
    }

    /// Run the `gradient` kernel: pixel (x, y) = (x / (w-1), y / (h-1), 0.25, 1).
    /// Row-major, `width * height` RGBA pixels.
    pub fn gradient(&self, width: usize, height: usize) -> Result<Vec<[f32; 4]>, String> {
        let pipeline = self.pipeline("gradient")?;
        let pixels = self.output_buffer::<[f32; 4]>(width * height)?;
        // Matches `uint2` in the shader: two packed u32, 8 bytes.
        let size = [width as u32, height as u32];
        self.submit(|encoder| {
            encoder.setComputePipelineState(&pipeline);
            // SAFETY: indices and types match the kernel's [[buffer(n)]]
            // arguments; setBytes copies `size` before returning.
            unsafe {
                encoder.setBuffer_offset_atIndex(Some(&pixels), 0, 0);
                encoder.setBytes_length_atIndex(NonNull::from(&size).cast(), size_of_val(&size), 1);
            }
            dispatch_2d(encoder, &pipeline, width, height);
        })?;
        Ok(read(&pixels, width * height))
    }

    /// Render `scene` with the `render` kernel: row-major RGBA pixels, one
    /// GPU thread each.
    pub fn render(&self, scene: &GpuScene) -> Result<Vec<[f32; 4]>, String> {
        let (width, height) = (scene.width(), scene.height());
        let pipeline = self.pipeline("render")?;
        let pixels = self.output_buffer::<[f32; 4]>(width * height)?;
        let (shapes, nodes, tables) = self.scene_buffers(scene)?;
        self.submit(|encoder| {
            encoder.setComputePipelineState(&pipeline);
            // SAFETY: indices and types match the kernel's [[buffer(n)]]
            // arguments; setBytes copies `info` before returning.
            unsafe {
                encoder.setBuffer_offset_atIndex(Some(&pixels), 0, 0);
                encoder.setBytes_length_atIndex(
                    NonNull::from(&scene.info).cast(),
                    size_of_val(&scene.info),
                    1,
                );
                encoder.setBuffer_offset_atIndex(Some(&shapes), 0, 2);
                encoder.setBuffer_offset_atIndex(Some(&nodes), 0, 3);
                encoder.setBuffer_offset_atIndex(Some(&tables), 0, 4);
            }
            dispatch_2d(encoder, &pipeline, width, height);
        })?;
        Ok(read(&pixels, width * height))
    }

    /// A progressive renderer for `scene`: its shapes and BVH are uploaded
    /// once and the pipelines compiled once, for any number of frames.
    pub fn accumulator(&self, scene: &GpuScene) -> Result<Accumulator, String> {
        let (width, height) = (scene.width(), scene.height());
        let (shapes, nodes, tables) = self.scene_buffers(scene)?;
        Ok(Accumulator {
            accum: self.output_buffer::<[f32; 4]>(width * height)?,
            shapes,
            nodes,
            tables,
            render: self.pipeline("render")?,
            resolve_buffer: self.pipeline("resolve_buffer")?,
            resolve_texture: self.pipeline("resolve_texture")?,
            width: width as u32,
            height: height as u32,
            frames: 0,
        })
    }

    pub fn device(&self) -> &ProtocolObject<dyn MTLDevice> {
        &self.device
    }

    pub fn queue(&self) -> &ProtocolObject<dyn MTLCommandQueue> {
        &self.queue
    }

    /// Run a 1D kernel with one thread per input element:
    ///
    /// `kernel void k(device const I *in [[buffer(0)]], device O *out [[buffer(1)]],
    ///                uint i [[thread_position_in_grid]])`
    ///
    /// `I` and `O` must be plain-old-data `#[repr(C)]` types matching the MSL
    /// structs byte for byte (see `types.rs`).
    pub fn run<I: Copy, O: Copy>(&self, kernel: &str, input: &[I]) -> Result<Vec<O>, String> {
        self.run_1d(kernel, input, None)
    }

    /// `run` for kernels that also read the scene, as
    /// `device const Shape *shapes [[buffer(2)]], device const BvhNode *nodes
    /// [[buffer(3)]], constant SceneCounts &counts [[buffer(4)]], device const
    /// uint *tables [[buffer(5)]]`.
    pub fn run_with_scene<I: Copy, O: Copy>(
        &self,
        kernel: &str,
        input: &[I],
        scene: &GpuScene,
    ) -> Result<Vec<O>, String> {
        let buffers = self.scene_buffers(scene)?;
        self.run_1d(kernel, input, Some((&buffers, scene)))
    }

    /// The shape, BVH node and table buffers of `scene`. Metal cannot make an
    /// empty buffer, so an empty array gets one unused element; the counts in
    /// `scene.info` stay 0 and the shaders never read it.
    fn scene_buffers(&self, scene: &GpuScene) -> Result<(Buffer, Buffer, Buffer), String> {
        let shapes = if scene.shapes.is_empty() {
            self.output_buffer::<GpuShape>(1)?
        } else {
            self.input_buffer(&scene.shapes)?
        };
        let nodes = if scene.nodes.is_empty() {
            self.output_buffer::<GpuBvhNode>(1)?
        } else {
            self.input_buffer(&scene.nodes)?
        };
        let tables = if scene.tables.is_empty() {
            self.output_buffer::<u32>(1)?
        } else {
            self.input_buffer(&scene.tables)?
        };
        Ok((shapes, nodes, tables))
    }

    fn run_1d<I: Copy, O: Copy>(
        &self,
        kernel: &str,
        input: &[I],
        scene: Option<(&(Buffer, Buffer, Buffer), &GpuScene)>,
    ) -> Result<Vec<O>, String> {
        if input.is_empty() {
            return Ok(Vec::new());
        }
        let pipeline = self.pipeline(kernel)?;
        let in_buf = self.input_buffer(input)?;
        let out_buf = self.output_buffer::<O>(input.len())?;
        self.submit(|encoder| {
            encoder.setComputePipelineState(&pipeline);
            // SAFETY: buffer indices match the kernel signatures documented
            // above; setBytes copies the counts before returning.
            unsafe {
                encoder.setBuffer_offset_atIndex(Some(&in_buf), 0, 0);
                encoder.setBuffer_offset_atIndex(Some(&out_buf), 0, 1);
                if let Some(((shapes, nodes, tables), scene)) = scene {
                    let counts = &scene.info.counts;
                    encoder.setBuffer_offset_atIndex(Some(shapes), 0, 2);
                    encoder.setBuffer_offset_atIndex(Some(nodes), 0, 3);
                    encoder.setBytes_length_atIndex(
                        NonNull::from(counts).cast(),
                        size_of_val(counts),
                        4,
                    );
                    encoder.setBuffer_offset_atIndex(Some(tables), 0, 5);
                }
            }
            // The grid is exactly `input.len()` threads, so kernels need no
            // bounds check.
            encoder.dispatchThreads_threadsPerThreadgroup(
                MTLSize {
                    width: input.len(),
                    height: 1,
                    depth: 1,
                },
                MTLSize {
                    width: pipeline.threadExecutionWidth(),
                    height: 1,
                    depth: 1,
                },
            );
        })?;
        Ok(read(&out_buf, input.len()))
    }

    /// Compile `kernel` from the library into a compute pipeline.
    fn pipeline(&self, kernel: &str) -> Result<Pipeline, String> {
        let function = self
            .library
            .newFunctionWithName(&NSString::from_str(kernel))
            .ok_or(format!("kernel `{kernel}` not found in library"))?;
        self.device
            .newComputePipelineStateWithFunction_error(&function)
            .map_err(|e| format!("pipeline `{kernel}` failed: {}", e.localizedDescription()))
    }

    /// A buffer holding a copy of `data`, which must not be empty.
    pub(crate) fn input_buffer<T: Copy>(&self, data: &[T]) -> Result<Buffer, String> {
        // SAFETY: `data` is a live slice of exactly this many bytes; Metal
        // copies it into the new buffer before returning.
        unsafe {
            self.device.newBufferWithBytes_length_options(
                NonNull::from(data).cast(),
                size_of_val(data),
                MTLResourceOptions::StorageModeShared,
            )
        }
        .ok_or("could not allocate input buffer".to_string())
    }

    /// Room for `count` values of `T`. Shared storage: CPU and GPU see the same
    /// memory (unified memory on Apple silicon), so reading back needs no copy.
    pub(crate) fn output_buffer<T>(&self, count: usize) -> Result<Buffer, String> {
        self.device
            .newBufferWithLength_options(
                count * size_of::<T>(),
                MTLResourceOptions::StorageModeShared,
            )
            .ok_or("could not allocate output buffer".to_string())
    }

    /// Record commands with `encode`, run them and wait for the GPU to finish.
    fn submit(&self, encode: impl FnOnce(&Encoder)) -> Result<(), String> {
        let commands = self.queue.commandBuffer().ok_or("no command buffer")?;
        let encoder = commands
            .computeCommandEncoder()
            .ok_or("no compute encoder")?;
        encode(&encoder);
        encoder.endEncoding();
        finish(&commands)
    }
}

/// Commit `commands`, wait for the GPU and report a failed command buffer.
pub(crate) fn finish(commands: &ProtocolObject<dyn MTLCommandBuffer>) -> Result<(), String> {
    {
        commands.commit();
        commands.waitUntilCompleted();
        if commands.status() == MTLCommandBufferStatus::Error {
            let reason = commands
                .error()
                .map(|e| e.localizedDescription().to_string());
            return Err(format!("GPU command failed: {reason:?}"));
        }
        Ok(())
    }
}

/// Sums sampled frames of one scene (`render` with `sampling = 1`) in a GPU
/// buffer; the image is the sum divided by the frame count. The f32 sum
/// loses precision at very high sample counts (~0.3% bias at 65k frames).
pub struct Accumulator {
    accum: Buffer,
    shapes: Buffer,
    nodes: Buffer,
    tables: Buffer,
    render: Pipeline,
    resolve_buffer: Pipeline,
    resolve_texture: Pipeline,
    width: u32,
    height: u32,
    frames: u32,
}

impl Accumulator {
    /// Samples accumulated so far.
    pub fn frames(&self) -> u32 {
        self.frames
    }

    /// Start over: the next frame overwrites the buffer.
    pub fn reset(&mut self) {
        self.frames = 0;
    }

    /// Record one sampled frame of `info`'s camera and light, with frame
    /// index `frames()` (frame 0 overwrites the buffer). The camera must keep
    /// the accumulator's size: the kernel writes one pixel per camera pixel.
    /// `info` must come from the scene the accumulator was built from: the
    /// accumulator holds that scene's shape and BVH buffers.
    pub fn encode_frame(&mut self, encoder: &Encoder, info: &GpuSceneInfo) {
        assert!(
            (info.camera.width, info.camera.height) == (self.width, self.height),
            "camera size differs from the accumulator's"
        );
        let info = GpuSceneInfo {
            frame: self.frames,
            sampling: 1,
            ..*info
        };
        encoder.setComputePipelineState(&self.render);
        // SAFETY: indices and types match the `render` kernel's
        // [[buffer(n)]] arguments; setBytes copies `info` before returning.
        unsafe {
            encoder.setBuffer_offset_atIndex(Some(&self.accum), 0, 0);
            encoder.setBytes_length_atIndex(NonNull::from(&info).cast(), size_of_val(&info), 1);
            encoder.setBuffer_offset_atIndex(Some(&self.shapes), 0, 2);
            encoder.setBuffer_offset_atIndex(Some(&self.nodes), 0, 3);
            encoder.setBuffer_offset_atIndex(Some(&self.tables), 0, 4);
        }
        dispatch_2d(encoder, &self.render, self.width(), self.height());
        self.frames += 1;
    }

    /// Record writing the clamped mean into `texture`, which must be
    /// `width x height` and writable. Needs at least one frame.
    pub fn encode_resolve_texture(
        &self,
        encoder: &Encoder,
        texture: &ProtocolObject<dyn MTLTexture>,
    ) {
        assert!(self.frames > 0, "no samples accumulated");
        let size_frames = self.size_frames();
        encoder.setComputePipelineState(&self.resolve_texture);
        // SAFETY: indices and types match the `resolve_texture` kernel's
        // arguments; setBytes copies `size_frames` before returning; the
        // caller keeps `texture` alive until the command buffer completes.
        unsafe {
            encoder.setBuffer_offset_atIndex(Some(&self.accum), 0, 0);
            encoder.setTexture_atIndex(Some(texture), 0);
            encoder.setBytes_length_atIndex(
                NonNull::from(&size_frames).cast(),
                size_of_val(&size_frames),
                2,
            );
        }
        dispatch_2d(encoder, &self.resolve_texture, self.width(), self.height());
    }

    /// Render one more sampled frame and wait for it.
    pub fn add_frame(&mut self, gpu: &Gpu, info: &GpuSceneInfo) -> Result<(), String> {
        gpu.submit(|encoder| self.encode_frame(encoder, info))
    }

    /// The mean of the frames so far: row-major RGBA pixels, unclamped.
    pub fn resolve(&self, gpu: &Gpu) -> Result<Vec<[f32; 4]>, String> {
        if self.frames == 0 {
            return Err("no samples accumulated".to_string());
        }
        let count = self.width() * self.height();
        let out = gpu.output_buffer::<[f32; 4]>(count)?;
        let size_frames = self.size_frames();
        gpu.submit(|encoder| {
            encoder.setComputePipelineState(&self.resolve_buffer);
            // SAFETY: indices and types match the `resolve_buffer` kernel's
            // arguments; setBytes copies `size_frames` before returning.
            unsafe {
                encoder.setBuffer_offset_atIndex(Some(&self.accum), 0, 0);
                encoder.setBuffer_offset_atIndex(Some(&out), 0, 1);
                encoder.setBytes_length_atIndex(
                    NonNull::from(&size_frames).cast(),
                    size_of_val(&size_frames),
                    2,
                );
            }
            dispatch_2d(encoder, &self.resolve_buffer, self.width(), self.height());
        })?;
        Ok(read(&out, count))
    }

    fn width(&self) -> usize {
        self.width as usize
    }

    fn height(&self) -> usize {
        self.height as usize
    }

    /// MSL `uint3`, which is 16 bytes: (width, height, frames, unused).
    fn size_frames(&self) -> [u32; 4] {
        [self.width, self.height, self.frames, 0]
    }
}

/// One thread per pixel. Threadgroup = one SIMD-group wide, as many rows as
/// the pipeline allows; dispatchThreads trims the edge groups to the image.
fn dispatch_2d(encoder: &Encoder, pipeline: &Pipeline, width: usize, height: usize) {
    let group_w = pipeline.threadExecutionWidth();
    let group_h = pipeline.maxTotalThreadsPerThreadgroup() / group_w;
    encoder.dispatchThreads_threadsPerThreadgroup(
        MTLSize {
            width,
            height,
            depth: 1,
        },
        MTLSize {
            width: group_w,
            height: group_h,
            depth: 1,
        },
    );
}

/// Copy `count` values of `T` out of a shared buffer the GPU has finished with.
fn read<T: Copy>(buffer: &Buffer, count: usize) -> Vec<T> {
    let data = buffer.contents().cast::<T>();
    // SAFETY: callers allocate `buffer` for `count` values of `T` and only read
    // after `submit` has waited for the GPU.
    unsafe { std::slice::from_raw_parts(data.as_ptr(), count) }.to_vec()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gradient_kernel_matches_cpu_formula() {
        let (w, h) = (5, 3);
        let gpu = Gpu::new().unwrap();
        let pixels = gpu.gradient(w, h).unwrap();
        assert_eq!(pixels.len(), w * h);
        for y in 0..h {
            for x in 0..w {
                let expected = [
                    x as f32 / (w - 1) as f32,
                    y as f32 / (h - 1) as f32,
                    0.25,
                    1.0,
                ];
                assert_eq!(pixels[y * w + x], expected, "pixel ({x}, {y})");
            }
        }
    }
}

#[cfg(test)]
mod accumulator_tests {
    use super::*;
    use crate::compare::chapter_7_scene;

    fn chapter_7() -> GpuScene {
        let (world, camera) = chapter_7_scene(20, 10);
        GpuScene::from_world(&world, &camera).unwrap()
    }

    #[test]
    fn accumulator_reset_overwrites() -> Result<(), String> {
        let (gpu, scene) = (Gpu::new()?, chapter_7());
        let mut acc = gpu.accumulator(&scene)?;
        for _ in 0..3 {
            acc.add_frame(&gpu, &scene.info)?;
        }
        acc.reset();
        assert_eq!(acc.frames(), 0);
        acc.add_frame(&gpu, &scene.info)?;
        assert_eq!(acc.frames(), 1);
        let mut fresh = gpu.accumulator(&scene)?;
        fresh.add_frame(&gpu, &scene.info)?;
        assert!(acc.resolve(&gpu)? == fresh.resolve(&gpu)?);
        Ok(())
    }

    /// Every sample of every pixel is the same colour, so the mean of 4
    /// frames must be that colour exactly: catches an off-by-one divisor.
    #[test]
    fn accumulator_resolve_is_mean_of_frames() -> Result<(), String> {
        use rtc_rust::light::Light;
        use rtc_rust::material::Material;
        use rtc_rust::matrix::Matrix4;
        use rtc_rust::shape::Shape;
        use rtc_rust::shapes::sphere::Sphere;
        use rtc_rust::tuple::Tuple;
        use rtc_rust::world::{World, view_transform};

        let (r, g, b) = (0.5, 0.25, 0.75);
        let mut world = World::new();
        world.light = Some(Light::point_light(
            Tuple::point(1.0, 2.0, -3.0),
            Tuple::color(1.0, 1.0, 1.0),
        ));
        let mut m = Material::new();
        m.color = Tuple::color(r, g, b);
        (m.ambient, m.diffuse, m.specular) = (1.0, 0.0, 0.0);
        (m.reflective, m.transparency) = (0.0, 0.0);
        let mut sky = Sphere::new();
        sky.set_transform(Matrix4::scaling(100.0, 100.0, 100.0));
        sky.set_material(m);
        world.add_shape(Box::new(sky));
        let camera = rtc_rust::camera::Camera::new(8, 6, std::f64::consts::PI / 3.0)
            .with_transform(view_transform(
                Tuple::point(0.0, 0.0, 0.0),
                Tuple::point(0.0, 0.0, 1.0),
                Tuple::vector(0.0, 1.0, 0.0),
            ));
        let scene = GpuScene::from_world(&world, &camera)?;

        let gpu = Gpu::new()?;
        let mut acc = gpu.accumulator(&scene)?;
        for _ in 0..4 {
            acc.add_frame(&gpu, &scene.info)?;
        }
        let expected = [r as f32, g as f32, b as f32, 1.0];
        for (i, pixel) in acc.resolve(&gpu)?.iter().enumerate() {
            for c in 0..4 {
                assert!(
                    (pixel[c] - expected[c]).abs() <= 1e-6,
                    "pixel {i}: {pixel:?} != {expected:?}"
                );
            }
        }
        Ok(())
    }

    #[test]
    fn accumulator_without_frames_errors() -> Result<(), String> {
        let (gpu, scene) = (Gpu::new()?, chapter_7());
        let acc = gpu.accumulator(&scene)?;
        assert_eq!(acc.resolve(&gpu), Err("no samples accumulated".to_string()));
        Ok(())
    }
}
