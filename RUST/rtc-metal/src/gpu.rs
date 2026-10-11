//! Metal host side: device, command queue and the compiled shader library.

use crate::scene::GpuScene;
use crate::types::GpuSphere;
use objc2::rc::Retained;
use objc2::runtime::ProtocolObject;
use objc2_foundation::NSString;
use objc2_metal::{
    MTLBuffer, MTLCommandBuffer, MTLCommandBufferStatus, MTLCommandEncoder, MTLCommandQueue,
    MTLComputeCommandEncoder, MTLComputePipelineState, MTLCreateSystemDefaultDevice, MTLDevice,
    MTLLibrary, MTLResourceOptions, MTLSize,
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
type Encoder = ProtocolObject<dyn MTLComputeCommandEncoder>;

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
        // Metal cannot make an empty buffer; an empty scene gets one unused
        // sphere (sphere_count stays 0).
        let spheres = if scene.spheres.is_empty() {
            self.output_buffer::<GpuSphere>(1)?
        } else {
            self.input_buffer(&scene.spheres)?
        };
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
                encoder.setBuffer_offset_atIndex(Some(&spheres), 0, 2);
            }
            dispatch_2d(encoder, &pipeline, width, height);
        })?;
        Ok(read(&pixels, width * height))
    }

    /// Run a 1D kernel with one thread per input element:
    ///
    /// `kernel void k(device const I *in [[buffer(0)]], device O *out [[buffer(1)]],
    ///                uint i [[thread_position_in_grid]])`
    ///
    /// `I` and `O` must be plain-old-data `#[repr(C)]` types matching the MSL
    /// structs byte for byte (see `types.rs`).
    pub fn run<I: Copy, O: Copy>(&self, kernel: &str, input: &[I]) -> Result<Vec<O>, String> {
        if input.is_empty() {
            return Ok(Vec::new());
        }
        let pipeline = self.pipeline(kernel)?;
        let in_buf = self.input_buffer(input)?;
        let out_buf = self.output_buffer::<O>(input.len())?;
        self.submit(|encoder| {
            encoder.setComputePipelineState(&pipeline);
            // SAFETY: buffer indices match the kernel signature documented above.
            unsafe {
                encoder.setBuffer_offset_atIndex(Some(&in_buf), 0, 0);
                encoder.setBuffer_offset_atIndex(Some(&out_buf), 0, 1);
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
    fn input_buffer<T: Copy>(&self, data: &[T]) -> Result<Buffer, String> {
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
    fn output_buffer<T>(&self, count: usize) -> Result<Buffer, String> {
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
