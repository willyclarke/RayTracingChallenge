//! Metal host side: device, command queue and one compute pipeline per kernel.

use objc2::rc::Retained;
use objc2::runtime::ProtocolObject;
use objc2_foundation::NSString;
use objc2_metal::{
    MTLBuffer, MTLCommandBuffer, MTLCommandBufferStatus, MTLCommandEncoder, MTLCommandQueue,
    MTLComputeCommandEncoder, MTLComputePipelineState, MTLCreateSystemDefaultDevice, MTLDevice,
    MTLLibrary, MTLResourceOptions, MTLSize,
};
use std::ptr::NonNull;

/// MSL source, compiled by the Metal driver when `Gpu::new` runs.
const GRADIENT_SRC: &str = include_str!("../shaders/gradient.metal");

/// Owns the Metal objects that live for the whole program.
pub struct Gpu {
    device: Retained<ProtocolObject<dyn MTLDevice>>,
    queue: Retained<ProtocolObject<dyn MTLCommandQueue>>,
    gradient: Retained<ProtocolObject<dyn MTLComputePipelineState>>,
}

impl Gpu {
    pub fn new() -> Result<Self, String> {
        let device = MTLCreateSystemDefaultDevice().ok_or("no Metal device")?;
        let queue = device
            .newCommandQueue()
            .ok_or("could not create command queue")?;
        let library = device
            .newLibraryWithSource_options_error(&NSString::from_str(GRADIENT_SRC), None)
            .map_err(|e| format!("shader compile failed: {}", e.localizedDescription()))?;
        let function = library
            .newFunctionWithName(&NSString::from_str("gradient"))
            .ok_or("kernel `gradient` not found in library")?;
        let gradient = device
            .newComputePipelineStateWithFunction_error(&function)
            .map_err(|e| format!("pipeline creation failed: {}", e.localizedDescription()))?;
        Ok(Self {
            device,
            queue,
            gradient,
        })
    }

    /// Run the `gradient` kernel: pixel (x, y) = (x / (w-1), y / (h-1), 0.25, 1).
    /// Row-major, `width * height` RGBA pixels.
    pub fn gradient(&self, width: usize, height: usize) -> Result<Vec<[f32; 4]>, String> {
        let count = width * height;
        // Shared storage: CPU and GPU see the same memory (unified memory on
        // Apple silicon), so reading the result back needs no copy.
        let pixels = self
            .device
            .newBufferWithLength_options(
                count * size_of::<[f32; 4]>(),
                MTLResourceOptions::StorageModeShared,
            )
            .ok_or("could not allocate pixel buffer")?;
        // Matches `uint2` in the shader: two packed u32, 8 bytes.
        let size = [width as u32, height as u32];

        let commands = self.queue.commandBuffer().ok_or("no command buffer")?;
        let encoder = commands
            .computeCommandEncoder()
            .ok_or("no compute encoder")?;
        encoder.setComputePipelineState(&self.gradient);
        // SAFETY: indices and types match the kernel's [[buffer(n)]] arguments;
        // setBytes copies `size`, so the borrow only has to outlive this call.
        unsafe {
            encoder.setBuffer_offset_atIndex(Some(&pixels), 0, 0);
            encoder.setBytes_length_atIndex(NonNull::from(&size).cast(), size_of_val(&size), 1);
        }
        // Threadgroup = one SIMD-group wide, as many rows as the pipeline allows.
        // dispatchThreads lets Metal trim the edge groups to the grid size.
        let group_w = self.gradient.threadExecutionWidth();
        let group_h = self.gradient.maxTotalThreadsPerThreadgroup() / group_w;
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
        encoder.endEncoding();
        commands.commit();
        commands.waitUntilCompleted();
        if commands.status() == MTLCommandBufferStatus::Error {
            let reason = commands
                .error()
                .map(|e| e.localizedDescription().to_string());
            return Err(format!("GPU command failed: {reason:?}"));
        }

        // SAFETY: the GPU has finished (waitUntilCompleted) and the buffer holds
        // exactly `count` float4 values.
        let data = pixels.contents().cast::<[f32; 4]>();
        Ok(unsafe { std::slice::from_raw_parts(data.as_ptr(), count) }.to_vec())
    }
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
