//! `#[repr(C)]` mirrors of the structs in `shaders/rt.metal`.
//!
//! MSL aligns `float4` and `float4x4` to 16 bytes, so these do too. The size
//! asserts here and the `static_assert`s in the shaders keep both sides in step.

use rtc_rust::matrix::Matrix4;
use rtc_rust::ray::Ray;
use rtc_rust::tuple::Tuple;

/// MSL `float4`. Points have w = 1, vectors w = 0, as in rtc-rust's `Tuple`.
pub type Float4 = [f32; 4];

pub fn float4(t: Tuple) -> Float4 {
    [t.x as f32, t.y as f32, t.z as f32, t.w as f32]
}

/// MSL `Ray`.
#[repr(C, align(16))]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpuRay {
    pub origin: Float4,
    pub direction: Float4,
}

impl From<&Ray> for GpuRay {
    fn from(r: &Ray) -> Self {
        Self {
            origin: float4(r.origin),
            direction: float4(r.direction),
        }
    }
}

/// MSL `float4x4`, which is column-major: `cols[c][r]` is row `r`, column `c`.
/// rtc-rust's `Matrix4` is indexed `(row, col)`, so conversion transposes the
/// storage order (not the matrix).
#[repr(C, align(16))]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpuMat4 {
    pub cols: [Float4; 4],
}

impl From<&Matrix4> for GpuMat4 {
    fn from(m: &Matrix4) -> Self {
        let col = |c: usize| [0, 1, 2, 3].map(|r| m[(r, c)] as f32);
        Self {
            cols: [col(0), col(1), col(2), col(3)],
        }
    }
}

const _: () = assert!(size_of::<GpuRay>() == 32);
const _: () = assert!(size_of::<GpuMat4>() == 64);
