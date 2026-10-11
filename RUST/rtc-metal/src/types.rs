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

/// MSL `Material`: rtc-rust's `Material` without patterns and bump maps.
#[repr(C, align(16))]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpuMaterial {
    pub color: Float4,
    pub ambient: f32,
    pub diffuse: f32,
    pub specular: f32,
    pub shininess: f32,
    pub reflective: f32,
    pub transparency: f32,
    pub refractive_index: f32,
}

/// `GpuShape::kind` values; must match the `SHAPE_*` constants in rt.metal.
pub const SHAPE_SPHERE: u32 = 0;
pub const SHAPE_PLANE: u32 = 1;

/// MSL `Shape`: a primitive (`kind`) seen through its world-to-object
/// transform. One tagged array holds every kind, like rtc-rust's `shapes`.
#[repr(C, align(16))]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpuShape {
    pub inverse: GpuMat4,
    pub material: GpuMaterial,
    pub kind: u32,
}

/// MSL `Light`: a point light.
#[repr(C, align(16))]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpuLight {
    pub position: Float4,
    pub intensity: Float4,
}

/// MSL `Camera`: what `ray_for_pixel` needs.
#[repr(C, align(16))]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpuCamera {
    pub inverse: GpuMat4,
    pub half_width: f32,
    pub half_height: f32,
    pub pixel_size: f32,
    pub width: u32,
    pub height: u32,
}

/// MSL `SceneInfo`: the per-render constants, passed with `setBytes`.
#[repr(C, align(16))]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpuSceneInfo {
    pub camera: GpuCamera,
    pub light: GpuLight,
    pub shape_count: u32,
}

// 16 + 7 * 4 = 44 bytes of fields, padded to the 16-byte alignment.
const _: () = assert!(size_of::<GpuMaterial>() == 48);
// 64 + 48 + 4 = 116 bytes of fields, padded to the 16-byte alignment.
const _: () = assert!(size_of::<GpuShape>() == 128);
const _: () = assert!(size_of::<GpuLight>() == 32);
// 64 + 5 * 4 = 84 bytes of fields, padded to the 16-byte alignment.
const _: () = assert!(size_of::<GpuCamera>() == 96);
const _: () = assert!(size_of::<GpuSceneInfo>() == 144);
