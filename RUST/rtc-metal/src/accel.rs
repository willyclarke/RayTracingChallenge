//! Metal's ray tracing API (milestone 10): one primitive acceleration
//! structure per scene over the bounded `shapes` entries, for the shaders'
//! `TRAVERSAL_METAL` path. Triangles go in as Metal triangle geometry, every
//! other bounded entry (spheres, cubes, cylinders, cones, CSG roots) as its
//! BVH box, intersected by our own code in the query loop.

use crate::gpu::{Gpu, finish};
use crate::scene::GpuScene;
use objc2::rc::Retained;
use objc2::runtime::ProtocolObject;
use objc2_foundation::NSArray;
use objc2_metal::{
    MTLAccelerationStructure, MTLAccelerationStructureBoundingBoxGeometryDescriptor,
    MTLAccelerationStructureCommandEncoder, MTLAccelerationStructureGeometryDescriptor,
    MTLAccelerationStructureTriangleGeometryDescriptor, MTLAttributeFormat,
    MTLAxisAlignedBoundingBox, MTLBuffer, MTLCommandBuffer, MTLCommandEncoder, MTLCommandQueue,
    MTLDevice, MTLPackedFloat3, MTLPrimitiveAccelerationStructureDescriptor, MTLResourceOptions,
};

pub type Structure = Retained<ProtocolObject<dyn MTLAccelerationStructure>>;

/// A primitive structure holds one geometry type only (Metal reads every
/// geometry descriptor as the first one's type), so triangles and boxes get
/// a structure each and the shader runs one query per structure.
pub struct Accel {
    /// The bounded triangles, opaque: Metal intersects them.
    pub triangles: Structure,
    /// Every other bounded entry's box, not opaque: our code intersects them.
    pub boxes: Structure,
    /// `scene.rt_ids`, or one unused word when empty.
    pub rt_ids: Retained<ProtocolObject<dyn MTLBuffer>>,
}

/// Build `scene`'s structures and wait for them. The shader only reads the
/// per-structure `primitive_id`, mapped through `rt_ids`. An empty side gets
/// one zero-size box at the origin, so the kernels always have both bound;
/// they skip a query when its count (`triangle_count`, `box_count`) is 0.
pub fn build(gpu: &Gpu, scene: &GpuScene) -> Result<Accel, String> {
    let triangles = if scene.rt_triangle_count > 0 {
        let vertices = gpu.input_buffer(&scene.rt_vertices)?;
        let geometry = MTLAccelerationStructureTriangleGeometryDescriptor::descriptor();
        geometry.setVertexBuffer(Some(&vertices));
        geometry.setVertexFormat(MTLAttributeFormat::Float3);
        geometry.setVertexStride(size_of::<[f32; 3]>());
        geometry.setTriangleCount(scene.rt_triangle_count as usize);
        geometry.setOpaque(true);
        primitive(gpu, Retained::into_super(geometry))?
    } else {
        placeholder(gpu)?
    };
    let boxes = if scene.rt_boxes.is_empty() {
        placeholder(gpu)?
    } else {
        boxes(gpu, &scene.rt_boxes)?
    };
    let rt_ids = if scene.rt_ids.is_empty() {
        gpu.output_buffer::<u32>(1)?
    } else {
        gpu.input_buffer(&scene.rt_ids)?
    };
    Ok(Accel {
        triangles,
        boxes,
        rt_ids,
    })
}

fn placeholder(gpu: &Gpu) -> Result<Structure, String> {
    boxes(gpu, &[([0.0; 3], [0.0; 3])])
}

fn boxes(gpu: &Gpu, bounds: &[([f32; 3], [f32; 3])]) -> Result<Structure, String> {
    let packed = |p: [f32; 3]| MTLPackedFloat3 {
        x: p[0],
        y: p[1],
        z: p[2],
    };
    let boxes: Vec<MTLAxisAlignedBoundingBox> = bounds
        .iter()
        .map(|&(min, max)| MTLAxisAlignedBoundingBox {
            min: packed(min),
            max: packed(max),
        })
        .collect();
    let buffer = gpu.input_buffer(&boxes)?;
    let geometry = MTLAccelerationStructureBoundingBoxGeometryDescriptor::descriptor();
    geometry.setBoundingBoxBuffer(Some(&buffer));
    // SAFETY: the buffer holds exactly `boxes.len()` boxes of this stride.
    unsafe { geometry.setBoundingBoxStride(size_of::<MTLAxisAlignedBoundingBox>()) };
    geometry.setBoundingBoxCount(boxes.len());
    geometry.setOpaque(false);
    primitive(gpu, Retained::into_super(geometry))
}

/// Build a primitive structure over one geometry. The geometry's buffers
/// are only referenced; they outlive the build, which this waits for.
fn primitive(
    gpu: &Gpu,
    geometry: Retained<MTLAccelerationStructureGeometryDescriptor>,
) -> Result<Structure, String> {
    let descriptor = MTLPrimitiveAccelerationStructureDescriptor::descriptor();
    descriptor.setGeometryDescriptors(Some(&NSArray::from_retained_slice(&[geometry])));
    let device = gpu.device();
    let sizes = device.accelerationStructureSizesWithDescriptor(&descriptor);
    let structure = device
        .newAccelerationStructureWithSize(sizes.accelerationStructureSize)
        .ok_or("could not allocate acceleration structure")?;
    let scratch = device
        .newBufferWithLength_options(
            sizes.buildScratchBufferSize.max(1),
            MTLResourceOptions::StorageModePrivate,
        )
        .ok_or("could not allocate scratch buffer")?;
    let commands = gpu.queue().commandBuffer().ok_or("no command buffer")?;
    let encoder = commands
        .accelerationStructureCommandEncoder()
        .ok_or("no acceleration structure encoder")?;
    encoder.buildAccelerationStructure_descriptor_scratchBuffer_scratchBufferOffset(
        &structure,
        &descriptor,
        &scratch,
        0,
    );
    encoder.endEncoding();
    finish(&commands)?;
    Ok(structure)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compare::scene_file;
    use rtc_rust::camera::Camera;
    use rtc_rust::light::Light;
    use rtc_rust::shape::Shape;
    use rtc_rust::shapes::plane::Plane;
    use rtc_rust::shapes::triangle::Triangle;
    use rtc_rust::tuple::Tuple;
    use rtc_rust::world::World;

    fn triangle() -> Triangle {
        Triangle::new(
            Tuple::point(0.0, 1.0, 0.0),
            Tuple::point(-1.0, 0.0, 0.0),
            Tuple::point(1.0, 0.0, 0.0),
        )
    }

    fn lit() -> World {
        let mut world = World::new();
        world.light = Some(Light::point_light(
            Tuple::point(-10.0, 10.0, -10.0),
            Tuple::color(1.0, 1.0, 1.0),
        ));
        world
    }

    fn builds(world: &World, camera: &Camera) -> Result<(), String> {
        let scene = GpuScene::from_world(world, camera)?;
        let accel = build(&Gpu::new()?, &scene)?;
        if accel.triangles.size() == 0 || accel.boxes.size() == 0 {
            return Err("empty acceleration structure".into());
        }
        Ok(())
    }

    #[test]
    fn test_chap_15_accel_builds_triangles() -> Result<(), String> {
        let mut world = lit();
        world.add_shape(Box::new(triangle()));
        world.build_bounds();
        builds(&world, &Camera::default())
    }

    #[test]
    fn test_chap_12_accel_builds_boxes() -> Result<(), String> {
        let (world, camera) = scene_file(include_str!("../scenes/chapter12.json"), 40, 20);
        builds(&world, &camera)
    }

    #[test]
    fn test_chap_16_accel_builds_mixed() -> Result<(), String> {
        let (mut world, camera) = scene_file(include_str!("../scenes/chapter16.json"), 40, 25);
        let mut t = triangle();
        t.set_transform(rtc_rust::matrix::Matrix4::translation(0.0, 3.0, 0.0));
        world.add_shape(Box::new(t));
        world.build_bounds();
        builds(&world, &camera)
    }

    #[test]
    fn test_chap_9_accel_builds_planes_only() -> Result<(), String> {
        let mut world = lit();
        world.add_shape(Box::new(Plane::new()));
        builds(&world, &Camera::default())
    }
}
