//! Flatten an rtc-rust `World` + `Camera` into GPU buffers.
//!
//! Only what the shaders implement so far is accepted; anything else is an
//! error rather than being silently dropped from the image.

use crate::types::{
    GpuCamera, GpuLight, GpuMat4, GpuMaterial, GpuSceneInfo, GpuShape, SHAPE_PLANE, SHAPE_SPHERE,
    float4,
};
use rtc_rust::camera::Camera;
use rtc_rust::light::Light;
use rtc_rust::material::Material;
use rtc_rust::shape::Shape;
use rtc_rust::shapes::plane::Plane;
use rtc_rust::shapes::sphere::Sphere;
use rtc_rust::world::World;
use std::any::Any;

pub struct GpuScene {
    pub info: GpuSceneInfo,
    pub shapes: Vec<GpuShape>,
}

impl GpuScene {
    pub fn from_world(world: &World, camera: &Camera) -> Result<Self, String> {
        if world.gi_samples > 0 {
            return Err("path tracing is not supported yet".into());
        }
        let light = world.light.as_ref().ok_or("scene has no light")?;
        let shapes = world
            .shapes
            .iter()
            .map(|s| gpu_shape(s.as_ref()))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            info: GpuSceneInfo {
                camera: gpu_camera(camera)?,
                light: gpu_light(light)?,
                shape_count: shapes.len() as u32,
            },
            shapes,
        })
    }

    pub fn width(&self) -> usize {
        self.info.camera.width as usize
    }

    pub fn height(&self) -> usize {
        self.info.camera.height as usize
    }
}

pub fn gpu_camera(camera: &Camera) -> Result<GpuCamera, String> {
    if camera.aa_samples > 1 || camera.aperture > 0.0 || camera.motion_samples > 1 {
        return Err("anti-aliasing, focal blur and motion blur are not supported yet".into());
    }
    let inverse = camera
        .transform()
        .inverse()
        .ok_or("camera transform is not invertible")?;
    Ok(GpuCamera {
        inverse: (&inverse).into(),
        half_width: camera.half_width as f32,
        half_height: camera.half_height as f32,
        pixel_size: camera.pixel_size as f32,
        width: camera.hsize as u32,
        height: camera.vsize as u32,
    })
}

pub fn gpu_light(light: &Light) -> Result<GpuLight, String> {
    if light.samples > 1 || light.spot.is_some() {
        return Err("only point lights are supported so far".into());
    }
    Ok(GpuLight {
        position: float4(light.position),
        intensity: float4(light.intensity),
    })
}

pub fn gpu_material(m: &Material) -> Result<GpuMaterial, String> {
    if m.pattern.is_some() || m.bump.is_some() || m.reflective > 0.0 || m.transparency > 0.0 {
        return Err("patterns, bump maps, reflection and refraction are not supported yet".into());
    }
    Ok(GpuMaterial {
        color: float4(m.color),
        ambient: m.ambient as f32,
        diffuse: m.diffuse as f32,
        specular: m.specular as f32,
        shininess: m.shininess as f32,
    })
}

fn gpu_shape(shape: &dyn Shape) -> Result<GpuShape, String> {
    // Upcast to `&dyn Any` to learn the concrete type. This must be the
    // `&dyn Shape`, not the `Box`: a `&Box<dyn Shape>` would coerce to an
    // `Any` whose type is the Box itself, and `is::<Sphere>` would be false.
    let any: &dyn Any = shape;
    let kind = if any.is::<Sphere>() {
        SHAPE_SPHERE
    } else if any.is::<Plane>() {
        SHAPE_PLANE
    } else {
        return Err(format!(
            "shape {}: only spheres and planes are supported so far",
            shape.id()
        ));
    };
    let data = shape.data();
    if data.parent.is_some() || data.motion.is_some() {
        return Err(format!(
            "shape {}: groups and motion are not supported yet",
            shape.id()
        ));
    }
    let inverse: GpuMat4 = (&data.transform_inv).into();
    Ok(GpuShape {
        inverse,
        material: gpu_material(&data.material)?,
        kind,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use rtc_rust::shapes::cube::Cube;
    use rtc_rust::tuple::Tuple;

    #[test]
    fn unsupported_shapes_are_rejected() {
        let mut world = World::new();
        world.light = Some(Light::point_light(
            Tuple::point(0.0, 0.0, 0.0),
            Tuple::color(1.0, 1.0, 1.0),
        ));
        world.add_shape(Box::new(Sphere::new()));
        world.add_shape(Box::new(Plane::new()));
        assert!(GpuScene::from_world(&world, &Camera::default()).is_ok());
        world.add_shape(Box::new(Cube::new()));
        let err = GpuScene::from_world(&world, &Camera::default()).err();
        assert!(err.is_some_and(|e| e.contains("only spheres and planes")));
    }
}
