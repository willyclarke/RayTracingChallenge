//! Light
//!
//! Let there be light
//!

use crate::{material::Material, tuple::Tuple};

#[derive(Debug, Clone, Copy)]
pub struct Light {
    pub position: Tuple,
    pub intensity: Tuple,
}

impl Light {
    pub fn point_light(position: Tuple, intensity: Tuple) -> Self {
        Self {
            position,
            intensity,
        }
    }

    ///
    /// add together the material’s ambient, diffuse, and specular components, weighted by the angles between the different vec- tors.
    ///
    pub fn lighting(&self, material: Material, point: Tuple, eyev: Tuple, normalv: Tuple) -> Tuple {
        // combine the surface color with the light's color/intensity
        let effective_color = material.color.mul(self.intensity);

        // find the direction to the light source
        let lightv = (self.position - point).normalize();

        // compute the ambient contribution
        let ambient = effective_color.mul(material.ambient);

        // light_dot_normal represents the cosine of the angle between the
        // light vector and the normal vector. A negative number means the
        // light is on the other side of the surface.
        let light_dot_normal = lightv.dot(normalv);

        let black = Tuple::color(0.0, 0.0, 0.0);
        let mut diffuse = black;
        let mut specular = black;

        if light_dot_normal >= 0.0 {
            // compute the diffuse contribution
            diffuse = effective_color.mul(material.diffuse.mul(light_dot_normal));

            // reflect_dot_eye represents the cosine of the angle between the
            // reflection vector and the eye vector. A negative number means the
            // light reflects away from the eye.
            let reflectv = -lightv.reflect(normalv);
            let reflect_dot_eye = reflectv.dot(eyev);

            if reflect_dot_eye > 0.0 {
                // compute the specular contribution
                let factor = reflect_dot_eye.powf(material.shininess);
                specular = self.intensity.mul(material.specular).mul(factor);
            }
        }

        ambient + diffuse + specular
    }
}

use std::{fmt, ops::Mul};
impl fmt::Display for Light {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Light Pos:({:.12}, {:.12}, {:.12}) Intensity:({:.12}, {:.12}, {:.12})",
            self.position.x,
            self.position.y,
            self.position.z,
            self.intensity.x,
            self.intensity.y,
            self.intensity.z
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{light::Light, loge, tuple::Tuple};

    /// Chap 6 - Lighting with the eye between the light and the surface
    #[test]
    fn test_chap_6_13() -> Result<(), String> {
        let m = Material::new();
        let position = Tuple::point(0.0, 0.0, 0.0);
        let eyev = Tuple::vector(0.0, 0.0, -1.0);
        let normalv = Tuple::vector(0.0, 0.0, -1.0);
        let light = Light::point_light(Tuple::point(0.0, 0.0, -10.0), Tuple::color(1.0, 1.0, 1.0));
        let result = light.lighting(m, position, eyev, normalv);
        let chk = Tuple::color(1.9, 1.9, 1.9).approx_eq(result);
        if chk {
            Ok(())
        } else {
            Err("Lighting with the eye between the light and the surface".into())
        }
    }

    /// Chap 6 - Lighting with the eye between light and surface, eye offset 45°
    #[test]
    fn test_chap_6_14() -> Result<(), String> {
        let sqrt2_o_2 = 2_f64.sqrt() / 2.0;
        let m = Material::new();
        let position = Tuple::point(0.0, 0.0, 0.0);
        let eyev = Tuple::vector(0.0, sqrt2_o_2, -sqrt2_o_2);
        let normalv = Tuple::vector(0.0, 0.0, -1.0);
        let light = Light::point_light(Tuple::point(0.0, 0.0, -10.0), Tuple::color(1.0, 1.0, 1.0));
        let result = light.lighting(m, position, eyev, normalv);
        let chk = Tuple::color(1.0, 1.0, 1.0).approx_eq(result);
        if chk {
            Ok(())
        } else {
            Err("Lighting with the eye between light and surface, eye offset 45°".into())
        }
    }

    /// Chap 6 - Lighting with eye opposite surface, light offset 45°
    #[test]
    fn test_chap_6_15() -> Result<(), String> {
        let m = Material::new();
        let position = Tuple::point(0.0, 0.0, 0.0);
        let eyev = Tuple::vector(0.0, 0.0, -1.0);
        let normalv = Tuple::vector(0.0, 0.0, -1.0);
        let light = Light::point_light(Tuple::point(0.0, 10.0, -10.0), Tuple::color(1.0, 1.0, 1.0));
        let result = light.lighting(m, position, eyev, normalv);
        let chk = Tuple::color(0.73639610306, 0.73639610306, 0.73639610306).approx_eq(result);
        if chk {
            Ok(())
        } else {
            loge!("test_chap_6_15", "result:{}", result);
            Err("Lighting with eye opposite surface, light offset 45°".into())
        }
    }

    /// Chap 6 -Lighting with eye in the path of the reflection vector
    #[test]
    fn test_chap_6_16() -> Result<(), String> {
        let sqrt2_o_2 = 2_f64.sqrt() / 2.0;
        let m = Material::new();
        let position = Tuple::point(0.0, 0.0, 0.0);
        let eyev = Tuple::vector(0.0, -sqrt2_o_2, -sqrt2_o_2);
        let normalv = Tuple::vector(0.0, 0.0, -1.0);
        let light = Light::point_light(Tuple::point(0.0, 10.0, -10.0), Tuple::color(1.0, 1.0, 1.0));
        let result = light.lighting(m, position, eyev, normalv);
        let chk = Tuple::color(1.636396103068, 1.636396103068, 1.636396103068).approx_eq(result);

        if chk {
            Ok(())
        } else {
            loge!("test_chap_6_16", "result:{}", result);
            Err("Lighting with eye in the path of the reflection vector".into())
        }
    }

    /// Chap 6 -Lighting with the light behind the surface
    #[test]
    fn test_chap_6_17() -> Result<(), String> {
        let m = Material::new();
        let position = Tuple::point(0.0, 0.0, 0.0);
        let eyev = Tuple::vector(0.0, 0.0, -1.0);
        let normalv = Tuple::vector(0.0, 0.0, -1.0);
        let light = Light::point_light(Tuple::point(0.0, 0.0, 10.0), Tuple::color(1.0, 1.0, 1.0));
        let result = light.lighting(m, position, eyev, normalv);
        let chk = Tuple::color(0.1, 0.1, 0.1).approx_eq(result);

        if chk {
            Ok(())
        } else {
            loge!("test_chap_6_17", "result:{}", result);
            Err("Lighting with the light behind the surface".into())
        }
    }
}
