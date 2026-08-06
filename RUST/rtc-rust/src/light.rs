//! Light
//!
//! Let there be light
//!

use crate::{shape::Shape, tuple::Tuple};

#[derive(Debug, Clone, Copy)]
pub struct Light {
    pub position: Tuple,
    pub intensity: Tuple,
}

impl Light {
    pub fn approx_eq(&self, other: Light) -> bool {
        self.position.approx_eq(other.position) && self.intensity.approx_eq(other.intensity)
    }

    pub fn point_light(position: Tuple, intensity: Tuple) -> Self {
        Self {
            position,
            intensity,
        }
    }

    ///
    /// Add together the material’s ambient, diffuse, and specular components, weighted by the angles between the different vec- tors.
    ///
    pub fn lighting(
        &self,
        shape: &dyn Shape,
        point: Tuple,        // WORLD point — still needed below
        object_point: Tuple, // NEW — group-aware object-space point, for the pattern
        eyev: Tuple,
        normalv: Tuple,
        in_shadow: bool,
    ) -> Tuple {
        let material = shape.material();

        // Check if pattern is borrowed and use that color as input, othewise use the material color.
        let color = match &material.pattern {
            Some(pattern) => pattern.color_at_local(object_point),
            None => material.color,
        };

        // combine the surface color with the light's color/intensity
        let effective_color = color.mul(self.intensity);

        // find the direction to the light source
        let lightv = (self.position - point).normalize(); // ← still uses the WORLD point

        // compute the ambient contribution
        let ambient = effective_color.mul(material.ambient);

        if in_shadow {
            return ambient;
        }

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

    ///
    /// Add together the material’s ambient, diffuse, and specular components, weighted by the angles between the different vec- tors.
    ///
    pub fn lighting_old(
        &self,
        shape: &dyn Shape,
        point: Tuple,
        eyev: Tuple,
        normalv: Tuple,
        in_shadow: bool,
    ) -> Tuple {
        let material = shape.material();

        // Check if pattern is borrowed and use that color as input, othewise use the material color.
        let color = match &material.pattern {
            Some(pattern) => pattern.color_at_shape(shape, point),
            None => material.color,
        };

        // combine the surface color with the light's color/intensity
        let effective_color = color.mul(self.intensity);

        // find the direction to the light source
        let lightv = (self.position - point).normalize();

        // compute the ambient contribution
        let ambient = effective_color.mul(material.ambient);

        if in_shadow {
            return ambient;
        }

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
    use crate::material::Material;
    use crate::patterns::stripepattern::*;
    use crate::shapes::sphere::Sphere;
    use crate::tuple::colors::*;
    use crate::{light::Light, loge, tuple::Tuple};

    /// Chap 6 - Lighting with the eye between the light and the surface
    #[test]
    fn test_chap_6_13() -> Result<(), String> {
        let m = Material::new();
        let mut s = Sphere::new();
        s.set_material(m);
        let position = Tuple::point(0.0, 0.0, 0.0);
        let eyev = Tuple::vector(0.0, 0.0, -1.0);
        let normalv = Tuple::vector(0.0, 0.0, -1.0);
        let light = Light::point_light(Tuple::point(0.0, 0.0, -10.0), Tuple::color(1.0, 1.0, 1.0));
        let in_shadow = false;
        let result = light.lighting(&s, position, position, eyev, normalv, in_shadow);
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
        let mut s = Sphere::new();
        s.set_material(m);
        let position = Tuple::point(0.0, 0.0, 0.0);
        let eyev = Tuple::vector(0.0, sqrt2_o_2, -sqrt2_o_2);
        let normalv = Tuple::vector(0.0, 0.0, -1.0);
        let light = Light::point_light(Tuple::point(0.0, 0.0, -10.0), Tuple::color(1.0, 1.0, 1.0));
        let in_shadow = false;
        let result = light.lighting(&s, position, position, eyev, normalv, in_shadow);
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
        let mut s = Sphere::new();
        s.set_material(m);
        let position = Tuple::point(0.0, 0.0, 0.0);
        let eyev = Tuple::vector(0.0, 0.0, -1.0);
        let normalv = Tuple::vector(0.0, 0.0, -1.0);
        let light = Light::point_light(Tuple::point(0.0, 10.0, -10.0), Tuple::color(1.0, 1.0, 1.0));
        let in_shadow = false;
        let result = light.lighting(&s, position, position, eyev, normalv, in_shadow);
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
        let mut s = Sphere::new();
        s.set_material(m);
        let position = Tuple::point(0.0, 0.0, 0.0);
        let eyev = Tuple::vector(0.0, -sqrt2_o_2, -sqrt2_o_2);
        let normalv = Tuple::vector(0.0, 0.0, -1.0);
        let light = Light::point_light(Tuple::point(0.0, 10.0, -10.0), Tuple::color(1.0, 1.0, 1.0));
        let in_shadow = false;
        let result = light.lighting(&s, position, position, eyev, normalv, in_shadow);
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
        let mut s = Sphere::new();
        s.set_material(m);
        let position = Tuple::point(0.0, 0.0, 0.0);
        let eyev = Tuple::vector(0.0, 0.0, -1.0);
        let normalv = Tuple::vector(0.0, 0.0, -1.0);
        let light = Light::point_light(Tuple::point(0.0, 0.0, 10.0), Tuple::color(1.0, 1.0, 1.0));
        let in_shadow = false;
        let result = light.lighting(&s, position, position, eyev, normalv, in_shadow);
        let chk = Tuple::color(0.1, 0.1, 0.1).approx_eq(result);

        if chk {
            Ok(())
        } else {
            loge!("test_chap_6_17", "result:{}", result);
            Err("Lighting with the light behind the surface".into())
        }
    }

    /// Chap 8 - Lighting with the surface in shadow
    #[test]
    fn test_chap_8_1() -> Result<(), String> {
        let eyev = Tuple::vector(0.0, 0.0, -1.0);
        let normalv = Tuple::vector(0.0, 0.0, -1.0);
        let position = Tuple::point(0.0, 0.1, 10.0);
        let intensity = Tuple::color(1.0, 1.0, 1.0);
        let light = Light::point_light(position, intensity);
        let in_shadow = true;
        let m = Material::new();
        let mut s = Sphere::new();
        s.set_material(m);
        let result = light.lighting(&s, position, position, eyev, normalv, in_shadow);
        let chk = result.approx_eq(Tuple::color(0.1, 0.1, 0.1));
        if chk {
            Ok(())
        } else {
            Err("Lighting with the surface in shadow".into())
        }
    }

    /// Chap 10 - Lighting with a pattern applied
    #[test]
    fn test_chap_10_5() -> Result<(), String> {
        let mut m = Material::new();
        m.pattern = Some(Box::new(StripePattern::new(WHITE, BLACK)));
        m.ambient = 1.0;
        m.diffuse = 0.0;
        m.specular = 0.0;

        let eyev = Tuple::vector(0.0, 0.0, -1.0);
        let normalv = Tuple::vector(0.0, 0.0, -1.0);
        let position_light = Tuple::point(0.0, 0.0, -10.0);
        let intensity = Tuple::color(1.0, 1.0, 1.0);
        let light = Light::point_light(position_light, intensity);
        let in_shadow = false;
        let position_c1 = Tuple::point(0.9, 0.0, 0.0);
        let position_c2 = Tuple::point(1.1, 0.0, 0.0);
        let mut s = Sphere::new();
        s.set_material(m);
        let c1 = light.lighting(&s, position_c1, position_c1, eyev, normalv, in_shadow);
        let c2 = light.lighting(&s, position_c2, position_c2, eyev, normalv, in_shadow);

        let chk = c1.approx_eq(WHITE);
        let chk = chk && c2.approx_eq(BLACK);
        if chk {
            Ok(())
        } else {
            loge!("test_chap_10_5", "color c1:{}", c1);
            loge!("test_chap_10_5", "color c2:{}", c2);
            Err("Lighting with a pattern applied".into())
        }
    }
}
