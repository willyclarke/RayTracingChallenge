//! Material
//!
//!

use crate::log::*;
use crate::math::approx_eq;
use crate::noise::Perlin;
use crate::pattern::Pattern;
use crate::tuple::Tuple;
use std::fmt;

/// Normal perturbation (book chapter 17): the object-space normal is tilted
/// by Perlin noise sampled at the object-space point, giving a bumpy or
/// rippled surface without changing the geometry.
#[derive(Debug, Clone)]
pub struct Bump {
    pub noise: Perlin,
    /// How far the normal is tilted, relative to its unit length.
    pub amplitude: f64,
    /// Spatial frequency of the bumps.
    pub frequency: f64,
}

impl Bump {
    pub fn new(amplitude: f64, frequency: f64) -> Self {
        Self {
            noise: Perlin::default(),
            amplitude,
            frequency,
        }
    }

    pub fn with_seed(mut self, seed: u64) -> Self {
        self.noise = Perlin::new(seed);
        self
    }

    /// Tilt `normal` by the noise at `point`; the result is unit length.
    pub fn perturb(&self, point: Tuple, normal: Tuple) -> Tuple {
        (normal + self.noise.vector(point * self.frequency) * self.amplitude).normalize()
    }
}

#[derive(Debug, Clone)]
pub struct Material {
    pub color: Tuple,
    pub ambient: f64,
    pub diffuse: f64,
    pub specular: f64,
    pub shininess: f64,
    pub reflective: f64,
    pub transparency: f64,
    pub refractive_index: f64,
    pub pattern: Option<Box<dyn Pattern>>,
    pub bump: Option<Bump>,
}

impl Material {
    pub fn new() -> Self {
        Self {
            color: Tuple::color(1.0, 1.0, 1.0),
            ambient: 0.1,
            diffuse: 0.9,
            specular: 0.9,
            shininess: 200.0,
            reflective: 0.0,
            transparency: 0.0,
            refractive_index: 1.0,
            pattern: None,
            bump: None,
        }
    }

    pub fn approx_eq(&self, other: &Material) -> bool {
        self.color.approx_eq(other.color)
            && approx_eq(self.ambient, other.ambient)
            && approx_eq(self.diffuse, other.diffuse)
            && approx_eq(self.specular, other.specular)
            && approx_eq(self.shininess, other.shininess)
            && approx_eq(self.reflective, other.reflective)
            && approx_eq(self.transparency, other.transparency)
            && approx_eq(self.refractive_index, other.refractive_index)
    }
}

impl Default for Material {
    fn default() -> Self {
        Self::new()
    }
}

impl fmt::Display for Material {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "{} Material:{} \nColor: {}{:.3}{} {}{:.3}{} {}{:.3}{} {}\n{:16}:{:.3} \n{:16}:{:.3} \n{:16}:{:.3} \n{:16}:{:.3} \n{:16}:{:.3} \n{:16}:{:.3} \n{:16}:{:.3}{}",
            Color::Yellow,
            Color::Reset,
            Color::Red,
            self.color.x,
            Color::Reset,
            Color::Green,
            self.color.y,
            Color::Reset,
            Color::Blue,
            self.color.z,
            Color::Reset,
            Color::Yellow,
            ".ambient         ",
            self.ambient,
            ".diffuse         ",
            self.diffuse,
            ".specular        ",
            self.specular,
            ".shininess       ",
            self.shininess,
            ".reflective      ",
            self.reflective,
            ".transparency    ",
            self.transparency,
            ".refractive_index",
            self.refractive_index,
            Color::Reset,
        )?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{logd, loge, logi, tuple::Tuple};

    /// Chap 6 - The default material
    #[test]
    fn test_chap_6_11() -> Result<(), String> {
        let m = Material::new();
        logi!("test_chap_6_11", "m:{}", m);

        let chk = Tuple::color(1.0, 1.0, 1.0).approx_eq(m.color);
        let chk = chk && approx_eq(m.ambient, 0.1);
        let chk = chk && approx_eq(m.diffuse, 0.9);
        let chk = chk && approx_eq(m.specular, 0.9);
        let chk = chk && approx_eq(m.shininess, 200.0);
        let chk = chk && m.approx_eq(&Material::new());
        if chk {
            Ok(())
        } else {
            logi!("test_chap_6_11", "m:{}", m);
            logd!("test_chap_6_11", "m:{:?}", m);
            loge!("test_chap_6_11", "chk:{:?}", chk);
            Err("The default material".into())
        }
    }

    /// Chap 11 - Reflectivity for the default material
    #[test]
    fn test_chap_11_1() -> Result<(), String> {
        let m = Material::new();

        let chk = approx_eq(m.reflective, 0.0);
        if chk {
            Ok(())
        } else {
            Err("Reflectivity for the default material".into())
        }
    }
}
