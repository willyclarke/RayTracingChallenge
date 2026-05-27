//! Material
//!
//!

use crate::log::*;
use crate::math::approx_eq;
use crate::tuple::Tuple;
use std::fmt;

#[derive(Debug, Clone, Copy)]
pub struct Material {
    pub color: Tuple,
    pub ambient: f64,
    pub diffuse: f64,
    pub specular: f64,
    pub shininess: f64,
}

impl Material {
    pub fn new() -> Self {
        Self {
            color: Tuple::color(1.0, 1.0, 1.0),
            ambient: 0.1,
            diffuse: 0.9,
            specular: 0.9,
            shininess: 200.0,
        }
    }

    pub fn approx_eq(&self, other: Material) -> bool {
        self.color.approx_eq(other.color)
            && approx_eq(self.ambient, other.ambient)
            && approx_eq(self.diffuse, other.diffuse)
            && approx_eq(self.specular, other.specular)
            && approx_eq(self.shininess, other.shininess)
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
            "{} Material:{} Color: {}{:.3}{} {}{:.3}{} {}{:.3}{} {}ambient:{:.3} diffuse:{:.3} specular:{:.3} shininess:{:.3} {}",
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
            self.ambient,
            self.diffuse,
            self.specular,
            self.shininess,
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
        let chk = chk && m.approx_eq(Material::new());
        if chk {
            Ok(())
        } else {
            logi!("test_chap_6_11", "m:{}", m);
            logd!("test_chap_6_11", "m:{:?}", m);
            loge!("test_chap_6_11", "chk:{:?}", chk);
            Err("The default material".into())
        }
    }
}
