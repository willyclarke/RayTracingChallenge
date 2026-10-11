//! Texture mapping (bonus chapter "Texture mapping")
//!
//! `TextureMap` maps a 3D pattern-space point to `(u, v)` with a spherical,
//! planar or cylindrical projection and looks the colour up in a
//! `UvPattern`. `CubeMap` does the same with one `UvPattern` per cube face.

use std::f64::consts::PI;

use crate::matrix::Matrix4;
use crate::pattern::{Pattern, PatternData};
use crate::patterns::uvpattern::UvPattern;
use crate::tuple::Tuple;

/// How a 3D point is projected onto the unit UV square.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UvMap {
    Spherical,
    Planar,
    Cylindrical,
}

impl UvMap {
    pub fn uv(self, p: Tuple) -> (f64, f64) {
        match self {
            UvMap::Spherical => spherical_map(p),
            UvMap::Planar => planar_map(p),
            UvMap::Cylindrical => cylindrical_map(p),
        }
    }
}

/// Azimuth around the y axis as `u` (0 at -z, increasing clockwise seen
/// from above), polar angle as `v` (0 at the south pole).
pub fn spherical_map(p: Tuple) -> (f64, f64) {
    let theta = p.x.atan2(p.z);
    let radius = Tuple::vector(p.x, p.y, p.z).magnitude();
    let phi = (p.y / radius).acos();
    let raw_u = theta / (2.0 * PI);
    let u = 1.0 - (raw_u + 0.5);
    let v = 1.0 - phi / PI;
    (u, v)
}

/// `x` and `z` modulo 1, repeating every unit.
pub fn planar_map(p: Tuple) -> (f64, f64) {
    (p.x.rem_euclid(1.0), p.z.rem_euclid(1.0))
}

/// Azimuth as `u` like the spherical map; `y` modulo 1 as `v`.
pub fn cylindrical_map(p: Tuple) -> (f64, f64) {
    let theta = p.x.atan2(p.z);
    let raw_u = theta / (2.0 * PI);
    let u = 1.0 - (raw_u + 0.5);
    let v = p.y.rem_euclid(1.0);
    (u, v)
}

#[derive(Debug, Clone)]
pub struct TextureMap {
    pub data: PatternData,
    pub uv_pattern: Box<dyn UvPattern>,
    pub uv_map: UvMap,
}

impl TextureMap {
    pub fn new(uv_pattern: Box<dyn UvPattern>, uv_map: UvMap) -> Self {
        Self {
            data: PatternData::new(),
            uv_pattern,
            uv_map,
        }
    }
}

impl Pattern for TextureMap {
    fn data(&self) -> &PatternData {
        &self.data
    }

    fn data_mut(&mut self) -> &mut PatternData {
        &mut self.data
    }

    fn color_at(&self, point: Tuple) -> Tuple {
        let (u, v) = self.uv_map.uv(point);
        self.uv_pattern.uv_pattern_at(u, v)
    }

    fn clone_box(&self) -> Box<dyn Pattern> {
        Box::new(self.clone())
    }

    fn set_transform(&mut self, m: Matrix4) {
        self.data.set_transform(m);
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CubeFace {
    Left,
    Right,
    Front,
    Back,
    Up,
    Down,
}

/// Which face of the unit cube a point lies on (or is nearest to): the axis
/// with the largest absolute coordinate.
pub fn face_from_point(p: Tuple) -> CubeFace {
    let coord = p.x.abs().max(p.y.abs()).max(p.z.abs());
    if coord == p.x {
        CubeFace::Right
    } else if coord == -p.x {
        CubeFace::Left
    } else if coord == p.y {
        CubeFace::Up
    } else if coord == -p.y {
        CubeFace::Down
    } else if coord == p.z {
        CubeFace::Front
    } else {
        CubeFace::Back
    }
}

pub fn cube_uv_front(p: Tuple) -> (f64, f64) {
    (
        (p.x + 1.0).rem_euclid(2.0) / 2.0,
        (p.y + 1.0).rem_euclid(2.0) / 2.0,
    )
}

pub fn cube_uv_back(p: Tuple) -> (f64, f64) {
    (
        (1.0 - p.x).rem_euclid(2.0) / 2.0,
        (p.y + 1.0).rem_euclid(2.0) / 2.0,
    )
}

pub fn cube_uv_left(p: Tuple) -> (f64, f64) {
    (
        (p.z - 1.0).rem_euclid(2.0) / 2.0,
        (p.y + 1.0).rem_euclid(2.0) / 2.0,
    )
}

pub fn cube_uv_right(p: Tuple) -> (f64, f64) {
    (
        (1.0 - p.z).rem_euclid(2.0) / 2.0,
        (p.y + 1.0).rem_euclid(2.0) / 2.0,
    )
}

pub fn cube_uv_up(p: Tuple) -> (f64, f64) {
    (
        (p.x + 1.0).rem_euclid(2.0) / 2.0,
        (1.0 - p.z).rem_euclid(2.0) / 2.0,
    )
}

pub fn cube_uv_down(p: Tuple) -> (f64, f64) {
    (
        (p.x + 1.0).rem_euclid(2.0) / 2.0,
        (p.z + 1.0).rem_euclid(2.0) / 2.0,
    )
}

/// One `UvPattern` per face of the unit cube.
#[derive(Debug, Clone)]
pub struct CubeMap {
    pub data: PatternData,
    pub left: Box<dyn UvPattern>,
    pub front: Box<dyn UvPattern>,
    pub right: Box<dyn UvPattern>,
    pub back: Box<dyn UvPattern>,
    pub up: Box<dyn UvPattern>,
    pub down: Box<dyn UvPattern>,
}

impl CubeMap {
    pub fn new(
        left: Box<dyn UvPattern>,
        front: Box<dyn UvPattern>,
        right: Box<dyn UvPattern>,
        back: Box<dyn UvPattern>,
        up: Box<dyn UvPattern>,
        down: Box<dyn UvPattern>,
    ) -> Self {
        Self {
            data: PatternData::new(),
            left,
            front,
            right,
            back,
            up,
            down,
        }
    }
}

impl Pattern for CubeMap {
    fn data(&self) -> &PatternData {
        &self.data
    }

    fn data_mut(&mut self) -> &mut PatternData {
        &mut self.data
    }

    fn color_at(&self, point: Tuple) -> Tuple {
        let (pattern, (u, v)) = match face_from_point(point) {
            CubeFace::Left => (&self.left, cube_uv_left(point)),
            CubeFace::Right => (&self.right, cube_uv_right(point)),
            CubeFace::Front => (&self.front, cube_uv_front(point)),
            CubeFace::Back => (&self.back, cube_uv_back(point)),
            CubeFace::Up => (&self.up, cube_uv_up(point)),
            CubeFace::Down => (&self.down, cube_uv_down(point)),
        };
        pattern.uv_pattern_at(u, v)
    }

    fn clone_box(&self) -> Box<dyn Pattern> {
        Box::new(self.clone())
    }

    fn set_transform(&mut self, m: Matrix4) {
        self.data.set_transform(m);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::loge;
    use crate::math::approx_eq;
    use crate::patterns::uvpattern::{UvAlignCheck, UvCheckers};
    use crate::tuple::colors::*;

    fn check_uv(
        name: &str,
        cases: &[(Tuple, f64, f64)],
        map: fn(Tuple) -> (f64, f64),
    ) -> Result<(), String> {
        for &(p, eu, ev) in cases {
            let (u, v) = map(p);
            if !(approx_eq(u, eu) && approx_eq(v, ev)) {
                loge!(name, "p:{} got:({}, {}) want:({}, {})", p, u, v, eu, ev);
                return Err(name.into());
            }
        }
        Ok(())
    }

    /// Bonus (texture mapping) - Using a spherical mapping on a 3D point
    #[test]
    fn test_bonus_texture_2() -> Result<(), String> {
        let s = std::f64::consts::FRAC_1_SQRT_2;
        let cases = [
            (Tuple::point(0.0, 0.0, -1.0), 0.0, 0.5),
            (Tuple::point(1.0, 0.0, 0.0), 0.25, 0.5),
            (Tuple::point(0.0, 0.0, 1.0), 0.5, 0.5),
            (Tuple::point(-1.0, 0.0, 0.0), 0.75, 0.5),
            (Tuple::point(0.0, 1.0, 0.0), 0.5, 1.0),
            (Tuple::point(0.0, -1.0, 0.0), 0.5, 0.0),
            (Tuple::point(s, s, 0.0), 0.25, 0.75),
        ];
        check_uv("test_bonus_texture_2", &cases, spherical_map)
    }

    /// Bonus (texture mapping) - Using a texture map pattern with a spherical map
    #[test]
    fn test_bonus_texture_3() -> Result<(), String> {
        let checkers = UvCheckers::new(16, 8, BLACK, WHITE);
        let pattern = TextureMap::new(Box::new(checkers), UvMap::Spherical);
        let cases = [
            (Tuple::point(0.4315, 0.4670, 0.7719), WHITE),
            (Tuple::point(-0.9654, 0.2552, -0.0534), BLACK),
            (Tuple::point(0.1039, 0.7090, 0.6975), WHITE),
            (Tuple::point(-0.4986, -0.7856, -0.3663), BLACK),
            (Tuple::point(-0.0317, -0.9395, 0.3411), BLACK),
            (Tuple::point(0.4809, -0.7721, 0.4154), BLACK),
            (Tuple::point(0.0285, -0.9612, -0.2745), BLACK),
            (Tuple::point(-0.5734, -0.2162, -0.7903), WHITE),
            (Tuple::point(0.7688, -0.1470, 0.6223), BLACK),
            (Tuple::point(-0.7652, 0.2175, 0.6060), BLACK),
        ];
        for (p, expected) in cases {
            let c = pattern.color_at(p);
            if !c.approx_eq(expected) {
                loge!("test_bonus_texture_3", "p:{} got:{}", p, c);
                return Err("Using a texture map pattern with a spherical map".into());
            }
        }
        Ok(())
    }

    /// Bonus (texture mapping) - Using a planar mapping on a 3D point
    #[test]
    fn test_bonus_texture_4() -> Result<(), String> {
        let cases = [
            (Tuple::point(0.25, 0.0, 0.5), 0.25, 0.5),
            (Tuple::point(0.25, 0.0, -0.25), 0.25, 0.75),
            (Tuple::point(0.25, 0.5, -0.25), 0.25, 0.75),
            (Tuple::point(1.25, 0.0, 0.5), 0.25, 0.5),
            (Tuple::point(0.25, 0.0, -1.75), 0.25, 0.25),
            (Tuple::point(1.0, 0.0, -1.0), 0.0, 0.0),
            (Tuple::point(0.0, 0.0, 0.0), 0.0, 0.0),
        ];
        check_uv("test_bonus_texture_4", &cases, planar_map)
    }

    /// Bonus (texture mapping) - Using a cylindrical mapping on a 3D point
    #[test]
    fn test_bonus_texture_5() -> Result<(), String> {
        let s = std::f64::consts::FRAC_1_SQRT_2;
        let cases = [
            (Tuple::point(0.0, 0.0, -1.0), 0.0, 0.0),
            (Tuple::point(0.0, 0.5, -1.0), 0.0, 0.5),
            (Tuple::point(0.0, 1.0, -1.0), 0.0, 0.0),
            (Tuple::point(s, 0.5, -s), 0.125, 0.5),
            (Tuple::point(1.0, 0.5, 0.0), 0.25, 0.5),
            (Tuple::point(s, 0.5, s), 0.375, 0.5),
            (Tuple::point(0.0, -0.25, 1.0), 0.5, 0.75),
            (Tuple::point(-s, 0.5, s), 0.625, 0.5),
            (Tuple::point(-1.0, 1.25, 0.0), 0.75, 0.25),
            (Tuple::point(-s, 0.5, -s), 0.875, 0.5),
        ];
        check_uv("test_bonus_texture_5", &cases, cylindrical_map)
    }

    /// Bonus (texture mapping) - Identifying the face of a cube from a point
    #[test]
    fn test_bonus_texture_10() -> Result<(), String> {
        let cases = [
            (Tuple::point(-1.0, 0.5, -0.25), CubeFace::Left),
            (Tuple::point(1.1, -0.5, 0.8), CubeFace::Right),
            (Tuple::point(0.1, 0.6, 0.9), CubeFace::Front),
            (Tuple::point(-0.7, 0.0, -2.0), CubeFace::Back),
            (Tuple::point(0.5, 1.0, 0.9), CubeFace::Up),
            (Tuple::point(-0.2, -1.3, 1.1), CubeFace::Down),
        ];
        for (p, expected) in cases {
            let face = face_from_point(p);
            if face != expected {
                loge!("test_bonus_texture_10", "p:{} got:{:?}", p, face);
                return Err("Identifying the face of a cube from a point".into());
            }
        }
        Ok(())
    }

    /// Bonus (texture mapping) - UV mapping the front face of a cube
    #[test]
    fn test_bonus_texture_11() -> Result<(), String> {
        let cases = [
            (Tuple::point(-0.5, 0.5, 1.0), 0.25, 0.75),
            (Tuple::point(0.5, -0.5, 1.0), 0.75, 0.25),
        ];
        check_uv("test_bonus_texture_11", &cases, cube_uv_front)
    }

    /// Bonus (texture mapping) - UV mapping the back face of a cube
    #[test]
    fn test_bonus_texture_12() -> Result<(), String> {
        let cases = [
            (Tuple::point(0.5, 0.5, -1.0), 0.25, 0.75),
            (Tuple::point(-0.5, -0.5, -1.0), 0.75, 0.25),
        ];
        check_uv("test_bonus_texture_12", &cases, cube_uv_back)
    }

    /// Bonus (texture mapping) - UV mapping the left face of a cube
    #[test]
    fn test_bonus_texture_13() -> Result<(), String> {
        let cases = [
            (Tuple::point(-1.0, 0.5, -0.5), 0.25, 0.75),
            (Tuple::point(-1.0, -0.5, 0.5), 0.75, 0.25),
        ];
        check_uv("test_bonus_texture_13", &cases, cube_uv_left)
    }

    /// Bonus (texture mapping) - UV mapping the right face of a cube
    #[test]
    fn test_bonus_texture_14() -> Result<(), String> {
        let cases = [
            (Tuple::point(1.0, 0.5, 0.5), 0.25, 0.75),
            (Tuple::point(1.0, -0.5, -0.5), 0.75, 0.25),
        ];
        check_uv("test_bonus_texture_14", &cases, cube_uv_right)
    }

    /// Bonus (texture mapping) - UV mapping the upper face of a cube
    #[test]
    fn test_bonus_texture_15() -> Result<(), String> {
        let cases = [
            (Tuple::point(-0.5, 1.0, -0.5), 0.25, 0.75),
            (Tuple::point(0.5, 1.0, 0.5), 0.75, 0.25),
        ];
        check_uv("test_bonus_texture_15", &cases, cube_uv_up)
    }

    /// Bonus (texture mapping) - UV mapping the lower face of a cube
    #[test]
    fn test_bonus_texture_16() -> Result<(), String> {
        let cases = [
            (Tuple::point(-0.5, -1.0, 0.5), 0.25, 0.75),
            (Tuple::point(0.5, -1.0, -0.5), 0.75, 0.25),
        ];
        check_uv("test_bonus_texture_16", &cases, cube_uv_down)
    }

    /// Bonus (texture mapping) - Finding the colors on a mapped cube
    #[test]
    fn test_bonus_texture_17() -> Result<(), String> {
        let red = Tuple::color(1.0, 0.0, 0.0);
        let yellow = Tuple::color(1.0, 1.0, 0.0);
        let brown = Tuple::color(1.0, 0.5, 0.0);
        let green = Tuple::color(0.0, 1.0, 0.0);
        let cyan = Tuple::color(0.0, 1.0, 1.0);
        let blue = Tuple::color(0.0, 0.0, 1.0);
        let purple = Tuple::color(1.0, 0.0, 1.0);
        let white = Tuple::color(1.0, 1.0, 1.0);

        let align = |main, ul, ur, bl, br| -> Box<dyn UvPattern> {
            Box::new(UvAlignCheck::new(main, ul, ur, bl, br))
        };
        let pattern = CubeMap::new(
            align(yellow, cyan, red, blue, brown),
            align(cyan, red, yellow, brown, green),
            align(red, yellow, purple, green, white),
            align(green, purple, cyan, white, blue),
            align(brown, cyan, purple, red, yellow),
            align(purple, brown, green, blue, white),
        );

        let cases = [
            // left
            (Tuple::point(-1.0, 0.0, 0.0), yellow),
            (Tuple::point(-1.0, 0.9, -0.9), cyan),
            (Tuple::point(-1.0, 0.9, 0.9), red),
            (Tuple::point(-1.0, -0.9, -0.9), blue),
            (Tuple::point(-1.0, -0.9, 0.9), brown),
            // front
            (Tuple::point(0.0, 0.0, 1.0), cyan),
            (Tuple::point(-0.9, 0.9, 1.0), red),
            (Tuple::point(0.9, 0.9, 1.0), yellow),
            (Tuple::point(-0.9, -0.9, 1.0), brown),
            (Tuple::point(0.9, -0.9, 1.0), green),
            // right
            (Tuple::point(1.0, 0.0, 0.0), red),
            (Tuple::point(1.0, 0.9, 0.9), yellow),
            (Tuple::point(1.0, 0.9, -0.9), purple),
            (Tuple::point(1.0, -0.9, 0.9), green),
            (Tuple::point(1.0, -0.9, -0.9), white),
            // back
            (Tuple::point(0.0, 0.0, -1.0), green),
            (Tuple::point(0.9, 0.9, -1.0), purple),
            (Tuple::point(-0.9, 0.9, -1.0), cyan),
            (Tuple::point(0.9, -0.9, -1.0), white),
            (Tuple::point(-0.9, -0.9, -1.0), blue),
            // up
            (Tuple::point(0.0, 1.0, 0.0), brown),
            (Tuple::point(-0.9, 1.0, -0.9), cyan),
            (Tuple::point(0.9, 1.0, -0.9), purple),
            (Tuple::point(-0.9, 1.0, 0.9), red),
            (Tuple::point(0.9, 1.0, 0.9), yellow),
            // down
            (Tuple::point(0.0, -1.0, 0.0), purple),
            (Tuple::point(-0.9, -1.0, 0.9), brown),
            (Tuple::point(0.9, -1.0, 0.9), green),
            (Tuple::point(-0.9, -1.0, -0.9), blue),
            (Tuple::point(0.9, -1.0, -0.9), white),
        ];
        for (p, expected) in cases {
            let c = pattern.color_at(p);
            if !c.approx_eq(expected) {
                loge!("test_bonus_texture_17", "p:{} got:{}", p, c);
                return Err("Finding the colors on a mapped cube".into());
            }
        }
        Ok(())
    }
}
