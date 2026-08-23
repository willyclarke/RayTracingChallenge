//! UV patterns (bonus chapter "Texture mapping")
//!
//! A `UvPattern` is a 2D pattern addressed by `(u, v)` in `[0, 1]`. A
//! `TextureMap` or `CubeMap` (see `texturemap`) maps a 3D point to `(u, v)`
//! and looks the colour up here.

use crate::canvas::Canvas;
use crate::tuple::Tuple;

pub trait UvPattern: std::fmt::Debug + Send + Sync {
    fn uv_pattern_at(&self, u: f64, v: f64) -> Tuple;
    fn clone_box(&self) -> Box<dyn UvPattern>;
}

impl Clone for Box<dyn UvPattern> {
    fn clone(&self) -> Self {
        self.clone_box()
    }
}

/// A `width` x `height` checkerboard over the unit UV square.
#[derive(Debug, Clone, Copy)]
pub struct UvCheckers {
    pub width: usize,
    pub height: usize,
    pub a: Tuple,
    pub b: Tuple,
}

impl UvCheckers {
    pub fn new(width: usize, height: usize, a: Tuple, b: Tuple) -> Self {
        Self {
            width,
            height,
            a,
            b,
        }
    }
}

impl UvPattern for UvCheckers {
    fn uv_pattern_at(&self, u: f64, v: f64) -> Tuple {
        let u2 = (u * self.width as f64).floor() as i64;
        let v2 = (v * self.height as f64).floor() as i64;
        if (u2 + v2) % 2 == 0 { self.a } else { self.b }
    }

    fn clone_box(&self) -> Box<dyn UvPattern> {
        Box::new(*self)
    }
}

/// A `main` colour with a distinct colour in each corner, for checking
/// which way a face is oriented when building cube maps.
#[derive(Debug, Clone, Copy)]
pub struct UvAlignCheck {
    pub main: Tuple,
    pub ul: Tuple,
    pub ur: Tuple,
    pub bl: Tuple,
    pub br: Tuple,
}

impl UvAlignCheck {
    pub fn new(main: Tuple, ul: Tuple, ur: Tuple, bl: Tuple, br: Tuple) -> Self {
        Self {
            main,
            ul,
            ur,
            bl,
            br,
        }
    }
}

impl UvPattern for UvAlignCheck {
    fn uv_pattern_at(&self, u: f64, v: f64) -> Tuple {
        if v > 0.8 {
            if u < 0.2 {
                return self.ul;
            }
            if u > 0.8 {
                return self.ur;
            }
        } else if v < 0.2 {
            if u < 0.2 {
                return self.bl;
            }
            if u > 0.8 {
                return self.br;
            }
        }
        self.main
    }

    fn clone_box(&self) -> Box<dyn UvPattern> {
        Box::new(*self)
    }
}

/// An image (a `Canvas`, typically read from a PPM) stretched over the unit
/// UV square. `v = 0` is the bottom row of the image.
#[derive(Debug, Clone)]
pub struct UvImage {
    pub canvas: Canvas,
}

impl UvImage {
    pub fn new(canvas: Canvas) -> Self {
        Self { canvas }
    }
}

impl UvPattern for UvImage {
    fn uv_pattern_at(&self, u: f64, v: f64) -> Tuple {
        // flip v so that v = 0 is at the bottom of the image
        let v = 1.0 - v;
        let x = u * (self.canvas.width() - 1) as f64;
        let y = v * (self.canvas.height() - 1) as f64;
        self.canvas.pixel_at(x.round() as usize, y.round() as usize)
    }

    fn clone_box(&self) -> Box<dyn UvPattern> {
        Box::new(self.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::loge;
    use crate::tuple::colors::*;

    /// Bonus (texture mapping) - Checker pattern in 2D
    #[test]
    fn test_bonus_texture_1() -> Result<(), String> {
        let checkers = UvCheckers::new(2, 2, BLACK, WHITE);
        let cases = [
            (0.0, 0.0, BLACK),
            (0.5, 0.0, WHITE),
            (0.0, 0.5, WHITE),
            (0.5, 0.5, BLACK),
            (1.0, 1.0, BLACK),
        ];
        for (u, v, expected) in cases {
            let c = checkers.uv_pattern_at(u, v);
            if !c.approx_eq(expected) {
                loge!("test_bonus_texture_1", "u:{} v:{} got:{}", u, v, c);
                return Err("Checker pattern in 2D".into());
            }
        }
        Ok(())
    }

    /// Bonus (texture mapping) - Layout of the "align check" pattern
    #[test]
    fn test_bonus_texture_9() -> Result<(), String> {
        let main = Tuple::color(1.0, 1.0, 1.0);
        let ul = Tuple::color(1.0, 0.0, 0.0);
        let ur = Tuple::color(1.0, 1.0, 0.0);
        let bl = Tuple::color(0.0, 1.0, 0.0);
        let br = Tuple::color(0.0, 1.0, 1.0);
        let pattern = UvAlignCheck::new(main, ul, ur, bl, br);
        let cases = [
            (0.5, 0.5, main),
            (0.1, 0.9, ul),
            (0.9, 0.9, ur),
            (0.1, 0.1, bl),
            (0.9, 0.1, br),
        ];
        for (u, v, expected) in cases {
            let c = pattern.uv_pattern_at(u, v);
            if !c.approx_eq(expected) {
                loge!("test_bonus_texture_9", "u:{} v:{} got:{}", u, v, c);
                return Err("Layout of the align check pattern".into());
            }
        }
        Ok(())
    }

    /// Bonus (texture mapping) - Checker pattern in 2D from an image
    #[test]
    fn test_bonus_texture_20() -> Result<(), String> {
        let mut ppm = String::from("P3\n10 10\n10\n");
        for row in 0..10 {
            let line: Vec<String> = (0..10)
                .map(|col| {
                    let v = (row + col) % 10;
                    format!("{v} {v} {v}")
                })
                .collect();
            ppm.push_str(&line.join("  "));
            ppm.push('\n');
        }
        let canvas = Canvas::from_ppm(&ppm)?;
        let pattern = UvImage::new(canvas);
        let cases = [
            (0.0, 0.0, Tuple::color(0.9, 0.9, 0.9)),
            (0.3, 0.0, Tuple::color(0.2, 0.2, 0.2)),
            (0.6, 0.3, Tuple::color(0.1, 0.1, 0.1)),
            (1.0, 1.0, Tuple::color(0.9, 0.9, 0.9)),
        ];
        for (u, v, expected) in cases {
            let c = pattern.uv_pattern_at(u, v);
            if !c.approx_eq(expected) {
                loge!("test_bonus_texture_20", "u:{} v:{} got:{}", u, v, c);
                return Err("Checker pattern in 2D from an image".into());
            }
        }
        Ok(())
    }
}
