//! Perturbed pattern (book chapter 17, "Next Steps")
//!
//! Wraps any pattern and looks it up at a point jittered by Perlin noise:
//! stripes become marble, checkers go wobbly.

use crate::matrix::Matrix4;
use crate::noise::Perlin;
use crate::pattern::{Pattern, PatternData};
use crate::tuple::Tuple;

#[derive(Debug, Clone)]
pub struct PerturbedPattern {
    pub data: PatternData,
    pub pattern: Box<dyn Pattern>,
    pub noise: Perlin,
    /// How far (in pattern units) a point can be displaced.
    pub amplitude: f64,
    /// Spatial frequency of the noise: higher = finer wobbles.
    pub frequency: f64,
}

impl PerturbedPattern {
    pub fn new(pattern: Box<dyn Pattern>, amplitude: f64, frequency: f64) -> Self {
        Self {
            data: PatternData::new(),
            pattern,
            noise: Perlin::default(),
            amplitude,
            frequency,
        }
    }

    pub fn with_seed(mut self, seed: u64) -> Self {
        self.noise = Perlin::new(seed);
        self
    }

    pub fn perturb(&self, point: Tuple) -> Tuple {
        point + self.noise.vector(point * self.frequency) * self.amplitude
    }
}

impl Pattern for PerturbedPattern {
    fn data(&self) -> &PatternData {
        &self.data
    }

    fn data_mut(&mut self) -> &mut PatternData {
        &mut self.data
    }

    fn color_at(&self, point: Tuple) -> Tuple {
        self.pattern.color_at_local(self.perturb(point))
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
    use crate::patterns::stripepattern::StripePattern;
    use crate::tuple::colors::*;

    /// Chap 17 - With zero amplitude a perturbed pattern is the pattern
    #[test]
    fn test_chap_17_8() -> Result<(), String> {
        let stripes = StripePattern::new(WHITE, BLACK);
        let p = PerturbedPattern::new(Box::new(stripes), 0.0, 1.0);
        for i in 0..50 {
            let pt = Tuple::point(i as f64 * 0.123, 0.3, -0.7);
            if !p.color_at(pt).approx_eq(stripes.color_at(pt)) {
                return Err("Zero amplitude must not change the pattern".into());
            }
        }
        Ok(())
    }

    /// Chap 17 - A perturbed pattern displaces the lookup point by at most
    /// the amplitude, and actually displaces it somewhere
    #[test]
    fn test_chap_17_9() -> Result<(), String> {
        let stripes = StripePattern::new(WHITE, BLACK);
        let p = PerturbedPattern::new(Box::new(stripes), 0.2, 1.0);
        let mut moved = false;
        for i in 0..200 {
            let t = i as f64 * 0.0713;
            let pt = Tuple::point(t, t * 0.4, -t * 0.9);
            let q = p.perturb(pt);
            let d = q - pt;
            if d.x.abs() > 0.2 || d.y.abs() > 0.2 || d.z.abs() > 0.2 {
                return Err("Displacement exceeds amplitude".into());
            }
            if d.magnitude() > 1e-3 {
                moved = true;
            }
        }
        if moved {
            Ok(())
        } else {
            Err("Perturbation never displaced the point".into())
        }
    }
}
