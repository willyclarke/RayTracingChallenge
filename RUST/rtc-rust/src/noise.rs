//! Perlin noise (book chapter 17, "Next Steps")
//!
//! Ken Perlin's 2002 "improved noise": a smooth, band-limited random
//! function of 3D space, zero at every integer lattice point and in `[-1, 1]`
//! everywhere. Used to jitter patterns (`PerturbedPattern`) and normals
//! (`Material::bump`).

use crate::tuple::Tuple;

#[derive(Debug, Clone)]
pub struct Perlin {
    /// Permutation table doubled to 512 entries so indices never wrap.
    perm: Box<[u8; 512]>,
}

impl Perlin {
    /// Deterministic for a given seed: the same seed always yields the same
    /// noise field, so renders are reproducible.
    pub fn new(seed: u64) -> Self {
        let mut table: [u8; 256] = std::array::from_fn(|i| i as u8);
        // xorshift64* to shuffle (Fisher-Yates); no external RNG dependency.
        let mut state = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
        let mut next = || {
            state ^= state >> 12;
            state ^= state << 25;
            state ^= state >> 27;
            state.wrapping_mul(0x2545_F491_4F6C_DD1D)
        };
        for i in (1..256).rev() {
            let j = (next() % (i as u64 + 1)) as usize;
            table.swap(i, j);
        }
        let mut perm = Box::new([0u8; 512]);
        for i in 0..512 {
            perm[i] = table[i % 256];
        }
        Self { perm }
    }

    /// Noise value at `p`, in `[-1, 1]`; exactly 0 at integer lattice points.
    pub fn noise(&self, p: Tuple) -> f64 {
        let (xf, yf, zf) = (p.x.floor(), p.y.floor(), p.z.floor());
        let (xi, yi, zi) = (
            (xf as i64 & 255) as usize,
            (yf as i64 & 255) as usize,
            (zf as i64 & 255) as usize,
        );
        let (x, y, z) = (p.x - xf, p.y - yf, p.z - zf);
        let (u, v, w) = (fade(x), fade(y), fade(z));

        let perm = &self.perm;
        let a = perm[xi] as usize + yi;
        let aa = perm[a] as usize + zi;
        let ab = perm[a + 1] as usize + zi;
        let b = perm[xi + 1] as usize + yi;
        let ba = perm[b] as usize + zi;
        let bb = perm[b + 1] as usize + zi;

        lerp(
            w,
            lerp(
                v,
                lerp(u, grad(perm[aa], x, y, z), grad(perm[ba], x - 1.0, y, z)),
                lerp(
                    u,
                    grad(perm[ab], x, y - 1.0, z),
                    grad(perm[bb], x - 1.0, y - 1.0, z),
                ),
            ),
            lerp(
                v,
                lerp(
                    u,
                    grad(perm[aa + 1], x, y, z - 1.0),
                    grad(perm[ba + 1], x - 1.0, y, z - 1.0),
                ),
                lerp(
                    u,
                    grad(perm[ab + 1], x, y - 1.0, z - 1.0),
                    grad(perm[bb + 1], x - 1.0, y - 1.0, z - 1.0),
                ),
            ),
        )
    }

    /// Three decorrelated noise samples as a vector, for jittering points
    /// and normals.
    pub fn vector(&self, p: Tuple) -> Tuple {
        Tuple::vector(
            self.noise(p),
            self.noise(p + Tuple::vector(31.4, 15.9, 26.5)),
            self.noise(p + Tuple::vector(-27.1, 82.8, -18.2)),
        )
    }

    /// Fractal sum of `octaves` noise layers, each at twice the frequency
    /// and half the amplitude of the last. Still roughly in `[-1, 1]`.
    pub fn octaves(&self, p: Tuple, octaves: u32) -> f64 {
        let mut sum = 0.0;
        let mut amplitude = 1.0;
        let mut frequency = 1.0;
        let mut max = 0.0;
        for _ in 0..octaves {
            sum += amplitude * self.noise(p * frequency);
            max += amplitude;
            amplitude *= 0.5;
            frequency *= 2.0;
        }
        sum / max
    }
}

impl Default for Perlin {
    fn default() -> Self {
        Self::new(0)
    }
}

/// 6t^5 - 15t^4 + 10t^3: zero first and second derivatives at 0 and 1, so
/// lattice cells join without visible seams.
fn fade(t: f64) -> f64 {
    t * t * t * (t * (t * 6.0 - 15.0) + 10.0)
}

fn lerp(t: f64, a: f64, b: f64) -> f64 {
    a + t * (b - a)
}

/// Dot product of the offset with one of 12 gradient directions (the cube
/// edge midpoints), picked by the low bits of the hash.
fn grad(hash: u8, x: f64, y: f64, z: f64) -> f64 {
    let h = hash & 15;
    let u = if h < 8 { x } else { y };
    let v = if h < 4 {
        y
    } else if h == 12 || h == 14 {
        x
    } else {
        z
    };
    (if h & 1 == 0 { u } else { -u }) + (if h & 2 == 0 { v } else { -v })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::loge;
    use crate::math::approx_eq;

    /// Chap 17 - Perlin noise is zero at every integer lattice point
    #[test]
    fn test_chap_17_5() -> Result<(), String> {
        let perlin = Perlin::new(42);
        for x in -3..4 {
            for y in -3..4 {
                for z in -3..4 {
                    let n = perlin.noise(Tuple::point(x as f64, y as f64, z as f64));
                    if !approx_eq(n, 0.0) {
                        loge!("test_chap_17_5", "({x},{y},{z}) noise:{n}");
                        return Err("Perlin noise is zero at lattice points".into());
                    }
                }
            }
        }
        Ok(())
    }

    /// Chap 17 - Perlin noise is bounded, non-trivial, and deterministic per seed
    #[test]
    fn test_chap_17_6() -> Result<(), String> {
        let a = Perlin::new(7);
        let b = Perlin::new(7);
        let c = Perlin::new(8);
        let mut max_abs: f64 = 0.0;
        let mut differs_by_seed = false;
        for i in 0..1000 {
            let t = i as f64 * 0.0371;
            let p = Tuple::point(t * 1.3, t * 0.7 + 0.5, -t * 1.1 + 0.25);
            let na = a.noise(p);
            if !approx_eq(na, b.noise(p)) {
                return Err("Same seed must give the same noise".into());
            }
            if !approx_eq(na, c.noise(p)) {
                differs_by_seed = true;
            }
            max_abs = max_abs.max(na.abs());
        }
        if max_abs <= 1.0 && max_abs > 0.1 && differs_by_seed {
            Ok(())
        } else {
            loge!(
                "test_chap_17_6",
                "max_abs:{max_abs} differs:{differs_by_seed}"
            );
            Err("Perlin noise is bounded, non-trivial and seed-dependent".into())
        }
    }

    /// Chap 17 - Perlin noise is continuous: a small step changes it a little
    #[test]
    fn test_chap_17_7() -> Result<(), String> {
        let perlin = Perlin::new(1);
        let step = 1e-4;
        for i in 0..200 {
            let t = i as f64 * 0.137;
            let p = Tuple::point(t, t * 0.5, -t);
            let d = (perlin.noise(p + Tuple::vector(step, 0.0, 0.0)) - perlin.noise(p)).abs();
            // |gradient| of improved noise is bounded well below 10
            if d > 10.0 * step {
                loge!("test_chap_17_7", "p:{} delta:{d}", p);
                return Err("Perlin noise is continuous".into());
            }
        }
        Ok(())
    }
}
