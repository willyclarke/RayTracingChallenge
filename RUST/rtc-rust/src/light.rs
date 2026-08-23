//! Light
//!
//! Let there be light
//!

use crate::{shape::Shape, tuple::Tuple, world::World};
use std::sync::atomic::{AtomicUsize, Ordering};

/// A cyclic number generator used to jitter the sample points on an area
/// light (bonus chapter "Rendering soft shadows"). `next()` returns the
/// values in order and wraps around.
///
/// ```
/// use rtc_rust::light::Sequence;
/// let seq = Sequence::new(vec![0.1, 0.5, 1.0]);
/// assert_eq!(seq.next(), 0.1);
/// assert_eq!(seq.next(), 0.5);
/// assert_eq!(seq.next(), 1.0);
/// assert_eq!(seq.next(), 0.1);
/// ```
#[derive(Debug)]
pub struct Sequence {
    values: Vec<f64>,
    // Atomic so a light shared across rayon render threads stays `Sync`; the
    // jitter only needs to be "not a grid", so the cross-thread order is fine.
    index: AtomicUsize,
}

impl Sequence {
    pub fn new(values: Vec<f64>) -> Self {
        assert!(!values.is_empty(), "Sequence needs at least one value");
        Self {
            values,
            index: AtomicUsize::new(0),
        }
    }

    #[allow(clippy::should_implement_trait)]
    pub fn next(&self) -> f64 {
        self.value(self.reserve(1))
    }

    /// Claim the next `n` positions with one atomic op; read them with
    /// `value(start + i)`. Equivalent to `n` calls of `next()` but without
    /// `n` contended atomics per pixel under `render_parallel`.
    pub fn reserve(&self, n: usize) -> usize {
        self.index.fetch_add(n, Ordering::Relaxed)
    }

    pub fn value(&self, i: usize) -> f64 {
        self.values[i % self.values.len()]
    }
}

impl Clone for Sequence {
    fn clone(&self) -> Self {
        Self {
            values: self.values.clone(),
            index: AtomicUsize::new(self.index.load(Ordering::Relaxed)),
        }
    }
}

/// A rectangular area light made of `usteps` x `vsteps` cells. A point light
/// is the 1x1 degenerate case with zero-length edge vectors, so every light
/// goes through the same sampling code.
#[derive(Debug, Clone)]
pub struct Light {
    /// Centre of the light (the single sample point for a point light).
    pub position: Tuple,
    pub intensity: Tuple,
    pub corner: Tuple,
    /// One cell's worth of the u edge (`full_uvec / usteps`).
    pub uvec: Tuple,
    /// One cell's worth of the v edge (`full_vvec / vsteps`).
    pub vvec: Tuple,
    pub usteps: usize,
    pub vsteps: usize,
    pub samples: usize,
    /// Offset within a cell for each sample, in `[0, 1)`. `0.5` = cell centre.
    pub jitter_by: Sequence,
    /// Restrict the light to a cone (spotlight); `None` shines everywhere.
    pub spot: Option<Spot>,
}

/// The cone of a spotlight (book chapter 17): full intensity within
/// `inner_angle` of `direction`, fading smoothly to zero at `outer_angle`.
#[derive(Debug, Clone, Copy)]
pub struct Spot {
    /// Unit vector the light points along.
    pub direction: Tuple,
    pub cos_inner: f64,
    pub cos_outer: f64,
}

impl Spot {
    /// 1.0 inside the inner cone, 0.0 outside the outer cone, smoothstep in
    /// between. `to_point` is the (unnormalised) vector from the light.
    pub fn factor(&self, to_point: Tuple) -> f64 {
        let cos_angle = to_point.normalize().dot(self.direction);
        if cos_angle >= self.cos_inner {
            1.0
        } else if cos_angle <= self.cos_outer {
            0.0
        } else {
            let t = (cos_angle - self.cos_outer) / (self.cos_inner - self.cos_outer);
            t * t * (3.0 - 2.0 * t)
        }
    }
}

impl Light {
    pub fn approx_eq(&self, other: &Light) -> bool {
        self.position.approx_eq(other.position) && self.intensity.approx_eq(other.intensity)
    }

    pub fn point_light(position: Tuple, intensity: Tuple) -> Self {
        Self {
            position,
            intensity,
            corner: position,
            uvec: Tuple::vector(0.0, 0.0, 0.0),
            vvec: Tuple::vector(0.0, 0.0, 0.0),
            usteps: 1,
            vsteps: 1,
            samples: 1,
            jitter_by: Sequence::new(vec![0.5]),
            spot: None,
        }
    }

    /// A point light that shines only within a cone pointed at `target`:
    /// full intensity within `inner_angle` (radians, from the axis), fading
    /// to nothing at `outer_angle`.
    pub fn spotlight(
        position: Tuple,
        target: Tuple,
        intensity: Tuple,
        inner_angle: f64,
        outer_angle: f64,
    ) -> Self {
        let mut light = Self::point_light(position, intensity);
        light.spot = Some(Spot {
            direction: (target - position).normalize(),
            cos_inner: inner_angle.cos(),
            cos_outer: outer_angle.cos(),
        });
        light
    }

    /// `corner` is one corner of the rectangle; `full_uvec` and `full_vvec`
    /// are its two edges, subdivided into `usteps` and `vsteps` cells.
    pub fn area_light(
        corner: Tuple,
        full_uvec: Tuple,
        usteps: usize,
        full_vvec: Tuple,
        vsteps: usize,
        intensity: Tuple,
    ) -> Self {
        Self {
            position: corner + full_uvec / 2.0 + full_vvec / 2.0,
            intensity,
            corner,
            uvec: full_uvec / usteps as f64,
            vvec: full_vvec / vsteps as f64,
            usteps,
            vsteps,
            samples: usteps * vsteps,
            jitter_by: Sequence::new(vec![0.5]),
            spot: None,
        }
    }

    /// Jittered world-space point inside cell `(u, v)` of the light.
    pub fn point_on_light(&self, u: usize, v: usize) -> Tuple {
        if self.samples == 1 {
            return self.position;
        }
        let start = self.jitter_by.reserve(2);
        self.cell_point(u, v, start)
    }

    /// Point in cell `(u, v)` using jitter values `start` and `start + 1`.
    fn cell_point(&self, u: usize, v: usize, start: usize) -> Tuple {
        self.corner
            + self.uvec * (u as f64 + self.jitter_by.value(start))
            + self.vvec * (v as f64 + self.jitter_by.value(start + 1))
    }

    /// Reserve jitter values for one full pass over the grid with a single
    /// atomic op (the per-value cursor was a contention hotspot across render
    /// threads) and return the sample point for cell `(u, v)` of that pass.
    /// A point light always yields its position.
    fn sample_grid(&self) -> impl Fn(usize, usize) -> Tuple + Copy + '_ {
        let start = if self.samples == 1 {
            0
        } else {
            self.jitter_by.reserve(2 * self.samples)
        };
        move |u, v| {
            if self.samples == 1 {
                self.position
            } else {
                self.cell_point(u, v, start + 2 * (v * self.usteps + u))
            }
        }
    }

    /// All sample points, row by row.
    fn sample_points(&self) -> impl Iterator<Item = Tuple> + '_ {
        let grid = self.sample_grid();
        (0..self.vsteps).flat_map(move |v| (0..self.usteps).map(move |u| grid(u, v)))
    }

    /// How much of this light reaches `point`: the fraction of its sample
    /// points not shadowed (0.0 fully shadowed, 1.0 fully lit, in between
    /// for the penumbra), scaled by the spotlight cone factor if any.
    ///
    /// Adaptive: the four corner cells are tested first, and the rest of the
    /// grid only when they disagree. Most pixels are entirely lit or entirely
    /// shadowed, so this skips nearly all shadow rays outside the penumbra.
    /// The sample positions are the same as a full pass, so a 2x2 light is
    /// unaffected.
    pub fn intensity_at(&self, point: Tuple, world: &World) -> f64 {
        self.intensity_at_time(point, world, 0.0)
    }

    /// `intensity_at` at shutter time `time` (motion blur).
    pub fn intensity_at_time(&self, point: Tuple, world: &World, time: f64) -> f64 {
        let spot = match &self.spot {
            Some(spot) => spot.factor(point - self.position),
            None => 1.0,
        };
        if spot <= 0.0 {
            return 0.0; // outside the cone: no shadow rays needed
        }
        spot * self.shadow_fraction(point, world, time)
    }

    /// Fraction of the light's sample points visible from `point`.
    fn shadow_fraction(&self, point: Tuple, world: &World, time: f64) -> f64 {
        let grid = self.sample_grid();
        let lit = |u: usize, v: usize| !world.is_shadowed_at(grid(u, v), point, time);

        if self.usteps < 2 || self.vsteps < 2 {
            let visible = (0..self.vsteps)
                .flat_map(|v| (0..self.usteps).map(move |u| (u, v)))
                .filter(|&(u, v)| lit(u, v))
                .count();
            return visible as f64 / self.samples as f64;
        }

        let (umax, vmax) = (self.usteps - 1, self.vsteps - 1);
        let corners = [(0, 0), (umax, 0), (0, vmax), (umax, vmax)];
        let corners_lit = corners.iter().filter(|&&(u, v)| lit(u, v)).count();
        if corners_lit == 0 {
            return 0.0;
        }
        if corners_lit == 4 {
            return 1.0;
        }

        let is_corner = |u: usize, v: usize| (u == 0 || u == umax) && (v == 0 || v == vmax);
        let rest_lit = (0..self.vsteps)
            .flat_map(|v| (0..self.usteps).map(move |u| (u, v)))
            .filter(|&(u, v)| !is_corner(u, v) && lit(u, v))
            .count();
        (corners_lit + rest_lit) as f64 / self.samples as f64
    }

    ///
    /// Add together the material's ambient, diffuse, and specular components,
    /// weighted by the angles between the different vectors. Diffuse and
    /// specular are averaged over the light's sample points and scaled by
    /// `intensity` (from `intensity_at`); ambient is unaffected by shadow.
    ///
    pub fn lighting(
        &self,
        shape: &dyn Shape,
        point: Tuple,        // WORLD point — still needed below
        object_point: Tuple, // NEW — group-aware object-space point, for the pattern
        eyev: Tuple,
        normalv: Tuple,
        intensity: f64,
    ) -> Tuple {
        self.lighting_with_ambient(shape, point, object_point, eyev, normalv, intensity, 1.0)
    }

    /// As `lighting`, with the ambient term scaled by `ambient_scale`.
    /// Path tracing (chapter 17) passes 0.0: the gathered indirect light
    /// replaces the constant ambient approximation of it.
    #[allow(clippy::too_many_arguments)]
    pub fn lighting_with_ambient(
        &self,
        shape: &dyn Shape,
        point: Tuple,
        object_point: Tuple,
        eyev: Tuple,
        normalv: Tuple,
        intensity: f64,
        ambient_scale: f64,
    ) -> Tuple {
        let material = shape.material();

        // Check if pattern is borrowed and use that color as input, othewise use the material color.
        let color = match &material.pattern {
            Some(pattern) => pattern.color_at_local(object_point),
            None => material.color,
        };

        // combine the surface color with the light's color/intensity
        let effective_color = color.mul(self.intensity);

        // compute the ambient contribution
        let ambient = effective_color.mul(material.ambient * ambient_scale);

        if intensity <= 0.0 {
            return ambient;
        }

        // Average the per-sample cosine factors as scalars and apply the
        // colours once at the end; it is the same sum, minus 64 tuple
        // multiplies per pixel.
        let mut diffuse_sum = 0.0;
        let mut specular_sum = 0.0;
        // powf dominates the loop, so skip it when specular can't contribute.
        let has_specular = material.specular > 0.0;

        for light_position in self.sample_points() {
            // find the direction to this sample point on the light
            let to_light = light_position - point;

            // light_dot_normal represents the cosine of the angle between the
            // light vector and the normal vector. A negative number means the
            // light is on the other side of the surface. (Cosine straight
            // from the unnormalised vector: one divide instead of four.)
            let light_dot_normal = to_light.dot(normalv) / to_light.magnitude();

            if light_dot_normal >= 0.0 {
                diffuse_sum += light_dot_normal;

                if has_specular {
                    let lightv = to_light.normalize();
                    // reflect_dot_eye represents the cosine of the angle between the
                    // reflection vector and the eye vector. A negative number means the
                    // light reflects away from the eye.
                    let reflectv = -lightv.reflect(normalv);
                    let reflect_dot_eye = reflectv.dot(eyev);

                    if reflect_dot_eye > 0.0 {
                        specular_sum += reflect_dot_eye.powf(material.shininess);
                    }
                }
            }
        }

        let n = self.samples as f64;
        let diffuse = effective_color.mul(material.diffuse * diffuse_sum / n);
        let specular = self.intensity.mul(material.specular * specular_sum / n);
        ambient + (diffuse + specular) * intensity
    }
}

use std::{fmt, ops::Mul};
impl fmt::Display for Light {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Light Pos:({:.12}, {:.12}, {:.12}) Intensity:({:.12}, {:.12}, {:.12}) Samples:{}",
            self.position.x,
            self.position.y,
            self.position.z,
            self.intensity.x,
            self.intensity.y,
            self.intensity.z,
            self.samples
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
        let result = light.lighting(&s, position, position, eyev, normalv, 1.0);
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
        let result = light.lighting(&s, position, position, eyev, normalv, 1.0);
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
        let result = light.lighting(&s, position, position, eyev, normalv, 1.0);
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
        let result = light.lighting(&s, position, position, eyev, normalv, 1.0);
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
        let result = light.lighting(&s, position, position, eyev, normalv, 1.0);
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
        let m = Material::new();
        let mut s = Sphere::new();
        s.set_material(m);
        let result = light.lighting(&s, position, position, eyev, normalv, 0.0);
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
        let position_c1 = Tuple::point(0.9, 0.0, 0.0);
        let position_c2 = Tuple::point(1.1, 0.0, 0.0);
        let mut s = Sphere::new();
        s.set_material(m);
        let c1 = light.lighting(&s, position_c1, position_c1, eyev, normalv, 1.0);
        let c2 = light.lighting(&s, position_c2, position_c2, eyev, normalv, 1.0);

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

    /// Bonus (soft shadows) - Point lights evaluate the light intensity at a given point
    #[test]
    fn test_bonus_soft_shadows_1() -> Result<(), String> {
        let w = World::default_world();
        let light = w.light.as_ref().unwrap();
        let cases = [
            (Tuple::point(0.0, 1.0001, 0.0), 1.0),
            (Tuple::point(-1.0001, 0.0, 0.0), 1.0),
            (Tuple::point(0.0, 0.0, -1.0001), 1.0),
            (Tuple::point(0.0, 0.0, 1.0001), 0.0),
            (Tuple::point(1.0001, 0.0, 0.0), 0.0),
            (Tuple::point(0.0, -1.0001, 0.0), 0.0),
            (Tuple::point(0.0, 0.0, 0.0), 0.0),
        ];
        for (pt, expected) in cases {
            let result = light.intensity_at(pt, &w);
            if result != expected {
                loge!("test_bonus_soft_shadows_1", "pt:{} got:{}", pt, result);
                return Err("Point lights evaluate the light intensity at a given point".into());
            }
        }
        Ok(())
    }

    /// Bonus (soft shadows) - lighting() uses light intensity to attenuate color
    #[test]
    fn test_bonus_soft_shadows_2() -> Result<(), String> {
        let light = Light::point_light(Tuple::point(0.0, 0.0, -10.0), WHITE);
        let mut m = Material::new();
        m.ambient = 0.1;
        m.diffuse = 0.9;
        m.specular = 0.0;
        m.color = WHITE;
        let mut s = Sphere::new();
        s.set_material(m);
        let pt = Tuple::point(0.0, 0.0, -1.0);
        let eyev = Tuple::vector(0.0, 0.0, -1.0);
        let normalv = Tuple::vector(0.0, 0.0, -1.0);
        let cases = [
            (1.0, Tuple::color(1.0, 1.0, 1.0)),
            (0.5, Tuple::color(0.55, 0.55, 0.55)),
            (0.0, Tuple::color(0.1, 0.1, 0.1)),
        ];
        for (intensity, expected) in cases {
            let result = light.lighting(&s, pt, pt, eyev, normalv, intensity);
            if !result.approx_eq(expected) {
                loge!(
                    "test_bonus_soft_shadows_2",
                    "intensity:{} got:{}",
                    intensity,
                    result
                );
                return Err("lighting() uses light intensity to attenuate color".into());
            }
        }
        Ok(())
    }

    /// Bonus (soft shadows) - Creating an area light
    #[test]
    fn test_bonus_soft_shadows_3() -> Result<(), String> {
        let corner = Tuple::point(0.0, 0.0, 0.0);
        let v1 = Tuple::vector(2.0, 0.0, 0.0);
        let v2 = Tuple::vector(0.0, 0.0, 1.0);
        let light = Light::area_light(corner, v1, 4, v2, 2, WHITE);
        let chk = light.corner.approx_eq(corner)
            && light.uvec.approx_eq(Tuple::vector(0.5, 0.0, 0.0))
            && light.usteps == 4
            && light.vvec.approx_eq(Tuple::vector(0.0, 0.0, 0.5))
            && light.vsteps == 2
            && light.samples == 8
            && light.position.approx_eq(Tuple::point(1.0, 0.0, 0.5));
        if chk {
            Ok(())
        } else {
            loge!("test_bonus_soft_shadows_3", "light:{}", light);
            Err("Creating an area light".into())
        }
    }

    /// Bonus (soft shadows) - Finding a single point on an area light
    #[test]
    fn test_bonus_soft_shadows_4() -> Result<(), String> {
        let corner = Tuple::point(0.0, 0.0, 0.0);
        let v1 = Tuple::vector(2.0, 0.0, 0.0);
        let v2 = Tuple::vector(0.0, 0.0, 1.0);
        let light = Light::area_light(corner, v1, 4, v2, 2, WHITE);
        let cases = [
            (0, 0, Tuple::point(0.25, 0.0, 0.25)),
            (1, 0, Tuple::point(0.75, 0.0, 0.25)),
            (0, 1, Tuple::point(0.25, 0.0, 0.75)),
            (2, 0, Tuple::point(1.25, 0.0, 0.25)),
            (3, 1, Tuple::point(1.75, 0.0, 0.75)),
        ];
        for (u, v, expected) in cases {
            let pt = light.point_on_light(u, v);
            if !pt.approx_eq(expected) {
                loge!("test_bonus_soft_shadows_4", "u:{} v:{} got:{}", u, v, pt);
                return Err("Finding a single point on an area light".into());
            }
        }
        Ok(())
    }

    /// Bonus (soft shadows) - The area light intensity function
    #[test]
    fn test_bonus_soft_shadows_5() -> Result<(), String> {
        let w = World::default_world();
        let corner = Tuple::point(-0.5, -0.5, -5.0);
        let v1 = Tuple::vector(1.0, 0.0, 0.0);
        let v2 = Tuple::vector(0.0, 1.0, 0.0);
        let light = Light::area_light(corner, v1, 2, v2, 2, WHITE);
        let cases = [
            (Tuple::point(0.0, 0.0, 2.0), 0.0),
            (Tuple::point(1.0, -1.0, 2.0), 0.25),
            (Tuple::point(1.5, 0.0, 2.0), 0.5),
            (Tuple::point(1.25, 1.25, 3.0), 0.75),
            (Tuple::point(0.0, 0.0, -2.0), 1.0),
        ];
        for (pt, expected) in cases {
            let result = light.intensity_at(pt, &w);
            if result != expected {
                loge!("test_bonus_soft_shadows_5", "pt:{} got:{}", pt, result);
                return Err("The area light intensity function".into());
            }
        }
        Ok(())
    }

    /// Bonus (soft shadows) - A number generator returns a cyclic sequence of numbers
    #[test]
    fn test_bonus_soft_shadows_6() -> Result<(), String> {
        let seq = Sequence::new(vec![0.1, 0.5, 1.0]);
        let chk = seq.next() == 0.1 && seq.next() == 0.5 && seq.next() == 1.0 && seq.next() == 0.1;
        if chk {
            Ok(())
        } else {
            Err("A number generator returns a cyclic sequence of numbers".into())
        }
    }

    /// Bonus (soft shadows) - Finding a single point on a jittered area light
    #[test]
    fn test_bonus_soft_shadows_7() -> Result<(), String> {
        let corner = Tuple::point(0.0, 0.0, 0.0);
        let v1 = Tuple::vector(2.0, 0.0, 0.0);
        let v2 = Tuple::vector(0.0, 0.0, 1.0);
        let mut light = Light::area_light(corner, v1, 4, v2, 2, WHITE);
        light.jitter_by = Sequence::new(vec![0.3, 0.7]);
        let cases = [
            (0, 0, Tuple::point(0.15, 0.0, 0.35)),
            (1, 0, Tuple::point(0.65, 0.0, 0.35)),
            (0, 1, Tuple::point(0.15, 0.0, 0.85)),
            (2, 0, Tuple::point(1.15, 0.0, 0.35)),
            (3, 1, Tuple::point(1.65, 0.0, 0.85)),
        ];
        for (u, v, expected) in cases {
            let pt = light.point_on_light(u, v);
            if !pt.approx_eq(expected) {
                loge!("test_bonus_soft_shadows_7", "u:{} v:{} got:{}", u, v, pt);
                return Err("Finding a single point on a jittered area light".into());
            }
        }
        Ok(())
    }

    /// Bonus (soft shadows) - The area light with jittered samples
    #[test]
    fn test_bonus_soft_shadows_8() -> Result<(), String> {
        let w = World::default_world();
        let corner = Tuple::point(-0.5, -0.5, -5.0);
        let v1 = Tuple::vector(1.0, 0.0, 0.0);
        let v2 = Tuple::vector(0.0, 1.0, 0.0);
        let cases = [
            (Tuple::point(0.0, 0.0, 2.0), 0.0),
            (Tuple::point(1.0, -1.0, 2.0), 0.5),
            (Tuple::point(1.5, 0.0, 2.0), 0.75),
            (Tuple::point(1.25, 1.25, 3.0), 0.75),
            (Tuple::point(0.0, 0.0, -2.0), 1.0),
        ];
        for (pt, expected) in cases {
            // The book resets the sequence for each example.
            let mut light = Light::area_light(corner, v1, 2, v2, 2, WHITE);
            light.jitter_by = Sequence::new(vec![0.7, 0.3, 0.9, 0.1, 0.5]);
            let result = light.intensity_at(pt, &w);
            if result != expected {
                loge!("test_bonus_soft_shadows_8", "pt:{} got:{}", pt, result);
                return Err("The area light with jittered samples".into());
            }
        }
        Ok(())
    }

    /// Bonus (soft shadows) - lighting() samples the area light
    #[test]
    fn test_bonus_soft_shadows_9() -> Result<(), String> {
        let corner = Tuple::point(-0.5, -0.5, -5.0);
        let v1 = Tuple::vector(1.0, 0.0, 0.0);
        let v2 = Tuple::vector(0.0, 1.0, 0.0);
        let mut light = Light::area_light(corner, v1, 2, v2, 2, WHITE);
        light.jitter_by = Sequence::new(vec![0.5]);
        let mut m = Material::new();
        m.ambient = 0.1;
        m.diffuse = 0.9;
        m.specular = 0.0;
        m.color = WHITE;
        let mut s = Sphere::new();
        s.set_material(m);
        let eye = Tuple::point(0.0, 0.0, -5.0);
        let cases = [
            (
                Tuple::point(0.0, 0.0, -1.0),
                Tuple::color(0.9965, 0.9965, 0.9965),
            ),
            (
                Tuple::point(
                    0.0,
                    std::f64::consts::FRAC_1_SQRT_2,
                    -std::f64::consts::FRAC_1_SQRT_2,
                ),
                Tuple::color(0.6232, 0.6232, 0.6232),
            ),
        ];
        for (pt, expected) in cases {
            let eyev = (eye - pt).normalize();
            let normalv = Tuple::vector(pt.x, pt.y, pt.z);
            let result = light.lighting(&s, pt, pt, eyev, normalv, 1.0);
            let close = |a: f64, b: f64| (a - b).abs() < 1e-4;
            if !(close(result.x, expected.x)
                && close(result.y, expected.y)
                && close(result.z, expected.z))
            {
                loge!("test_bonus_soft_shadows_9", "pt:{} got:{}", pt, result);
                return Err("lighting() samples the area light".into());
            }
        }
        Ok(())
    }

    /// Chap 17 - A spotlight is full strength inside its inner cone, off
    /// outside its outer cone, and fades smoothly in between
    #[test]
    fn test_chap_17_15() -> Result<(), String> {
        use std::f64::consts::PI;
        let light = Light::spotlight(
            Tuple::point(0.0, 5.0, 0.0),
            Tuple::point(0.0, 0.0, 0.0),
            WHITE,
            PI / 8.0, // 22.5° inner
            PI / 4.0, // 45° outer
        );
        let w = World::new(); // nothing to cast shadows
        // On the axis, and at 20° off it: full.
        let on_axis = light.intensity_at(Tuple::point(0.0, 0.0, 0.0), &w);
        let inside =
            light.intensity_at(Tuple::point(5.0 * (20f64).to_radians().tan(), 0.0, 0.0), &w);
        // At 60° off the axis: nothing.
        let outside =
            light.intensity_at(Tuple::point(5.0 * (60f64).to_radians().tan(), 0.0, 0.0), &w);
        // At 33.75° (halfway between the angles): strictly between 0 and 1.
        let edge = light.intensity_at(
            Tuple::point(5.0 * (33.75f64).to_radians().tan(), 0.0, 0.0),
            &w,
        );
        // Behind the light: nothing.
        let behind = light.intensity_at(Tuple::point(0.0, 10.0, 0.0), &w);

        let chk = on_axis == 1.0
            && inside == 1.0
            && outside == 0.0
            && edge > 0.0
            && edge < 1.0
            && behind == 0.0;
        if chk {
            Ok(())
        } else {
            loge!(
                "test_chap_17_15",
                "on_axis:{on_axis} inside:{inside} outside:{outside} edge:{edge} behind:{behind}"
            );
            Err("Spotlight cone factor".into())
        }
    }

    /// Chap 17 - A spotlight still casts shadows inside its cone
    #[test]
    fn test_chap_17_16() -> Result<(), String> {
        use std::f64::consts::PI;
        // Default world: unit sphere at the origin. Light above, pointing down.
        let mut w = World::default_world();
        w.set_light(Light::spotlight(
            Tuple::point(0.0, 5.0, 0.0),
            Tuple::point(0.0, 0.0, 0.0),
            WHITE,
            PI / 6.0,
            PI / 4.0,
        ));
        let light = w.light.as_ref().unwrap();
        let below_sphere = light.intensity_at(Tuple::point(0.0, -2.0, 0.0), &w);
        let above_sphere = light.intensity_at(Tuple::point(0.0, 1.5, 0.0), &w);
        if below_sphere == 0.0 && above_sphere == 1.0 {
            Ok(())
        } else {
            loge!(
                "test_chap_17_16",
                "below:{below_sphere} above:{above_sphere}"
            );
            Err("Spotlight must be shadowed by the sphere".into())
        }
    }
}
