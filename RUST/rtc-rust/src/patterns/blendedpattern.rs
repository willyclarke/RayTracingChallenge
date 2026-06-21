//! Blended pattern definition
//!
//! A composite pattern whose two "colors" are themselves patterns.
//! Evaluates BOTH sub-patterns `a` and `b` at the point and averages them
//! 50/50, recursing into each child so its own transform is applied.

use crate::{
    matrix::Matrix4,
    pattern::{Pattern, PatternData},
    tuple::Tuple,
};

// NOT Copy: a/b own heap-allocated trait objects.
#[derive(Debug, Clone)]
pub struct BlendedPattern {
    pub data: PatternData,
    pub a: Box<dyn Pattern>,
    pub b: Box<dyn Pattern>,
}

impl BlendedPattern {
    pub fn new(a: Box<dyn Pattern>, b: Box<dyn Pattern>) -> Self {
        Self {
            data: PatternData::new(),
            a,
            b,
        }
    }
}

impl Pattern for BlendedPattern {
    fn clone_box(&self) -> Box<dyn Pattern> {
        // self.clone() (not *self) — BlendedPattern is not Copy.
        Box::new(self.clone())
    }

    ///
    /// `point` is already in THIS pattern's space (color_at_local stripped our
    /// transform). Evaluate BOTH children via `color_at_local` (so each applies
    /// ITS OWN transform), then return the 50/50 average of the two colors.
    ///
    fn color_at(&self, point: Tuple) -> Tuple {
        self.a.color_at_local(point) * 0.5 + self.b.color_at_local(point) * 0.5
    }

    fn data(&self) -> &PatternData {
        &self.data
    }

    fn data_mut(&mut self) -> &mut PatternData {
        &mut self.data
    }

    fn set_transform(&mut self, m: Matrix4) {
        self.data_mut().set_transform(m);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::patterns::stripepattern::StripePattern;
    use crate::shape::Shape;
    use crate::shapes::sphere::Sphere;
    use crate::tuple::colors::{BLACK, WHITE};
    use crate::{loge, tuple::Tuple};

    /// Chap 10 (bonus) - A blended pattern averages both children 50/50.
    #[test]
    fn test_chap_10_blended_selection() -> Result<(), String> {
        let red = Tuple::color(1.0, 0.0, 0.0);
        let green = Tuple::color(0.0, 1.0, 0.0);
        let a = StripePattern::new(WHITE, BLACK);
        let b = StripePattern::new(red, green);
        let blended = BlendedPattern::new(Box::new(a), Box::new(b));

        // at x=0: a -> WHITE, b -> red; blend = half WHITE + half red.
        let c0 = blended.color_at(Tuple::point(0.0, 0.0, 0.0));
        let expect_at_0 = WHITE * 0.5 + red * 0.5;

        // at x=1: a -> BLACK, b -> green; blend = half BLACK + half green.
        let c1 = blended.color_at(Tuple::point(1.0, 0.0, 0.0));
        let expect_at_1 = BLACK * 0.5 + green * 0.5;

        let chk = c0.approx_eq(expect_at_0) && c1.approx_eq(expect_at_1);
        if chk {
            Ok(())
        } else {
            loge!("test_chap_10_blended_selection", "c0 : {}", c0);
            loge!("test_chap_10_blended_selection", "c1 : {}", c1);
            Err(format!("blended selection: c0={c0} c1={c1}"))
        }
    }

    /// Chap 10 (bonus) - Each child's OWN transform is applied during recursion
    /// (proves color_at recurses through color_at_local, not color_at).
    #[test]
    fn test_chap_10_blended_child_transform() -> Result<(), String> {
        // Child a is translated +0.5 in x. Its inverse shifts the sample -0.5,
        // so at x=0.4 child a sees x=-0.1 -> floor -1 -> BLACK.
        // Child b is untransformed, so at x=0.4 it sees x=0.4 -> WHITE.
        // Blend = half BLACK + half WHITE = WHITE * 0.5.
        // Without the child transform, a would also see WHITE -> blend WHITE.
        let mut a = StripePattern::new(WHITE, BLACK);
        a.set_transform(Matrix4::translation(0.5, 0.0, 0.0));
        let b = StripePattern::new(WHITE, BLACK);
        let blended = BlendedPattern::new(Box::new(a), Box::new(b));

        let c = blended.color_at(Tuple::point(0.4, 0.0, 0.0));
        if c.approx_eq(WHITE * 0.5) {
            Ok(())
        } else {
            loge!("test_chap_10_blended_child_transform", "c : {}", c);
            Err(format!(
                "child transform not applied: expected WHITE * 0.5, got {c}"
            ))
        }
    }

    /// Chap 10 (bonus) - End to end through color_at_shape: object transform +
    /// blended pattern still produces the correct averaged color.
    #[test]
    fn test_chap_10_blended_via_shape() -> Result<(), String> {
        let mut shape = Sphere::new();
        shape.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));
        let a = StripePattern::new(WHITE, BLACK);
        let b = StripePattern::new(WHITE, BLACK);
        let blended = BlendedPattern::new(Box::new(a), Box::new(b));

        // world 1.5 -> object 0.75 -> both children -> WHITE; blend = WHITE
        let c = blended.color_at_shape(&shape, Tuple::point(1.5, 0.0, 0.0));
        if c.approx_eq(WHITE) {
            Ok(())
        } else {
            loge!("test_chap_10_blended_via_shape", "c : {}", c);
            Err(format!("blended via shape: expected WHITE, got {c}"))
        }
    }
}
