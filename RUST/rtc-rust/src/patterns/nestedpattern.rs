//! Nested pattern definition
//!
//! A composite pattern whose two "colors" are themselves patterns.
//! Alternates between sub-pattern `a` and `b` modulo 2 in x (like a stripe),
//! but instead of returning a solid color it recurses into the chosen child.

use crate::{
    matrix::Matrix4,
    pattern::{Pattern, PatternData},
    tuple::Tuple,
};

// NOT Copy: a/b own heap-allocated trait objects.
#[derive(Debug, Clone)]
pub struct NestedPattern {
    pub data: PatternData,
    pub a: Box<dyn Pattern>,
    pub b: Box<dyn Pattern>,
}

impl NestedPattern {
    pub fn new(a: Box<dyn Pattern>, b: Box<dyn Pattern>) -> Self {
        Self {
            data: PatternData::new(),
            a,
            b,
        }
    }
}

impl Pattern for NestedPattern {
    fn clone_box(&self) -> Box<dyn Pattern> {
        // self.clone() (not *self) — NestedPattern is not Copy.
        Box::new(self.clone())
    }

    ///
    /// `point` is already in THIS pattern's space (color_at_local stripped our
    /// transform). Select a child by the stripe rule, then delegate via
    /// `color_at_local` so the child applies ITS OWN transform before evaluating.
    ///
    fn color_at(&self, point: Tuple) -> Tuple {
        if (point.x.floor() as i32).rem_euclid(2) == 0 {
            self.a.color_at_local(point)
        } else {
            self.b.color_at_local(point)
        }
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

    /// Chap 10 (bonus) - A nested pattern selects the correct child by x,
    /// then returns that child's color.
    #[test]
    fn test_chap_10_nested_selection() -> Result<(), String> {
        let red = Tuple::color(1.0, 0.0, 0.0);
        let green = Tuple::color(0.0, 1.0, 0.0);
        let a = StripePattern::new(WHITE, BLACK);
        let b = StripePattern::new(red, green);
        let nested = NestedPattern::new(Box::new(a), Box::new(b));

        // x in [0,1) -> child a; a at x=0 -> WHITE
        let c0 = nested.color_at(Tuple::point(0.0, 0.0, 0.0));
        // x in [1,2) -> child b; b at x=1 -> green (b's second color)
        let c1 = nested.color_at(Tuple::point(1.0, 0.0, 0.0));

        let chk = c0.approx_eq(WHITE) && c1.approx_eq(green);
        if chk {
            Ok(())
        } else {
            Err(format!("nested selection: c0={c0} c1={c1}"))
        }
    }

    /// Chap 10 (bonus) - The chosen child's OWN transform is applied during
    /// recursion (proves color_at recurses through color_at_local, not color_at).
    #[test]
    fn test_chap_10_nested_child_transform() -> Result<(), String> {
        // Child a is translated +0.5 in x. Its inverse shifts the sample -0.5,
        // so at parent x=0.4 the child sees x=-0.1 -> floor -1 -> BLACK.
        // Without the child transform it would see x=0.4 -> WHITE.
        let mut a = StripePattern::new(WHITE, BLACK);
        a.set_transform(Matrix4::translation(0.5, 0.0, 0.0));
        let b = StripePattern::new(WHITE, BLACK);
        let nested = NestedPattern::new(Box::new(a), Box::new(b));

        let c = nested.color_at(Tuple::point(0.4, 0.0, 0.0)); // routes to child a
        if c.approx_eq(BLACK) {
            Ok(())
        } else {
            Err(format!(
                "child transform not applied: expected BLACK, got {c}"
            ))
        }
    }

    /// Chap 10 (bonus) - End to end through color_at_shape: object transform +
    /// nested pattern still resolves to the right child color.
    #[test]
    fn test_chap_10_nested_via_shape() -> Result<(), String> {
        let mut shape = Sphere::new();
        shape.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));
        let a = StripePattern::new(WHITE, BLACK);
        let b = StripePattern::new(WHITE, BLACK);
        let nested = NestedPattern::new(Box::new(a), Box::new(b));

        // world 1.5 -> object 0.75 -> parent picks child a -> a at 0.75 -> WHITE
        let c = nested.color_at_shape(&shape, Tuple::point(1.5, 0.0, 0.0));
        if c.approx_eq(WHITE) {
            Ok(())
        } else {
            Err(format!("nested via shape: expected WHITE, got {c}"))
        }
    }
}
