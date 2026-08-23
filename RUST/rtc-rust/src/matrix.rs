//! Matrix defintion an implementation
//!

use crate::log::*;
use crate::math::approx_eq;
use crate::ray::Ray;
use crate::tuple::Tuple;
use std::fmt;

#[derive(Debug, Copy, Clone)]
pub struct Matrix2 {
    data: [[f64; 2]; 2],
}

#[derive(Debug, Copy, Clone)]
pub struct Matrix3 {
    data: [[f64; 3]; 3],
}

#[derive(Debug, Copy, Clone)]
pub struct Matrix4 {
    data: [[f64; 4]; 4],
}

impl Matrix2 {
    pub fn approx_eq(self, other: Self) -> bool {
        for r in 0..2 {
            for c in 0..2 {
                if !approx_eq(self[(r, c)], other[(r, c)]) {
                    return false;
                }
            }
        }
        true
    }

    pub fn determinant(&self) -> f64 {
        self[(0, 0)] * self[(1, 1)] - self[(0, 1)] * self[(1, 0)]
    }

    pub fn new(data: [[f64; 2]; 2]) -> Self {
        Self { data }
    }
}

impl Matrix3 {
    pub fn approx_eq(self, other: Self) -> bool {
        for r in 0..3 {
            for c in 0..3 {
                if !approx_eq(self[(r, c)], other[(r, c)]) {
                    return false;
                }
            }
        }
        true
    }

    /// Cofactor - like minor but sign may be flipped
    ///
    /// The sign is negated when (row + col) is odd.
    ///
    /// # Examples
    /// ```
    /// # use rtc_rust::matrix::Matrix3;
    ///
    /// let a = Matrix3::new([
    ///     [3.0,  5.0,  0.0],
    ///     [2.0, -1.0, -7.0],
    ///     [6.0, -1.0,  5.0],
    /// ]);
    ///
    /// assert!(a.cofactor(0, 0) == -12.0);  // minor is -12, sign stays
    /// assert!(a.cofactor(1, 0) == -25.0);  // minor is 25, sign flips
    /// ```
    pub fn cofactor(&self, row: usize, col: usize) -> f64 {
        let mi = self.minor(row, col);
        if (row + col).is_multiple_of(2) {
            mi
        } else {
            -mi
        }
    }

    /// Determinant via cofactor expansion along row 0
    ///
    /// # Examples
    /// ```
    /// # use rtc_rust::matrix::Matrix3;
    ///
    /// let a = Matrix3::new([
    ///     [1.0, 2.0,  6.0],
    ///     [-5.0, 8.0, -4.0],
    ///     [2.0, 6.0,  4.0],
    /// ]);
    ///
    /// assert!(a.cofactor(0, 0) == 56.0);
    /// assert!(a.cofactor(0, 1) == 12.0);
    /// assert!(a.cofactor(0, 2) == -46.0);
    /// assert!(a.determinant() == -196.0);
    /// ```
    pub fn determinant(&self) -> f64 {
        (0..3).map(|c| self.data[0][c] * self.cofactor(0, c)).sum()
    }

    /// A minor is the determinant of the submatrix after removing row and col
    ///
    ///
    /// # Examples
    /// ```
    /// # use rtc_rust::matrix::Matrix3;
    ///
    /// let a = Matrix3::new([[0.0, 0.0, 1.0], [0.0, 1.0, 1.0],[0.0, 0.0, 1.0]]);
    /// let minor = a.minor(0, 0);
    ///
    /// assert!(minor == 1.0);
    /// ```
    pub fn minor(&self, row: usize, col: usize) -> f64 {
        let sm = self.submatrix(row, col);
        sm.determinant()
    }

    pub fn new(data: [[f64; 3]; 3]) -> Self {
        Self { data }
    }

    /// Create a SubMatrix by deleting row and col
    ///
    /// # Examples
    /// ```
    /// # use rtc_rust::matrix::Matrix3;
    ///
    /// let result = Matrix3::new([[0.0, 0.0, 1.0], [0.0, 0.0, 1.0],[0.0, 0.0, 1.0]]);
    /// let sm = result.submatrix(0, 2);
    ///
    /// assert!(true);
    /// ```
    pub fn submatrix(&self, row: usize, col: usize) -> Matrix2 {
        debug_assert!(row < 3);
        debug_assert!(col < 3);

        let mut out = [[0.0; 2]; 2];
        let mut out_row = 0;

        for r in 0..3 {
            if r == row {
                continue;
            }

            let mut out_col = 0;

            for c in 0..3 {
                if c == col {
                    continue;
                }

                out[out_row][out_col] = self.data[r][c];
                out_col += 1;
            }
            out_row += 1;
        }

        Matrix2 { data: out }
    }
}

impl Matrix4 {
    /// CTOR for Matrix4
    ///
    /// # Examples
    /// ```
    /// # use rtc_rust::matrix::Matrix4;
    /// let m = Matrix4::new([
    ///     [1.0, 2.0, 3.0, 4.0],
    ///     [5.5, 6.5, 7.5, 8.5],
    ///     [9.0, 10.0, 11.0, 12.0],
    ///     [13.5, 14.5, 15.5, 16.5],
    /// ]);
    ///
    /// assert!(m[(0, 0)] == 1.0_f64);
    /// assert!(m[(3, 3)] == 16.5_f64);
    /// ```
    pub fn new(data: [[f64; 4]; 4]) -> Self {
        Self { data }
    }

    pub fn approx_eq(self, other: Self) -> bool {
        for r in 0..4 {
            for c in 0..4 {
                if !approx_eq(self[(r, c)], other[(r, c)]) {
                    return false;
                }
            }
        }
        true
    }

    pub fn minor(&self, row: usize, col: usize) -> f64 {
        let sm = self.submatrix(row, col);
        sm.determinant()
    }

    pub fn cofactor(&self, row: usize, col: usize) -> f64 {
        let mi = self.minor(row, col);
        if (row + col).is_multiple_of(2) {
            mi
        } else {
            -mi
        }
    }

    /// Determinant via cofactor expansion along row 0
    ///
    /// # Examples
    /// ```
    /// # use rtc_rust::matrix::Matrix4;
    ///
    /// let a = Matrix4::new([
    ///     [-2.0, -8.0,  3.0,  5.0],
    ///     [-3.0,  1.0,  7.0,  3.0],
    ///     [ 1.0,  2.0, -9.0,  6.0],
    ///     [-6.0,  7.0,  7.0, -9.0],
    /// ]);
    ///
    /// assert!(a.cofactor(0, 0) == 690.0);
    /// assert!(a.cofactor(0, 1) == 447.0);
    /// assert!(a.cofactor(0, 2) == 210.0);
    /// assert!(a.cofactor(0, 3) == 51.0);
    /// assert!(a.determinant() == -4071.0);
    /// ```
    pub fn determinant(&self) -> f64 {
        (0..4).map(|c| self.data[0][c] * self.cofactor(0, c)).sum()
    }

    /// Identity matrix
    ///
    /// # Examples
    /// ```
    /// # use rtc_rust::matrix::Matrix4;
    ///
    /// let i = Matrix4::identity();
    ///
    /// assert!(i[(0, 0)] == 1.0_f64);
    /// assert!(i[(1, 1)] == 1.0_f64);
    /// assert!(i[(2, 2)] == 1.0_f64);
    /// assert!(i[(3, 3)] == 1.0_f64);
    /// ```
    pub fn identity() -> Self {
        Self {
            data: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        }
    }

    pub fn is_invertible(&self) -> bool {
        !approx_eq(self.determinant(), 0.0)
    }

    /// Compute the inverse of the matrix.
    ///
    /// Returns None if the matrix is not invertible (determinant ≈ 0).
    ///
    /// # Examples
    /// ```
    /// # use rtc_rust::matrix::Matrix4;
    ///
    /// let a = Matrix4::new([
    ///     [-5.0,  2.0,  6.0, -8.0],
    ///     [ 1.0, -5.0,  1.0,  8.0],
    ///     [ 7.0,  7.0, -6.0, -7.0],
    ///     [ 1.0, -3.0,  7.0,  4.0],
    /// ]);
    ///
    /// let b = a.inverse().unwrap();
    /// // a * inverse(a) should equal identity
    /// assert!((a * b).approx_eq(Matrix4::identity()));
    /// ```
    pub fn inverse(&self) -> Option<Matrix4> {
        let det = self.determinant();
        if approx_eq(det, 0.0) {
            return None;
        }
        let mut out = [[0.0; 4]; 4];

        // Do not let clippy fool you, using a range based loop could invalidate the cache.
        #[allow(clippy::needless_range_loop)]
        for row in 0..4 {
            for col in 0..4 {
                // Transposed assignment: [col][row] handles the transpose
                out[col][row] = self.cofactor(row, col) / det;
            }
        }

        Some(Matrix4 { data: out })
    }

    pub fn submatrix(&self, row: usize, col: usize) -> Matrix3 {
        debug_assert!(row < 4);
        debug_assert!(col < 4);

        let mut out = [[0.0; 3]; 3];
        let mut out_row = 0;

        for r in 0..4 {
            if r == row {
                continue;
            }

            let mut out_col = 0;

            for c in 0..4 {
                if c == col {
                    continue;
                }

                out[out_row][out_col] = self.data[r][c];
                out_col += 1;
            }
            out_row += 1;
        }

        Matrix3 { data: out }
    }

    pub fn rotation_x(r: f64) -> Matrix4 {
        let mut out = Matrix4::identity().data;
        out[1][1] = r.cos();
        out[1][2] = -r.sin();
        out[2][1] = r.sin();
        out[2][2] = r.cos();
        Matrix4 { data: out }
    }

    pub fn rotation_y(r: f64) -> Matrix4 {
        let mut out = Matrix4::identity().data;
        out[0][0] = r.cos();
        out[2][0] = -r.sin();
        out[0][2] = r.sin();
        out[2][2] = r.cos();
        Matrix4 { data: out }
    }

    pub fn rotation_z(r: f64) -> Matrix4 {
        let mut out = Matrix4::identity().data;
        out[0][0] = r.cos();
        out[0][1] = -r.sin();
        out[1][0] = r.sin();
        out[1][1] = r.cos();
        Matrix4 { data: out }
    }

    pub fn shearing(x_y: f64, x_z: f64, y_x: f64, y_z: f64, z_x: f64, z_y: f64) -> Matrix4 {
        let mut out = Matrix4::identity().data;
        out[0][1] = x_y;
        out[0][2] = x_z;
        out[1][0] = y_x;
        out[1][2] = y_z;
        out[2][0] = z_x;
        out[2][1] = z_y;
        Matrix4 { data: out }
    }

    /// CTOR for scaling matrix
    ///
    /// # return identity() with translation in third column.
    ///
    /// # Examples
    /// ```
    /// # use rtc_rust::matrix::Matrix4;
    ///
    /// let transform = Matrix4::scaling(5.0, -3.0, 2.0);
    /// ```
    pub fn scaling(x: f64, y: f64, z: f64) -> Matrix4 {
        let mut out = Matrix4::identity().data;
        out[0][0] = x;
        out[1][1] = y;
        out[2][2] = z;
        Matrix4 { data: out }
    }

    /// CTOR for translation matrix
    ///
    /// # return identity() with translation in third column.
    ///
    /// # Examples
    /// ```
    /// # use rtc_rust::matrix::Matrix4;
    ///
    /// let transform = Matrix4::translation(5.0, -3.0, 2.0);
    /// ```
    pub fn translation(x: f64, y: f64, z: f64) -> Matrix4 {
        let mut out = Matrix4::identity().data;
        out[0][3] = x;
        out[1][3] = y;
        out[2][3] = z;
        Matrix4 { data: out }
    }

    pub fn transpose(&self) -> Matrix4 {
        // ✔ a is copied
        // ✔ cost is negligible
        // ✔ often helps optimization
        // ✔ idiomatic for small fixed-size math types
        let a = self.data;
        let mut out = [[0.0; 4]; 4];

        for row in 0..4 {
            for col in 0..4 {
                out[row][col] = a[col][row];
            }
        }

        Matrix4 { data: out }
    }
}

impl fmt::Display for Matrix4 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let width = 17;
        let prec = 15;
        let coldelta = width + 4;

        writeln!(
            f,
            "\n{}COL: {:3}{:coldelta$}{:coldelta$}{:coldelta$}{}",
            Color::Cyan,
            0,
            1,
            2,
            3,
            Color::Reset
        )?;

        let eps = crate::math::EPSILON;

        for row in 0..4 {
            write!(f, "{}{}{} ", Color::Cyan, row, Color::Reset)?;
            write!(f, "{}||{} ", Color::Cyan, Color::Reset)?;

            for col in 0..4 {
                let v = self[(row, col)];
                let space = if v >= 0.0 { "  " } else { " " };

                let color = if v.abs() < eps {
                    Color::Yellow
                } else if (v - 1.0).abs() < eps {
                    Color::Green
                } else if v < 0.0 {
                    Color::Red
                } else {
                    Color::Reset
                };

                f.write_fmt(core::format_args!(
                    "{color}{space}{value:>width$.precision$}{reset} |",
                    color = color,
                    space = space,
                    value = v,
                    width = width,
                    precision = prec,
                    reset = Color::Reset,
                ))?;
            }

            writeln!(f, "{}|{}", Color::Cyan, Color::Reset)?;
        }

        Ok(())
    }
}

// use std::ops::Add;
// use std::ops::Div;
use std::ops::Mul;
// use std::ops::Neg;
// use std::ops::Sub;
use std::ops::{Index, IndexMut};

impl Index<(usize, usize)> for Matrix2 {
    type Output = f64;

    fn index(&self, (row, col): (usize, usize)) -> &Self::Output {
        &self.data[row][col]
    }
}

impl IndexMut<(usize, usize)> for Matrix2 {
    fn index_mut(&mut self, (row, col): (usize, usize)) -> &mut Self::Output {
        &mut self.data[row][col]
    }
}

impl Index<(usize, usize)> for Matrix3 {
    type Output = f64;

    fn index(&self, (row, col): (usize, usize)) -> &Self::Output {
        &self.data[row][col]
    }
}

impl IndexMut<(usize, usize)> for Matrix3 {
    fn index_mut(&mut self, (row, col): (usize, usize)) -> &mut Self::Output {
        &mut self.data[row][col]
    }
}

impl Index<(usize, usize)> for Matrix4 {
    type Output = f64;

    fn index(&self, (row, col): (usize, usize)) -> &Self::Output {
        &self.data[row][col]
    }
}

impl IndexMut<(usize, usize)> for Matrix4 {
    fn index_mut(&mut self, (row, col): (usize, usize)) -> &mut Self::Output {
        &mut self.data[row][col]
    }
}

impl Mul<Tuple> for Matrix4 {
    type Output = Tuple;

    fn mul(self, rhs: Tuple) -> Tuple {
        &self * rhs
    }
}

/// By-reference variant so hot paths (ray transforms in the BVH walk) don't
/// copy the 128-byte matrix for every multiply.
impl Mul<Tuple> for &Matrix4 {
    type Output = Tuple;

    fn mul(self, rhs: Tuple) -> Tuple {
        let a = &self.data;
        Tuple::new(
            a[0][0] * rhs.x + a[0][1] * rhs.y + a[0][2] * rhs.z + a[0][3] * rhs.w,
            a[1][0] * rhs.x + a[1][1] * rhs.y + a[1][2] * rhs.z + a[1][3] * rhs.w,
            a[2][0] * rhs.x + a[2][1] * rhs.y + a[2][2] * rhs.z + a[2][3] * rhs.w,
            a[3][0] * rhs.x + a[3][1] * rhs.y + a[3][2] * rhs.z + a[3][3] * rhs.w,
        )
    }
}

impl Matrix4 {
    /// Transform a ray without the `w` work: an affine transform leaves the
    /// bottom row as `0 0 0 1`, so a point picks up the translation column and
    /// a vector ignores it.
    pub fn transform_ray(&self, ray: &Ray) -> Ray {
        let a = &self.data;
        let o = ray.origin;
        let d = ray.direction;
        Ray::new(
            Tuple::point(
                a[0][0] * o.x + a[0][1] * o.y + a[0][2] * o.z + a[0][3],
                a[1][0] * o.x + a[1][1] * o.y + a[1][2] * o.z + a[1][3],
                a[2][0] * o.x + a[2][1] * o.y + a[2][2] * o.z + a[2][3],
            ),
            Tuple::vector(
                a[0][0] * d.x + a[0][1] * d.y + a[0][2] * d.z,
                a[1][0] * d.x + a[1][1] * d.y + a[1][2] * d.z,
                a[2][0] * d.x + a[2][1] * d.y + a[2][2] * d.z,
            ),
        )
    }
}

impl Mul<Matrix4> for Matrix4 {
    type Output = Matrix4;

    fn mul(self, rhs: Matrix4) -> Matrix4 {
        // ✔ a and b are copied
        // ✔ cost is negligible
        // ✔ often helps optimization
        // ✔ idiomatic for small fixed-size math types
        let a = self.data;
        let b = rhs.data;
        let mut out = [[0.0; 4]; 4];

        for row in 0..4 {
            let [a0, a1, a2, a3] = a[row];

            for col in 0..4 {
                out[row][col] = a0 * b[0][col] + a1 * b[1][col] + a2 * b[2][col] + a3 * b[3][col];
            }
        }

        Matrix4 { data: out }
    }
}

#[cfg(test)]
mod tests {
    use core::f64;

    use super::*;
    use crate::canvas::Canvas;
    use crate::{logd, loge, logi, tuple::Tuple};

    /// Chap 3 - Constructing and inspecting a 4x4 matrix
    #[test]
    fn test_chap_3_1() -> Result<(), String> {
        let m = Matrix4::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.5, 6.5, 7.5, 8.5],
            [9.0, 10.0, 11.0, 12.0],
            [13.5, 14.5, 15.5, 16.5],
        ]);
        logi!("test_chap_3_1", "m:{:?}", m);
        logi!("test_chap_3_1", "\nm:\n{}\n", m);

        let chk = true;

        let chk = chk && approx_eq(m[(0, 0)], 1.0);
        let chk = chk && approx_eq(m[(0, 1)], 2.0);
        let chk = chk && approx_eq(m[(0, 2)], 3.0);
        let chk = chk && approx_eq(m[(0, 3)], 4.0);

        let chk = chk && approx_eq(m[(1, 0)], 5.5);
        let chk = chk && approx_eq(m[(1, 1)], 6.5);
        let chk = chk && approx_eq(m[(1, 2)], 7.5);
        let chk = chk && approx_eq(m[(1, 3)], 8.5);

        let chk = chk && approx_eq(m[(2, 0)], 9.0);
        let chk = chk && approx_eq(m[(2, 1)], 10.0);
        let chk = chk && approx_eq(m[(2, 2)], 11.0);
        let chk = chk && approx_eq(m[(2, 3)], 12.0);

        let chk = chk && approx_eq(m[(3, 0)], 13.5);
        let chk = chk && approx_eq(m[(3, 1)], 14.5);
        let chk = chk && approx_eq(m[(3, 2)], 15.5);
        let chk = chk && approx_eq(m[(3, 3)], 16.5);

        if chk {
            Ok(())
        } else {
            logi!("test_chap_3_1", "m:{:?}", m);
            logd!("test_chap_3_1", "m:{:?}", m);
            loge!("test_chap_3_1", "m:{:?}", m);
            Err("Constructing and inspecting a 4x4 matrix".into())
        }
    }

    /// Chap 3 - Matrix equality with identical matrices
    #[test]
    fn test_chap_3_3() -> Result<(), String> {
        let a = Matrix4::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 10.0, 11.0, 12.0],
            [13.0, 14.0, 15.0, 16.0],
        ]);
        let b = Matrix4::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 10.0, 11.0, 12.0],
            [13.0, 14.0, 15.0, 16.0],
        ]);
        logi!("test_chap_3_3", "\na:\n{}\n", a);
        logi!("test_chap_3_3", "\nb:\n{}\n", b);

        let chk = a.approx_eq(b);
        if chk {
            Ok(())
        } else {
            Err("Matrix equality with identical matrices".into())
        }
    }

    // Chap 3 - Matrix equality with different matrices
    #[test]
    fn test_chap_3_4() -> Result<(), String> {
        let a = Matrix4::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 8.0, 7.0, 6.0],
            [5.0, 4.0, 3.0, 2.0],
        ]);
        let b = Matrix4::new([
            [2.0, 3.0, 4.0, 5.0],
            [6.0, 7.0, 8.0, 9.0],
            [8.0, 7.0, 6.0, 5.0],
            [4.0, 3.0, 2.0, 1.0],
        ]);
        logi!("test_chap_3_4", "\na:\n{}\n", a);
        logi!("test_chap_3_4", "\nb:\n{}\n", b);

        let chk = !a.approx_eq(b);
        if chk {
            Ok(())
        } else {
            Err("Matrix equality with different matrices".into())
        }
    }

    // Chap 3 - Multiplying two matrices
    #[test]
    fn test_chap_3_5() -> Result<(), String> {
        let a = Matrix4::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 8.0, 7.0, 6.0],
            [5.0, 4.0, 3.0, 2.0],
        ]);
        let b = Matrix4::new([
            [-2.0, 1.0, 2.0, 3.0],
            [3.0, 2.0, 1.0, -1.0],
            [4.0, 3.0, 6.0, 5.0],
            [1.0, 2.0, 7.0, 8.0],
        ]);
        let expect = Matrix4::new([
            [20.0, 22.0, 50.0, 48.0],
            [44.0, 54.0, 114.0, 108.0],
            [40.0, 58.0, 110.0, 102.0],
            [16.0, 26.0, 46.0, 42.0],
        ]);
        let result = a * b;
        logi!("test_chap_3_5", "\na:\n{}\n", a);
        logi!("test_chap_3_5", "\nb:\n{}\n", b);
        logi!("test_chap_3_5", "\nexpect:\n{}\n", expect);
        logi!("test_chap_3_5", "\nresult:\n{}\n", result);

        let chk = a.mul(b).approx_eq(expect);
        if chk {
            Ok(())
        } else {
            Err("Multiplying two matrices".into())
        }
    }

    /// Chap 3 - A matrix multiplied by a tuple
    #[test]
    fn test_chap_3_6() -> Result<(), String> {
        let a = Matrix4::new([
            [1.0, 2.0, 3.0, 4.0],
            [2.0, 4.0, 4.0, 2.0],
            [8.0, 6.0, 4.0, 1.0],
            [0.0, 0.0, 0.0, 1.0],
        ]);
        let b = Tuple::point(1.0, 2.0, 3.0);
        let result = a * b;
        let chk = Tuple::point(18.0, 24.0, 33.0).approx_eq(result);
        if chk {
            Ok(())
        } else {
            Err("A matrix multiplied by a tuple".into())
        }
    }

    /// Chap 3 - Multiplying a matrix by the identity matrix
    #[test]
    fn test_chap_3_7() -> Result<(), String> {
        let a = Matrix4::new([
            [0.0, 1.0, 2.0, 4.0],
            [1.0, 2.0, 4.0, 8.0],
            [2.0, 4.0, 8.0, 16.0],
            [4.0, 8.0, 16.0, 32.0],
        ]);
        let b = a * Matrix4::identity();
        let chk = b.approx_eq(a);
        if chk {
            Ok(())
        } else {
            Err("Multiplying a matrix by the identity matrix".into())
        }
    }

    /// Chap 3 - Multiplying the identity matrix by a tuple
    #[test]
    fn test_chap_3_8() -> Result<(), String> {
        let a = Tuple::new(1.0, 2.0, 3.0, 4.0);
        let chk = Tuple::new(1.0, 2.0, 3.0, 4.0).approx_eq(Matrix4::identity() * a);
        if chk {
            Ok(())
        } else {
            Err("Multiplying the identity matrix by a tuple".into())
        }
    }

    /// Chap 3 - Transposing a matrix
    #[test]
    fn test_chap_3_9() -> Result<(), String> {
        let a = Matrix4::new([
            [0.0, 9.0, 3.0, 0.0],
            [9.0, 8.0, 0.0, 8.0],
            [1.0, 8.0, 5.0, 3.0],
            [0.0, 0.0, 5.0, 8.0],
        ]);
        let expect = Matrix4::new([
            [0.0, 9.0, 1.0, 0.0],
            [9.0, 8.0, 8.0, 0.0],
            [3.0, 0.0, 5.0, 5.0],
            [0.0, 8.0, 3.0, 8.0],
        ]);
        let b = a.transpose();
        let c = a.transpose();
        let chk = expect.transpose().approx_eq(a);
        let chk = chk && b.approx_eq(c);
        let chk = chk && a.transpose().approx_eq(expect);
        if chk {
            Ok(())
        } else {
            Err("Transposing a matrix".into())
        }
    }

    /// Chap 3 - Transposing the identity matrix
    #[test]
    fn test_chap_3_10() -> Result<(), String> {
        let i = Matrix4::identity();
        let i_t = i.transpose();
        let chk = i_t.approx_eq(Matrix4::identity());
        if chk {
            Ok(())
        } else {
            Err("Transposing the identity matrix".into())
        }
    }

    /// Chap x - Calculating the determinant of a 2x2 matrix
    #[test]
    fn test_chap_3_11() -> Result<(), String> {
        let a = Matrix2::new([[1.0, 5.0], [-3.0, 2.0]]);
        let det_a = a.determinant();
        let chk = approx_eq(det_a, 17.0);
        if chk {
            Ok(())
        } else {
            Err("Calculating the determinant of a 2x2 matrix".into())
        }
    }

    /// Chap 3 - A submatrix of a 3x3 matrix is a 2x2 matrix
    #[test]
    fn test_chap_3_12() -> Result<(), String> {
        let a = Matrix3::new([[1.0, 5.0, 0.0], [-3.0, 2.0, 7.0], [0.0, 6.0, -3.0]]);
        let sm_a_0_2 = a.submatrix(0, 2);
        let expect = Matrix2::new([[-3.0, 2.0], [0.0, 6.0]]);

        let chk = sm_a_0_2.approx_eq(expect);
        let chk = chk && approx_eq(sm_a_0_2[(0, 0)], -3.0);
        let chk = chk && approx_eq(sm_a_0_2[(0, 1)], 2.0);
        let chk = chk && approx_eq(sm_a_0_2[(1, 0)], 0.0);
        let chk = chk && approx_eq(sm_a_0_2[(1, 1)], 6.0);
        if chk {
            Ok(())
        } else {
            Err("A submatrix of a 3x3 matrix is a 2x2 matrix".into())
        }
    }

    /// Chap 3 -A submatrix of a 4x4 matrix is a 3x3 matrix
    #[test]
    fn test_chap_3_13() -> Result<(), String> {
        let a = Matrix4::new([
            [-6.0, 1.0, 1.0, 6.0],
            [-8.0, 5.0, 8.0, 6.0],
            [-1.0, 0.0, 8.0, 2.0],
            [-7.0, 1.0, -1.0, 1.0],
        ]);

        let sm_a_2_1 = a.submatrix(2, 1);
        let chk = Matrix3::new([[-6.0, 1.0, 6.0], [-8.0, 8.0, 6.0], [-7.0, -1.0, 1.0]])
            .approx_eq(sm_a_2_1);

        if chk {
            Ok(())
        } else {
            Err("A submatrix of a 4x4 matrix is a 3x3 matrix".into())
        }
    }

    /// Chap 3 -Calculating a minor of a 3x3 matrix
    #[test]
    fn test_chap_3_14() -> Result<(), String> {
        let a = Matrix3::new([
            [3.0, 5.0, 0.0],   //
            [2.0, -1.0, -7.0], //
            [6.0, -1.0, 5.0],  //
        ]);
        let b = a.submatrix(1, 0);
        let chk = approx_eq(b.determinant(), 25.0);
        let chk = chk && approx_eq(a.minor(1, 0), 25.0);

        if chk {
            Ok(())
        } else {
            Err("Calculating a minor of a 3x3 matrix".into())
        }
    }

    /// Chap 3 - Calculating a cofactor of a 3x3 matrix
    #[test]
    fn test_chap_3_15() -> Result<(), String> {
        let a = Matrix3::new([
            [3.0, 5.0, 0.0],   //
            [2.0, -1.0, -7.0], //
            [6.0, -1.0, 5.0],  //
        ]);

        let minor_a_0_0 = a.minor(0, 0);
        let cofactor_a_0_0 = a.cofactor(0, 0);
        let minor_a_1_0 = a.minor(1, 0);
        let cofactor_1_0 = a.cofactor(1, 0);
        let chk = approx_eq(minor_a_0_0, cofactor_a_0_0);
        let chk = chk && approx_eq(minor_a_1_0, -cofactor_1_0);

        if chk {
            Ok(())
        } else {
            Err("Calculating a cofactor of a 3x3 matrix".into())
        }
    }

    /// Chap x - Calculating the determinant of a 3x3 matrix
    #[test]
    fn test_chap_3_16() -> Result<(), String> {
        let a = Matrix3::new([
            [1.0, 2.0, 6.0],   //
            [-5.0, 8.0, -4.0], //
            [2.0, 6.0, 4.0],   //
        ]);

        let cofactor_a_0_0 = a.cofactor(0, 0);
        let cofactor_a_0_1 = a.cofactor(0, 1);
        let cofactor_a_0_2 = a.cofactor(0, 2);
        let a_determinant = a.determinant();

        let chk = approx_eq(cofactor_a_0_0, 56.0)
            && approx_eq(cofactor_a_0_1, 12.0)
            && approx_eq(cofactor_a_0_2, -46.0)
            && approx_eq(a_determinant, -196.0);
        if chk {
            Ok(())
        } else {
            Err("Calculating the determinant of a 3x3 matrix".into())
        }
    }

    /// Chap 3 - Calculating the determinant of a 4x4 matrix
    #[test]
    fn test_chap_3_17() -> Result<(), String> {
        let a = Matrix4::new([
            [-2.0, -8.0, 3.0, 5.0],
            [-3.0, 1.0, 7.0, 3.0],
            [1.0, 2.0, -9.0, 6.0],
            [-6.0, 7.0, 7.0, -9.0],
        ]);

        let cofactor_a_0_0 = a.cofactor(0, 0);
        let cofactor_a_0_1 = a.cofactor(0, 1);
        let cofactor_a_0_2 = a.cofactor(0, 2);
        let cofactor_a_0_3 = a.cofactor(0, 3);
        let a_determinant = a.determinant();

        let chk = approx_eq(cofactor_a_0_0, 690.0)
            && approx_eq(cofactor_a_0_1, 447.0)
            && approx_eq(cofactor_a_0_2, 210.0)
            && approx_eq(cofactor_a_0_3, 51.0)
            && approx_eq(a_determinant, -4071.0);
        if chk {
            Ok(())
        } else {
            Err("Calculating the determinant of a 4x4 matrix".into())
        }
    }

    /// Chap 3 - Testing an invertible matrix for invertibility
    #[test]
    fn test_chap_3_18() -> Result<(), String> {
        let a = Matrix4::new([
            [6.0, 4.0, 4.0, 4.0],
            [5.0, 5.0, 7.0, 6.0],
            [4.0, -9.0, 3.0, -7.0],
            [9.0, 1.0, 7.0, -6.0],
        ]);

        let chk = approx_eq(a.determinant(), -2120.0) && a.is_invertible();
        if chk {
            Ok(())
        } else {
            Err("Testing an invertible matrix for invertibility".into())
        }
    }

    /// Chap x - Testing a noninvertible matrix for invertibility
    #[test]
    fn test_chap_3_19() -> Result<(), String> {
        let a = Matrix4::new([
            [4.0, 2.0, -2.0, -3.0],
            [9.0, 6.0, 2.0, 6.0],
            [0.0, -5.0, 1.0, -5.0],
            [0.0, 0.0, 0.0, 0.0],
        ]);

        let chk = approx_eq(a.determinant(), 0.0) && !a.is_invertible();
        if chk {
            Ok(())
        } else {
            Err("Testing a noninvertible matrix for invertibility".into())
        }
    }

    /// Chap 3 - Calculating the inverse of a matrix
    #[test]
    fn test_chap_3_20() -> Result<(), String> {
        let a = Matrix4::new([
            [-5.0, 2.0, 6.0, -8.0],
            [1.0, -5.0, 1.0, 8.0],
            [7.0, 7.0, -6.0, -7.0],
            [1.0, -3.0, 7.0, 4.0],
        ]);

        // NOTE: Need to have quite a few decimal places to satisfy approx_eq...
        let expect = Matrix4::new([
            [
                0.2180451127820,
                0.4511278195489,
                0.2406015037594,
                -0.0451127819549,
            ],
            [
                -0.8082706766917,
                -1.4567669172932,
                -0.4436090225564,
                0.5206766917293,
            ],
            [
                -0.0789473684211,
                -0.2236842105263,
                -0.0526315789474,
                0.1973684210526,
            ],
            [
                -0.5225563909774,
                -0.8139097744361,
                -0.3007518796992,
                0.3063909774436,
            ],
        ]);

        let b = a.inverse().unwrap_or(Matrix4::identity());
        logi!("test_chap_3_20", "\na:\n{}\n", a);
        logi!("test_chap_3_20", "\nb:\n{}\n", b);
        logi!("test_chap_3_20", "\nexpect:\n{}\n", expect);
        let determinant_a = a.determinant();
        let cofactor_a_2_3 = a.cofactor(2, 3);
        let cofactor_a_3_2 = a.cofactor(3, 2);

        let chk = true;
        let chk = chk && approx_eq(determinant_a, 532.0);
        let chk = chk && approx_eq(cofactor_a_2_3, -160.0);
        let chk = chk && approx_eq(cofactor_a_3_2, 105.0);
        let chk = chk && approx_eq(b[(2, 3)], 105.0 / 532.0);
        let chk = chk && approx_eq(b[(3, 2)], -160.0 / 532.0);
        let chk = chk && b.approx_eq(expect);
        if chk {
            Ok(())
        } else {
            Err("Calculating the inverse of a matrix".into())
        }
    }

    /// Chap 3 - Calculating the inverse of another matrix
    #[test]
    fn test_chap_3_21() -> Result<(), String> {
        let a = Matrix4::new([
            [8.0, -5.0, 9.0, 2.0],
            [7.0, 5.0, 6.0, 1.0],
            [-6.0, 0.0, 9.0, 6.0],
            [-3.0, 0.0, -9.0, -4.0],
        ]);

        let expect = Matrix4::new([
            [
                -0.1538461538462,
                -0.1538461538462,
                -0.2820512820513,
                -0.5384615384615,
            ],
            [
                -0.0769230769231,
                0.1230769230769,
                0.0256410256410,
                0.0307692307692,
            ],
            [
                0.3589743589744,
                0.3589743589744,
                0.4358974358974,
                0.9230769230769,
            ],
            [
                -0.6923076923077,
                -0.6923076923077,
                -0.7692307692308,
                -1.9230769230769,
            ],
        ]);

        let b = a.inverse().unwrap_or(Matrix4::identity());
        let c = a * b;

        let chk = true;
        let chk = chk && c.approx_eq(Matrix4::identity());
        let chk = chk && b.approx_eq(expect);

        logi!("test_chap_3_21", "\na:\n{}\n", a);
        logi!("test_chap_3_21", "\nb:\n{}\n", b);
        logi!("test_chap_3_21", "\nexpect:\n{}\n", c);

        if chk {
            Ok(())
        } else {
            Err("Calculating the inverse of another matrix".into())
        }
    }

    /// Chap 3 - Calculating the inverse of a third matrix
    #[test]
    fn test_chap_3_22() -> Result<(), String> {
        let a = Matrix4::new([
            [9.0, 3.0, 0.0, 9.0],
            [-5.0, -2.0, -6.0, -3.0],
            [-4.0, 9.0, 6.0, 4.0],
            [-7.0, 6.0, 6.0, 2.0],
        ]);

        let expect = Matrix4::new([
            [
                -0.0407407407407,
                -0.0777777777778,
                0.1444444444444,
                -0.2222222222222,
            ],
            [
                -0.0777777777778,
                0.0333333333333,
                0.3666666666667,
                -0.3333333333333,
            ],
            [
                -0.0290123456790,
                -0.1462962962963,
                -0.1092592592593,
                0.1296296296296,
            ],
            [
                0.1777777777778,
                0.0666666666667,
                -0.2666666666667,
                0.3333333333333,
            ],
        ]);

        let b = a.inverse().unwrap_or(Matrix4::identity());
        let c = a * b;

        let chk = true;
        let chk = chk && c.approx_eq(Matrix4::identity());
        let chk = chk && b.approx_eq(expect);

        logi!("test_chap_3_22", "\na:\n{}\n", a);
        logi!("test_chap_3_22", "\nb:\n{}\n", b);
        logi!("test_chap_3_22", "\nexpect:\n{}\n", c);

        if chk {
            Ok(())
        } else {
            Err("Calculating the inverse of a third matrix".into())
        }
    }

    /// Chap 3 - Multiplying a product by its inverse
    #[test]
    fn test_chap_3_23() -> Result<(), String> {
        let a = Matrix4::new([
            [3.0, -9.0, 7.0, 3.0],
            [3.0, -8.0, 2.0, -9.0],
            [-4.0, 4.0, 4.0, 1.0],
            [-6.0, 5.0, -1.0, 1.0],
        ]);
        let b = Matrix4::new([
            [8.0, 2.0, 2.0, 2.0],
            [3.0, -1.0, 7.0, 0.0],
            [7.0, 0.0, 5.0, 4.0],
            [6.0, -2.0, 0.0, 5.0],
        ]);
        let c = a * b;
        let a2 = c * b.inverse().unwrap_or(Matrix4::identity());
        let chk = a.approx_eq(a2);
        if chk {
            Ok(())
        } else {
            Err("Multiplying a product by its inverse".into())
        }
    }

    /// Chap 4 - Multiplying by a translation matrix
    #[test]
    fn test_chap_4_1() -> Result<(), String> {
        let transform = Matrix4::translation(5.0, -3.0, 2.0);
        let p = Tuple::point(-3.0, 4.0, 5.0);
        let p_translated = transform * p;
        let p_expected = Tuple::point(2.0, 1.0, 7.0);
        let chk = p_translated.approx_eq(p_expected);
        if chk {
            Ok(())
        } else {
            Err("Multiplying by a translation matrix".into())
        }
    }

    /// Chap 4 - Multiplying by the inverse of a translation matrix
    #[test]
    fn test_chap_4_2() -> Result<(), String> {
        let transform = Matrix4::translation(5.0, -3.0, 2.0);
        let inv = transform.inverse().unwrap_or(Matrix4::identity());
        let p = Tuple::point(-3.0, 4.0, 5.0);
        let result = inv * p;
        let chk = result.approx_eq(Tuple::point(-8.0, 7.0, 3.0));

        if chk {
            Ok(())
        } else {
            Err("Multiplying by the inverse of a translation matrix".into())
        }
    }

    /// Chap 4 - Translation does not affect vectors
    #[test]
    fn test_chap_4_3() -> Result<(), String> {
        let transform = Matrix4::translation(5.0, -3.0, 2.0);
        let v = Tuple::vector(-3.0, 4.0, 5.0);
        let result = transform * v;
        let chk = result.approx_eq(v);
        if chk {
            Ok(())
        } else {
            Err("Translation does not affect vectors".into())
        }
    }

    /// Chap 4 - A scaling matrix applied to a point
    #[test]
    fn test_chap_4_4() -> Result<(), String> {
        let transform = Matrix4::scaling(2.0, 3.0, 4.0);
        let p = Tuple::point(-4.0, 6.0, 8.0);
        let result = transform * p;
        let chk = result.approx_eq(Tuple::point(-8.0, 18.0, 32.0));
        if chk {
            Ok(())
        } else {
            Err("A scaling matrix applied to a point".into())
        }
    }

    /// Chap x - A scaling matrix applied to a vector
    #[test]
    fn test_chap_4_5() -> Result<(), String> {
        let transform = Matrix4::scaling(2.0, 3.0, 4.0);
        let v = Tuple::vector(-4.0, 6.0, 8.0);
        let result = transform * v;
        let chk = result.approx_eq(Tuple::vector(-8.0, 18.0, 32.0));
        if chk {
            Ok(())
        } else {
            Err("A scaling matrix applied to a vector".into())
        }
    }

    /// Chap 4 - Multiplying by the inverse of a scaling matrix
    #[test]
    fn test_chap_4_6() -> Result<(), String> {
        let transform = Matrix4::scaling(2.0, 3.0, 4.0);
        let inv = transform.inverse().unwrap_or(Matrix4::identity());
        let v = Tuple::vector(-4.0, 6.0, 8.0);
        let result = inv * v;
        let chk = result.approx_eq(Tuple::vector(-2.0, 2.0, 2.0));
        if chk {
            Ok(())
        } else {
            Err("Multiplying by the inverse of a scaling matrix".into())
        }
    }

    /// Chap 4 - Reflection is scaling by a negative value
    #[test]
    fn test_chap_4_7() -> Result<(), String> {
        let transform = Matrix4::scaling(-1.0, 1.0, 1.0);
        let p = Tuple::point(2.0, 3.0, 4.0);
        let result = transform * p;
        let chk = result.approx_eq(Tuple::point(-2.0, 3.0, 4.0));
        if chk {
            Ok(())
        } else {
            Err("Reflection is scaling by a negative value".into())
        }
    }

    /// Chap 4 - Rotating a point around the x axis
    #[test]
    fn test_chap_4_8() -> Result<(), String> {
        let p = Tuple::point(0.0, 1.0, 0.0);
        let half_quarter = Matrix4::rotation_x(std::f64::consts::PI / 4.0);
        let full_quarter = Matrix4::rotation_x(std::f64::consts::PI / 2.0);
        let phq = half_quarter * p;
        let pfq = full_quarter * p;

        let sqrt2 = std::f64::consts::SQRT_2;
        let chk = true;
        let chk = chk && phq.approx_eq(Tuple::point(0.0, sqrt2 / 2.0, sqrt2 / 2.0));
        let chk = chk && pfq.approx_eq(Tuple::point(0.0, 0.0, 1.0));

        if chk {
            Ok(())
        } else {
            Err("Rotating a point around the x axis".into())
        }
    }

    /// Chap 4 - The inverse of an x-rotation rotates in the opposite direction
    #[test]
    fn test_chap_4_9() -> Result<(), String> {
        let p = Tuple::point(0.0, 1.0, 0.0);
        let half_quarter = Matrix4::rotation_x(std::f64::consts::PI / 4.0);
        let inv = half_quarter.inverse().unwrap_or(Matrix4::identity());
        let phq = inv * p;

        let sqrt2 = std::f64::consts::SQRT_2;
        let chk = true;
        let chk = chk && phq.approx_eq(Tuple::point(0.0, sqrt2 / 2.0, -sqrt2 / 2.0));
        if chk {
            Ok(())
        } else {
            Err("The inverse of an x-rotation rotates in the opposite direction".into())
        }
    }

    /// Chap 4 - Rotating a point around the y axis
    #[test]
    fn test_chap_4_10() -> Result<(), String> {
        let p = Tuple::point(0.0, 0.0, 1.0);
        let half_quarter = Matrix4::rotation_y(std::f64::consts::PI / 4.0);
        let full_quarter = Matrix4::rotation_y(std::f64::consts::PI / 2.0);
        let phq = half_quarter * p;
        let pfq = full_quarter * p;

        let sqrt2 = std::f64::consts::SQRT_2;
        let chk = true;
        let chk = chk && phq.approx_eq(Tuple::point(sqrt2 / 2.0, 0.0, sqrt2 / 2.0));
        let chk = chk && pfq.approx_eq(Tuple::point(1.0, 0.0, 0.0));
        if chk {
            Ok(())
        } else {
            Err("Rotating a point around the z axis".into())
        }
    }

    /// Chap 4 - Rotating a point around the z axis
    #[test]
    fn test_chap_4_11() -> Result<(), String> {
        let p = Tuple::point(0.0, 1.0, 0.0);
        let half_quarter = Matrix4::rotation_z(std::f64::consts::PI / 4.0);
        let full_quarter = Matrix4::rotation_z(std::f64::consts::PI / 2.0);
        let phq = half_quarter * p;
        let pfq = full_quarter * p;

        let sqrt2 = std::f64::consts::SQRT_2;
        let chk = true;
        let chk = chk && phq.approx_eq(Tuple::point(-sqrt2 / 2.0, sqrt2 / 2.0, 0.0));
        let chk = chk && pfq.approx_eq(Tuple::point(-1.0, 0.0, 0.0));
        if chk {
            Ok(())
        } else {
            Err("Rotating a point around the z axis".into())
        }
    }

    /// Chap x - A shearing transformation moves x in proportion to y
    #[test]
    fn test_chap_4_12() -> Result<(), String> {
        let transform = Matrix4::shearing(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let p = Tuple::point(2.0, 3.0, 4.0);
        let result = transform * p;
        let chk = result.approx_eq(Tuple::point(5.0, 3.0, 4.0));
        if chk {
            Ok(())
        } else {
            loge!("test_chap_4_12", "result:{}", result);
            loge!("test_chap_4_12", "result:{}", transform);
            Err("A shearing transformation moves x in proportion to y".into())
        }
    }

    /// Chap x - Individual transformations are applied in sequence
    #[test]
    fn test_chap_4_13() -> Result<(), String> {
        let p = Tuple::point(1.0, 0.0, 1.0);
        let a = Matrix4::rotation_x(std::f64::consts::PI / 2.0);
        let b = Matrix4::scaling(5.0, 5.0, 5.0);
        let c = Matrix4::translation(10.0, 5.0, 7.0);

        // apply rotation first
        let p2 = a * p;
        let chk = p2.approx_eq(Tuple::point(1.0, -1.0, 0.0));

        // then apply scaling
        let p3 = b * p2;
        let chk = chk && p3.approx_eq(Tuple::point(5.0, -5.0, 0.0));

        // then apply translation
        let p4 = c * p3;
        let chk = chk && p4.approx_eq(Tuple::point(15.0, 0.0, 7.0));

        if chk {
            Ok(())
        } else {
            Err("Individual transformations are applied in sequence".into())
        }
    }

    /// Chap 4 - Chained transformations must be applied in reverse order
    #[test]
    fn test_chap_4_14() -> Result<(), String> {
        let p = Tuple::point(1.0, 0.0, 1.0);
        let a = Matrix4::rotation_x(std::f64::consts::PI / 2.0);
        let b = Matrix4::scaling(5.0, 5.0, 5.0);
        let c = Matrix4::translation(10.0, 5.0, 7.0);
        let t = c * b * a;
        let expect = t * p;

        let chk = expect.approx_eq(Tuple::point(15.0, 0.0, 7.0));
        if chk {
            Ok(())
        } else {
            Err("Chained transformations must be applied in reverse order".into())
        }
    }

    /// Chap 4 - Putting it together
    #[test]
    fn test_chap_4_15() -> std::io::Result<()> {
        let mut canvas = Canvas::new(500, 500);
        let p = Tuple::point(1.0, 0.0, 0.0);

        // canvas.width() as f64 - 10.0,
        // canvas.height() as f64 - 10.0,
        let b = Matrix4::scaling(
            canvas.width() as f64 / 2.5,
            canvas.width() as f64 / 2.5,
            4.0,
        );
        let c = Matrix4::translation(
            canvas.width() as f64 / 2.0,
            canvas.height() as f64 / 2.0,
            5.0,
        );

        let alpha_increment = std::f64::consts::PI / 6.0;
        let alpha_max = std::f64::consts::PI * 2.0;

        let mut alpha = 0.0;

        while alpha < alpha_max {
            let a = Matrix4::rotation_z(alpha);

            let p_clock = c * b * a * p;

            // logi!(
            //     "test_chap_4_15",
            //     "\nalpha: {}. clock point:\n{}\n",
            //     alpha * 180.0 / std::f64::consts::PI,
            //     p_clock
            // );

            alpha += alpha_increment;

            let dot_radius = 4;
            for dy in -dot_radius..=dot_radius {
                for dx in -dot_radius..=dot_radius {
                    let px = (p_clock.x.floor() as i32 + dx) as usize;
                    let py = (p_clock.y.floor() as i32 + dy) as usize;
                    canvas[(px, py)] = Tuple::color(1.0, 1.0, 0.0);
                }
            }
        }

        canvas.write_ppm("image_4_15.ppm")?;

        Ok(())
    }
}
