
pub const EPSILON: f64 = 1e-9;


/// Return true when delta between to floating point numbers are less than epsilon.
pub fn approx_eq(a: f64, b: f64) -> bool {
    
    (a - b).abs() < EPSILON
}


