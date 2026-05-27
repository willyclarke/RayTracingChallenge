//! Intersections
//! Defintion of the intersections that occurs between rays and shapes.

use std::ops::Index;

#[derive(Debug, Copy, Clone)]
pub struct Intersection {
    pub t: f64,
    pub object_id: usize,
}

impl Intersection {
    pub fn new(t: f64, object_id: usize) -> Self {
        Self { t, object_id }
    }
}

#[derive(Debug, Clone)]
pub struct Intersections {
    data: Vec<Intersection>,
}

impl Intersections {
    pub fn hit(&self) -> Option<Intersection> {
        self.data.iter().find(|i| i.t >= 0.0).copied()
    }

    pub fn new() -> Self {
        Self { data: vec![] }
    }

    /// Push the Intersection by first finding the partition_point of where to insert
    pub fn push(&mut self, i: Intersection) {
        let pos = self.data.partition_point(|x| x.t < i.t);
        self.data.insert(pos, i);
    }

    pub fn count(&self) -> usize {
        self.data.len()
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = Intersection> + '_ {
        self.data.iter().copied()
    }
}

impl Default for Intersections {
    fn default() -> Self {
        Self::new()
    }
}

impl Index<usize> for Intersections {
    type Output = Intersection;
    fn index(&self, i: usize) -> &Self::Output {
        &self.data[i]
    }
}
