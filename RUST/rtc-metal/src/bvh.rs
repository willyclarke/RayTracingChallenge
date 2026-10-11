//! Median-split bounding volume hierarchy over the flattened shapes (book
//! chapter 15's bounding boxes, built on the CPU for the GPU to walk).

use crate::types::{GpuBvhNode, GpuShape};

pub const MAX_DEPTH: usize = 32;
pub const LEAF_SIZE: usize = 4;

/// A flattened leaf. `bounds` is its f32 world box, None when unbounded.
pub struct Leaf {
    pub shape: GpuShape,
    pub id: usize,
    pub bounds: Option<([f32; 3], [f32; 3])>,
}

pub struct Bvh {
    /// Unbounded shapes first, then the bounded ones ordered so that every
    /// BVH leaf covers a contiguous range.
    pub shapes: Vec<GpuShape>,
    /// rtc-rust shape id per entry of `shapes`.
    pub ids: Vec<usize>,
    /// Depth-first; empty when nothing is bounded.
    pub nodes: Vec<GpuBvhNode>,
    pub unbounded_count: u32,
}

type Bounds = ([f32; 3], [f32; 3]);

fn union(leaves: &[(Leaf, Bounds)]) -> Bounds {
    let mut lo = [f32::INFINITY; 3];
    let mut hi = [f32::NEG_INFINITY; 3];
    for (_, (min, max)) in leaves {
        for a in 0..3 {
            lo[a] = lo[a].min(min[a]);
            hi[a] = hi[a].max(max[a]);
        }
    }
    (lo, hi)
}

fn centroid(b: &Bounds, axis: usize) -> f32 {
    (b.0[axis] + b.1[axis]) * 0.5
}

fn node(
    range: &mut [(Leaf, Bounds)],
    start: usize,
    depth: usize,
    nodes: &mut Vec<GpuBvhNode>,
) -> Result<(), String> {
    if depth > MAX_DEPTH {
        return Err(format!("BVH deeper than {MAX_DEPTH} levels"));
    }
    let (min, max) = union(range);
    let index = nodes.len();
    nodes.push(GpuBvhNode {
        min: [min[0], min[1], min[2], 0.0],
        max: [max[0], max[1], max[2], 0.0],
        start: 0,
        count: 0,
        right: 0,
        _pad: 0,
    });
    if range.len() <= LEAF_SIZE {
        nodes[index].start = start as u32;
        nodes[index].count = range.len() as u32;
        return Ok(());
    }
    // Split the longest axis of the centroids' box, which separates clusters
    // better than the shapes' own box when they vary in size.
    let mut clo = [f32::INFINITY; 3];
    let mut chi = [f32::NEG_INFINITY; 3];
    for (_, b) in range.iter() {
        for a in 0..3 {
            clo[a] = clo[a].min(centroid(b, a));
            chi[a] = chi[a].max(centroid(b, a));
        }
    }
    let axis = (0..3)
        .max_by(|&a, &b| (chi[a] - clo[a]).total_cmp(&(chi[b] - clo[b])))
        .unwrap();
    range.sort_by(|x, y| centroid(&x.1, axis).total_cmp(&centroid(&y.1, axis)));
    // Splitting by count, not position, terminates even when every centroid
    // coincides.
    let mid = range.len() / 2;
    let (left, right) = range.split_at_mut(mid);
    node(left, start, depth + 1, nodes)?;
    nodes[index].right = nodes.len() as u32;
    node(right, start + mid, depth + 1, nodes)
}

/// Err if the tree would be deeper than `MAX_DEPTH`, the traversal stack size.
pub fn build(leaves: Vec<Leaf>) -> Result<Bvh, String> {
    let mut shapes = Vec::new();
    let mut ids = Vec::new();
    let mut bounded = Vec::new();
    for leaf in leaves {
        match leaf.bounds {
            None => {
                shapes.push(leaf.shape);
                ids.push(leaf.id);
            }
            Some(b) => bounded.push((leaf, b)),
        }
    }
    let unbounded_count = shapes.len();
    let mut nodes = Vec::new();
    if !bounded.is_empty() {
        node(&mut bounded, unbounded_count, 0, &mut nodes)?;
    }
    for (leaf, _) in bounded {
        shapes.push(leaf.shape);
        ids.push(leaf.id);
    }
    Ok(Bvh {
        shapes,
        ids,
        nodes,
        unbounded_count: unbounded_count as u32,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn leaf(i: usize, bounds: Option<Bounds>) -> Leaf {
        Leaf {
            shape: unsafe { std::mem::zeroed() },
            id: i,
            bounds,
        }
    }

    fn unit_box(x: usize, y: usize, z: usize) -> Option<Bounds> {
        let p = [x as f32, y as f32, z as f32];
        Some((p, [p[0] + 1.0, p[1] + 1.0, p[2] + 1.0]))
    }

    fn sample() -> (Vec<Bounds>, Bvh) {
        let mut leaves = vec![leaf(0, None), leaf(1, None)];
        let mut boxes = vec![];
        for i in 2..39 {
            let b = unit_box((i * 7) % 11, (i * 3) % 5, i % 13);
            boxes.push(b.unwrap());
            leaves.push(leaf(i, b));
        }
        (boxes, build(leaves).unwrap())
    }

    fn leaf_nodes(bvh: &Bvh) -> Vec<&GpuBvhNode> {
        bvh.nodes.iter().filter(|n| n.count > 0).collect()
    }

    #[test]
    fn test_chap_15_bvh_covers_every_leaf_once() {
        let (_, bvh) = sample();
        assert_eq!(bvh.unbounded_count, 2);
        assert_eq!(bvh.ids[..2], [0, 1]);
        let mut seen = [0; 39];
        let mut next = 2;
        let mut leaves = leaf_nodes(&bvh);
        leaves.sort_by_key(|n| n.start);
        for n in leaves {
            // Disjoint and gap-free: each starts where the last ended.
            assert_eq!(n.start as usize, next);
            for k in n.start..n.start + n.count {
                seen[bvh.ids[k as usize]] += 1;
            }
            next += n.count as usize;
        }
        assert_eq!(next, bvh.shapes.len());
        assert_eq!(seen[..2], [0, 0]);
        assert!(seen[2..].iter().all(|&c| c == 1));
    }

    #[test]
    fn test_chap_15_bvh_boxes_nest() {
        let (boxes, bvh) = sample();
        let contains = |n: &GpuBvhNode, lo: [f32; 3], hi: [f32; 3]| {
            (0..3).all(|a| n.min[a] <= lo[a] && n.max[a] >= hi[a])
        };
        for (i, n) in bvh.nodes.iter().enumerate() {
            assert_eq!((n.min[3], n.max[3]), (0.0, 0.0));
            if n.count > 0 {
                assert!(n.count as usize <= LEAF_SIZE);
                for k in n.start..n.start + n.count {
                    let b = boxes[bvh.ids[k as usize] - 2];
                    assert!(contains(n, b.0, b.1));
                }
            } else {
                for c in [&bvh.nodes[i + 1], &bvh.nodes[n.right as usize]] {
                    assert!(contains(
                        n,
                        [c.min[0], c.min[1], c.min[2]],
                        [c.max[0], c.max[1], c.max[2]]
                    ));
                }
            }
        }
    }

    #[test]
    fn test_chap_15_bvh_identical_boxes_terminate() {
        let leaves = (0..100).map(|i| leaf(i, unit_box(0, 0, 0))).collect();
        let bvh = build(leaves).unwrap();
        assert!(leaf_nodes(&bvh).iter().all(|n| n.count as usize <= 4));
        let total: u32 = leaf_nodes(&bvh).iter().map(|n| n.count).sum();
        assert_eq!(total, 100);
        // Depth-first with the left child at index + 1: track depth by walking.
        fn depth(bvh: &Bvh, i: usize) -> usize {
            let n = &bvh.nodes[i];
            if n.count > 0 {
                0
            } else {
                1 + depth(bvh, i + 1).max(depth(bvh, n.right as usize))
            }
        }
        assert!(depth(&bvh, 0) <= MAX_DEPTH);
    }

    #[test]
    fn test_chap_15_bvh_empty_and_unbounded_only() {
        let bvh = build(vec![]).unwrap();
        assert!(bvh.nodes.is_empty() && bvh.shapes.is_empty());
        let bvh = build((0..3).map(|i| leaf(i, None)).collect()).unwrap();
        assert!(bvh.nodes.is_empty());
        assert_eq!(bvh.unbounded_count, 3);
        assert_eq!(bvh.shapes.len(), 3);
    }
}
