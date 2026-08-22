//! Read data from OBJ files.

use crate::{shapes::triangleuv::TriangleUV, tuple::Tuple};
use std::collections::HashMap;
use std::fmt;
use std::fs;
use std::io;
use std::path::Path;

pub struct Parser {
    pub vertices: Vec<Tuple>,
    pub normals: Vec<Tuple>,
    pub default_group: Vec<TriangleUV>,
    pub named_groups: HashMap<String, Vec<TriangleUV>>,
    pub ignored: Vec<(usize, String)>, // (1-based line number, the line text)
    current_group: Option<String>,
}

impl Default for Parser {
    fn default() -> Self {
        Self::new()
    }
}

impl fmt::Display for Parser {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let _ = writeln!(f, " ");
        for idx in 1..self.vertices.len() {
            let _ = writeln!(f, "Vertice[{}] = {}", idx, self.vertices[idx]);
        }

        for (name, triangles) in &self.named_groups {
            let _ = writeln!(f, "group \"{name}\": {} triangle(s)", triangles.len());
            for t in triangles {
                let _ = writeln!(f, "  {:?} {:?} {:?}", t.p1, t.p2, t.p3);
            }
        }

        write!(f, " ")
    }
}

impl Parser {
    /// The Vec the current group's triangles go into.
    fn current_group_faces(&mut self) -> &mut Vec<TriangleUV> {
        match self.current_group.clone() {
            Some(name) => self.named_groups.entry(name).or_default(),
            None => &mut self.default_group,
        }
    }

    pub fn new() -> Self {
        Self {
            vertices: vec![Tuple::point(0.0, 0.0, 0.0)], // dummy at index 0
            normals: vec![Tuple::vector(0.0, 0.0, 0.0)], // dummy at index 0
            default_group: Vec::new(),
            ignored: Vec::new(),
            named_groups: HashMap::new(),
            current_group: None,
        }
    }

    pub fn parse(text: &str) -> Self {
        let mut p = Self::new();
        for (i, line) in text.lines().enumerate() {
            let line_no = i + 1; // 1-based, matches an editor
            let mut tok = line.split_whitespace();
            match tok.next() {
                None | Some("#") => {} // blank line or comment → valid, skip silently
                Some("v") => {
                    let coords: Option<Vec<f64>> = tok.map(|s| s.parse::<f64>().ok()).collect();
                    match coords {
                        Some(c) if c.len() >= 3 => {
                            p.vertices.push(Tuple::point(c[0], c[1], c[2]));
                        }
                        _ => p.ignored.push((line_no, line.to_string())), // malformed vertex
                    }
                }
                Some("vn") => {
                    let coords: Option<Vec<f64>> = tok.map(|s| s.parse::<f64>().ok()).collect();
                    match coords {
                        Some(c) if c.len() >= 3 => {
                            p.normals.push(Tuple::vector(c[0], c[1], c[2]));
                        }
                        _ => p.ignored.push((line_no, line.to_string())), // malformed normal
                    }
                }
                Some("f") => {
                    // Handle face elements. Each token is `v`, `v/vt`, `v//vn`, or
                    // `v/vt/vn`; the texture index (vt) is ignored.
                    // Indices parse SIGNED so negative (relative) indices are allowed.
                    let raw: Option<Vec<(isize, Option<isize>)>> = tok
                        .map(|s| {
                            let mut parts = s.split('/');
                            let v = parts.next()?.parse::<isize>().ok()?;
                            let _vt = parts.next(); // texture index, unused
                            let vn = match parts.next() {
                                None | Some("") => None,
                                Some(n) => Some(n.parse::<isize>().ok()?),
                            };
                            Some((v, vn))
                        })
                        .collect();
                    match raw {
                        Some(raw) if raw.len() >= 3 => {
                            // resolve EACH index independently to an absolute vec index
                            let resolve = |r: isize, len: usize| {
                                if r < 0 {
                                    (len as isize + r) as usize // relative: from the end
                                } else {
                                    r as usize // absolute (1-based; dummy at 0 handles it)
                                }
                            };
                            let resolved: Vec<(usize, Option<usize>)> = raw
                                .iter()
                                .map(|&(v, vn)| {
                                    (
                                        resolve(v, p.vertices.len()),
                                        vn.map(|n| resolve(n, p.normals.len())),
                                    )
                                })
                                .collect();

                            // validate the RESOLVED indices point at real vertices assumption: the
                            // incoming data always describes convex polygons (interior angles are
                            // all less than or equal to 180°) When this is the case you can break
                            // them into triangles using a fan triangulation.
                            // This means that the triangle starting index is always found in
                            // resolved[0].
                            let verts_ok = resolved
                                .iter()
                                .all(|&(v, _)| v >= 1 && v < p.vertices.len());
                            let norms_ok = resolved
                                .iter()
                                .all(|&(_, vn)| vn.is_none_or(|n| n >= 1 && n < p.normals.len()));
                            // smooth shading only when EVERY vertex brings a normal
                            let smooth = resolved.iter().all(|&(_, vn)| vn.is_some());
                            if verts_ok && norms_ok {
                                // 1. build (this READS p.vertices / p.normals)
                                let mut tris = Vec::new();

                                for k in 1..resolved.len() - 1 {
                                    let p1 = p.vertices[resolved[0].0];
                                    let p2 = p.vertices[resolved[k].0];
                                    let p3 = p.vertices[resolved[k + 1].0];
                                    if smooth {
                                        let n1 = p.normals[resolved[0].1.unwrap()];
                                        let n2 = p.normals[resolved[k].1.unwrap()];
                                        let n3 = p.normals[resolved[k + 1].1.unwrap()];
                                        tris.push(TriangleUV::new(p1, p2, p3, n1, n2, n3));
                                    } else {
                                        tris.push(TriangleUV::flat(p1, p2, p3));
                                    }
                                }

                                // 2. route (this WRITES to a group)
                                p.current_group_faces().extend(tris);
                            } else {
                                p.ignored.push((line_no, line.to_string()));
                            }
                        }
                        _ => p.ignored.push((line_no, line.to_string())),
                    }
                }
                // Handle named groups
                Some("g") => {
                    // `g Name` just switches the current group; following faces route into it
                    match tok.next() {
                        Some(name) => p.current_group = Some(name.to_string()),
                        None => p.ignored.push((line_no, line.to_string())), // `g` with no name
                    }
                }
                _ => p.ignored.push((line_no, line.to_string())),
            }
        }
        p
    }

    /// Parse an OBJ file from disk.
    pub fn parse_obj_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let text = fs::read_to_string(path)?;
        Ok(Self::parse(&text))
    }

    pub fn report(&self) -> String {
        if self.ignored.is_empty() {
            return "OBJ parsed cleanly (no ignored lines).".to_string();
        }
        let mut s = format!("{} ignored line(s):\n", self.ignored.len());
        for (line_no, text) in &self.ignored {
            s.push_str(&format!("  line {line_no}: {text}\n"));
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use crate::loge;

    use super::*;

    /// Chap 16 - Ignoring unrecognized lines
    #[test]
    fn test_chap_16_1() -> Result<(), String> {
        let gibberish = "\"There was a young lady named Bright\n who traveled much faster than light.\nShe set out one day\n in a relative way,\nand came back the previous night.\"";

        let p = Parser::parse(gibberish);

        let chk = p.ignored.len() == 5;
        if chk {
            Ok(())
        } else {
            loge!("test_chap_16_1", "str:{}", gibberish);
            Err("Ignoring unrecognized lines".into())
        }
    }

    /// Chap 16 - Vertex records
    #[test]
    fn test_chap_16_2() -> std::io::Result<()> {
        let path = std::env::temp_dir().join("test_chap_16_2.obj");
        let contents = "v -1 1 0\nv -1.0000 0.5000 0.0000\nv 1 0 0\nv 1 1 0";

        std::fs::write(&path, contents)?; // ? → I/O errors become the Err
        let p = Parser::parse_obj_file(&path)?;
        std::fs::remove_file(&path)?; // cleaned up BEFORE the asserts

        println!("{} {}", path.display(), p.report());

        assert!(p.vertices[1].approx_eq(Tuple::point(-1.0, 1.0, 0.0)));
        assert!(p.vertices[2].approx_eq(Tuple::point(-1.0, 0.5, 0.0)));
        assert!(p.vertices[3].approx_eq(Tuple::point(1.0, 0.0, 0.0)));
        assert!(p.vertices[4].approx_eq(Tuple::point(1.0, 1.0, 0.0)));
        Ok(())
    }

    /// Chap 16 - Parsing triangle faces
    #[test]
    fn test_chap_16_3() -> std::io::Result<()> {
        let path = std::env::temp_dir().join("test_chap_16_3.obj");
        let contents = "v -1 1 0\nv -1 0 0\nv 1 0 0\nv 1 1 0\n \nf 1 2 3\nf 1 3 4\nv 2 3 4\nv 5 6 7\nv 8 9 0\nf -1 -2 -3";

        std::fs::write(&path, contents)?; // ? → I/O errors become the Err
        let p = Parser::parse_obj_file(&path)?;
        std::fs::remove_file(&path)?; // cleaned up BEFORE the asserts

        println!("{} {}", path.display(), p.report());

        assert!(p.vertices[1].approx_eq(Tuple::point(-1.0, 1.0, 0.0)));
        assert!(p.vertices[2].approx_eq(Tuple::point(-1.0, 0.0, 0.0)));
        assert!(p.vertices[3].approx_eq(Tuple::point(1.0, 0.0, 0.0)));
        assert!(p.vertices[4].approx_eq(Tuple::point(1.0, 1.0, 0.0)));

        let g = &p.default_group;
        assert!(!g.is_empty());
        assert!(g.len() == 3);

        let t1 = &g[0];
        let t2 = &g[1];
        assert!(t1.p1.approx_eq(p.vertices[1]));
        assert!(t1.p2.approx_eq(p.vertices[2]));
        assert!(t1.p3.approx_eq(p.vertices[3]));
        assert!(t2.p1.approx_eq(p.vertices[1]));
        assert!(t2.p2.approx_eq(p.vertices[3]));
        assert!(t2.p3.approx_eq(p.vertices[4]));

        Ok(())
    }

    /// Chap 16 - Triangulating polygons
    #[test]
    fn test_chap_16_4() -> std::io::Result<()> {
        let path = std::env::temp_dir().join("test_chap_16_4.obj");
        let contents = "v -1 1 0\nv -1 0 0\nv 1 0 0\nv 1 1 0\n \nv 0 2 0\nf 1 2 3 4 5";

        std::fs::write(&path, contents)?; // ? → I/O errors become the Err
        let p = Parser::parse_obj_file(&path)?;
        std::fs::remove_file(&path)?; // cleaned up BEFORE the asserts

        println!("{} {}", path.display(), p.report());

        let g = &p.default_group;
        assert!(!g.is_empty());
        assert!(g.len() == 3);

        let t1 = &g[0];
        let t2 = &g[1];
        let t3 = &g[2];
        assert!(p.vertices[1].approx_eq(Tuple::point(-1.0, 1.0, 0.0)));
        assert!(t1.p1.approx_eq(p.vertices[1]));
        assert!(t1.p2.approx_eq(p.vertices[2]));
        assert!(t1.p3.approx_eq(p.vertices[3]));
        assert!(t2.p1.approx_eq(p.vertices[1]));
        assert!(t2.p2.approx_eq(p.vertices[3]));
        assert!(t2.p3.approx_eq(p.vertices[4]));
        assert!(t3.p1.approx_eq(p.vertices[1]));
        assert!(t3.p2.approx_eq(p.vertices[4]));
        assert!(t3.p3.approx_eq(p.vertices[5]));

        Ok(())
    }

    /// Chap 16 - Triangles in groups
    #[test]
    fn test_chap_16_5() -> std::io::Result<()> {
        let path = std::env::temp_dir().join("test_chap_16_5.obj");
        let contents = "v -1 1 0\nv -1 0 0\nv 1 0 0\nv 1 1 0\n \ng FirstGroup\nf 1 2 3\ng SecondGroup\nf 1 3 4";

        std::fs::write(&path, contents)?; // ? → I/O errors become the Err
        let p = Parser::parse_obj_file(&path)?;
        std::fs::remove_file(&path)?; // cleaned up BEFORE the asserts

        println!("{} {}", path.display(), p.report());

        // ...parse...
        let t1 = &p.named_groups["FirstGroup"][0];
        let t2 = &p.named_groups["SecondGroup"][0];
        assert!(t1.p1.approx_eq(p.vertices[1]));
        assert!(t1.p2.approx_eq(p.vertices[2]));
        assert!(t1.p3.approx_eq(p.vertices[3]));
        assert!(t2.p1.approx_eq(p.vertices[1]));
        assert!(t2.p2.approx_eq(p.vertices[3]));
        assert!(t2.p3.approx_eq(p.vertices[4]));

        Ok(())
    }

    /// Chap 16 - Vertex normal records
    #[test]
    fn test_chap_16_6() -> std::io::Result<()> {
        let path = std::env::temp_dir().join("test_chap_16_6.obj");
        let contents = "vn 0 0 1\nvn 0.707 0 -0.707\nvn 1 2 3";

        std::fs::write(&path, contents)?; // ? → I/O errors become the Err
        let p = Parser::parse_obj_file(&path)?;
        std::fs::remove_file(&path)?; // cleaned up BEFORE the asserts

        println!("{} {}", path.display(), p.report());

        assert!(p.normals[1].approx_eq(Tuple::vector(0.0, 0.0, 1.0)));
        assert!(p.normals[2].approx_eq(Tuple::vector(0.707, 0.0, -0.707)));
        assert!(p.normals[3].approx_eq(Tuple::vector(1.0, 2.0, 3.0)));
        Ok(())
    }

    /// Chap 16 - Faces with normals
    #[test]
    fn test_chap_16_7() -> std::io::Result<()> {
        let path = std::env::temp_dir().join("test_chap_16_7.obj");
        let contents = "v 0 1 0\nv -1 0 0\nv 1 0 0\n \nvn -1 0 0\nvn 1 0 0\nvn 0 1 0\n \nf 1//3 2//1 3//2\nf 1/0/3 2/102/1 3/14/2";

        std::fs::write(&path, contents)?; // ? → I/O errors become the Err
        let p = Parser::parse_obj_file(&path)?;
        std::fs::remove_file(&path)?; // cleaned up BEFORE the asserts

        println!("{} {}", path.display(), p.report());

        let g = &p.default_group;
        assert!(g.len() == 2);

        let t1 = &g[0];
        let t2 = &g[1];
        assert!(t1.p1.approx_eq(p.vertices[1]));
        assert!(t1.p2.approx_eq(p.vertices[2]));
        assert!(t1.p3.approx_eq(p.vertices[3]));
        assert!(t1.n1.approx_eq(p.normals[3]));
        assert!(t1.n2.approx_eq(p.normals[1]));
        assert!(t1.n3.approx_eq(p.normals[2]));

        // the `v/vt/vn` form parses to the same triangle as `v//vn`
        assert!(t2.p1.approx_eq(t1.p1));
        assert!(t2.p2.approx_eq(t1.p2));
        assert!(t2.p3.approx_eq(t1.p3));
        assert!(t2.n1.approx_eq(t1.n1));
        assert!(t2.n2.approx_eq(t1.n2));
        assert!(t2.n3.approx_eq(t1.n3));

        Ok(())
    }
}
