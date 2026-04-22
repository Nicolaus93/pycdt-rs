use std::collections::HashSet;

use crate::types::{Point, TriangleNeighbors, TriangleVertices, NO_NEIGHBOR};

pub struct Triangulation {
    pub points: Vec<Point>,
    pub triangle_vertices: Vec<TriangleVertices>,
    pub triangle_neighbors: Vec<TriangleNeighbors>,
    pub constrained_edges: HashSet<(usize, usize)>,
    pub num_super_triangle_points: usize,
}

impl Triangulation {
    pub fn new() -> Self {
        Self {
            points: Vec::new(),
            triangle_vertices: Vec::new(),
            triangle_neighbors: Vec::new(),
            constrained_edges: HashSet::new(),
            num_super_triangle_points: 0,
        }
    }

    pub fn num_triangles(&self) -> usize {
        self.triangle_vertices.len()
    }

    pub fn num_points(&self) -> usize {
        self.points.len()
    }

    pub fn find_adjacent_triangle(&self, tri: usize, v1: usize, v2: usize) -> Option<usize> {
        let neighbors = self.triangle_neighbors.get(tri)?;
        let edge = Self::edge_key(v1, v2);

        for &neighbor_idx in neighbors {
            if neighbor_idx == NO_NEIGHBOR {
                continue;
            }

            let Some(vertices) = self.triangle_vertices.get(neighbor_idx) else {
                continue;
            };

            let mut found = false;
            for i in 0..3 {
                for j in (i + 1)..3 {
                    if Self::edge_key(vertices[i], vertices[j]) == edge {
                        found = true;
                        break;
                    }
                }
                if found {
                    break;
                }
            }

            if found {
                return Some(neighbor_idx);
            }
        }

        None
    }

    pub fn edge_key(a: usize, b: usize) -> (usize, usize) {
        (a.min(b), a.max(b))
    }
}

impl Default for Triangulation {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_creates_empty_triangulation() {
        let tri = Triangulation::new();

        assert!(tri.points.is_empty());
        assert!(tri.triangle_vertices.is_empty());
        assert!(tri.triangle_neighbors.is_empty());
        assert!(tri.constrained_edges.is_empty());
        assert_eq!(tri.num_super_triangle_points, 0);
    }

    #[test]
    fn num_counts_reflect_vectors() {
        let mut tri = Triangulation::new();
        tri.points.push([0.0, 0.0]);
        tri.points.push([1.0, 0.0]);
        tri.triangle_vertices.push([0, 1, 0]);

        assert_eq!(tri.num_points(), 2);
        assert_eq!(tri.num_triangles(), 1);
    }

    #[test]
    fn edge_key_normalizes_vertex_order() {
        assert_eq!(Triangulation::edge_key(3, 1), (1, 3));
    }
}
