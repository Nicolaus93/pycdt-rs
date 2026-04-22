use pycdt_rs::build::triangulate;
use pycdt_rs::constrained::{add_constraints, find_intersecting_edges, segments_intersect};
use pycdt_rs::types::NO_NEIGHBOR;
use pycdt_rs::Triangulation;

fn assert_neighbors_consistent(t: &Triangulation) {
    let n = t.num_triangles();
    for i in 0..n {
        for j in 0..3 {
            let nb = t.triangle_neighbors[i][j];
            if nb != NO_NEIGHBOR {
                assert!(nb < n, "triangle {} neighbor {} out of range", i, nb);
                assert!(
                    t.triangle_neighbors[nb].contains(&i),
                    "triangle {} has neighbor {} but {} does not list {} as neighbor",
                    i,
                    nb,
                    nb,
                    i
                );
            }
        }
    }
}

fn find_two_adjacent_triangles(t: &Triangulation) -> Option<(usize, usize)> {
    for i in 0..t.num_triangles() {
        for j in (i + 1)..t.num_triangles() {
            let vi: std::collections::HashSet<usize> =
                t.triangle_vertices[i].iter().copied().collect();
            let vj: std::collections::HashSet<usize> =
                t.triangle_vertices[j].iter().copied().collect();
            if vi.intersection(&vj).count() == 2 {
                return Some((i, j));
            }
        }
    }
    None
}

#[test]
fn test_segments_intersect_parallel_no_intersection() {
    let p1 = [0.0, 0.0];
    let p2 = [1.0, 0.0];
    let q1 = [0.0, 1.0];
    let q2 = [1.0, 1.0];
    assert!(!segments_intersect(&p1, &p2, &q1, &q2));
}

#[test]
fn test_segments_intersect_crossing() {
    let p1 = [0.0, 0.0];
    let p2 = [1.0, 1.0];
    let q1 = [0.0, 1.0];
    let q2 = [1.0, 0.0];
    assert!(segments_intersect(&p1, &p2, &q1, &q2));
}

#[test]
fn test_segments_intersect_collinear_overlapping_no_crossing() {
    let p1 = [0.0, 0.0];
    let p2 = [2.0, 0.0];
    let q1 = [1.0, 0.0];
    let q2 = [3.0, 0.0];
    assert!(!segments_intersect(&p1, &p2, &q1, &q2));
}

#[test]
fn test_segments_intersect_diagonal() {
    let p1 = [-1.0, -1.0];
    let p2 = [1.0, 1.0];
    let q1 = [-1.0, 1.0];
    let q2 = [1.0, -1.0];
    assert!(segments_intersect(&p1, &p2, &q1, &q2));
}

#[test]
fn test_find_intersecting_edges_existing_edge_returns_empty() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]];
    let t = triangulate(points);
    let edges = find_intersecting_edges(&t, 0, 1);
    assert_eq!(edges.len(), 0);
}

#[test]
fn test_find_intersecting_edges_returns_valid_entries() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
    let t = triangulate(points);
    let n = t.num_points();
    let mut found = false;
    'outer: for v1 in 0..n {
        for v2 in (v1 + 1)..n {
            let edges = find_intersecting_edges(&t, v1, v2);
            if edges.len() > 0 {
                for &(t1, t2) in &edges {
                    assert!(t1 < t.num_triangles());
                    assert!(t2 < t.num_triangles());
                }
                found = true;
                break 'outer;
            }
        }
    }
    assert!(
        found,
        "at least one vertex pair should yield intersecting edges in a square triangulation"
    );
}

#[test]
fn test_find_intersecting_edges_returns_vec() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0], [1.0, 0.8]];
    let t = triangulate(points);
    let n = t.num_points();
    let mut found_intersecting = false;
    for v1 in 0..n {
        for v2 in (v1 + 1)..n {
            let edges = find_intersecting_edges(&t, v1, v2);
            if edges.len() > 0 {
                found_intersecting = true;
                for &(t1, t2) in &edges {
                    assert!(t1 < t.num_triangles());
                    assert!(t2 < t.num_triangles());
                }
            }
        }
    }
    assert!(
        found_intersecting,
        "5-point square should have some non-existing edges"
    );
}

#[test]
fn test_add_constraints_single() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
    let mut t = triangulate(points);
    let num_triangles_before = t.num_triangles();

    let success = add_constraints(&mut t, &[(0, 2)]);

    assert!(success);
    assert_eq!(t.num_triangles(), num_triangles_before);
}

#[test]
fn test_add_constraints_existing_edge() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
    let mut t = triangulate(points);
    let num_triangles_before = t.num_triangles();

    let success = add_constraints(&mut t, &[(0, 1)]);

    assert!(success);
    assert_eq!(t.num_triangles(), num_triangles_before);
}

#[test]
fn test_add_constraints_preserves_connectivity() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [3.0, 0.0], [3.0, 3.0], [0.0, 3.0], [1.5, 1.5]];
    let mut t = triangulate(points);

    let success = add_constraints(&mut t, &[(0, 2)]);
    assert!(success);

    for i in 0..t.num_triangles() {
        for &nb in &t.triangle_neighbors[i] {
            if nb != NO_NEIGHBOR {
                assert!(nb < t.num_triangles());
                assert!(t.triangle_neighbors[nb].contains(&i));
            }
        }
    }
}

#[test]
fn test_add_constraints_maintains_point_count() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0], [1.0, 1.0]];
    let mut t = triangulate(points);
    let num_points_before = t.num_points();

    add_constraints(&mut t, &[(0, 2)]);

    assert_eq!(t.num_points(), num_points_before);
}

#[test]
fn test_add_constraints_multiple() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
    let mut t = triangulate(points);

    let success = add_constraints(&mut t, &[(0, 2), (1, 3)]);
    assert!(success);

    assert!(t.num_triangles() > 0);
    for &[a, b, c] in &t.triangle_vertices {
        assert!(a < t.num_points());
        assert!(b < t.num_points());
        assert!(c < t.num_points());
    }
}

#[test]
fn test_add_constraints_diagonal_grid() {
    let mut points: Vec<[f64; 2]> = Vec::new();
    for i in 0..5 {
        for j in 0..5 {
            points.push([i as f64, j as f64]);
        }
    }
    let n = points.len();
    let mut t = triangulate(&points);
    let num_triangles_before = t.num_triangles();

    let success = add_constraints(&mut t, &[(0, n - 1)]);
    assert!(success);

    assert_eq!(t.num_triangles(), num_triangles_before);

    for &[a, b, c] in &t.triangle_vertices {
        assert_eq!(
            {
                let mut v = [a, b, c];
                v.sort();
                v.len()
            },
            3
        );
        let mut verts = [a, b, c];
        verts.sort();
        assert_eq!(verts[0], verts[0]);
        assert_ne!(a, b);
        assert_ne!(b, c);
        assert_ne!(a, c);
    }
}

#[test]
fn test_add_constraints_neighbor_consistency() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [3.0, 0.0], [3.0, 3.0], [0.0, 3.0], [1.5, 1.5]];
    let mut t = triangulate(points);
    add_constraints(&mut t, &[(0, 2)]);
    assert_neighbors_consistent(&t);
}

#[test]
fn test_add_constraints_marks_constrained_edges() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0], [1.0, 1.0]];
    let mut t = triangulate(points);
    add_constraints(&mut t, &[(0, 2)]);

    let edge = pycdt_rs::Triangulation::edge_key(0, 2);
    assert!(t.constrained_edges.contains(&edge));
}

#[test]
fn test_add_constraints_two_adjacent_triangles_consistent() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [1.0, 0.0], [0.5, 1.0], [1.5, 1.0]];
    let mut t = triangulate(points);
    if let Some((i, j)) = find_two_adjacent_triangles(&t) {
        let vi: Vec<usize> = t.triangle_vertices[i].iter().copied().collect();
        let vj: Vec<usize> = t.triangle_vertices[j].iter().copied().collect();
        let shared: Vec<usize> = vi.iter().copied().filter(|v| vj.contains(v)).collect();
        if shared.len() == 2 {
            add_constraints(&mut t, &[(shared[0], shared[1])]);
            assert_neighbors_consistent(&t);
        }
    }
}
