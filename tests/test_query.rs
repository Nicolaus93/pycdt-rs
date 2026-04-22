use pycdt_rs::build::{find_containing_triangle, triangulate, update_triangulation};
use pycdt_rs::constrained::segments_intersect;
use pycdt_rs::types::{PointLocation, NO_NEIGHBOR};
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
                    "triangle {} has neighbor {} but neighbor does not list it",
                    i,
                    nb
                );
            }
        }
    }
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
fn test_segments_intersect_touching_endpoint_not_proper_crossing() {
    let p1 = [0.0, 0.0];
    let p2 = [1.0, 1.0];
    let q1 = [1.0, 1.0];
    let q2 = [2.0, 0.0];
    assert!(!segments_intersect(&p1, &p2, &q1, &q2));
}

#[test]
fn test_segments_intersect_collinear_overlapping_not_crossing() {
    let p1 = [0.0, 0.0];
    let p2 = [2.0, 0.0];
    let q1 = [1.0, 0.0];
    let q2 = [3.0, 0.0];
    assert!(!segments_intersect(&p1, &p2, &q1, &q2));
}

#[test]
fn test_segments_intersect_collinear_non_overlapping() {
    let p1 = [0.0, 0.0];
    let p2 = [1.0, 0.0];
    let q1 = [2.0, 0.0];
    let q2 = [3.0, 0.0];
    assert!(!segments_intersect(&p1, &p2, &q1, &q2));
}

#[test]
fn test_segments_intersect_t_shape_not_proper_crossing() {
    let p1 = [0.0, 0.0];
    let p2 = [2.0, 0.0];
    let q1 = [1.0, 0.0];
    let q2 = [1.0, 1.0];
    assert!(!segments_intersect(&p1, &p2, &q1, &q2));
}

#[test]
fn test_segments_intersect_perpendicular_non_intersecting() {
    let p1 = [0.0, 0.0];
    let p2 = [1.0, 0.0];
    let q1 = [2.0, -1.0];
    let q2 = [2.0, 1.0];
    assert!(!segments_intersect(&p1, &p2, &q1, &q2));
}

#[test]
fn test_segments_intersect_diagonal_crossing() {
    let p1 = [-1.0, -1.0];
    let p2 = [1.0, 1.0];
    let q1 = [-1.0, 1.0];
    let q2 = [1.0, -1.0];
    assert!(segments_intersect(&p1, &p2, &q1, &q2));
}

#[test]
fn test_find_containing_triangle_center() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let t = triangulate(points);

    let point = [0.5, 0.5];
    let result = find_containing_triangle(&t, &point);
    match result {
        PointLocation::Interior(idx) | PointLocation::OnEdge(idx, _) => {
            assert!(idx < t.num_triangles())
        }
        PointLocation::NotFound => panic!("center point not found"),
    }
}

#[test]
fn test_find_containing_triangle_vertex() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let t = triangulate(points);

    let point = [0.0, 0.0];
    let result = find_containing_triangle(&t, &point);
    match result {
        PointLocation::Interior(idx) | PointLocation::OnEdge(idx, _) => {
            assert!(idx < t.num_triangles())
        }
        PointLocation::NotFound => panic!("vertex point not found"),
    }
}

#[test]
fn test_find_containing_triangle_outside_returns_not_found() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let t = triangulate(points);

    let point = [10.0, 10.0];
    let result = find_containing_triangle(&t, &point);
    assert_eq!(
        result,
        PointLocation::NotFound,
        "point far outside should return NotFound"
    );
}

#[test]
fn test_find_containing_triangle_near_edge() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
    let t = triangulate(points);

    let point = [1.0, 0.1];
    let result = find_containing_triangle(&t, &point);
    match result {
        PointLocation::Interior(idx) | PointLocation::OnEdge(idx, _) => {
            assert!(idx < t.num_triangles())
        }
        PointLocation::NotFound => panic!("near-edge point not found"),
    }
}

#[test]
fn test_update_triangulation_adds_points() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [3.0, 0.0], [3.0, 3.0], [0.0, 3.0], [1.5, 1.5]];
    let mut t = triangulate(points);
    let new_points: &[[f64; 2]] = &[[0.8, 0.8], [2.2, 0.8], [1.5, 2.2]];
    update_triangulation(&mut t, new_points);
    assert_eq!(t.num_points(), 8);
}

#[test]
fn test_update_triangulation_preserves_neighbor_consistency() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0], [2.0, 2.0]];
    let mut t = triangulate(points);
    let new_points: &[[f64; 2]] = &[[1.0, 1.0], [3.0, 1.0], [3.0, 3.0], [1.0, 3.0], [2.0, 1.5]];
    update_triangulation(&mut t, new_points);
    assert_neighbors_consistent(&t);
}

#[test]
fn test_update_triangulation_no_out_of_range_vertices() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [3.0, 0.0], [3.0, 3.0], [0.0, 3.0], [1.5, 1.5]];
    let mut t = triangulate(points);
    let new_points: &[[f64; 2]] = &[[0.8, 0.8], [2.2, 0.8], [1.5, 2.2]];
    update_triangulation(&mut t, new_points);

    let n = t.num_points();
    for &[a, b, c] in &t.triangle_vertices {
        assert!(a < n, "vertex {} out of range", a);
        assert!(b < n, "vertex {} out of range", b);
        assert!(c < n, "vertex {} out of range", c);
    }
}
