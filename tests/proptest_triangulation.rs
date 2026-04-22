use proptest::prelude::*;
use pycdt_rs::build::triangulate;
use pycdt_rs::types::NO_NEIGHBOR;

fn assert_neighbors_consistent_result(t: &pycdt_rs::Triangulation) -> Result<(), TestCaseError> {
    let n = t.num_triangles();
    for i in 0..n {
        for j in 0..3 {
            let nb = t.triangle_neighbors[i][j];
            if nb != NO_NEIGHBOR {
                prop_assert!(nb < n, "neighbor index out of range: {} >= {}", nb, n);
                prop_assert!(
                    t.triangle_neighbors[nb].contains(&i),
                    "neighbor symmetry violated: tri {} lists {} as neighbor, but not vice versa",
                    i,
                    nb
                );
            }
        }
    }
    Ok(())
}

fn make_non_collinear_points(raw: Vec<(f64, f64)>) -> Option<Vec<[f64; 2]>> {
    if raw.len() < 3 {
        return None;
    }
    let mut pts: Vec<[f64; 2]> = raw.into_iter().map(|(x, y)| [x, y]).collect();
    pts.dedup_by(|a, b| (a[0] - b[0]).abs() < 1e-9 && (a[1] - b[1]).abs() < 1e-9);
    if pts.len() < 3 {
        return None;
    }
    let (p0, p1, p2) = (pts[0], pts[1], pts[2]);
    let cross = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p1[1] - p0[1]) * (p2[0] - p0[0]);
    if cross.abs() < 1e-9 {
        return None;
    }
    Some(pts)
}

proptest! {
    #[test]
    fn prop_triangulate_vertex_indices_in_range(
        raw in proptest::collection::vec(((-100.0f64..100.0), (-100.0f64..100.0)), 3..20)
    ) {
        if let Some(pts) = make_non_collinear_points(raw) {
            let t = triangulate(&pts);
            let n = t.num_points();
            for &[a, b, c] in &t.triangle_vertices {
                prop_assert!(a < n, "vertex {} out of range (n={})", a, n);
                prop_assert!(b < n, "vertex {} out of range (n={})", b, n);
                prop_assert!(c < n, "vertex {} out of range (n={})", c, n);
            }
        }
    }

    #[test]
    fn prop_triangulate_neighbor_symmetry(
        raw in proptest::collection::vec(((-100.0f64..100.0), (-100.0f64..100.0)), 3..20)
    ) {
        if let Some(pts) = make_non_collinear_points(raw) {
            let t = triangulate(&pts);
            assert_neighbors_consistent_result(&t)?;
        }
    }

    #[test]
    fn prop_triangulate_no_duplicate_triangles(
        raw in proptest::collection::vec(((-100.0f64..100.0), (-100.0f64..100.0)), 3..15)
    ) {
        if let Some(pts) = make_non_collinear_points(raw) {
            let t = triangulate(&pts);
            let mut seen = std::collections::HashSet::new();
            for &[a, b, c] in &t.triangle_vertices {
                let mut key = [a, b, c];
                key.sort();
                prop_assert!(seen.insert(key), "duplicate triangle {:?}", key);
            }
        }
    }

    #[test]
    fn prop_triangulate_no_super_triangle_vertices(
        raw in proptest::collection::vec(((-100.0f64..100.0), (-100.0f64..100.0)), 3..20)
    ) {
        if let Some(pts) = make_non_collinear_points(raw) {
            let n = pts.len();
            let t = triangulate(&pts);
            prop_assert_eq!(t.num_super_triangle_points, 0);
            prop_assert_eq!(t.num_points(), n);
            for &[a, b, c] in &t.triangle_vertices {
                prop_assert!(a < n);
                prop_assert!(b < n);
                prop_assert!(c < n);
            }
        }
    }

    #[test]
    fn prop_triangulate_at_least_one_triangle(
        raw in proptest::collection::vec(((-100.0f64..100.0), (-100.0f64..100.0)), 3..20)
    ) {
        if let Some(pts) = make_non_collinear_points(raw) {
            let t = triangulate(&pts);
            prop_assert!(t.num_triangles() > 0);
        }
    }
}

#[test]
fn regression_three_nearly_collinear_points_one_triangle() {
    let pts = [
        [-36.753982043234636f64, 42.84112520673667],
        [-69.64673440752645, -62.08972773575242],
        [-31.059590703171967, 62.63123281964592],
    ];
    let t = triangulate(&pts);
    assert_eq!(t.num_triangles(), 1);
}
