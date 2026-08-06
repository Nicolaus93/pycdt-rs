use pycdt_rs::build::{find_containing_triangle_from, triangulate};
use pycdt_rs::types::{PointLocation, NO_NEIGHBOR};
use pycdt_rs::Triangulation;
use rand::{rngs::StdRng, Rng, SeedableRng};

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

#[test]
fn test_neighbor_indices_valid_after_triangulation() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
    let t = triangulate(points);
    let n = t.num_triangles();

    for i in 0..n {
        for &nb in &t.triangle_neighbors[i] {
            assert!(
                nb == NO_NEIGHBOR || nb < n,
                "triangle {} has invalid neighbor {} (n={})",
                i,
                nb,
                n
            );
        }
    }
}

#[test]
fn test_neighbor_symmetry_after_triangulation() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
    let t = triangulate(points);
    assert_neighbors_consistent(&t);
}

#[test]
fn test_square_gives_two_triangles() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
    let t = triangulate(points);

    assert_eq!(t.num_triangles(), 2);

    for &verts in &t.triangle_vertices {
        assert_eq!(verts.len(), 3);
        assert_ne!(verts[0], verts[1]);
        assert_ne!(verts[1], verts[2]);
        assert_ne!(verts[0], verts[2]);
    }

    for &neighbors in &t.triangle_neighbors {
        assert_eq!(neighbors.len(), 3);
    }
}

#[test]
fn test_no_super_triangle_vertices_after_triangulation() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let t = triangulate(points);

    assert_eq!(t.num_super_triangle_points, 0);
    assert_eq!(t.num_points(), points.len());

    let n = points.len();
    for &[a, b, c] in &t.triangle_vertices {
        assert!(a < n, "vertex {} out of range", a);
        assert!(b < n, "vertex {} out of range", b);
        assert!(c < n, "vertex {} out of range", c);
    }
}

#[test]
fn test_triangulation_with_grid_neighbor_validity() {
    let mut points: Vec<[f64; 2]> = Vec::new();
    for i in 0..5 {
        for j in 0..5 {
            points.push([i as f64 * 0.75, j as f64 * 0.75]);
        }
    }
    let t = triangulate(&points);
    let n = t.num_triangles();

    for i in 0..n {
        for &nb in &t.triangle_neighbors[i] {
            assert!(
                nb == NO_NEIGHBOR || nb < n,
                "invalid neighbor index {} for triangle {} (n={})",
                nb,
                i,
                n
            );
        }
    }
}

#[test]
fn test_find_containing_triangle_center() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let t = triangulate(points);

    let point = [0.5, 0.5];
    let result = find_containing_triangle_from(&t, &point, 0);
    match result {
        PointLocation::Interior(idx) => assert!(idx < t.num_triangles()),
        PointLocation::OnEdge(idx, _) => assert!(idx < t.num_triangles()),
        PointLocation::NotFound => panic!("point in center not found in any triangle"),
    }
}

#[test]
fn test_find_containing_triangle_multiple_interior_points() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
    let t = triangulate(points);

    let test_pts = [[0.5, 0.5], [1.5, 0.5], [0.5, 1.5], [1.0, 1.0]];
    for p in &test_pts {
        let result = find_containing_triangle_from(&t, p, 0);
        match result {
            PointLocation::Interior(idx) => assert!(idx < t.num_triangles()),
            PointLocation::OnEdge(idx, _) => assert!(idx < t.num_triangles()),
            PointLocation::NotFound => panic!("point {:?} not found in any triangle", p),
        }
    }
}

#[test]
fn test_find_containing_triangle_near_edge() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
    let t = triangulate(points);

    let p = [1.0, 0.1];
    let result = find_containing_triangle_from(&t, &p, 0);
    match result {
        PointLocation::Interior(idx) | PointLocation::OnEdge(idx, _) => {
            assert!(idx < t.num_triangles())
        }
        PointLocation::NotFound => panic!("point near edge not found"),
    }
}

#[test]
fn test_triangulation_preserves_structure() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
    let t = triangulate(points);

    assert_eq!(t.num_triangles(), 2);

    for &verts in &t.triangle_vertices {
        assert_eq!(verts.len(), 3);
    }
    for &neighbors in &t.triangle_neighbors {
        assert_eq!(neighbors.len(), 3);
    }
}

#[test]
fn test_uniform_10k_preserves_convex_hull_triangle_count() {
    let mut rng = StdRng::seed_from_u64(0x1234_5678);
    let points: Vec<[f64; 2]> = (0..10_000)
        .map(|_| [rng.gen::<f64>() * 1000.0, rng.gen::<f64>() * 1000.0])
        .collect();

    let triangulation = triangulate(&points);
    let boundary_edges = triangulation
        .triangle_neighbors
        .iter()
        .flatten()
        .filter(|&&neighbor| neighbor == NO_NEIGHBOR)
        .count();

    // This deterministic set has 29 convex-hull vertices, so Euler gives
    // 2n - 2 - h = 19,969 triangles. A nearby finite super triangle used to
    // leave eight hull triangles attached to its vertices and remove them.
    assert_eq!(boundary_edges, 29);
    assert_eq!(triangulation.num_triangles(), 19_969);
}
