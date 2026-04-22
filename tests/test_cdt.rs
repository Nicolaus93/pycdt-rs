use pycdt_rs::build::triangulate;
use pycdt_rs::types::NO_NEIGHBOR;
use pycdt_rs::Triangulation;

fn assert_neighbors_consistent(t: &Triangulation) {
    let n = t.num_triangles();
    for i in 0..n {
        for j in 0..3 {
            let nb = t.triangle_neighbors[i][j];
            if nb != NO_NEIGHBOR {
                assert!(nb < n, "neighbor index {} out of range (n={})", nb, n);
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

fn assert_all_vertex_indices_in_range(t: &Triangulation) {
    let n = t.num_points();
    for &[a, b, c] in &t.triangle_vertices {
        assert!(a < n, "vertex {} out of range (n={})", a, n);
        assert!(b < n, "vertex {} out of range (n={})", b, n);
        assert!(c < n, "vertex {} out of range (n={})", c, n);
    }
}

#[test]
fn test_triangulate_circle_points() {
    let points: &[[f64; 2]] = &[
        [24.311, -7.358],
        [23.574, -9.456],
        [22.657, -11.481],
        [21.566, -13.419],
        [20.31, -15.253],
        [18.899, -16.971],
        [17.342, -18.558],
        [15.653, -20.004],
        [13.844, -21.296],
        [11.929, -22.425],
        [6.53, -24.546],
        [0.791, -25.388],
        [-4.989, -24.905],
        [-10.509, -23.124],
        [-15.48, -20.137],
        [-19.645, -16.1],
        [-22.786, -11.224],
        [-24.738, -5.762],
        [-25.4, 0.0],
        [-23.868, 8.687],
        [-19.458, 16.327],
        [-12.7, 21.997],
        [-4.411, 25.014],
        [4.411, 25.014],
        [12.7, 21.997],
        [19.458, 16.327],
        [23.868, 8.687],
        [25.4, 0.0],
        [25.386, -0.829],
        [25.346, -1.658],
        [25.278, -2.485],
        [25.184, -3.309],
        [25.062, -4.129],
        [24.914, -4.945],
        [24.739, -5.756],
        [24.538, -6.561],
    ];

    let t = triangulate(points);

    assert!(t.num_triangles() > 0, "expected at least one triangle");
    assert_eq!(t.num_points(), points.len());
    assert_eq!(t.num_super_triangle_points, 0);
    assert_all_vertex_indices_in_range(&t);
    assert_neighbors_consistent(&t);
}

#[test]
fn test_triangulate_square() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
    let t = triangulate(points);
    assert_eq!(t.num_triangles(), 2);
    assert_eq!(t.num_points(), 4);
    assert_eq!(t.num_super_triangle_points, 0);
    assert_all_vertex_indices_in_range(&t);
    assert_neighbors_consistent(&t);
}

#[test]
fn test_triangulate_five_points() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.5, 0.5]];
    let t = triangulate(points);
    assert!(t.num_triangles() > 0);
    assert_eq!(t.num_points(), 5);
    assert_all_vertex_indices_in_range(&t);
    assert_neighbors_consistent(&t);
}

#[test]
fn test_triangulate_grid() {
    let mut points: Vec<[f64; 2]> = Vec::new();
    for i in 0..5 {
        for j in 0..5 {
            points.push([i as f64, j as f64]);
        }
    }
    let t = triangulate(&points);
    assert!(t.num_triangles() > 0);
    assert_eq!(t.num_points(), 25);
    assert_all_vertex_indices_in_range(&t);
    assert_neighbors_consistent(&t);
}

#[test]
fn test_triangulate_no_duplicate_triangles() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [3.0, 0.0], [3.0, 3.0], [0.0, 3.0], [1.5, 1.5]];
    let t = triangulate(points);

    let mut seen = std::collections::HashSet::new();
    for &[a, b, c] in &t.triangle_vertices {
        let mut key = [a, b, c];
        key.sort();
        assert!(seen.insert(key), "duplicate triangle {:?}", key);
    }
}

#[test]
fn test_triangulate_euler_formula() {
    let points: &[[f64; 2]] = &[
        [0.0, 0.0],
        [2.0, 0.0],
        [4.0, 0.0],
        [4.0, 2.0],
        [4.0, 4.0],
        [2.0, 4.0],
        [0.0, 4.0],
        [0.0, 2.0],
        [2.0, 2.0],
    ];
    let t = triangulate(points);
    let n = points.len();
    assert!(t.num_triangles() >= n - 2);
    assert_all_vertex_indices_in_range(&t);
    assert_neighbors_consistent(&t);
}
