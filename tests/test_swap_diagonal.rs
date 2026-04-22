use pycdt_rs::build::triangulate;
use pycdt_rs::topology::swap_diagonal;
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
fn test_swap_diagonal_basic_vertices_change() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0], [1.0, 1.0]];
    let mut t = triangulate(points);
    let (tri1, tri2) = find_two_adjacent_triangles(&t).expect("no adjacent triangles found");

    let verts1_before: std::collections::HashSet<usize> =
        t.triangle_vertices[tri1].iter().copied().collect();
    let verts2_before: std::collections::HashSet<usize> =
        t.triangle_vertices[tri2].iter().copied().collect();
    let shared_before: std::collections::HashSet<usize> = verts1_before
        .intersection(&verts2_before)
        .copied()
        .collect();
    assert_eq!(shared_before.len(), 2);

    swap_diagonal(&mut t, tri1, tri2);

    let verts1_after: std::collections::HashSet<usize> =
        t.triangle_vertices[tri1].iter().copied().collect();
    let verts2_after: std::collections::HashSet<usize> =
        t.triangle_vertices[tri2].iter().copied().collect();
    let shared_after: std::collections::HashSet<usize> =
        verts1_after.intersection(&verts2_after).copied().collect();

    assert_eq!(shared_after.len(), 2);
    assert_ne!(
        shared_before, shared_after,
        "diagonal should have changed after swap"
    );
}

#[test]
fn test_swap_diagonal_triangles_still_neighbors() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0], [1.0, 1.0]];
    let mut t = triangulate(points);
    let (tri1, tri2) = find_two_adjacent_triangles(&t).expect("no adjacent triangles found");

    swap_diagonal(&mut t, tri1, tri2);

    assert!(
        t.triangle_neighbors[tri1].contains(&tri2),
        "tri1 should have tri2 as neighbor after swap"
    );
    assert!(
        t.triangle_neighbors[tri2].contains(&tri1),
        "tri2 should have tri1 as neighbor after swap"
    );
}

#[test]
fn test_swap_diagonal_neighbor_consistency() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0], [1.0, 1.0]];
    let mut t = triangulate(points);
    let (tri1, tri2) = find_two_adjacent_triangles(&t).expect("no adjacent triangles found");

    swap_diagonal(&mut t, tri1, tri2);

    assert_neighbors_consistent(&t);
}

#[test]
fn test_swap_diagonal_new_diagonal_in_both_triangles() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0], [1.0, 1.0]];
    let mut t = triangulate(points);
    let (tri1, tri2) = find_two_adjacent_triangles(&t).expect("no adjacent triangles found");

    let opp1 = {
        let vi: std::collections::HashSet<usize> =
            t.triangle_vertices[tri1].iter().copied().collect();
        let vj: std::collections::HashSet<usize> =
            t.triangle_vertices[tri2].iter().copied().collect();
        let shared: std::collections::HashSet<usize> = vi.intersection(&vj).copied().collect();
        *vi.difference(&shared).next().unwrap()
    };
    let opp2 = {
        let vi: std::collections::HashSet<usize> =
            t.triangle_vertices[tri1].iter().copied().collect();
        let vj: std::collections::HashSet<usize> =
            t.triangle_vertices[tri2].iter().copied().collect();
        let shared: std::collections::HashSet<usize> = vi.intersection(&vj).copied().collect();
        *vj.difference(&shared).next().unwrap()
    };

    swap_diagonal(&mut t, tri1, tri2);

    let verts1: std::collections::HashSet<usize> =
        t.triangle_vertices[tri1].iter().copied().collect();
    let verts2: std::collections::HashSet<usize> =
        t.triangle_vertices[tri2].iter().copied().collect();

    assert!(
        verts1.contains(&opp1),
        "new diagonal vertex opp1 not in tri1"
    );
    assert!(
        verts1.contains(&opp2),
        "new diagonal vertex opp2 not in tri1"
    );
    assert!(
        verts2.contains(&opp1),
        "new diagonal vertex opp1 not in tri2"
    );
    assert!(
        verts2.contains(&opp2),
        "new diagonal vertex opp2 not in tri2"
    );
}

#[test]
fn test_swap_diagonal_multiple_consecutive_swaps() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [3.0, 0.0], [3.0, 3.0], [0.0, 3.0], [1.5, 1.5]];
    let mut t = triangulate(points);

    for _ in 0..3 {
        if let Some((tri1, tri2)) = find_two_adjacent_triangles(&t) {
            swap_diagonal(&mut t, tri1, tri2);
        }
    }

    assert_neighbors_consistent(&t);
    assert!(t.num_triangles() > 0);
}

#[test]
fn test_swap_diagonal_explicit_two_triangle_setup() {
    let mut t = Triangulation {
        points: vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
        triangle_vertices: vec![[0, 1, 2], [3, 2, 1]],
        triangle_neighbors: vec![[1, NO_NEIGHBOR, NO_NEIGHBOR], [0, NO_NEIGHBOR, NO_NEIGHBOR]],
        constrained_edges: Default::default(),
        num_super_triangle_points: 0,
    };

    swap_diagonal(&mut t, 0, 1);

    let v0: std::collections::HashSet<usize> = t.triangle_vertices[0].iter().copied().collect();
    let v1: std::collections::HashSet<usize> = t.triangle_vertices[1].iter().copied().collect();

    assert!(
        v0.contains(&0) && v0.contains(&3),
        "tri0 must contain new diagonal vertices 0 and 3"
    );
    assert!(
        v1.contains(&0) && v1.contains(&3),
        "tri1 must contain new diagonal vertices 0 and 3"
    );

    assert!(t.triangle_neighbors[0].contains(&1));
    assert!(t.triangle_neighbors[1].contains(&0));
}

#[test]
fn test_swap_diagonal_preserves_vertex_union() {
    let points: &[[f64; 2]] = &[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0], [1.0, 1.0]];
    let mut t = triangulate(points);
    let (tri1, tri2) = find_two_adjacent_triangles(&t).expect("no adjacent triangles found");

    let verts_before: std::collections::HashSet<usize> = t.triangle_vertices[tri1]
        .iter()
        .chain(t.triangle_vertices[tri2].iter())
        .copied()
        .collect();

    swap_diagonal(&mut t, tri1, tri2);

    let verts_after: std::collections::HashSet<usize> = t.triangle_vertices[tri1]
        .iter()
        .chain(t.triangle_vertices[tri2].iter())
        .copied()
        .collect();

    assert_eq!(
        verts_before, verts_after,
        "vertex union of the two triangles must be preserved after swap"
    );
}

#[test]
fn test_swap_diagonal_with_outer_neighbors_updates_back_refs() {
    let mut t = Triangulation {
        points: vec![[0.0, 1.0], [1.0, 2.0], [2.0, 1.0], [1.0, 0.0], [1.0, 1.0]],
        triangle_vertices: vec![[0, 4, 1], [1, 4, 2], [0, 3, 4], [3, 2, 4]],
        triangle_neighbors: vec![
            [NO_NEIGHBOR, NO_NEIGHBOR, NO_NEIGHBOR],
            [NO_NEIGHBOR, NO_NEIGHBOR, NO_NEIGHBOR],
            [NO_NEIGHBOR, 3, 0],
            [NO_NEIGHBOR, NO_NEIGHBOR, 2],
        ],
        constrained_edges: Default::default(),
        num_super_triangle_points: 0,
    };

    swap_diagonal(&mut t, 2, 3);

    assert!(t.triangle_neighbors[2].contains(&3));
    assert!(t.triangle_neighbors[3].contains(&2));
}
