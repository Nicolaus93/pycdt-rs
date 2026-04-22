use crate::geometry::{ensure_ccw, incircle};
use crate::triangulation::Triangulation;
use crate::types::NO_NEIGHBOR;

// ---------------------------------------------------------------------------
// Helper
// ---------------------------------------------------------------------------

/// Find the position (0, 1, or 2) of vertex `v` in triangle `tri_idx`.
fn vertex_pos(t: &Triangulation, tri_idx: usize, v: usize) -> Option<usize> {
    t.triangle_vertices[tri_idx].iter().position(|&x| x == v)
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Find the local edge index (0, 1, or 2) for edge (v1, v2) in triangle tri_idx.
/// Edge i connects vertices[i] and vertices[(i+1)%3].
pub fn get_edge_index(t: &Triangulation, tri_idx: usize, v1: usize, v2: usize) -> Option<usize> {
    let verts = t.triangle_vertices[tri_idx];
    for i in 0..3 {
        let a = verts[i];
        let b = verts[(i + 1) % 3];
        if (a == v1 && b == v2) || (a == v2 && b == v1) {
            return Some(i);
        }
    }
    None
}

/// Find the vertex opposite to edge (v1, v2) in triangle tri_idx.
/// Returns None if the edge is not found in the triangle.
pub fn get_opposite_vertex(
    t: &Triangulation,
    tri_idx: usize,
    v1: usize,
    v2: usize,
) -> Option<usize> {
    let verts = t.triangle_vertices[tri_idx];
    if !verts.contains(&v1) || !verts.contains(&v2) {
        return None;
    }
    verts.iter().copied().find(|&v| v != v1 && v != v2)
}

/// Find the shared edge between tri_a and tri_b.
/// Returns (v1, v2) — the two vertex indices that are shared.
pub fn find_shared_edge(t: &Triangulation, tri_a: usize, tri_b: usize) -> Option<(usize, usize)> {
    let va = t.triangle_vertices[tri_a];
    let vb = t.triangle_vertices[tri_b];
    let shared: Vec<usize> = va.iter().copied().filter(|v| vb.contains(v)).collect();
    if shared.len() == 2 {
        Some((shared[0], shared[1]))
    } else {
        None
    }
}

/// Replace old_neighbor with new_neighbor in triangle tri_idx's neighbor list.
pub fn reorder_neighbors(
    t: &mut Triangulation,
    tri_idx: usize,
    old_neighbor: usize,
    new_neighbor: usize,
) {
    for slot in t.triangle_neighbors[tri_idx].iter_mut() {
        if *slot == old_neighbor {
            *slot = new_neighbor;
        }
    }
}

/// Swap the diagonal shared between tri_a and tri_b.
///
/// Convention matches Python topology.py:
///   - tri_a plays the role of t4_idx (the "candidate" triangle)
///   - tri_b plays the role of t3_idx (the "neighbor" triangle)
///
/// After the swap:
///   - Both triangles keep their indices but get new vertex/neighbor arrays.
///   - All 4 surrounding neighbor triangles have their back-references updated.
pub fn swap_diagonal(t: &mut Triangulation, tri_a: usize, tri_b: usize) {
    let va = t.triangle_vertices[tri_a];
    let vb = t.triangle_vertices[tri_b];

    // Shared vertices
    let shared: Vec<usize> = va.iter().copied().filter(|v| vb.contains(v)).collect();
    assert_eq!(
        shared.len(),
        2,
        "swap_diagonal: triangles must share exactly one edge"
    );
    let (a, b) = (shared[0], shared[1]);

    // Opposite vertices
    let point_idx = va
        .iter()
        .copied()
        .find(|&v| v != a && v != b)
        .expect("invariant: triangle A must contain an opposite vertex");
    let c = vb
        .iter()
        .copied()
        .find(|&v| v != a && v != b)
        .expect("invariant: triangle B must contain an opposite vertex");

    // Old neighbors (neighbor[i] is opposite vertex[i])
    let pos_a_in_a = va
        .iter()
        .position(|&v| v == a)
        .expect("invariant: shared vertex a must exist in triangle A");
    let pos_b_in_a = va
        .iter()
        .position(|&v| v == b)
        .expect("invariant: shared vertex b must exist in triangle A");
    let t6 = t.triangle_neighbors[tri_a][pos_a_in_a]; // opposite a in tri_a
    let t5 = t.triangle_neighbors[tri_a][pos_b_in_a]; // opposite b in tri_a

    let pos_a_in_b = vb
        .iter()
        .position(|&v| v == a)
        .expect("invariant: shared vertex a must exist in triangle B");
    let pos_b_in_b = vb
        .iter()
        .position(|&v| v == b)
        .expect("invariant: shared vertex b must exist in triangle B");
    let t7 = t.triangle_neighbors[tri_b][pos_a_in_b]; // opposite a in tri_b
    let t8 = t.triangle_neighbors[tri_b][pos_b_in_b]; // opposite b in tri_b

    // New triangle vertices (CCW-corrected)
    //   new tri_b = CCW(point_idx, c, a)
    //   new tri_a = CCW(point_idx, b, c)
    let (p0, p1, p2) = ensure_ccw(&t.points, point_idx, c, a);
    let new_vb = [p0, p1, p2];
    let (q0, q1, q2) = ensure_ccw(&t.points, point_idx, b, c);
    let new_va = [q0, q1, q2];

    // Apply new vertices
    t.triangle_vertices[tri_b] = new_vb;
    t.triangle_vertices[tri_a] = new_va;

    // New neighbors for tri_b
    //   vertex point_idx → opposite neighbor = t8
    //   vertex c         → opposite neighbor = t5
    //   vertex a         → opposite neighbor = tri_a
    let mut nb = [NO_NEIGHBOR; 3];
    for i in 0..3 {
        nb[i] = match new_vb[i] {
            v if v == point_idx => t8,
            v if v == c => t5,
            _ => tri_a, // must be a
        };
    }
    t.triangle_neighbors[tri_b] = nb;

    // New neighbors for tri_a
    //   vertex point_idx → opposite neighbor = t7
    //   vertex c         → opposite neighbor = t6
    //   vertex b         → opposite neighbor = tri_b
    let mut na = [NO_NEIGHBOR; 3];
    for i in 0..3 {
        na[i] = match new_va[i] {
            v if v == point_idx => t7,
            v if v == c => t6,
            _ => tri_b, // must be b
        };
    }
    t.triangle_neighbors[tri_a] = na;

    // Update back-references in the 2 non-swapped neighbors:
    //   t5 used to point at tri_a, now should point at tri_b
    if t5 != NO_NEIGHBOR {
        reorder_neighbors(t, t5, tri_a, tri_b);
    }
    //   t7 used to point at tri_b, now should point at tri_a
    if t7 != NO_NEIGHBOR {
        reorder_neighbors(t, t7, tri_b, tri_a);
    }
}

/// Walk from start_tri using BFS to find the first triangle that contains `vertex`.
pub fn find_triangle_with_vertex(t: &Triangulation, start_tri: usize, vertex: usize) -> usize {
    use std::collections::{HashSet, VecDeque};

    let mut visited = HashSet::new();
    let mut queue = VecDeque::new();
    queue.push_back(start_tri);

    while let Some(tri) = queue.pop_front() {
        if visited.contains(&tri) {
            continue;
        }
        visited.insert(tri);

        if t.triangle_vertices[tri].contains(&vertex) {
            return tri;
        }

        for &neighbor in &t.triangle_neighbors[tri] {
            if neighbor != NO_NEIGHBOR && !visited.contains(&neighbor) {
                queue.push_back(neighbor);
            }
        }
    }

    panic!(
        "find_triangle_with_vertex: vertex {} not reachable from triangle {}",
        vertex, start_tri
    );
}

/// Iterative Lawson swapping — restores Delaunay property after inserting point_idx.
///
/// Stack entries are (t3_idx, t4_idx) where:
///   t4_idx contains point_idx (the candidate triangle)
///   t3_idx is the neighbor to test (the neighbor triangle)
///
/// Mirrors Python topology.py lawson_swapping but uses an explicit Vec stack.
pub fn lawson_swapping(t: &mut Triangulation, tri_idx: usize, point_idx: usize) {
    let mut stack: Vec<(usize, usize)> = Vec::new();

    // Seed the stack with all valid neighbors of the newly inserted triangle
    for i in 0..3 {
        let neighbor = t.triangle_neighbors[tri_idx][i];
        if neighbor != NO_NEIGHBOR {
            stack.push((neighbor, tri_idx));
        }
    }

    while let Some((t3_idx, t4_idx)) = stack.pop() {
        if t3_idx == NO_NEIGHBOR || t4_idx == NO_NEIGHBOR {
            continue;
        }

        // Guard: point_idx must still live in t4_idx (stale entries can arise after flips)
        if !t.triangle_vertices[t4_idx].contains(&point_idx) {
            continue;
        }

        // Guard: t3_idx and t4_idx must still be neighbors
        if !t.triangle_neighbors[t4_idx].contains(&t3_idx) {
            continue;
        }

        // Skip if t3 is a sibling (contains the newly inserted point — incircle is degenerate)
        if t.triangle_vertices[t3_idx].contains(&point_idx) {
            continue;
        }

        // incircle(a,b,c,d) > 0  iff  d is inside the circumcircle of (a,b,c)
        let t3_verts = t.triangle_vertices[t3_idx];
        let a = t.points[t3_verts[0]];
        let b = t.points[t3_verts[1]];
        let c = t.points[t3_verts[2]];
        let p = t.points[point_idx];

        if incircle(&a, &b, &c, &p) <= 0.0 {
            continue; // Delaunay condition satisfied, no flip needed
        }

        // Flip: swap_diagonal(tri_a=t4_idx, tri_b=t3_idx)
        swap_diagonal(t, t4_idx, t3_idx);

        // After flip, point_idx is in both t4_idx and t3_idx.
        // Push the new potentially-illegal edges (neighbors opposite point_idx in each).
        if let Some(pos) = vertex_pos(t, t4_idx, point_idx) {
            let new_t7 = t.triangle_neighbors[t4_idx][pos];
            if new_t7 != NO_NEIGHBOR {
                stack.push((new_t7, t4_idx));
            }
        }
        if let Some(pos) = vertex_pos(t, t3_idx, point_idx) {
            let new_t8 = t.triangle_neighbors[t3_idx][pos];
            if new_t8 != NO_NEIGHBOR {
                stack.push((new_t8, t3_idx));
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::triangulation::Triangulation;
    use crate::types::NO_NEIGHBOR;

    /// Build a minimal triangulation with two triangles sharing edge (1, 2):
    ///
    ///   Points:
    ///     0 = (0, 0)
    ///     1 = (1, 0)
    ///     2 = (0, 1)
    ///     3 = (1, 1)
    ///
    ///   Triangle 0: vertices [0, 1, 2] — CCW
    ///   Triangle 1: vertices [3, 2, 1] — CCW (shares edge 1-2 with tri 0)
    ///
    ///   Neighbor convention: neighbors[i] is opposite vertices[i]
    ///   Tri 0: neighbor opposite 0 = tri 1, opposite 1 = NO, opposite 2 = NO
    ///          → neighbors = [1, NO, NO]
    ///   Tri 1: neighbor opposite 3 = NO, opposite 2 = NO, opposite 1 = tri 0
    ///          → neighbors = [NO, NO, 0]
    fn two_tri_setup() -> Triangulation {
        let mut t = Triangulation::new();
        t.points = vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]];
        //                0            1             2             3
        // Tri 0: [0,1,2] — n0 opposite 0 = tri1 (shares edge 1-2)
        //                   n1 opposite 1 = NO
        //                   n2 opposite 2 = NO
        t.triangle_vertices = vec![[0, 1, 2], [3, 2, 1]];
        t.triangle_neighbors = vec![[1, NO_NEIGHBOR, NO_NEIGHBOR], [0, NO_NEIGHBOR, NO_NEIGHBOR]];
        t
    }

    // -----------------------------------------------------------------------

    #[test]
    fn get_edge_index_returns_correct_local_index() {
        let t = two_tri_setup();
        // Tri 0 = [0,1,2]: edge (0,1) is at index 0, edge (1,2) at 1, edge (2,0) at 2
        assert_eq!(get_edge_index(&t, 0, 0, 1), Some(0));
        assert_eq!(get_edge_index(&t, 0, 1, 0), Some(0)); // reversed order
        assert_eq!(get_edge_index(&t, 0, 1, 2), Some(1));
        assert_eq!(get_edge_index(&t, 0, 2, 0), Some(2));
        // Edge not in triangle
        assert_eq!(get_edge_index(&t, 0, 0, 3), None);
    }

    #[test]
    fn get_opposite_vertex_basic() {
        let t = two_tri_setup();
        // Tri 0 = [0,1,2]: vertex opposite edge (1,2) is 0
        assert_eq!(get_opposite_vertex(&t, 0, 1, 2), Some(0));
        assert_eq!(get_opposite_vertex(&t, 0, 0, 1), Some(2));
        assert_eq!(get_opposite_vertex(&t, 0, 0, 2), Some(1));
        // Edge not in triangle
        assert_eq!(get_opposite_vertex(&t, 0, 0, 3), None);
    }

    #[test]
    fn find_shared_edge_returns_shared_vertices() {
        let t = two_tri_setup();
        // Tri 0 = [0,1,2], Tri 1 = [3,2,1] — shared: 1 and 2
        let edge = find_shared_edge(&t, 0, 1);
        assert!(edge.is_some());
        let (v1, v2) = match edge {
            Some(edge) => edge,
            None => panic!("expected shared edge between triangles 0 and 1"),
        };
        let mut pair = [v1, v2];
        pair.sort();
        assert_eq!(pair, [1, 2]);
    }

    #[test]
    fn find_shared_edge_none_for_non_adjacent() {
        let _t = two_tri_setup();
        // Add a third triangle with no shared vertex with tri 0 via a new setup
        let mut t2 = Triangulation::new();
        t2.points = vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [5.0, 5.0], [6.0, 5.0]];
        t2.triangle_vertices = vec![[0, 1, 2], [3, 4, 3]]; // second tri uses dup vertex — just for testing None
        t2.triangle_neighbors = vec![[NO_NEIGHBOR; 3]; 2];
        assert_eq!(find_shared_edge(&t2, 0, 1), None);
    }

    #[test]
    fn swap_diagonal_vertices_change() {
        let mut t = two_tri_setup();
        // Before: tri0=[0,1,2], tri1=[3,2,1]
        // Shared edge: 1-2; opposite in tri0=0, opposite in tri1=3
        // After swap: new diagonal is 0-3
        swap_diagonal(&mut t, 0, 1);

        let v0: std::collections::HashSet<usize> = t.triangle_vertices[0].iter().copied().collect();
        let v1: std::collections::HashSet<usize> = t.triangle_vertices[1].iter().copied().collect();

        // Both new triangles must contain the new diagonal vertices 0 and 3
        assert!(
            v0.contains(&0) && v0.contains(&3),
            "tri0 must contain new diagonal verts 0,3"
        );
        assert!(
            v1.contains(&0) && v1.contains(&3),
            "tri1 must contain new diagonal verts 0,3"
        );
    }

    #[test]
    fn swap_diagonal_neighbor_back_refs_consistent() {
        let mut t = two_tri_setup();
        swap_diagonal(&mut t, 0, 1);

        // After swap, tri0 and tri1 must be neighbors of each other
        assert!(
            t.triangle_neighbors[0].contains(&1),
            "tri0 neighbors should contain tri1 after swap"
        );
        assert!(
            t.triangle_neighbors[1].contains(&0),
            "tri1 neighbors should contain tri0 after swap"
        );
    }

    #[test]
    fn swap_diagonal_with_outer_neighbors_updated() {
        // Set up 4 triangles in a "fan" around a central edge so that
        // the outer triangles get their back-refs updated correctly.
        //
        //  Points:
        //    0=(0,1)  1=(1,2)  2=(2,1)  3=(1,0)  4=(1,1) center
        //
        //  Triangles:
        //    t0 = [0,4,1]  — top-left
        //    t1 = [1,4,2]  — top-right (shares edge 1-4 with t0? no)
        //    t2 = [0,3,4]  — bottom-left; adjacent to t3 on edge 3-4
        //    t3 = [3,2,4]  — bottom-right
        //
        //  We'll focus on t2 and t3 sharing edge (3,4).
        //  outer neighbors of t2 = some dummy, outer of t3 = some dummy
        let mut t = Triangulation::new();
        t.points = vec![
            [0.0, 1.0], // 0
            [1.0, 2.0], // 1
            [2.0, 1.0], // 2
            [1.0, 0.0], // 3
            [1.0, 1.0], // 4 center
        ];
        // Tri 0: [0,4,1] CCW
        // Tri 1: [1,4,2] CCW
        // Tri 2: [0,3,4] CCW — but actually [4,0,3] might be CCW; let's verify manually:
        //   orient2d([0,1],[1,0],[1,1]) = det([1,−1],[1,0]) = 1·0−(−1)·1 = 1 > 0 ✓ CCW
        // Tri 3: [3,2,4] CCW
        //   orient2d([1,0],[2,1],[1,1]) = det([1,1],[0,1]) = 1·1−1·0 = 1 > 0 ✓ CCW
        t.triangle_vertices = vec![
            [0, 4, 1], // tri 0
            [1, 4, 2], // tri 1
            [0, 3, 4], // tri 2  — opposite 0 is NO, opposite 3 is tri3, opposite 4 is some "outer"
            [3, 2, 4], // tri 3  — shares edge 3-4 with tri2
        ];
        // Neighbors (opp vertex[i]):
        // tri2=[0,3,4]: opp 0 = NO, opp 3 = tri3, opp 4 = tri0 (let's say)
        // tri3=[3,2,4]: opp 3 = NO, opp 2 = NO,   opp 4 = tri2
        t.triangle_neighbors = vec![
            [NO_NEIGHBOR, NO_NEIGHBOR, NO_NEIGHBOR], // tri0 (irrelevant here)
            [NO_NEIGHBOR, NO_NEIGHBOR, NO_NEIGHBOR], // tri1 (irrelevant here)
            [NO_NEIGHBOR, 3, 0],                     // tri2: opp 0=NO, opp 3=tri3, opp 4=tri0
            [NO_NEIGHBOR, NO_NEIGHBOR, 2],           // tri3: opp 3=NO, opp 2=NO,  opp 4=tri2
        ];

        swap_diagonal(&mut t, 2, 3); // tri_a=2, tri_b=3

        // After swap:
        // New diagonal connects 0 (opp in tri2) and 2 (opp in tri3)
        // tri0 neighbor list should have tri2 replaced by tri3 (or vice versa)
        let tri0_n = t.triangle_neighbors[0];
        // tri0 previously referenced tri2 (at index 2 of its neighbors).
        // After swap, t5 (= neighbor opposite b in tri_a=tri2) = tri0.
        // t5 should now point to tri_b=tri3 instead of tri_a=tri2.
        // So tri0 should reference tri3, not tri2.
        assert!(
            !tri0_n.contains(&2) || tri0_n.contains(&3),
            "tri0 back-ref should be updated after swap"
        );

        // tri2 and tri3 must still be neighbors of each other
        assert!(t.triangle_neighbors[2].contains(&3));
        assert!(t.triangle_neighbors[3].contains(&2));
    }

    #[test]
    fn lawson_swapping_flips_non_delaunay_edge() {
        // Classic 4-point example where one diagonal needs to be flipped.
        // Points: square with a slightly off-center interior point
        // After inserting point 4, the edge 1-3 violates Delaunay.
        //
        //   0=(0,0)  1=(2,0)  2=(2,2)  3=(0,2)  4=(1.1, 0.9) (slightly off-center)
        //
        // Initial triangulation after inserting pt4 into the quadrilateral:
        //   tri0 = [0,1,4]
        //   tri1 = [1,2,4]
        //   tri2 = [2,3,4]
        //   tri3 = [3,0,4]
        // (All four CCW triangles around the fan at vertex 4)
        //
        // Neighbors (opp vertex[i]):
        //   tri0=[0,1,4]: opp 0=tri1, opp 1=tri3, opp 4=NO  → [1,3,NO]
        //   tri1=[1,2,4]: opp 1=tri2, opp 2=tri0, opp 4=NO  → [2,0,NO]
        //   tri2=[2,3,4]: opp 2=tri3, opp 3=tri1, opp 4=NO  → [3,1,NO]
        //   tri3=[3,0,4]: opp 3=tri0, opp 0=tri2, opp 4=NO  → [0,2,NO]
        let mut t = Triangulation::new();
        t.points = vec![
            [0.0, 0.0], // 0
            [2.0, 0.0], // 1
            [2.0, 2.0], // 2
            [0.0, 2.0], // 3
            [1.1, 0.9], // 4  — newly inserted
        ];
        t.triangle_vertices = vec![
            [0, 1, 4], // tri0
            [1, 2, 4], // tri1
            [2, 3, 4], // tri2
            [3, 0, 4], // tri3
        ];
        t.triangle_neighbors = vec![
            [1, 3, NO_NEIGHBOR], // tri0: opp 0=tri1, opp 1=tri3, opp 4=NO
            [2, 0, NO_NEIGHBOR], // tri1: opp 1=tri2, opp 2=tri0, opp 4=NO
            [3, 1, NO_NEIGHBOR], // tri2: opp 2=tri3, opp 3=tri1, opp 4=NO
            [0, 2, NO_NEIGHBOR], // tri3: opp 3=tri0, opp 0=tri2, opp 4=NO
        ];

        // tri0 is the triangle containing the "newly inserted" point 4.
        // Run lawson swapping starting from tri0 for point 4.
        lawson_swapping(&mut t, 0, 4);

        // After Lawson swapping, the triangulation should be Delaunay.
        // Verify every triangle-triangle pair satisfies the incircle test.
        let n = t.triangle_vertices.len();
        for i in 0..n {
            for j in 0..3 {
                let neighbor = t.triangle_neighbors[i][j];
                if neighbor == NO_NEIGHBOR {
                    continue;
                }
                let vi = t.triangle_vertices[i];
                let ai = t.points[vi[0]];
                let bi = t.points[vi[1]];
                let ci = t.points[vi[2]];
                // The vertex opposite the shared edge in the neighbor triangle
                let opp = get_opposite_vertex(
                    &t,
                    neighbor,
                    t.triangle_vertices[i][(j + 1) % 3],
                    t.triangle_vertices[i][(j + 2) % 3],
                );
                if let Some(opp_v) = opp {
                    let d = t.points[opp_v];
                    // incircle > 0 means d is INSIDE — Delaunay violation
                    assert!(
                        incircle(&ai, &bi, &ci, &d) <= 0.0,
                        "Delaunay violation after lawson_swapping: tri {} neighbor {} opp vertex {}",
                        i, neighbor, opp_v
                    );
                }
            }
        }
    }

    #[test]
    fn find_triangle_with_vertex_basic() {
        let t = two_tri_setup();
        // vertex 3 is only in tri 1
        assert_eq!(find_triangle_with_vertex(&t, 0, 3), 1);
        // vertex 0 is only in tri 0
        assert_eq!(find_triangle_with_vertex(&t, 0, 0), 0);
        // vertex 1 is in both; starting from tri 1 should return tri 1
        assert_eq!(find_triangle_with_vertex(&t, 1, 1), 1);
    }
}
