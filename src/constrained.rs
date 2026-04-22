use std::collections::VecDeque;

use crate::geometry::{incircle, orient2d};
use crate::topology::{find_shared_edge, swap_diagonal};
use crate::triangulation::Triangulation;
use crate::types::{Point, NO_NEIGHBOR};

/// Test if segments (p1,p2) and (p3,p4) properly intersect (not touching at endpoints).
/// Port from Python constrained.py:segments_intersect
///
/// Returns true only for proper interior crossings — collinear or endpoint-touching
/// cases return false (matching Python behavior where o==0 means not a proper intersection).
pub fn segments_intersect(p1: &Point, p2: &Point, p3: &Point, p4: &Point) -> bool {
    let o1 = orient2d(p3, p4, p1);
    let o2 = orient2d(p3, p4, p2);
    let o3 = orient2d(p1, p2, p3);
    let o4 = orient2d(p1, p2, p4);

    // If any are exactly collinear, this is not a proper intersection
    if o1 == 0.0 || o2 == 0.0 || o3 == 0.0 || o4 == 0.0 {
        return false;
    }

    // General case: segments intersect if orientations differ on both sides
    if o1 * o2 < 0.0 && o3 * o4 < 0.0 {
        return true;
    }

    false
}

/// Check if the quadrilateral formed by two triangles sharing an edge is strictly convex.
/// tri_a and tri_b must share exactly one edge.
/// Port from Python constrained.py:is_quadrilateral_convex (4-vertex version).
pub fn is_quadrilateral_convex(t: &Triangulation, tri_a: usize, tri_b: usize) -> bool {
    let va = t.triangle_vertices[tri_a];
    let vb = t.triangle_vertices[tri_b];

    // Find shared edge
    let shared: Vec<usize> = va.iter().copied().filter(|v| vb.contains(v)).collect();
    if shared.len() != 2 {
        return false;
    }
    let (vk, vl) = (shared[0], shared[1]);
    let vm = va
        .iter()
        .copied()
        .find(|&v| v != vk && v != vl)
        .expect("invariant: triangle A must have one opposite vertex");
    let vn = vb
        .iter()
        .copied()
        .find(|&v| v != vk && v != vl)
        .expect("invariant: triangle B must have one opposite vertex");

    let pk = &t.points[vk];
    let pl = &t.points[vl];
    let pm = &t.points[vm];
    let pn = &t.points[vn];

    // vm and vn must be on opposite sides of edge vk-vl
    let o_vm = orient2d(pk, pl, pm);
    let o_vn = orient2d(pk, pl, pn);
    if o_vm * o_vn >= 0.0 {
        return false;
    }

    // vk and vl must be on opposite sides of edge vm-vn
    let o_vk = orient2d(pm, pn, pk);
    let o_vl = orient2d(pm, pn, pl);
    if o_vk * o_vl >= 0.0 {
        return false;
    }

    true
}

/// Walk the triangulation from v1 toward v2, collecting all triangle edges that
/// properly intersect segment (v1,v2). Returns Vec<(tri_idx, neighbor_tri_idx)> pairs.
/// Port from Python constrained.py:find_intersecting_edges
pub fn find_intersecting_edges(t: &Triangulation, v1: usize, v2: usize) -> Vec<(usize, usize)> {
    let p = t.points[v1];
    let q = t.points[v2];

    // Find all triangles containing v1
    let tris_with_v1: Vec<usize> = t
        .triangle_vertices
        .iter()
        .enumerate()
        .filter_map(|(i, verts)| if verts.contains(&v1) { Some(i) } else { None })
        .collect();

    let tris_with_v2: Vec<usize> = t
        .triangle_vertices
        .iter()
        .enumerate()
        .filter_map(|(i, verts)| if verts.contains(&v2) { Some(i) } else { None })
        .collect();

    // If any triangle contains both v1 and v2, the edge is already in the triangulation
    for &tri in &tris_with_v1 {
        if tris_with_v2.contains(&tri) {
            return vec![];
        }
    }

    let tp = tris_with_v1[0];

    let mut intersecting: Vec<(usize, usize)> = Vec::new();
    let mut current = tp;
    let mut visited: std::collections::HashSet<usize> = std::collections::HashSet::new();
    visited.insert(current);

    let max_iterations = t.triangle_vertices.len() + 1;

    for _iteration in 0..max_iterations {
        // Check if we've reached a triangle containing v2
        if t.triangle_vertices[current].contains(&v2) {
            break;
        }

        // Check if q is inside/on the current triangle
        let verts = t.triangle_vertices[current];
        let pa = &t.points[verts[0]];
        let pb = &t.points[verts[1]];
        let pc = &t.points[verts[2]];
        use crate::geometry::point_in_triangle;
        use crate::geometry::PointInTriangle;
        let loc = point_in_triangle(&q, pa, pb, pc);
        if loc != PointInTriangle::Outside {
            break;
        }

        // Try Case C: proper crossing
        let result = check_proper_crossing(t, current, &p, &q, v1, v2, &visited);
        if let Some((next_tri, Some(edge))) = result {
            intersecting.push(edge);
            current = next_tri;
            visited.insert(current);
            continue;
        }
        if let Some((next_tri, None)) = result {
            current = next_tri;
            visited.insert(current);
            continue;
        }

        // Try Case A: collinear overlap
        let result = check_collinear_overlap(t, current, &p, &q, &visited);
        if let Some(next_tri) = result {
            current = next_tri;
            visited.insert(current);
            continue;
        }

        // Try Case B: one endpoint collinear
        let result = check_one_endpoint_collinear(t, current, &p, &q, &visited);
        if let Some(next_tri) = result {
            current = next_tri;
            visited.insert(current);
            continue;
        }

        panic!(
            "find_intersecting_edges: failed to advance from triangle {} toward v2={}",
            current, v2
        );
    }

    intersecting
}

/// Case C: check for proper crossing in triangle current_tri
/// Returns Some((next_tri, Option<edge>)) or None if not found
fn check_proper_crossing(
    t: &Triangulation,
    current_tri: usize,
    p: &Point,
    q: &Point,
    _v1: usize,
    _v2: usize,
    visited: &std::collections::HashSet<usize>,
) -> Option<(usize, Option<(usize, usize)>)> {
    let tri_verts = t.triangle_vertices[current_tri];
    let a = &t.points[tri_verts[0]];
    let b = &t.points[tri_verts[1]];
    let c = &t.points[tri_verts[2]];

    // edges: (start, end, opposite_vertex_local_idx, v_start_local, v_end_local)
    let edges: [(&Point, &Point, usize, usize, usize); 3] = [
        (a, b, 2, 0, 1), // edge v0-v1, opposite v2
        (b, c, 0, 1, 2), // edge v1-v2, opposite v0
        (c, a, 1, 2, 0), // edge v2-v0, opposite v1
    ];

    for (edge_start, edge_end, opp_local, _v_start_local, _v_end_local) in edges {
        let pqs_orient = orient2d(p, q, edge_start);
        let pqe_orient = orient2d(p, q, edge_end);

        // Skip if any endpoint is collinear (handled by other cases)
        if pqs_orient == 0.0 || pqe_orient == 0.0 {
            continue;
        }

        if !segments_intersect(p, q, edge_start, edge_end) {
            continue;
        }

        // Check if q is on opposite side from the opposite vertex
        let opp_point = &t.points[tri_verts[opp_local]];
        let o_q = orient2d(edge_start, edge_end, q);
        let o_opp = orient2d(edge_start, edge_end, opp_point);

        if o_q * o_opp < 0.0 {
            let neighbor_idx = t.triangle_neighbors[current_tri][opp_local];
            if neighbor_idx == NO_NEIGHBOR {
                panic!("find_intersecting_edges: ended up outside triangulation");
            }
            if visited.contains(&neighbor_idx) {
                continue;
            }

            return Some((neighbor_idx, Some((current_tri, neighbor_idx))));
        }
    }

    None
}

/// Case A: both edge endpoints are collinear with pq → walk across without recording
fn check_collinear_overlap(
    t: &Triangulation,
    current_tri: usize,
    p: &Point,
    q: &Point,
    visited: &std::collections::HashSet<usize>,
) -> Option<usize> {
    let tri_verts = t.triangle_vertices[current_tri];
    let pts = [
        &t.points[tri_verts[0]],
        &t.points[tri_verts[1]],
        &t.points[tri_verts[2]],
    ];

    // edges: (local_start, local_end, opposite_local)
    let edges = [(0usize, 1usize, 2usize), (1, 2, 0), (2, 0, 1)];

    for (ls, le, opp_local) in edges {
        let pqs = orient2d(p, q, pts[ls]);
        let pqe = orient2d(p, q, pts[le]);

        if pqs != 0.0 || pqe != 0.0 {
            continue;
        }

        // Both endpoints collinear — check overlap
        if !collinear_overlap(p, q, pts[ls], pts[le]) {
            continue;
        }

        // Find an unvisited neighbor sharing this edge
        let neighbor = t.triangle_neighbors[current_tri][opp_local];
        if neighbor != NO_NEIGHBOR && !visited.contains(&neighbor) {
            return Some(neighbor);
        }
        // Also check other neighbors that share a vertex
        for i in 0..3 {
            let nb = t.triangle_neighbors[current_tri][i];
            if nb == NO_NEIGHBOR || visited.contains(&nb) {
                continue;
            }
            let nb_verts = t.triangle_vertices[nb];
            if nb_verts.contains(&tri_verts[ls]) || nb_verts.contains(&tri_verts[le]) {
                return Some(nb);
            }
        }
    }

    None
}

/// Case B: exactly one endpoint collinear with pq → walk to adjacent triangle through that vertex
fn check_one_endpoint_collinear(
    t: &Triangulation,
    current_tri: usize,
    p: &Point,
    q: &Point,
    visited: &std::collections::HashSet<usize>,
) -> Option<usize> {
    let tri_verts = t.triangle_vertices[current_tri];
    let pts = [
        &t.points[tri_verts[0]],
        &t.points[tri_verts[1]],
        &t.points[tri_verts[2]],
    ];

    let edges = [(0usize, 1usize), (1, 2), (2, 0)];

    for (ls, le) in edges {
        let pqs = orient2d(p, q, pts[ls]);
        let pqe = orient2d(p, q, pts[le]);

        // Only handle exactly one endpoint collinear
        if (pqs == 0.0) == (pqe == 0.0) {
            continue;
        }

        let (collinear_pt, collinear_v) = if pqs == 0.0 {
            (pts[ls], tri_verts[ls])
        } else {
            (pts[le], tri_verts[le])
        };

        // Check if the collinear point is within segment pq bounding box
        if !point_in_segment_bbox(p, q, collinear_pt) {
            continue;
        }

        // Walk into an unvisited neighbor containing this vertex
        for i in 0..3 {
            let nb = t.triangle_neighbors[current_tri][i];
            if nb == NO_NEIGHBOR || visited.contains(&nb) {
                continue;
            }
            if t.triangle_vertices[nb].contains(&collinear_v) {
                return Some(nb);
            }
        }
    }

    None
}

/// Check if segment [a,b] and [c,d] overlap given they are all collinear.
fn collinear_overlap(p: &Point, q: &Point, a: &Point, b: &Point) -> bool {
    point_in_segment_bbox(p, q, a)
        || point_in_segment_bbox(p, q, b)
        || point_in_segment_bbox(a, b, p)
        || point_in_segment_bbox(a, b, q)
}

/// Check if point p is in the bounding box of segment [a, b].
fn point_in_segment_bbox(a: &Point, b: &Point, p: &Point) -> bool {
    use crate::types::EPS;
    let min_x = a[0].min(b[0]);
    let max_x = a[0].max(b[0]);
    let min_y = a[1].min(b[1]);
    let max_y = a[1].max(b[1]);
    p[0] >= min_x - EPS && p[0] <= max_x + EPS && p[1] >= min_y - EPS && p[1] <= max_y + EPS
}

fn collect_intersecting_edge_pairs(t: &Triangulation, v1: usize, v2: usize) -> Vec<(usize, usize)> {
    let p = t.points[v1];
    let q = t.points[v2];
    let mut result = Vec::new();

    for tri_a in 0..t.triangle_vertices.len() {
        for &tri_b in &t.triangle_neighbors[tri_a] {
            if tri_b == NO_NEIGHBOR || tri_a >= tri_b {
                continue;
            }

            let Some((edge_v1, edge_v2)) = find_shared_edge(t, tri_a, tri_b) else {
                continue;
            };

            if segments_intersect(&p, &q, &t.points[edge_v1], &t.points[edge_v2]) {
                result.push((tri_a, tri_b));
            }
        }
    }

    result
}

pub fn remove_intersecting_edges(
    t: &mut Triangulation,
    v1: usize,
    v2: usize,
    edges: Vec<(usize, usize)>,
) -> Vec<(usize, usize)> {
    if edges.is_empty() {
        return vec![];
    }

    let p = t.points[v1];
    let q = t.points[v2];
    let constraint_edge = Triangulation::edge_key(v1, v2);
    let mut newly_created = Vec::new();
    let max_iterations = edges.len().max(1) * 10;

    for _ in 0..max_iterations {
        let mut intersecting: VecDeque<(usize, usize)> =
            collect_intersecting_edge_pairs(t, v1, v2).into();
        let Some((tri_a, tri_b)) = intersecting.pop_front() else {
            break;
        };

        let mut candidate = Some((tri_a, tri_b));
        while let Some((cand_a, cand_b)) = candidate {
            if is_quadrilateral_convex(t, cand_a, cand_b) {
                swap_diagonal(t, cand_a, cand_b);

                let Some((new_v1, new_v2)) = find_shared_edge(t, cand_a, cand_b) else {
                    panic!(
                        "remove_intersecting_edges: swapped triangles {} and {} no longer share an edge",
                        cand_a, cand_b
                    );
                };

                let new_edge = Triangulation::edge_key(new_v1, new_v2);
                if new_edge == constraint_edge {
                    return newly_created;
                }

                if !segments_intersect(&p, &q, &t.points[new_edge.0], &t.points[new_edge.1]) {
                    newly_created.push(new_edge);
                }

                candidate = None;
                break;
            }

            candidate = intersecting.pop_front();
        }

        if candidate.is_some() {
            panic!("remove_intersecting_edges: failed to find a convex intersecting edge to swap");
        }
    }

    assert!(
        collect_intersecting_edge_pairs(t, v1, v2).is_empty(),
        "remove_intersecting_edges: failed to remove all intersecting edges within safety limit"
    );

    newly_created
}

pub fn find_triangles_sharing_edge(t: &Triangulation, v1: usize, v2: usize) -> (usize, usize) {
    let mut first = NO_NEIGHBOR;
    let mut second = NO_NEIGHBOR;

    for (tri_idx, tri_verts) in t.triangle_vertices.iter().enumerate() {
        if tri_verts.contains(&v1) && tri_verts.contains(&v2) {
            if first == NO_NEIGHBOR {
                first = tri_idx;
            } else {
                second = tri_idx;
                break;
            }
        }
    }

    (first, second)
}

fn restore_delaunay_edges(
    t: &mut Triangulation,
    mut edges: Vec<(usize, usize)>,
    constraint_edge: (usize, usize),
) -> bool {
    if edges.is_empty() {
        return true;
    }

    let max_iterations = edges.len().max(1) * 10;

    for _ in 0..max_iterations {
        let mut swapped = false;
        let mut next_edges = Vec::with_capacity(edges.len());

        for (vk, vl) in edges {
            let edge = Triangulation::edge_key(vk, vl);
            if edge == constraint_edge {
                next_edges.push(edge);
                continue;
            }

            let (tri_a, tri_b) = find_triangles_sharing_edge(t, edge.0, edge.1);
            if tri_a == NO_NEIGHBOR || tri_b == NO_NEIGHBOR {
                next_edges.push(edge);
                continue;
            }

            if !is_quadrilateral_convex(t, tri_a, tri_b) {
                next_edges.push(edge);
                continue;
            }

            let tri_a_verts = t.triangle_vertices[tri_a];
            let tri_b_verts = t.triangle_vertices[tri_b];
            let vn = tri_b_verts
                .iter()
                .copied()
                .find(|&vertex| vertex != edge.0 && vertex != edge.1)
                .expect("neighbor triangle must have one opposite vertex");

            let a = t.points[tri_a_verts[0]];
            let b = t.points[tri_a_verts[1]];
            let c = t.points[tri_a_verts[2]];
            let d = t.points[vn];

            if incircle(&a, &b, &c, &d) > 0.0 {
                swap_diagonal(t, tri_a, tri_b);
                let Some((new_v1, new_v2)) = find_shared_edge(t, tri_a, tri_b) else {
                    panic!(
                        "restore_delaunay_edges: swapped triangles {} and {} no longer share an edge",
                        tri_a, tri_b
                    );
                };
                next_edges.push(Triangulation::edge_key(new_v1, new_v2));
                swapped = true;
            } else {
                next_edges.push(edge);
            }
        }

        edges = next_edges;
        if !swapped {
            return true;
        }
    }

    false
}

pub fn add_constraints(t: &mut Triangulation, constraints: &[(usize, usize)]) -> bool {
    for &(v1, v2) in constraints {
        let constraint_edge = Triangulation::edge_key(v1, v2);
        let intersecting = collect_intersecting_edge_pairs(t, v1, v2);

        let newly_created = if intersecting.is_empty() {
            Vec::new()
        } else {
            remove_intersecting_edges(t, v1, v2, intersecting)
        };

        t.constrained_edges.insert(constraint_edge);

        if !restore_delaunay_edges(t, newly_created, constraint_edge) {
            return false;
        }
    }

    true
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::build::triangulate;
    use crate::geometry::incircle;
    use crate::triangulation::Triangulation;
    use crate::types::NO_NEIGHBOR;
    use std::collections::HashSet;

    #[test]
    fn segments_intersect_crossing() {
        // Two crossing segments: (0,0)-(1,1) and (0,1)-(1,0)
        let p1 = [0.0, 0.0];
        let p2 = [1.0, 1.0];
        let p3 = [0.0, 1.0];
        let p4 = [1.0, 0.0];
        assert!(segments_intersect(&p1, &p2, &p3, &p4));
    }

    #[test]
    fn segments_intersect_parallel() {
        // Two parallel horizontal segments
        let p1 = [0.0, 0.0];
        let p2 = [1.0, 0.0];
        let p3 = [0.0, 1.0];
        let p4 = [1.0, 1.0];
        assert!(!segments_intersect(&p1, &p2, &p3, &p4));
    }

    #[test]
    fn segments_intersect_t_junction() {
        // T-junction: endpoint of one segment on the other (not a proper interior crossing)
        let p1 = [0.0, 0.0];
        let p2 = [1.0, 0.0];
        let p3 = [0.5, 0.0]; // on segment p1-p2
        let p4 = [0.5, 1.0];
        // p3 is collinear with p1-p2, so orient2d returns 0 → not a proper intersection
        assert!(!segments_intersect(&p1, &p2, &p3, &p4));
    }

    #[test]
    fn segments_intersect_touching_endpoint() {
        // Endpoint touching: p2 == p3
        let p1 = [0.0, 0.0];
        let p2 = [1.0, 1.0];
        let p3 = [1.0, 1.0];
        let p4 = [2.0, 0.0];
        // p2 == p3, so orient2d(p1,p2,p3) == 0 → not a proper intersection
        assert!(!segments_intersect(&p1, &p2, &p3, &p4));
    }

    #[test]
    fn segments_intersect_non_crossing() {
        // Segments that don't cross and don't touch
        let p1 = [0.0, 0.0];
        let p2 = [1.0, 0.0];
        let p3 = [2.0, 1.0];
        let p4 = [3.0, 2.0];
        assert!(!segments_intersect(&p1, &p2, &p3, &p4));
    }

    fn two_tri_quad() -> Triangulation {
        // Square quad split by diagonal (0,2):
        //   0=(0,0)  1=(1,0)  2=(1,1)  3=(0,1)
        //   tri0 = [0,1,2] CCW
        //   tri1 = [0,2,3] CCW
        //   Shared edge: 0-2; opposite in tri0=1(local 2), opposite in tri1=3(local 2? let's check)
        //   tri0=[0,1,2]: neighbors opp 0=? opp 1=? opp 2=tri1 (edge 0-1 opp to 2)
        //   Wait: neighbor[i] is opposite to vertex[i]
        //   tri0=[0,1,2]: edge 1-2 is opposite to 0 → nb[0]; edge 0-2 opposite to 1 → nb[1]; edge 0-1 opposite to 2 → nb[2]
        //   Shared edge is 0-2, which is edge 0-2, opposite to vertex 1 in tri0 → nb[1] = tri1
        //   tri1=[0,2,3]: edge 2-3 opposite to 0 → nb[0]; edge 0-3 opposite to 2 → nb[1]; edge 0-2 opposite to 3 → nb[2] = tri0
        let mut t = Triangulation::new();
        t.points = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        t.triangle_vertices = vec![[0, 1, 2], [0, 2, 3]];
        t.triangle_neighbors = vec![
            [NO_NEIGHBOR, 1, NO_NEIGHBOR], // tri0: opp 0=NO, opp 1=tri1, opp 2=NO
            [NO_NEIGHBOR, NO_NEIGHBOR, 0], // tri1: opp 0=NO, opp 2=NO, opp 3=tri0
        ];
        t
    }

    fn assert_neighbors_consistent(t: &Triangulation) {
        for (tri_idx, neighbors) in t.triangle_neighbors.iter().enumerate() {
            for &neighbor in neighbors {
                if neighbor != NO_NEIGHBOR {
                    assert!(
                        t.triangle_neighbors[neighbor].contains(&tri_idx),
                        "triangle {} references neighbor {} but not vice versa",
                        tri_idx,
                        neighbor
                    );
                }
            }
        }
    }

    fn assert_unconstrained_edges_delaunay(t: &Triangulation) {
        let mut visited_edges = HashSet::new();

        for (tri_idx, neighbors) in t.triangle_neighbors.iter().enumerate() {
            let tri_verts = t.triangle_vertices[tri_idx];
            let a = t.points[tri_verts[0]];
            let b = t.points[tri_verts[1]];
            let c = t.points[tri_verts[2]];

            for &neighbor in neighbors {
                if neighbor == NO_NEIGHBOR {
                    continue;
                }

                let (shared_a, shared_b) = find_shared_edge(t, tri_idx, neighbor)
                    .expect("neighboring triangles must share an edge");
                let edge = Triangulation::edge_key(shared_a, shared_b);

                if !visited_edges.insert(edge) {
                    continue;
                }

                if t.constrained_edges.contains(&edge) {
                    continue;
                }

                let opposite = t.triangle_vertices[neighbor]
                    .iter()
                    .copied()
                    .find(|vertex| *vertex != edge.0 && *vertex != edge.1)
                    .expect("neighbor triangle must have opposite vertex");

                assert!(
                    incircle(&a, &b, &c, &t.points[opposite]) <= 0.0,
                    "unconstrained edge {:?} violates Delaunay condition",
                    edge
                );
            }
        }
    }

    fn point_index(t: &Triangulation, point: [f64; 2]) -> usize {
        t.points
            .iter()
            .position(|&candidate| candidate == point)
            .expect("point must exist in triangulation")
    }

    #[test]
    fn is_quadrilateral_convex_square_quad() {
        let t = two_tri_quad();
        // The square is convex
        assert!(is_quadrilateral_convex(&t, 0, 1));
    }

    #[test]
    fn is_quadrilateral_convex_degenerate() {
        // Non-convex: make a concave quad
        let mut t = Triangulation::new();
        t.points = vec![
            [0.0, 0.0], // 0
            [2.0, 0.0], // 1
            [1.0, 0.5], // 2 — inside the square, making it concave
            [1.0, 2.0], // 3
        ];
        t.triangle_vertices = vec![[0, 1, 2], [0, 2, 3]];
        t.triangle_neighbors = vec![[NO_NEIGHBOR, 1, NO_NEIGHBOR], [NO_NEIGHBOR, NO_NEIGHBOR, 0]];
        // 2 is inside the triangle 0-1-3, so the quad is concave
        // is_quadrilateral_convex should return false
        assert!(!is_quadrilateral_convex(&t, 0, 1));
    }

    #[test]
    fn find_intersecting_edges_simple_grid() {
        // 4-point square triangulation:
        //   0=(0,0)  1=(1,0)  2=(1,1)  3=(0,1)
        //   tri0=[0,1,2], tri1=[0,2,3]
        //   Constraint from vertex 1 to vertex 3 — should cross edge 0-2 (shared diagonal)
        let t = two_tri_quad();
        let edges = find_intersecting_edges(&t, 1, 3);
        assert_eq!(edges.len(), 1, "Should find exactly 1 intersecting edge");
        // The edge returned should be (tri0, tri1) or (tri1, tri0)
        let (a, b) = edges[0];
        let pair = (a.min(b), a.max(b));
        assert_eq!(
            pair,
            (0, 1),
            "Intersecting edge should be between tri0 and tri1"
        );
    }

    #[test]
    fn find_intersecting_edges_already_in_triangulation() {
        // Constraint from 0 to 2 — already the shared edge of both triangles
        let t = two_tri_quad();
        let edges = find_intersecting_edges(&t, 0, 2);
        assert_eq!(
            edges.len(),
            0,
            "Edge already in triangulation → no intersections"
        );
    }

    #[test]
    fn constraint_edge_already_present() {
        let points = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let mut t = triangulate(&points);
        let shared =
            find_shared_edge(&t, 0, 1).expect("square triangulation should have one shared edge");
        let existing_edge = Triangulation::edge_key(shared.0, shared.1);

        assert!(add_constraints(&mut t, &[existing_edge]));
        assert!(t.constrained_edges.contains(&existing_edge));
        assert_neighbors_consistent(&t);
    }

    #[test]
    fn single_constraint_crosses_triangulation() {
        let points = [
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.1],
            [2.0, 1.0],
            [0.0, 2.0],
            [1.0, 2.0],
            [2.0, 2.0],
        ];
        let mut t = triangulate(&points);
        let v1 = point_index(&t, [0.0, 0.0]);
        let v2 = point_index(&t, [2.0, 2.0]);

        assert!(add_constraints(&mut t, &[(v1, v2)]));
        assert!(t
            .constrained_edges
            .contains(&Triangulation::edge_key(v1, v2)));

        let (tri_a, tri_b) = find_triangles_sharing_edge(&t, v1, v2);
        assert_ne!(
            tri_a, NO_NEIGHBOR,
            "constraint edge should exist in triangulation"
        );
        assert_ne!(
            tri_b, NO_NEIGHBOR,
            "constraint edge should be internal on this fixture"
        );
        assert_neighbors_consistent(&t);
    }

    #[test]
    fn multiple_constraints_are_tracked() {
        let points = [
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.1],
            [2.0, 1.0],
            [0.0, 2.0],
            [1.0, 2.0],
            [2.0, 2.0],
        ];
        let mut t = triangulate(&points);
        let v00 = point_index(&t, [0.0, 0.0]);
        let v21 = point_index(&t, [2.0, 1.0]);
        let v02 = point_index(&t, [0.0, 2.0]);
        let constraints = [(v00, v21), (v21, v02)];

        assert!(add_constraints(&mut t, &constraints));
        for edge in constraints {
            assert!(t
                .constrained_edges
                .contains(&Triangulation::edge_key(edge.0, edge.1)));
        }
        assert_neighbors_consistent(&t);
    }

    #[test]
    fn delaunay_restoration_after_constraint() {
        let points = [
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.1],
            [2.0, 1.0],
            [0.0, 2.0],
            [1.0, 2.0],
            [2.0, 2.0],
        ];
        let mut t = triangulate(&points);
        let v1 = point_index(&t, [0.0, 0.0]);
        let v2 = point_index(&t, [2.0, 2.0]);

        assert!(add_constraints(&mut t, &[(v1, v2)]));
        assert_neighbors_consistent(&t);
        assert_unconstrained_edges_delaunay(&t);
    }

    #[test]
    fn constraint_neighbors_consistent() {
        let points = [
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.1],
            [2.0, 1.0],
            [0.0, 2.0],
            [1.0, 2.0],
            [2.0, 2.0],
        ];
        let mut t = triangulate(&points);
        let v1 = point_index(&t, [0.0, 0.0]);
        let v2 = point_index(&t, [2.0, 2.0]);

        assert!(add_constraints(&mut t, &[(v1, v2)]));
        assert_neighbors_consistent(&t);
    }
}
