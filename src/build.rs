use crate::geometry::{ensure_ccw, is_point_inside_polygon, point_in_triangle, PointInTriangle};
use crate::topology::{lawson_swapping, reorder_neighbors};
use crate::triangulation::Triangulation;
use crate::types::{Point, PointLocation, NO_NEIGHBOR};
use std::collections::{HashMap, HashSet};

fn initialize_triangulation(points: &[[f64; 2]]) -> Triangulation {
    assert!(
        !points.is_empty(),
        "initialize_triangulation requires at least one point"
    );

    let mut min_x = points[0][0];
    let mut min_y = points[0][1];
    let mut max_x = points[0][0];
    let mut max_y = points[0][1];

    for point in &points[1..] {
        min_x = min_x.min(point[0]);
        min_y = min_y.min(point[1]);
        max_x = max_x.max(point[0]);
        max_y = max_y.max(point[1]);
    }

    let dx = max_x - min_x;
    let dy = max_y - min_y;
    let delta_max = dx.max(dy);
    let mid_x = (min_x + max_x) / 2.0;
    let mid_y = (min_y + max_y) / 2.0;

    let p0 = [mid_x - 100.0 * delta_max, mid_y - delta_max];
    let p1 = [mid_x, mid_y + 100.0 * delta_max];
    let p2 = [mid_x + 100.0 * delta_max, mid_y - delta_max];

    Triangulation {
        points: vec![p0, p1, p2],
        triangle_vertices: vec![[0, 1, 2]],
        triangle_neighbors: vec![[NO_NEIGHBOR, NO_NEIGHBOR, NO_NEIGHBOR]],
        constrained_edges: Default::default(),
        num_super_triangle_points: 3,
    }
}

pub fn remove_super_triangle(t: &mut Triangulation) {
    let num_super = t.num_super_triangle_points;
    if num_super == 0 {
        return;
    }

    let keep: Vec<usize> = t
        .triangle_vertices
        .iter()
        .enumerate()
        .filter_map(|(idx, vertices)| {
            if vertices.iter().all(|&v| v >= num_super) {
                Some(idx)
            } else {
                None
            }
        })
        .collect();

    let mut old_to_new = vec![NO_NEIGHBOR; t.triangle_vertices.len()];
    for (new_idx, &old_idx) in keep.iter().enumerate() {
        old_to_new[old_idx] = new_idx;
    }

    let new_triangle_vertices: Vec<[usize; 3]> = keep
        .iter()
        .map(|&old_idx| {
            let [a, b, c] = t.triangle_vertices[old_idx];
            [a - num_super, b - num_super, c - num_super]
        })
        .collect();

    let new_triangle_neighbors: Vec<[usize; 3]> = keep
        .iter()
        .map(|&old_idx| {
            let neighbors = t.triangle_neighbors[old_idx];
            neighbors.map(|neighbor| {
                if neighbor == NO_NEIGHBOR {
                    NO_NEIGHBOR
                } else {
                    old_to_new[neighbor]
                }
            })
        })
        .collect();

    t.triangle_vertices = new_triangle_vertices;
    t.triangle_neighbors = new_triangle_neighbors;
    t.points.drain(0..num_super);
    t.num_super_triangle_points = 0;
}

pub fn triangulate(input_points: &[[f64; 2]]) -> Triangulation {
    let mut t = initialize_triangulation(input_points);

    for &point in input_points {
        t.points.push(point);
        let point_idx = t.points.len() - 1;
        insert_point(&mut t, point_idx);
    }

    remove_super_triangle(&mut t);
    t
}

pub fn find_containing_triangle(t: &Triangulation, point: &Point) -> PointLocation {
    for (tri_idx, vertices) in t.triangle_vertices.iter().enumerate() {
        let [v0, v1, v2] = *vertices;
        let position = point_in_triangle(point, &t.points[v0], &t.points[v1], &t.points[v2]);

        match position {
            PointInTriangle::Inside => return PointLocation::Interior(tri_idx),
            PointInTriangle::OnEdge0 => return PointLocation::OnEdge(tri_idx, 0),
            PointInTriangle::OnEdge1 => return PointLocation::OnEdge(tri_idx, 1),
            PointInTriangle::OnEdge2 => return PointLocation::OnEdge(tri_idx, 2),
            PointInTriangle::Outside => continue,
        }
    }

    PointLocation::NotFound
}

fn reorder_triangle_neighbors(
    original_vertices: [usize; 3],
    final_vertices: [usize; 3],
    original_neighbors: [usize; 3],
) -> [usize; 3] {
    let mut final_neighbors = [NO_NEIGHBOR; 3];

    for (i, vertex) in final_vertices.iter().enumerate() {
        let original_idx = original_vertices
            .iter()
            .position(|&v| v == *vertex)
            .expect("vertex must exist in original triangle ordering");
        final_neighbors[i] = original_neighbors[original_idx];
    }

    final_neighbors
}

fn update_external_neighbor(
    t: &mut Triangulation,
    neighbor_idx: usize,
    old_idx: usize,
    new_idx: usize,
) {
    if neighbor_idx == NO_NEIGHBOR || old_idx == new_idx {
        return;
    }

    reorder_neighbors(t, neighbor_idx, old_idx, new_idx);
}

fn insert_point_interior(t: &mut Triangulation, tri_idx: usize, point_idx: usize) {
    let [v0, v1, v2] = t.triangle_vertices[tri_idx];
    let [n0, n1, n2] = t.triangle_neighbors[tri_idx];

    let tri_b = t.triangle_vertices.len();
    let tri_c = tri_b + 1;

    let tri_a_vertices = {
        let (a, b, c) = ensure_ccw(&t.points, point_idx, v0, v1);
        [a, b, c]
    };
    let tri_b_vertices = {
        let (a, b, c) = ensure_ccw(&t.points, point_idx, v1, v2);
        [a, b, c]
    };
    let tri_c_vertices = {
        let (a, b, c) = ensure_ccw(&t.points, point_idx, v2, v0);
        [a, b, c]
    };

    t.triangle_vertices[tri_idx] = tri_a_vertices;
    t.triangle_vertices.push(tri_b_vertices);
    t.triangle_vertices.push(tri_c_vertices);

    t.triangle_neighbors[tri_idx] =
        reorder_triangle_neighbors([point_idx, v0, v1], tri_a_vertices, [n2, tri_b, tri_c]);
    t.triangle_neighbors.push(reorder_triangle_neighbors(
        [point_idx, v1, v2],
        tri_b_vertices,
        [n0, tri_c, tri_idx],
    ));
    t.triangle_neighbors.push(reorder_triangle_neighbors(
        [point_idx, v2, v0],
        tri_c_vertices,
        [n1, tri_idx, tri_b],
    ));

    update_external_neighbor(t, n0, tri_idx, tri_b);
    update_external_neighbor(t, n1, tri_idx, tri_c);
    update_external_neighbor(t, n2, tri_idx, tri_idx);

    lawson_swapping(t, tri_idx, point_idx);
    lawson_swapping(t, tri_b, point_idx);
    lawson_swapping(t, tri_c, point_idx);
}

fn insert_point_on_edge(t: &mut Triangulation, tri_idx: usize, edge_idx: usize, point_idx: usize) {
    let containing_vertices = t.triangle_vertices[tri_idx];
    let containing_neighbors = t.triangle_neighbors[tri_idx];

    let (shared_a_pos, shared_b_pos, opposite_pos) = match edge_idx {
        0 => (0, 1, 2),
        1 => (1, 2, 0),
        2 => (2, 0, 1),
        _ => panic!("invalid edge index {edge_idx}"),
    };

    let v1 = containing_vertices[shared_a_pos];
    let v2 = containing_vertices[shared_b_pos];
    let v3 = containing_vertices[opposite_pos];

    let ty = containing_neighbors[shared_a_pos];
    let tx = containing_neighbors[shared_b_pos];
    let adjacent_tri = containing_neighbors[opposite_pos];
    assert_ne!(
        adjacent_tri, NO_NEIGHBOR,
        "edge insertion requires adjacent triangle"
    );

    let adjacent_vertices = t.triangle_vertices[adjacent_tri];
    let adjacent_neighbors = t.triangle_neighbors[adjacent_tri];
    let v4 = adjacent_vertices
        .iter()
        .copied()
        .find(|&v| v != v1 && v != v2)
        .expect("adjacent triangle must contain opposite vertex");

    let v1_pos_in_adjacent = adjacent_vertices
        .iter()
        .position(|&v| v == v1)
        .expect("shared vertex must exist in adjacent triangle");
    let v2_pos_in_adjacent = adjacent_vertices
        .iter()
        .position(|&v| v == v2)
        .expect("shared vertex must exist in adjacent triangle");

    let tw = adjacent_neighbors[v1_pos_in_adjacent];
    let tz = adjacent_neighbors[v2_pos_in_adjacent];

    let tri_b = t.triangle_vertices.len();
    let tri_d = tri_b + 1;

    let tri_a_vertices = {
        let (a, b, c) = ensure_ccw(&t.points, point_idx, v3, v1);
        [a, b, c]
    };
    let tri_b_vertices = {
        let (a, b, c) = ensure_ccw(&t.points, point_idx, v3, v2);
        [a, b, c]
    };
    let tri_c_vertices = {
        let (a, b, c) = ensure_ccw(&t.points, point_idx, v1, v4);
        [a, b, c]
    };
    let tri_d_vertices = {
        let (a, b, c) = ensure_ccw(&t.points, point_idx, v2, v4);
        [a, b, c]
    };

    t.triangle_vertices[tri_idx] = tri_a_vertices;
    t.triangle_vertices[adjacent_tri] = tri_c_vertices;
    t.triangle_vertices.push(tri_b_vertices);
    t.triangle_vertices.push(tri_d_vertices);

    t.triangle_neighbors[tri_idx] = reorder_triangle_neighbors(
        [point_idx, v3, v1],
        tri_a_vertices,
        [tx, adjacent_tri, tri_b],
    );
    t.triangle_neighbors[adjacent_tri] =
        reorder_triangle_neighbors([point_idx, v1, v4], tri_c_vertices, [tz, tri_d, tri_idx]);
    t.triangle_neighbors.push(reorder_triangle_neighbors(
        [point_idx, v3, v2],
        tri_b_vertices,
        [ty, tri_d, tri_idx],
    ));
    t.triangle_neighbors.push(reorder_triangle_neighbors(
        [point_idx, v2, v4],
        tri_d_vertices,
        [tw, adjacent_tri, tri_b],
    ));

    update_external_neighbor(t, tx, tri_idx, tri_idx);
    update_external_neighbor(t, ty, tri_idx, tri_b);
    update_external_neighbor(t, tz, adjacent_tri, adjacent_tri);
    update_external_neighbor(t, tw, adjacent_tri, tri_d);

    lawson_swapping(t, tri_idx, point_idx);
    lawson_swapping(t, tri_b, point_idx);
    lawson_swapping(t, adjacent_tri, point_idx);
    lawson_swapping(t, tri_d, point_idx);
}

/// Incrementally add new_points to an existing triangulation.
/// Does NOT rebuild from scratch — inserts each new point incrementally.
pub fn update_triangulation(t: &mut Triangulation, new_points: &[[f64; 2]]) {
    for &point in new_points {
        t.points.push(point);
        let point_idx = t.points.len() - 1;
        insert_point(t, point_idx);
    }
}

pub fn insert_point(t: &mut Triangulation, point_idx: usize) {
    let point = t.points[point_idx];

    match find_containing_triangle(t, &point) {
        PointLocation::Interior(tri_idx) => insert_point_interior(t, tri_idx, point_idx),
        PointLocation::OnEdge(tri_idx, edge_idx) => {
            insert_point_on_edge(t, tri_idx, edge_idx, point_idx)
        }
        PointLocation::NotFound => panic!("point {} not found in any triangle", point_idx),
    }
}

pub fn build_polygons_from_edges(edges: &[(usize, usize)]) -> Vec<Vec<usize>> {
    let mut adj: HashMap<usize, Vec<usize>> = HashMap::new();
    for &(v1, v2) in edges {
        adj.entry(v1).or_default().push(v2);
        adj.entry(v2).or_default().push(v1);
    }

    let mut visited: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let mut polygons: Vec<Vec<usize>> = Vec::new();

    let mut starts: Vec<usize> = adj.keys().copied().collect();
    starts.sort();

    for start in starts {
        if visited.contains(&start) {
            continue;
        }
        let mut polygon = vec![start];
        let mut current = start;
        let mut prev: Option<usize> = None;
        loop {
            visited.insert(current);
            let neighbors = adj
                .get(&current)
                .expect("invariant: current vertex must exist in adjacency map");
            assert_eq!(
                neighbors.len(),
                2,
                "Vertex {} does not have degree 2 after cleanup",
                current
            );
            let next_v = if Some(neighbors[0]) != prev {
                neighbors[0]
            } else {
                neighbors[1]
            };
            if next_v == start {
                break;
            }
            polygon.push(next_v);
            prev = Some(current);
            current = next_v;
        }
        polygons.push(polygon);
    }

    polygons
}

pub fn polygon_area(points: &[Point], polygon: &[usize]) -> f64 {
    let n = polygon.len();
    let mut area = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        let xi = points[polygon[i]][0];
        let yi = points[polygon[i]][1];
        let xj = points[polygon[j]][0];
        let yj = points[polygon[j]][1];
        area += xi * yj - xj * yi;
    }
    0.5 * area
}

fn extract_planar_faces(points: &[Point], edges: &[(usize, usize)]) -> Vec<Vec<usize>> {
    let mut adjacency: HashMap<usize, HashSet<usize>> = HashMap::new();
    for &(a, b) in edges {
        adjacency.entry(a).or_default().insert(b);
        adjacency.entry(b).or_default().insert(a);
    }

    let mut sorted_neighbors: HashMap<usize, Vec<usize>> = HashMap::new();
    for (&vertex, neighbors) in &adjacency {
        let [px, py] = points[vertex];
        let mut neighbors: Vec<usize> = neighbors.iter().copied().collect();
        neighbors.sort_by(|&lhs, &rhs| {
            let lhs_angle = (points[lhs][1] - py).atan2(points[lhs][0] - px);
            let rhs_angle = (points[rhs][1] - py).atan2(points[rhs][0] - px);
            lhs_angle
                .partial_cmp(&rhs_angle)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        sorted_neighbors.insert(vertex, neighbors);
    }

    let mut next_edge: HashMap<(usize, usize), (usize, usize)> = HashMap::new();
    for (&vertex, neighbors) in &sorted_neighbors {
        let count = neighbors.len();
        for (index, &neighbor) in neighbors.iter().enumerate() {
            let next_neighbor = neighbors[(index + count - 1) % count];
            next_edge.insert((neighbor, vertex), (vertex, next_neighbor));
        }
    }

    let mut visited: HashSet<(usize, usize)> = HashSet::new();
    let mut faces = Vec::new();
    for &(a, b) in edges {
        for start_edge in [(a, b), (b, a)] {
            if visited.contains(&start_edge) {
                continue;
            }

            let mut face = Vec::new();
            let mut edge = start_edge;
            while !visited.contains(&edge) {
                visited.insert(edge);
                face.push(edge.0);
                edge = next_edge[&edge];
            }

            if face.len() >= 3 {
                faces.push(face);
            }
        }
    }

    faces
}

fn sanitize_constrained_edges(
    point_count: usize,
    constrained_edges: &[(usize, usize)],
) -> Vec<(usize, usize)> {
    let mut seen = HashSet::new();
    let mut sanitized = Vec::new();

    for &(a, b) in constrained_edges {
        if a >= point_count || b >= point_count || a == b {
            continue;
        }

        let edge = Triangulation::edge_key(a, b);
        if seen.insert(edge) {
            sanitized.push((a, b));
        }
    }

    sanitized
}

fn rebuild_triangle_neighbors(triangle_vertices: &[[usize; 3]]) -> Vec<[usize; 3]> {
    let mut edge_to_tri: HashMap<(usize, usize), usize> = HashMap::new();
    for (tri_idx, &[a, b, c]) in triangle_vertices.iter().enumerate() {
        edge_to_tri.insert((b, c), tri_idx);
        edge_to_tri.insert((a, c), tri_idx);
        edge_to_tri.insert((a, b), tri_idx);
    }

    triangle_vertices
        .iter()
        .map(|&[a, b, c]| {
            let n0 = edge_to_tri.get(&(c, b)).copied().unwrap_or(NO_NEIGHBOR);
            let n1 = edge_to_tri.get(&(c, a)).copied().unwrap_or(NO_NEIGHBOR);
            let n2 = edge_to_tri.get(&(b, a)).copied().unwrap_or(NO_NEIGHBOR);
            [n0, n1, n2]
        })
        .collect()
}

fn retain_triangles(t: &mut Triangulation, keep: &[usize]) {
    let new_triangle_vertices: Vec<[usize; 3]> = keep
        .iter()
        .map(|&tri_idx| t.triangle_vertices[tri_idx])
        .collect();
    let new_triangle_neighbors = rebuild_triangle_neighbors(&new_triangle_vertices);
    t.triangle_vertices = new_triangle_vertices;
    t.triangle_neighbors = new_triangle_neighbors;
}

fn remove_holes_by_edges_with_constraint_inserter<F>(
    t: &mut Triangulation,
    constrained_edges: &[(usize, usize)],
    add_constraints: F,
) where
    F: FnOnce(&mut Triangulation, &[(usize, usize)]) -> bool,
{
    if constrained_edges.is_empty() {
        return;
    }

    let constrained_edges = sanitize_constrained_edges(t.points.len(), constrained_edges);
    if constrained_edges.is_empty() {
        return;
    }

    let original_triangulation = t.clone();
    let constraint_result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        add_constraints(t, &constrained_edges)
    }));
    if !matches!(constraint_result, Ok(true)) {
        *t = original_triangulation;
        return;
    }

    let mut degree: HashMap<usize, usize> = HashMap::new();
    for &(a, b) in &constrained_edges {
        *degree.entry(a).or_insert(0) += 1;
        *degree.entry(b).or_insert(0) += 1;
    }
    let all_vertices_are_degree_two = degree.values().all(|&vertex_degree| vertex_degree == 2);

    let polygons = if all_vertices_are_degree_two {
        build_polygons_from_edges(&constrained_edges)
    } else {
        extract_planar_faces(&t.points, &constrained_edges)
            .into_iter()
            .filter(|face| polygon_area(&t.points, face) > 0.0)
            .collect()
    };

    if polygons.is_empty() {
        return;
    }

    let outer_idx = polygons
        .iter()
        .enumerate()
        .max_by(|(_, lhs), (_, rhs)| {
            polygon_area(&t.points, lhs)
                .abs()
                .partial_cmp(&polygon_area(&t.points, rhs).abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(idx, _)| idx)
        .expect("polygons is not empty");

    let outer_polygon: Vec<Point> = polygons[outer_idx]
        .iter()
        .map(|&vertex| t.points[vertex])
        .collect();
    let hole_polygons: Vec<Vec<Point>> = polygons
        .iter()
        .enumerate()
        .filter(|(idx, _)| *idx != outer_idx)
        .map(|(_, polygon)| polygon.iter().map(|&vertex| t.points[vertex]).collect())
        .collect();

    let mut to_delete = HashSet::new();
    for (tri_idx, vertices) in t.triangle_vertices.iter().enumerate() {
        let centroid = [
            (t.points[vertices[0]][0] + t.points[vertices[1]][0] + t.points[vertices[2]][0]) / 3.0,
            (t.points[vertices[0]][1] + t.points[vertices[1]][1] + t.points[vertices[2]][1]) / 3.0,
        ];

        if !is_point_inside_polygon(&centroid, &outer_polygon)
            || hole_polygons
                .iter()
                .any(|hole_polygon| is_point_inside_polygon(&centroid, hole_polygon))
        {
            to_delete.insert(tri_idx);
        }
    }

    if to_delete.is_empty() {
        return;
    }

    let keep: Vec<usize> = (0..t.triangle_vertices.len())
        .filter(|tri_idx| !to_delete.contains(tri_idx))
        .collect();
    retain_triangles(t, &keep);
}

pub fn remove_holes_by_edges(t: &mut Triangulation, constrained_edges: &[(usize, usize)]) {
    remove_holes_by_edges_with_constraint_inserter(
        t,
        constrained_edges,
        crate::constrained::add_constraints,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_neighbors_consistent(t: &Triangulation) {
        for (tri_idx, &neighbors) in t.triangle_neighbors.iter().enumerate() {
            let vertices = t.triangle_vertices[tri_idx];
            for (neighbor_slot, &neighbor_idx) in neighbors.iter().enumerate() {
                if neighbor_idx == NO_NEIGHBOR {
                    continue;
                }

                assert!(
                    t.triangle_neighbors[neighbor_idx].contains(&tri_idx),
                    "triangle {} has neighbor {} but {} does not have {} as neighbor",
                    tri_idx,
                    neighbor_idx,
                    neighbor_idx,
                    tri_idx
                );

                let shared_edge = [
                    vertices[(neighbor_slot + 1) % 3],
                    vertices[(neighbor_slot + 2) % 3],
                ];
                let shared_vertices = t.triangle_vertices[neighbor_idx]
                    .iter()
                    .copied()
                    .filter(|vertex| shared_edge.contains(vertex))
                    .count();
                assert_eq!(
                    shared_vertices, 2,
                    "triangle {} and neighbor {} do not share the expected edge",
                    tri_idx, neighbor_idx
                );
            }
        }
    }

    fn triangle_centroid(t: &Triangulation, tri_idx: usize) -> Point {
        let [a, b, c] = t.triangle_vertices[tri_idx];
        [
            (t.points[a][0] + t.points[b][0] + t.points[c][0]) / 3.0,
            (t.points[a][1] + t.points[b][1] + t.points[c][1]) / 3.0,
        ]
    }

    fn two_triangle_triangulation() -> Triangulation {
        Triangulation {
            points: vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
            triangle_vertices: vec![[0, 1, 2], [0, 2, 3]],
            triangle_neighbors: vec![[NO_NEIGHBOR, NO_NEIGHBOR, 1], [0, NO_NEIGHBOR, NO_NEIGHBOR]],
            constrained_edges: Default::default(),
            num_super_triangle_points: 0,
        }
    }

    fn single_triangle_triangulation() -> Triangulation {
        Triangulation {
            points: vec![[0.0, 0.0], [2.0, 0.0], [0.0, 2.0], [0.5, 0.5]],
            triangle_vertices: vec![[0, 1, 2]],
            triangle_neighbors: vec![[NO_NEIGHBOR, NO_NEIGHBOR, NO_NEIGHBOR]],
            constrained_edges: Default::default(),
            num_super_triangle_points: 0,
        }
    }

    fn edge_split_triangulation() -> Triangulation {
        Triangulation {
            points: vec![[0.0, 0.0], [2.0, 0.0], [0.0, 2.0], [2.0, 2.0], [1.0, 1.0]],
            triangle_vertices: vec![[0, 1, 2], [3, 2, 1]],
            triangle_neighbors: vec![[1, NO_NEIGHBOR, NO_NEIGHBOR], [0, NO_NEIGHBOR, NO_NEIGHBOR]],
            constrained_edges: Default::default(),
            num_super_triangle_points: 0,
        }
    }

    fn assert_all_vertex_indices_in_range(t: &Triangulation, num_points: usize) {
        for vertices in &t.triangle_vertices {
            for &vertex in vertices {
                assert!(
                    vertex < num_points,
                    "vertex index {} out of range for {} points",
                    vertex,
                    num_points
                );
            }
        }
    }

    fn assert_delaunay(t: &Triangulation) {
        for (tri_idx, neighbors) in t.triangle_neighbors.iter().enumerate() {
            let [a_idx, b_idx, c_idx] = t.triangle_vertices[tri_idx];
            let a = t.points[a_idx];
            let b = t.points[b_idx];
            let c = t.points[c_idx];

            for (opposite_vertex_idx, &neighbor_idx) in neighbors.iter().enumerate() {
                if neighbor_idx == NO_NEIGHBOR {
                    continue;
                }

                let edge = [
                    t.triangle_vertices[tri_idx][(opposite_vertex_idx + 1) % 3],
                    t.triangle_vertices[tri_idx][(opposite_vertex_idx + 2) % 3],
                ];
                let opposite = t.triangle_vertices[neighbor_idx]
                    .iter()
                    .copied()
                    .find(|vertex| !edge.contains(vertex))
                    .expect("neighbor must have one opposite vertex");

                assert!(
                    crate::geometry::incircle(&a, &b, &c, &t.points[opposite]) <= 0.0,
                    "triangle {} violates Delaunay condition with neighbor {}",
                    tri_idx,
                    neighbor_idx
                );
            }
        }
    }

    #[test]
    fn finds_point_inside_first_triangle() {
        let triangulation = two_triangle_triangulation();
        let point = [0.75, 0.25];

        assert_eq!(
            find_containing_triangle(&triangulation, &point),
            PointLocation::Interior(0)
        );
    }

    #[test]
    fn finds_point_inside_second_triangle() {
        let triangulation = two_triangle_triangulation();
        let point = [0.25, 0.75];

        assert_eq!(
            find_containing_triangle(&triangulation, &point),
            PointLocation::Interior(1)
        );
    }

    #[test]
    fn finds_point_on_shared_edge() {
        let triangulation = two_triangle_triangulation();
        let point = [0.5, 0.5];

        assert_eq!(
            find_containing_triangle(&triangulation, &point),
            PointLocation::OnEdge(0, 2)
        );
    }

    #[test]
    fn returns_not_found_for_point_outside_all_triangles() {
        let triangulation = two_triangle_triangulation();
        let point = [1.5, 0.5];

        assert_eq!(
            find_containing_triangle(&triangulation, &point),
            PointLocation::NotFound
        );
    }

    #[test]
    fn insert_point_interior_splits_one_triangle_into_three() {
        let mut triangulation = single_triangle_triangulation();

        insert_point_interior(&mut triangulation, 0, 3);

        assert_eq!(triangulation.num_triangles(), 3);
        assert_neighbors_consistent(&triangulation);
    }

    #[test]
    fn insert_point_on_edge_splits_two_triangles_into_four() {
        let mut triangulation = edge_split_triangulation();

        insert_point_on_edge(&mut triangulation, 0, 1, 4);

        assert_eq!(triangulation.num_triangles(), 4);
        assert_neighbors_consistent(&triangulation);
    }

    #[test]
    fn triangulate_five_points() {
        let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.5, 0.5]];

        let triangulation = triangulate(&points);

        assert_eq!(triangulation.num_super_triangle_points, 0);
        assert!(triangulation.num_triangles() > 0);
        assert_eq!(triangulation.num_points(), 5);
        assert_all_vertex_indices_in_range(&triangulation, points.len());
        assert_neighbors_consistent(&triangulation);
        assert_delaunay(&triangulation);
    }

    #[test]
    fn triangulate_no_super_triangle_vertices() {
        let points = [
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [2.0, 1.0],
            [0.5, 0.5],
            [1.5, 0.5],
            [0.5, 1.5],
            [1.5, 1.5],
        ];

        let triangulation = triangulate(&points);

        assert_eq!(triangulation.num_super_triangle_points, 0);
        assert_eq!(triangulation.num_points(), points.len());
        assert_all_vertex_indices_in_range(&triangulation, points.len());
    }

    #[test]
    fn update_adds_points() {
        let points = [[0.0, 0.0], [3.0, 0.0], [3.0, 3.0], [0.0, 3.0], [1.5, 1.5]];
        let mut t = triangulate(&points);
        let new_points = [[0.8, 0.8], [2.2, 0.8], [1.5, 2.2]];
        update_triangulation(&mut t, &new_points);
        assert_eq!(t.num_points(), 8);
    }

    #[test]
    fn update_preserves_neighbor_consistency() {
        let points = [[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0], [2.0, 2.0]];
        let mut t = triangulate(&points);
        let new_points = [[1.0, 1.0], [3.0, 1.0], [3.0, 3.0], [1.0, 3.0], [2.0, 1.5]];
        update_triangulation(&mut t, &new_points);
        assert_neighbors_consistent(&t);
    }

    #[test]
    fn update_no_super_triangle_vertices() {
        let points = [[0.0, 0.0], [3.0, 0.0], [3.0, 3.0], [0.0, 3.0], [1.5, 1.5]];
        let mut t = triangulate(&points);
        let new_points = [[0.8, 0.8], [2.2, 0.8], [1.5, 2.2]];
        update_triangulation(&mut t, &new_points);
        let n = t.num_points();
        for &[a, b, c] in &t.triangle_vertices {
            assert!(a < n, "vertex {} out of range", a);
            assert!(b < n, "vertex {} out of range", b);
            assert!(c < n, "vertex {} out of range", c);
        }
    }

    #[test]
    fn polygon_square_edges() {
        let edges = vec![(0usize, 1usize), (1, 2), (2, 3), (3, 0)];
        let polygons = build_polygons_from_edges(&edges);
        assert_eq!(polygons.len(), 1);
        let poly = &polygons[0];
        assert_eq!(poly.len(), 4);
        // All vertices 0..3 present
        let mut sorted = poly.clone();
        sorted.sort();
        assert_eq!(sorted, vec![0, 1, 2, 3]);
    }

    #[test]
    fn polygon_area_ccw() {
        // CCW unit square: (0,0)->(1,0)->(1,1)->(0,1)
        let points: Vec<Point> = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let polygon = vec![0, 1, 2, 3];
        let area = polygon_area(&points, &polygon);
        assert!((area - 1.0).abs() < 1e-12);
    }

    #[test]
    fn polygon_area_cw() {
        // CW unit square: (0,0)->(0,1)->(1,1)->(1,0)
        let points: Vec<Point> = vec![[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0]];
        let polygon = vec![0, 1, 2, 3];
        let area = polygon_area(&points, &polygon);
        assert!((area + 1.0).abs() < 1e-12);
    }

    #[test]
    fn two_separate_polygons() {
        // Triangle 0-1-2 and triangle 3-4-5
        let edges = vec![(0usize, 1usize), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)];
        let polygons = build_polygons_from_edges(&edges);
        assert_eq!(polygons.len(), 2);
        assert_eq!(polygons[0].len(), 3);
        assert_eq!(polygons[1].len(), 3);
    }

    #[test]
    fn triangulate_neighbor_consistency() {
        let points = [
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [2.0, 1.0],
            [0.5, 0.5],
            [1.5, 0.5],
            [0.5, 1.5],
            [1.5, 1.5],
        ];

        let triangulation = triangulate(&points);

        assert_neighbors_consistent(&triangulation);
    }

    #[test]
    fn remove_holes_by_edges_trims_to_simple_outer_polygon() {
        let points = [
            [0.0, 0.0],
            [4.0, 0.0],
            [4.0, 4.0],
            [0.0, 4.0],
            [1.0, 1.0],
            [3.0, 1.0],
            [3.0, 3.0],
            [1.0, 3.0],
            [2.0, 2.0],
        ];
        let mut t = triangulate(&points);
        let before_count = t.num_triangles();

        remove_holes_by_edges(&mut t, &[(4, 5), (5, 6), (6, 7), (7, 4)]);

        assert!(t.num_triangles() < before_count);
        for tri_idx in 0..t.num_triangles() {
            let centroid = triangle_centroid(&t, tri_idx);
            assert!(is_point_inside_polygon(
                &centroid,
                &[[1.0, 1.0], [3.0, 1.0], [3.0, 3.0], [1.0, 3.0]]
            ));
        }
        assert_neighbors_consistent(&t);
    }

    #[test]
    fn remove_holes_by_edges_trims_outer_polygon_with_hole() {
        let points = [
            [0.0, 0.0],
            [6.0, 0.0],
            [6.0, 6.0],
            [0.0, 6.0],
            [2.0, 2.0],
            [4.0, 2.0],
            [4.0, 4.0],
            [2.0, 4.0],
            [3.0, 1.0],
            [5.0, 3.0],
            [3.0, 5.0],
            [1.0, 3.0],
            [3.0, 3.0],
        ];
        let mut t = triangulate(&points);
        let before_count = t.num_triangles();
        let edges = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),
        ];

        remove_holes_by_edges(&mut t, &edges);

        assert!(t.num_triangles() < before_count);
        let outer = [[0.0, 0.0], [6.0, 0.0], [6.0, 6.0], [0.0, 6.0]];
        let hole = [[2.0, 2.0], [4.0, 2.0], [4.0, 4.0], [2.0, 4.0]];
        for tri_idx in 0..t.num_triangles() {
            let centroid = triangle_centroid(&t, tri_idx);
            assert!(is_point_inside_polygon(&centroid, &outer));
            assert!(!is_point_inside_polygon(&centroid, &hole));
        }
        assert_neighbors_consistent(&t);
    }

    #[test]
    fn remove_holes_by_edges_uses_face_extraction_for_branching_edges() {
        let points = [[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0], [2.0, 2.0]];
        let mut t = triangulate(&points);
        let before_count = t.num_triangles();
        let edges = [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4)];

        remove_holes_by_edges(&mut t, &edges);

        assert_eq!(t.num_triangles(), before_count);
        assert_neighbors_consistent(&t);
    }

    #[test]
    fn remove_holes_by_edges_noop_when_no_triangles_deleted() {
        let points = [[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0], [2.0, 2.0]];
        let mut t = triangulate(&points);
        let before_vertices = t.triangle_vertices.clone();
        let before_neighbors = t.triangle_neighbors.clone();

        remove_holes_by_edges(&mut t, &[(0, 1), (1, 2), (2, 3), (3, 0)]);

        assert_eq!(t.triangle_vertices, before_vertices);
        assert_eq!(t.triangle_neighbors, before_neighbors);
        assert_neighbors_consistent(&t);
    }

    #[test]
    fn remove_holes_by_edges_empty_edges_is_noop() {
        let points = [[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0], [2.0, 2.0]];
        let mut t = triangulate(&points);
        let before_vertices = t.triangle_vertices.clone();
        let before_neighbors = t.triangle_neighbors.clone();
        let before_constraints = t.constrained_edges.clone();

        remove_holes_by_edges(&mut t, &[]);

        assert_eq!(t.triangle_vertices, before_vertices);
        assert_eq!(t.triangle_neighbors, before_neighbors);
        assert_eq!(t.constrained_edges, before_constraints);
        assert_neighbors_consistent(&t);
    }

    #[test]
    fn remove_holes_by_edges_invalid_edges_are_ignored() {
        let points = [[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0], [2.0, 2.0]];
        let mut t = triangulate(&points);
        let before_vertices = t.triangle_vertices.clone();
        let before_neighbors = t.triangle_neighbors.clone();
        let before_constraints = t.constrained_edges.clone();

        remove_holes_by_edges(&mut t, &[(0, 0), (0, 99), (99, 100)]);

        assert_eq!(t.triangle_vertices, before_vertices);
        assert_eq!(t.triangle_neighbors, before_neighbors);
        assert_eq!(t.constrained_edges, before_constraints);
        assert_neighbors_consistent(&t);
    }

    #[test]
    fn remove_holes_by_edges_open_chain_does_not_panic() {
        let points = [[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0], [2.0, 2.0]];
        let mut t = triangulate(&points);
        let before_vertices = t.triangle_vertices.clone();
        let before_neighbors = t.triangle_neighbors.clone();

        remove_holes_by_edges(&mut t, &[(0, 1), (1, 2)]);

        assert_eq!(t.triangle_vertices, before_vertices);
        assert_eq!(t.triangle_neighbors, before_neighbors);
        assert_neighbors_consistent(&t);
    }

    #[test]
    fn remove_holes_by_edges_restores_state_when_constraint_insertion_fails() {
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
        let before_vertices = t.triangle_vertices.clone();
        let before_neighbors = t.triangle_neighbors.clone();
        let before_constraints = t.constrained_edges.clone();

        remove_holes_by_edges_with_constraint_inserter(
            &mut t,
            &[(0, 8), (2, 6)],
            |triangulation, edges| {
                triangulation
                    .constrained_edges
                    .insert(Triangulation::edge_key(edges[0].0, edges[0].1));
                triangulation.triangle_vertices.reverse();
                triangulation.triangle_neighbors.reverse();
                false
            },
        );

        assert_eq!(t.triangle_vertices, before_vertices);
        assert_eq!(t.triangle_neighbors, before_neighbors);
        assert_eq!(t.constrained_edges, before_constraints);
        assert_neighbors_consistent(&t);
    }

    #[test]
    fn remove_holes_by_edges_restores_state_when_constraint_insertion_panics() {
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
        let before_vertices = t.triangle_vertices.clone();
        let before_neighbors = t.triangle_neighbors.clone();
        let before_constraints = t.constrained_edges.clone();

        remove_holes_by_edges_with_constraint_inserter(
            &mut t,
            &[(0, 8), (2, 6)],
            |triangulation, edges| {
                triangulation
                    .constrained_edges
                    .insert(Triangulation::edge_key(edges[0].0, edges[0].1));
                triangulation.triangle_vertices.reverse();
                triangulation.triangle_neighbors.reverse();
                panic!("simulated constraint insertion panic");
            },
        );

        assert_eq!(t.triangle_vertices, before_vertices);
        assert_eq!(t.triangle_neighbors, before_neighbors);
        assert_eq!(t.constrained_edges, before_constraints);
        assert_neighbors_consistent(&t);
    }
}
