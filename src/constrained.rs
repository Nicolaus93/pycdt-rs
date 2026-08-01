use std::collections::{HashMap, HashSet, VecDeque};

use crate::geometry::{incircle, orient2d};
#[cfg(test)]
use crate::topology::find_shared_edge;
use crate::triangulation::Triangulation;
use crate::types::{Point, NO_NEIGHBOR};

/// One current incident triangle for every vertex in a constraint batch.
///
/// Construction is one linear pass. Flips keep triangle indices stable, so
/// refreshing a rewritten pair preserves O(1) walk starts.
struct IncidentMap {
    triangles: Vec<usize>,
}

impl IncidentMap {
    /// Builds the batch incidence map with one scan of triangle storage.
    fn new(t: &Triangulation) -> Self {
        let mut triangles = vec![NO_NEIGHBOR; t.points.len()];
        for (triangle, vertices) in t.triangle_vertices.iter().enumerate() {
            for &vertex in vertices {
                triangles[vertex] = triangle;
            }
        }
        Self { triangles }
    }

    /// Refreshes incidence for the constant-size neighborhood of a flip.
    fn refresh(&mut self, t: &Triangulation, changed: [usize; 2]) {
        for triangle in changed {
            for &vertex in &t.triangle_vertices[triangle] {
                self.triangles[vertex] = triangle;
            }
        }
    }
}

/// A directed view of one triangle edge, identified by its opposite vertex.
///
/// `local` is the local vertex index in `triangle`; consequently the edge is
/// formed by the other two local vertices and its adjacent triangle is stored
/// in `triangle_neighbors[triangle][local]`. This representation makes all
/// topology access constant-time and avoids constructing temporary vertex
/// sets when two adjacent triangles are inspected.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct LocalEdge {
    triangle: usize,
    local: usize,
}

impl LocalEdge {
    /// Returns the edge endpoints in the triangle's cyclic (CCW) order.
    fn vertices(self, t: &Triangulation) -> (usize, usize) {
        let next = self.next();
        (next.opposite(t), next.next().opposite(t))
    }

    /// Returns the vertex opposite the represented edge.
    fn opposite(self, t: &Triangulation) -> usize {
        t.triangle_vertices[self.triangle][self.local]
    }

    /// Returns the adjacent triangle, if this is not a hull edge.
    fn neighbor(self, t: &Triangulation) -> Option<usize> {
        let neighbor = t.triangle_neighbors[self.triangle][self.local];
        (neighbor != NO_NEIGHBOR).then_some(neighbor)
    }

    /// Rotates the handle counter-clockwise to the next edge of its triangle.
    fn next(self) -> Self {
        Self {
            triangle: self.triangle,
            local: (self.local + 1) % 3,
        }
    }

    /// Finds the reciprocal local handle in the adjacent triangle.
    fn across(self, t: &Triangulation) -> Option<Self> {
        let neighbor = self.neighbor(t)?;
        let local = t.triangle_neighbors[neighbor]
            .iter()
            .position(|&candidate| candidate == self.triangle)?;
        Some(Self {
            triangle: neighbor,
            local,
        })
    }
}

/// Result of flipping the edge designated by a [`LocalEdge`].
#[derive(Clone, Copy, Debug)]
struct FlipResult {
    /// The new diagonal as owned by the first rewritten triangle.
    diagonal: LocalEdge,
    /// The same diagonal as owned by the second rewritten triangle.
    diagonal_twin: LocalEdge,
}

/// Replaces a designated shared edge by the other diagonal of its quadrilateral.
///
/// Rotate the first triangle to `[c,u,v]` and its neighbor to `[d,v,u]`, where
/// `u-v` is the designated edge. The result is `[c,d,v]` and `[c,u,d]`.
/// Triangle indices are stable, and only the two outside neighbors whose owner
/// changes need reciprocal-link repairs. Both returned handles designate the
/// new `c-d` diagonal.
fn flip_designated_edge(t: &mut Triangulation, edge: LocalEdge) -> Option<FlipResult> {
    let twin = edge.across(t)?;
    let a = edge.triangle;
    let b = twin.triangle;
    let av = t.triangle_vertices[a];
    let an = t.triangle_neighbors[a];
    let bv = t.triangle_vertices[b];
    let bn = t.triangle_neighbors[b];
    let ia = edge.local;
    let ib = twin.local;

    let c = av[ia];
    let u = av[(ia + 1) % 3];
    let v = av[(ia + 2) % 3];
    let d = bv[ib];
    // The adjacent CCW triangle sees the shared edge in reverse.
    debug_assert_eq!(bv[(ib + 1) % 3], v);
    debug_assert_eq!(bv[(ib + 2) % 3], u);

    let a_across_u = an[(ia + 1) % 3]; // old edge v-c
    let a_across_v = an[(ia + 2) % 3]; // old edge c-u
    let b_across_v = bn[(ib + 1) % 3]; // old edge u-d
    let b_across_u = bn[(ib + 2) % 3]; // old edge d-v

    t.triangle_vertices[a] = [c, d, v];
    t.triangle_neighbors[a] = [b_across_u, a_across_u, b];
    t.triangle_vertices[b] = [c, u, d];
    t.triangle_neighbors[b] = [b_across_v, a, a_across_v];

    // These are the only boundary edges that changed triangle ownership.
    replace_neighbor(t, b_across_u, b, a);
    replace_neighbor(t, a_across_v, a, b);

    Some(FlipResult {
        diagonal: LocalEdge {
            triangle: a,
            local: 2,
        },
        diagonal_twin: LocalEdge {
            triangle: b,
            local: 1,
        },
    })
}

/// Updates one reciprocal neighbor reference, ignoring a hull sentinel.
fn replace_neighbor(t: &mut Triangulation, triangle: usize, old: usize, new: usize) {
    if triangle == NO_NEIGHBOR {
        return;
    }
    let slot = t.triangle_neighbors[triangle]
        .iter_mut()
        .find(|neighbor| **neighbor == old)
        .expect("topology invariant: outside neighbor must reference old owner");
    *slot = new;
}

/// Resolves an adjacent triangle pair to the edge owned by `tri_a`.
fn local_edge_between(t: &Triangulation, tri_a: usize, tri_b: usize) -> Option<LocalEdge> {
    let local = t
        .triangle_neighbors
        .get(tri_a)?
        .iter()
        .position(|&n| n == tri_b)?;
    Some(LocalEdge {
        triangle: tri_a,
        local,
    })
}

/// Tests strict convexity using the four vertices around a local edge.
fn local_quadrilateral_is_convex(t: &Triangulation, edge: LocalEdge) -> bool {
    let Some(twin) = edge.across(t) else {
        return false;
    };
    let (u, v) = edge.vertices(t);
    let c = edge.opposite(t);
    let d = twin.opposite(t);

    orient2d(&t.points[u], &t.points[v], &t.points[c])
        * orient2d(&t.points[u], &t.points[v], &t.points[d])
        < 0.0
        && orient2d(&t.points[c], &t.points[d], &t.points[u])
            * orient2d(&t.points[c], &t.points[d], &t.points[v])
            < 0.0
}

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
    local_edge_between(t, tri_a, tri_b).is_some_and(|edge| local_quadrilateral_is_convex(t, edge))
}

/// Walk the triangulation from v1 toward v2, collecting all triangle edges that
/// properly intersect segment (v1,v2). Returns None if the walk fails.
pub fn find_intersecting_edges(
    t: &Triangulation,
    v1: usize,
    v2: usize,
) -> Option<Vec<(usize, usize)>> {
    find_intersecting_edges_from(t, v1, v2, &IncidentMap::new(t))
}

/// Internal segment walk seeded from the batch incidence map.
fn find_intersecting_edges_from(
    t: &Triangulation,
    v1: usize,
    v2: usize,
    incidence: &IncidentMap,
) -> Option<Vec<(usize, usize)>> {
    if v1 == v2 || v1 >= t.points.len() || v2 >= t.points.len() {
        return None;
    }
    let start = *incidence.triangles.get(v1)?;
    if start == NO_NEIGHBOR {
        return None;
    }
    let fan = vertex_fan_triangles(t, start, v1);
    if fan
        .iter()
        .any(|&tri| t.triangle_vertices[tri].contains(&v2))
    {
        return Some(Vec::new());
    }

    let mut current = triangle_in_ray_wedge(t, &fan, v1, v2, &HashSet::new())?;
    let mut entry_local = None;
    let mut crossed = Vec::new();
    let mut visited = HashSet::new();

    while visited.insert(current) {
        if t.triangle_vertices[current].contains(&v2) {
            return Some(crossed);
        }

        let mut next = None;
        for local in 0..3 {
            // The segment entered through this edge, so only the other two
            // edges can be exits from a non-degenerate triangle.
            if entry_local == Some(local) {
                continue;
            }
            let handle = LocalEdge {
                triangle: current,
                local,
            };
            let (a, b) = handle.vertices(t);
            let oa = orient2d(&t.points[v1], &t.points[v2], &t.points[a]);
            let ob = orient2d(&t.points[v1], &t.points[v2], &t.points[b]);

            if oa != 0.0 && ob != 0.0 && oa.signum() != ob.signum() {
                let twin = handle.across(t)?;
                crossed.push((current, twin.triangle));
                next = Some((twin.triangle, Some(twin.local)));
                break;
            }

            // Passing through a vertex does not cross an edge. Continue in the
            // unique outgoing wedge of that vertex's triangle fan.
            let through = if oa == 0.0
                && point_strictly_between(&t.points[v1], &t.points[v2], &t.points[a])
            {
                Some(a)
            } else if ob == 0.0
                && point_strictly_between(&t.points[v1], &t.points[v2], &t.points[b])
            {
                Some(b)
            } else {
                None
            };
            if let Some(vertex) = through {
                let fan_start = incidence.triangles[vertex];
                let vertex_fan = vertex_fan_triangles(t, fan_start, vertex);
                if let Some(triangle) = triangle_in_ray_wedge(t, &vertex_fan, vertex, v2, &visited)
                {
                    next = Some((triangle, None));
                    break;
                }
            }
        }

        let (triangle, entry) = next?;
        current = triangle;
        entry_local = entry;
    }
    None
}

/// Collects the connected triangle fan around one vertex without global scans.
fn vertex_fan_triangles(t: &Triangulation, start: usize, vertex: usize) -> Vec<usize> {
    let mut result = Vec::new();
    let mut queue = VecDeque::from([start]);
    let mut visited = HashSet::new();
    while let Some(triangle) = queue.pop_front() {
        if triangle == NO_NEIGHBOR
            || !visited.insert(triangle)
            || !t.triangle_vertices[triangle].contains(&vertex)
        {
            continue;
        }
        result.push(triangle);
        queue.extend(t.triangle_neighbors[triangle].iter().copied());
    }
    result
}

/// Selects the fan wedge entered by the ray `vertex -> target`.
///
/// Rotate a CCW triangle to `[vertex,a,b]`. The ray is inside its closed wedge
/// exactly when it is left of `vertex-a` and right of `vertex-b`. This sign
/// test is scale-independent and therefore replaces the former epsilon probe.
fn triangle_in_ray_wedge(
    t: &Triangulation,
    fan: &[usize],
    vertex: usize,
    target: usize,
    excluded: &HashSet<usize>,
) -> Option<usize> {
    let p = &t.points[vertex];
    let q = &t.points[target];
    fan.iter().copied().find(|&triangle| {
        if excluded.contains(&triangle) {
            return false;
        }
        let vertices = t.triangle_vertices[triangle];
        let Some(local) = vertices.iter().position(|&v| v == vertex) else {
            return false;
        };
        let a = &t.points[vertices[(local + 1) % 3]];
        let b = &t.points[vertices[(local + 2) % 3]];
        orient2d(p, a, q) >= 0.0 && orient2d(p, q, b) >= 0.0
    })
}

/// True when `p` is a non-endpoint point of the collinear segment `a-b`.
fn point_strictly_between(a: &Point, b: &Point, p: &Point) -> bool {
    point_in_segment_bbox(a, b, p) && *p != *a && *p != *b
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

/// Diagonals created while removing a constraint corridor.
struct RemovedEdges {
    /// Current local owners of the new physical diagonals.
    newly_created: Vec<LocalEdge>,
}

/// Removes crossed edges with an incrementally updated queue.
///
/// `Q` owns physical edge keys, so a queued edge remains identifiable when a
/// neighboring flip moves it to another triangle slot. `handles` maps those
/// keys back to their current local owner. A flip refreshes this map by looking
/// only at the six slots of its two rewritten triangles. `queued` guarantees
/// that every physical edge has at most one live queue entry.
fn remove_intersecting_edges_local(
    t: &mut Triangulation,
    v1: usize,
    v2: usize,
    edges: Vec<(usize, usize)>,
    incidence: &mut IncidentMap,
) -> Option<RemovedEdges> {
    if edges.is_empty() {
        return Some(RemovedEdges {
            newly_created: Vec::new(),
        });
    }

    let p = t.points[v1];
    let q = t.points[v2];
    let constraint_edge = Triangulation::edge_key(v1, v2);
    let mut newly_created = Vec::new();
    let mut queue = VecDeque::with_capacity(edges.len());
    let mut queued = HashSet::with_capacity(edges.len());
    let mut handles = HashMap::with_capacity(edges.len() * 2);
    let mut sides = HashMap::with_capacity(edges.len() * 2 + 2);
    sides.insert(v1, 0);
    sides.insert(v2, 0);

    // The segment walk orders the crossed edges, but the flip algorithm only
    // needs queue ownership. Seed each physical edge exactly once.
    for (tri_a, tri_b) in edges {
        let handle = local_edge_between(t, tri_a, tri_b)?;
        let key = edge_key(t, handle);
        handles.insert(key, handle);
        sides
            .entry(key.0)
            .or_insert_with(|| orientation_side(&p, &q, &t.points[key.0]));
        sides
            .entry(key.1)
            .or_insert_with(|| orientation_side(&p, &q, &t.points[key.1]));
        enqueue_edge(key, &mut queue, &mut queued);
    }

    let mut flips = 0usize;
    let mut deferred = 0usize;
    let flip_limit = t.triangle_vertices.len().saturating_mul(32).max(32);

    while let Some(key) = queue.pop_front() {
        queued.remove(&key);
        let handle = *handles.get(&key)?;
        if edge_key(t, handle) != key {
            return None;
        }

        if !local_quadrilateral_is_convex(t, handle) {
            enqueue_edge(key, &mut queue, &mut queued);
            deferred += 1;
            // A complete pass without a flip means no queued diagonal can be
            // removed by the combinatorial flip algorithm.
            if deferred >= queue.len() {
                return None;
            }
            continue;
        }

        let first_triangle = handle.triangle;
        let second_triangle = handle.neighbor(t)?;
        let flipped = flip_designated_edge(t, handle)?;
        incidence.refresh(t, [first_triangle, second_triangle]);
        debug_assert_eq!(flipped.diagonal.across(t), Some(flipped.diagonal_twin));
        flips += 1;
        deferred = 0;
        if flips > flip_limit {
            return None;
        }

        // Boundary edges survive the flip but may change triangle ownership.
        // Refresh only tracked edges visible in the rewritten pair.
        refresh_local_handles(t, [first_triangle, second_triangle], &mut handles);
        let new_key = edge_key(t, flipped.diagonal);
        handles.insert(new_key, flipped.diagonal);

        if new_key == constraint_edge {
            return Some(RemovedEdges {
                newly_created: newly_created
                    .into_iter()
                    .filter_map(|key| handles.get(&key).copied())
                    .collect(),
            });
        }

        if flipped_edge_crosses_constraint(t, v1, v2, new_key, &mut sides) {
            enqueue_edge(new_key, &mut queue, &mut queued);
        } else {
            newly_created.push(new_key);
        }
    }

    Some(RemovedEdges {
        newly_created: newly_created
            .into_iter()
            .filter_map(|key| handles.get(&key).copied())
            .collect(),
    })
}

/// Removes intersecting edges while preserving the historical public result.
pub fn remove_intersecting_edges(
    t: &mut Triangulation,
    v1: usize,
    v2: usize,
    edges: Vec<(usize, usize)>,
) -> Option<Vec<(usize, usize)>> {
    let mut incidence = IncidentMap::new(t);
    remove_intersecting_edges_local(t, v1, v2, edges, &mut incidence).map(|removed| {
        removed
            .newly_created
            .into_iter()
            .map(|edge| edge_key(t, edge))
            .collect()
    })
}

/// Classifies a post-flip diagonal with at most one new orientation predicate.
///
/// The initial corridor labels every crossed-edge endpoint as left (`+`) or
/// right (`-`) of the directed constraint. A convex flip exposes one previously
/// labelled endpoint of the new diagonal and normally one new endpoint. The
/// four non-degenerate configurations are therefore:
///
/// | first side | second side | owner |
/// |------------|-------------|-------|
/// | left       | left        | `N`   |
/// | left       | right       | `Q`   |
/// | right      | left        | `Q`   |
/// | right      | right       | `N`   |
///
/// Only the missing side is evaluated. Zero sides and a layout where neither
/// endpoint was labelled are outside this combinatorial truth table and use the
/// safe full segment-intersection predicate.
fn flipped_edge_crosses_constraint(
    t: &Triangulation,
    v1: usize,
    v2: usize,
    diagonal: (usize, usize),
    sides: &mut HashMap<usize, i8>,
) -> bool {
    let p = &t.points[v1];
    let q = &t.points[v2];
    let known_a = sides.get(&diagonal.0).copied();
    let known_b = sides.get(&diagonal.1).copied();
    let optimized = match (known_a, known_b) {
        (Some(a), Some(b)) if a != 0 && b != 0 => a != b,
        (Some(a), None) if a != 0 => {
            let b = orientation_side(p, q, &t.points[diagonal.1]);
            sides.insert(diagonal.1, b);
            if b == 0 {
                segments_intersect(p, q, &t.points[diagonal.0], &t.points[diagonal.1])
            } else {
                a != b
            }
        }
        (None, Some(b)) if b != 0 => {
            let a = orientation_side(p, q, &t.points[diagonal.0]);
            sides.insert(diagonal.0, a);
            if a == 0 {
                segments_intersect(p, q, &t.points[diagonal.0], &t.points[diagonal.1])
            } else {
                a != b
            }
        }
        _ => segments_intersect(p, q, &t.points[diagonal.0], &t.points[diagonal.1]),
    };
    debug_assert_eq!(
        optimized,
        segments_intersect(p, q, &t.points[diagonal.0], &t.points[diagonal.1]),
        "one-orientation flip classification disagrees with full predicate"
    );
    optimized
}

/// Converts an exact orientation determinant to its combinatorial side label.
fn orientation_side(a: &Point, b: &Point, p: &Point) -> i8 {
    let determinant = orient2d(a, b, p);
    if determinant > 0.0 {
        1
    } else if determinant < 0.0 {
        -1
    } else {
        0
    }
}

/// Returns a local edge's normalized physical vertex pair.
fn edge_key(t: &Triangulation, edge: LocalEdge) -> (usize, usize) {
    let (a, b) = edge.vertices(t);
    Triangulation::edge_key(a, b)
}

/// Adds an edge to `Q` if no live entry already owns it.
fn enqueue_edge(
    edge: (usize, usize),
    queue: &mut VecDeque<(usize, usize)>,
    queued: &mut HashSet<(usize, usize)>,
) {
    if queued.insert(edge) {
        queue.push_back(edge);
    }
}

/// Refreshes tracked local owners from a constant-size rewritten neighborhood.
fn refresh_local_handles(
    t: &Triangulation,
    triangles: [usize; 2],
    handles: &mut HashMap<(usize, usize), LocalEdge>,
) {
    for triangle in triangles {
        for local in 0..3 {
            let handle = LocalEdge { triangle, local };
            let key = edge_key(t, handle);
            if handles.contains_key(&key) {
                handles.insert(key, handle);
            }
        }
    }
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

/// Restores Delaunay legality using only handles produced by insertion.
///
/// `N` owns stable physical keys while `handles` records their current slots.
/// A flip refreshes tracked boundary owners from its two rewritten triangles;
/// no global scan rediscovers an edge. Constrained edges are never flipped.
fn restore_delaunay_edges(
    t: &mut Triangulation,
    edges: Vec<LocalEdge>,
    incidence: &mut IncidentMap,
) -> bool {
    let mut queue = VecDeque::with_capacity(edges.len());
    let mut queued = HashSet::with_capacity(edges.len());
    let mut handles = HashMap::with_capacity(edges.len() * 2);
    for handle in edges {
        let key = edge_key(t, handle);
        handles.insert(key, handle);
        enqueue_edge(key, &mut queue, &mut queued);
    }

    let mut flips = 0usize;
    let flip_limit = t.triangle_vertices.len().saturating_mul(32).max(32);
    while let Some(key) = queue.pop_front() {
        queued.remove(&key);
        if t.constrained_edges.contains(&key) {
            continue;
        }
        let Some(handle) = handles.get(&key).copied() else {
            return false;
        };
        if edge_key(t, handle) != key {
            return false;
        }
        let Some(twin) = handle.across(t) else {
            continue;
        };
        if !local_quadrilateral_is_convex(t, handle) {
            continue;
        }

        let vertices = t.triangle_vertices[handle.triangle];
        if incircle(
            &t.points[vertices[0]],
            &t.points[vertices[1]],
            &t.points[vertices[2]],
            &t.points[twin.opposite(t)],
        ) <= 0.0
        {
            continue;
        }

        let triangles = [handle.triangle, twin.triangle];
        let Some(flipped) = flip_designated_edge(t, handle) else {
            return false;
        };
        incidence.refresh(t, triangles);
        debug_assert_eq!(flipped.diagonal.across(t), Some(flipped.diagonal_twin));
        refresh_local_handles(t, triangles, &mut handles);
        let new_key = edge_key(t, flipped.diagonal);
        handles.insert(new_key, flipped.diagonal);
        if !t.constrained_edges.contains(&new_key) {
            enqueue_edge(new_key, &mut queue, &mut queued);
        }
        flips += 1;
        if flips > flip_limit {
            return false;
        }
    }
    true
}

/// Adds constraints, representing every logical segment by physical subedges.
pub fn add_constraints(t: &mut Triangulation, constraints: &[(usize, usize)]) -> bool {
    let mut incidence = IncidentMap::new(t);
    for &(v1, v2) in constraints {
        let Some(chain) = constraint_vertex_chain(t, v1, v2) else {
            return false;
        };
        for endpoints in chain.windows(2) {
            if !insert_constraint_subedge(t, endpoints[0], endpoints[1], &mut incidence) {
                return false;
            }
        }
    }
    true
}

/// Inserts one subedge whose open segment contains no triangulation vertex.
fn insert_constraint_subedge(
    t: &mut Triangulation,
    v1: usize,
    v2: usize,
    incidence: &mut IncidentMap,
) -> bool {
    let constraint_edge = Triangulation::edge_key(v1, v2);
    let Some(intersecting) = find_intersecting_edges_from(t, v1, v2, incidence) else {
        return false;
    };
    let newly_created = if intersecting.is_empty() {
        RemovedEdges {
            newly_created: Vec::new(),
        }
    } else {
        let Some(removed) = remove_intersecting_edges_local(t, v1, v2, intersecting, incidence)
        else {
            return false;
        };
        removed
    };
    t.constrained_edges.insert(constraint_edge);
    restore_delaunay_edges(t, newly_created.newly_created, incidence)
}

/// Returns all existing vertices on a segment in endpoint-to-endpoint order.
///
/// This point pass is performed before topology mutation. Consecutive entries
/// define the real constrained edges, avoiding a non-physical logical edge when
/// the segment passes through an existing vertex.
fn constraint_vertex_chain(t: &Triangulation, v1: usize, v2: usize) -> Option<Vec<usize>> {
    let (p, q) = (*t.points.get(v1)?, *t.points.get(v2)?);
    if v1 == v2 {
        return None;
    }
    let dx = q[0] - p[0];
    let dy = q[1] - p[1];
    let denominator = dx * dx + dy * dy;
    if denominator == 0.0 {
        return None;
    }
    let mut vertices: Vec<(f64, usize)> = t
        .points
        .iter()
        .enumerate()
        .filter_map(|(vertex, point)| {
            if orient2d(&p, &q, point) != 0.0 || !point_in_segment_bbox(&p, &q, point) {
                return None;
            }
            let parameter = ((point[0] - p[0]) * dx + (point[1] - p[1]) * dy) / denominator;
            Some((parameter, vertex))
        })
        .collect();
    vertices.sort_by(|a, b| a.0.total_cmp(&b.0));
    Some(vertices.into_iter().map(|(_, vertex)| vertex).collect())
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
    fn local_edge_rotation_mapping_and_flip() {
        let mut t = two_tri_quad();
        let edge = local_edge_between(&t, 0, 1).expect("triangles are adjacent");
        assert_eq!(edge.vertices(&t), (2, 0));
        assert_eq!(edge.opposite(&t), 1);
        assert_eq!(edge.next().opposite(&t), 2);
        assert_eq!(edge.across(&t).unwrap().across(&t), Some(edge));

        let flipped = flip_designated_edge(&mut t, edge).expect("square is flippable");
        assert_eq!(
            Triangulation::edge_key(
                flipped.diagonal.vertices(&t).0,
                flipped.diagonal.vertices(&t).1
            ),
            (1, 3)
        );
        assert_eq!(flipped.diagonal.across(&t), Some(flipped.diagonal_twin));
        assert_neighbors_consistent(&t);
    }

    #[test]
    fn queue_membership_has_single_owner() {
        let mut queue = VecDeque::new();
        let mut queued = HashSet::new();
        enqueue_edge((2, 5), &mut queue, &mut queued);
        enqueue_edge((2, 5), &mut queue, &mut queued);
        assert_eq!(queue.into_iter().collect::<Vec<_>>(), vec![(2, 5)]);
        assert_eq!(queued.len(), 1);
    }

    #[test]
    fn restoration_never_flips_a_constrained_handle() {
        let mut t = two_tri_quad();
        let handle = local_edge_between(&t, 0, 1).unwrap();
        let key = edge_key(&t, handle);
        t.constrained_edges.insert(key);
        let before = t.triangle_vertices.clone();
        let mut incidence = IncidentMap::new(&t);
        assert!(restore_delaunay_edges(&mut t, vec![handle], &mut incidence));
        assert_eq!(t.triangle_vertices, before);
    }

    #[test]
    fn incidence_map_is_maintained_after_flip() {
        let mut t = two_tri_quad();
        let mut incidence = IncidentMap::new(&t);
        let edge = local_edge_between(&t, 0, 1).unwrap();
        flip_designated_edge(&mut t, edge).unwrap();
        incidence.refresh(&t, [0, 1]);
        for vertex in 0..t.points.len() {
            assert!(t.triangle_vertices[incidence.triangles[vertex]].contains(&vertex));
        }
    }

    #[test]
    fn flipped_edge_four_configuration_truth_table() {
        let mut t = Triangulation::new();
        t.points = vec![
            [0.0, 0.0],
            [4.0, 0.0],
            [1.0, 1.0],
            [3.0, 1.0],
            [1.0, -1.0],
            [3.0, -1.0],
        ];
        let cases = [
            ((2, 3), false),
            ((2, 4), true),
            ((4, 2), true),
            ((4, 5), false),
        ];
        for (diagonal, expected) in cases {
            let mut sides = HashMap::from([
                (
                    diagonal.0,
                    orientation_side(&t.points[0], &t.points[1], &t.points[diagonal.0]),
                ),
                (
                    diagonal.1,
                    orientation_side(&t.points[0], &t.points[1], &t.points[diagonal.1]),
                ),
            ]);
            assert_eq!(
                flipped_edge_crosses_constraint(&t, 0, 1, diagonal, &mut sides),
                expected
            );
        }
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
        let edges = find_intersecting_edges(&t, 1, 3).expect("walk should succeed");
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
        let edges = find_intersecting_edges(&t, 0, 2).expect("walk should succeed");
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
