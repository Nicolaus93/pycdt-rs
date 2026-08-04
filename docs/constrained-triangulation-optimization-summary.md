# Constrained-triangulation optimization summary

## Algorithm changes

| Concern | Previous implementation | Optimized implementation |
|---|---|---|
| Constraint flip queue | Walked the complete segment again after every flip | Walks once to seed `Q`, then classifies and queues only the new diagonal |
| Topology lookup | Built temporary shared-vertex vectors and rediscovered triangle pairs from vertex pairs | Uses private triangle/local-edge handles and designated-edge flips |
| Walk/restoration starts | Scanned all triangles for endpoints and for every Delaunay edge | Builds one batch incidence map and carries local owners through flips |
| Segment predicates | Repeated full segment predicates and epsilon probe points | Uses orientation-defined fan wedges and a four-case, one-orientation flip classifier with a degenerate fallback |
| Allocation | Allocated queues, sets, and maps for each insertion/restoration call | Reuses batch buffers and O(1) generation-cleared membership marks |
| Unsupported input | Could partially mutate or fail unpredictably on crossing constraints | Preflights crossings, T-junctions, and partial overlap before mutation; exact duplicate chains are idempotent |

Constraints through existing vertices are split into their real consecutive
physical subedges. The public `Triangulation`, `add_constraints() -> bool`,
`find_intersecting_edges()`, and existing public helper signatures remain
unchanged; local handles, incidence, and workspace state are private.

## Complexity

Let `T` be the triangle count, `V` the point count, `k` the crossed-edge
corridor length, `f` the number of constraint-removal flips, `n` the number of
new diagonals restored, and `C` the number of existing constrained edges.

Previously, endpoint discovery and Delaunay edge recovery performed `O(T)`
scans, and removal repeated an `O(T + k)` segment walk after flips, giving an
`O(f(T + k) + nT)` dominant insertion cost. Temporary set/vector allocation
was repeated during these operations.

The optimized batch pays `O(T)` once for incidence construction. Each logical
constraint pays `O(V)` to identify exact on-segment vertices and `O(Ck)` for the
unsupported-interaction preflight. The topology portion is `O(k + f + n)`
expected time: queue membership and edge ownership use hash tables, and every
flip inspects a constant-size two-triangle neighborhood. Fan traversal is local
to incident triangles. Workspace buffers retain capacity across the batch.

## Release-mode measurement

The benchmark timed only constraint insertion: an already triangulated regular
5×5 grid was cloned and its established corner-to-corner long constraint was
inserted 10,000 times. Both revisions were built with `cargo run --release` on
the same machine.

| Revision | Elapsed | Relative |
|---|---:|---:|
| `3f7b383` (before) | 192.372 ms | 1.00× |
| `d59024e` (optimized code) | 86.283 ms | 2.23× faster |

Attempted 60×60 long/multiple-constraint variants were not reportable as a
before/after comparison because the pre-optimization inserter returned failure
on those fixtures. No failed run is included in the speedup figure.

## Verification and limitations

Focused tests cover local-handle rotation/flips, reciprocal neighbors,
incremental queue ownership, constrained-edge protection during Delaunay
restoration, maintained incidence, orientation fan walks, intermediate
vertices, all four post-flip configurations, generation reuse, duplicate
constraints, partial overlap, and topology-preserving crossing rejection. The
full all-target suite is run for every series commit.

Properly intersecting and partially overlapping constraints remain unsupported
and are rejected. Symbolic perturbation and general predicate caching remain
deferred; their required data model, API, constructions, and tests are described
in `constrained-triangulation-future-work.md`.

## Commit series

- `3242e06` document cdt optimization plan
- `0d7edd0` use triangle-local edge handles
- `0969f54` update constraint flip queues incrementally
- `11c2410` restore delaunay edges from local handles
- `7ccc97b` remove global scans from constraint walks
- `ed5cb51` streamline constraint segment walks
- `4a9247f` classify flipped edges with one orientation
- `05a2973` reuse constraint insertion work buffers
- `37a76a5` reject unsupported constraint intersections safely
- `d59024e` document advanced constrained triangulation features

The hash of this summary commit is intentionally not self-referential; it is
available from the branch log as `summarize cdt constraint insertion optimizations`.
