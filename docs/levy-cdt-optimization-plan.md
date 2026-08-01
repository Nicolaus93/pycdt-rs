# Lévy CDT constraint-insertion optimization plan

This branch applies the constraint-insertion optimizations described by Lévy
where they fit the crate's current two-dimensional CDT API. Public signatures
remain unchanged. Intersecting and partially overlapping constraints,
symbolic perturbation, and predicate caching are outside this implementation.

## Implementation sequence

1. Introduce private triangle/local-edge handles for constant-time adjacency,
   convexity tests, and designated-edge flips.
2. Walk a constraint once to seed the flip queue, then update that queue
   locally after each flip with constant-time membership marks.
3. Track newly created diagonals by local handles and restore their Delaunay
   legality without global triangle scans.
4. Build a private vertex-to-incident-triangle map once per constraint batch,
   maintain it during flips, and use it to seed segment walks.
5. Replace epsilon probes with orientation-based traversal of the triangle fan,
   splitting a constraint at every existing vertex on the segment.
6. Classify post-flip diagonals using Lévy's one-orientation four-case test,
   retaining a full-intersection fallback for degenerate orientations.
7. Reuse queues, fan storage, incidence data, and generation markers across the
   whole `add_constraints` batch.
8. Preflight each requested segment and reject proper crossings or partial
   overlaps with existing constrained edges before mutating topology. Exact
   duplicate physical constraints remain idempotent.

## Invariants and verification

The implementation will document the designated-edge layout, queue ownership,
fan traversal, marker generations, and post-flip classification truth table.
Tests will cover handle rotation and mapping, incremental queue maintenance,
Delaunay restoration, intermediate on-segment vertices, all four combinatorial
classifications, workspace reuse, duplicates, unsupported-input rejection,
dense grids, reciprocal neighbors, physical constrained-edge existence, and
Delaunay legality of unconstrained edges.

Every implementation commit is verified with
`cargo test --all-targets --verbose`. The completed series is additionally
checked with `just check` when a `justfile` is available, and with the complete
test suite. Release-mode timings for identical dense workloads are recorded in
the final optimization summary.
