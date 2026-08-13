# Triangulation optimization results

This report measures the optimization sequence proposed in
[`delaunator-performance-comparison.md`](delaunator-performance-comparison.md).
Each optimization lives on a separate cumulative branch: branch `N` contains
one new optimization on top of branch `N - 1`.

## Method

- Hardware: AMD Ryzen 7 7730U, 8 cores / 16 threads
- Build: `cargo run --release`, LTO enabled, one codegen unit
- Input: 100,000 deterministic `StdRng` points uniformly distributed in a
  1000 × 1000 square, seed `0x5EED`
- Timing boundary: `triangulate(&points)` only; point generation is excluded
- Samples: two warmups followed by nine measured calls in one process
- Statistic: median wall-clock time
- Correctness: triangle count and boundary-edge count are checked against an
  independent robust monotone-chain convex hull in
  [`profile_uniform.rs`](../examples/profile_uniform.rs)

Run the profile with:

```bash
cargo run --release --example profile_uniform -- 100000 9
```

Hardware performance counters were unavailable because the host has
`perf_event_paranoid=4`, so “profile” here means repeatable wall-clock and
topology measurements rather than CPU event sampling.

## Branches and results

| Step | Branch | Commit | Median | Change from parent | Versus baseline | Topology |
|---:|---|---|---:|---:|---:|---|
| 0 | `opt/00-profile-baseline` | `eef6a4f` | 1310.451 ms | — | 1.00× | 199,970 / 28 |
| 1 | `opt/01-direct-neighbor` | `485e7fe` | 1403.764 ms | +7.1% | 0.93× | 199,970 / 28 |
| 2 | `opt/02-no-swap-allocation` | `44efef3` | 1326.533 ms | −5.5% | 0.99× | 199,970 / 28 |
| 3 | `opt/03-reuse-legalization-stack` | `2b3ce7f` | 1363.290 ms | +2.8% | 0.96× | 199,970 / 28 |
| 4 | `opt/04-preallocate-storage` | `b42bba3` | 1217.410 ms | −10.7% | 1.08× | 199,970 / 28 |
| 5 | `opt/05-spatial-sort` | `90e9a14` | 83.027 ms | −93.2% | 15.78× | 199,970 / 28 |
| 6 | `opt/06-u32-halfedges` | `94e20de` | 113.600 ms | +36.8% | 11.54× | 199,970 / 28 |
| 7 | `opt/07-sweep-hull` | `85e2840` | 61.552 ms | −45.8% | 21.29× | 199,971 / 27 |

“Topology” is `triangle count / boundary-edge count`. The independent hull
oracle gives 27 hull vertices and therefore `2n - 2 - h = 199,971` triangles.
Only step 7 matches that result. Steps 0–6 retain one missing boundary triangle
from the finite-super-triangle incremental construction.

These branches were measured sequentially on a frequency-scaling laptop, so
small differences should not be over-interpreted. The order-of-magnitude result
at step 5 and the further step-7 improvement are much larger than run-to-run
noise.

## Step-by-step interpretation

### 1. Direct neighbor indexing

Point location now uses the neighbor slot opposite the selected exit vertex
directly. It no longer scans all three neighbor triangles and compares their
vertices. This also exposed and corrected one hand-built test fixture whose
neighbor slots did not follow the production opposite-vertex invariant.

The measured result was 7.1% slower, which means the removed three-element
search was not a meaningful part of total time on this workload. Robust
orientation predicates and the number of triangles crossed dominate point
location. The change remains valuable because it makes the topology invariant
explicit and removes unnecessary work.

### 2. Allocation-free diagonal swap

`swap_diagonal` no longer collects shared vertices into a temporary `Vec`.
It identifies the one non-shared vertex and derives the two shared vertices and
their local positions directly from the fixed triangle arrays.

This recovered 5.5% relative to step 1. The cumulative result remains within
about 1% of baseline, consistent with allocator traffic being real but not the
dominant cost.

### 3. Reused legalization stack

One stack is reused across a complete batch or incremental update. An insertion
seeds all cavity-boundary edges into a single legalization call instead of
starting three or four separate Lawson walks.

The result was 2.8% slower than step 2. Consolidating the calls changes stack
and edge-flip order; on this input, that secondary effect outweighs the saved
small allocations. The branch is useful evidence that stack allocation was not
the main bottleneck.

### 4. Preallocated storage

Point storage reserves `n + 3`, while vertex and neighbor topology reserve the
planar upper bound `2n + 1`. This avoids geometric `Vec` growth and copies while
the mesh is built.

The median improved by 10.7% relative to step 3 and by 7.1% relative to the
baseline.

### 5. Spatially sorted insertion

Batch input indices are sorted by a 16-bit-per-axis Morton key. Points remain in
their original order in the public point array, so triangle and constraint
indices still refer to the caller's input. Only insertion order changes.

This is the decisive optimization: median time falls from 1217 ms to 83 ms.
Successive point-location walks are local instead of crossing a random fraction
of the mesh, changing the practical scaling behavior of incremental insertion.

The finite super triangle still causes one missing boundary face for this
particular 100k set. Spatial sorting improves speed, not that underlying
finalization limitation.

### 6. Exact `u32` twin half-edges

The triangulation gains an exact twin-half-edge array encoded as
`triangle * 3 + opposite_vertex_slot`. Point location and Lawson legalization
use it for direct navigation, and a new test verifies twin symmetry.

For compatibility, this experimental branch retains the existing public
triangle-neighbor array and synchronizes both representations after local
topology changes. That transitional design is 36.8% slower than step 5 because
it duplicates storage and bookkeeping. It demonstrates that half-edges are not
automatically faster: they need to replace the old representation, not shadow
it.

### 7. Dedicated sweep-hull batch path

Batch triangulation uses the Rust `delaunator` sweep-hull implementation. Its
native exact half-edges are converted once into this crate's CCW,
opposite-vertex topology convention. Incremental updates, constraints, and hole
removal continue to use the existing mutable topology after construction.

The median falls to 61.6 ms, a 21.29× speedup over baseline and 25.9% faster
than Morton-sorted incremental insertion. It is also the only branch that
matches the independent hull oracle: 27 boundary edges and 199,971 triangles.

The sweep-hull dependency uses robust orientation but a fast floating-point
in-circle test. This is the same performance/robustness trade-off discussed in
the comparison document and should be considered when selecting the final
production branch.

## Conclusions

The measurements support three conclusions:

1. General point location on randomly ordered input—not Rust's
   `Vec<[usize; 3]>` layout—is the dominant original cost.
2. Spatial locality alone closes most of the performance gap without changing
   the public index convention.
3. Exact half-edges pay off when they are the native construction topology;
   maintaining them as a compatibility shadow is slower.

For the smallest change with most of the speedup, use
`opt/05-spatial-sort`. For maximum batch speed and correct convex-hull
finalization on this workload, use `opt/07-sweep-hull` after accepting the
in-circle predicate trade-off and the new dependency.
