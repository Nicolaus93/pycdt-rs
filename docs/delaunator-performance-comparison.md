# Why Delaunator Is Faster

The algorithm is the primary reason. The data structure contributes, but
`Vec<[usize; 3]>` is not inherently inefficient.

## The representation difference

This crate stores:

```rust
triangle_vertices: Vec<[usize; 3]>
triangle_neighbors: Vec<[usize; 3]>
```

This is already contiguous memory. It is not like `Vec<Vec<usize>>`; there are
no per-triangle allocations or pointer chasing.

On a 64-bit machine, however:

| Per triangle | This crate | Delaunator |
|---|---:|---:|
| Vertex indices | 3 × 8 = 24 bytes | 3 × 4 = 12 bytes |
| Adjacency | 3 × 8 = 24 bytes | 3 × 4 = 12 bytes |
| Total topology | 48 bytes | 24 bytes |

At roughly 200k triangles, that is about:

- This crate: 9.6 MB of topology
- Delaunator: 4.8 MB
- Coordinates add another 1.6 MB to both

Delaunator therefore gets about twice as many topology entries into each cache
line. Its active data is also more likely to fit in the CPU's shared cache.

But the more important distinction is what the adjacency indices mean.

### Triangle neighbors versus half-edges

This crate stores the adjacent triangle:

```text
triangle_neighbors[triangle][opposite_vertex] = adjacent_triangle
```

Delaunator stores the exact opposite half-edge:

```text
halfedges[edge] = opposite_edge
```

Given half-edge `e`, Delaunator immediately knows:

- Its triangle: `e / 3`
- Its next edge: simple index arithmetic
- Its previous edge: simple index arithmetic
- Its precise twin: `halfedges[e]`
- The twin's local position: already encoded in the index

The [Delaunator data-structure guide](https://mapbox.github.io/delaunator/)
describes this mapping.

By contrast, this implementation frequently recovers information that its
neighbor representation does not directly encode. For example:

- Point location identifies an exit edge, but then searches all three neighbors
  and inspects their vertices to find the matching one
  ([`build.rs`](../src/build.rs#L153)).
- Diagonal swapping constructs a temporary `Vec` to discover the two shared
  vertices ([`topology.rs`](../src/topology.rs#L83)).
- It repeatedly uses `.contains()` and `.position()` to recover local
  vertex/edge positions.
- Updating a back-reference scans the neighboring triangle's three slots
  ([`topology.rs`](../src/topology.rs#L60)).

Each scan is only three elements, but these operations occur millions of times.
The allocation in `swap_diagonal` is particularly expensive.

A Rust half-edge representation using `u32` could therefore be faster for two
separate reasons:

1. Half the memory bandwidth.
2. Direct topology navigation without repeated searches.

Simply changing this to `Vec<u32>` or flattening `[usize; 3]` would help cache
density, but would not eliminate the topology searches.

## The larger difference: point insertion order

This crate processes points in caller-provided order:

```rust
for &point in input_points {
    // locate point in current triangulation
    // split containing triangle
    // legalize surrounding edges
}
```

For uniformly random input, consecutive points are spatially unrelated. Point
location begins at the triangle used by the previous insertion and walks across
the mesh to the next point ([`build.rs`](../src/build.rs#L106)).

In a roughly square triangulation of `n` uniform points, crossing from one random
location to another can traverse approximately `O(√n)` triangles. That makes
this particular walking strategy behave closer to:

```text
O(n√n)
```

on the uniform workload, despite incremental Delaunay construction having
better bounds when paired with an effective point-location structure.

The measured scaling supports this:

```text
10k points:   ~25 ms
100k points: ~974 ms
```

Ten times as many points takes roughly 39 times as long—an effective exponent
around 1.6, close to `n^1.5`.

The point-location walk also performs substantial work per crossed triangle:

- Three robust orientation tests in `point_in_triangle`
- Up to three more orientation tests to choose the exit edge
- Neighbor scanning
- Neighbor-vertex membership tests

So this is probably the largest hot spot.

## What Delaunator does instead

Delaunator is also incremental, but it carefully controls insertion order:

1. Select a seed triangle near the center.
2. Sort points by distance from the seed circumcenter.
3. Maintain an advancing convex hull.
4. Find a likely visible hull edge using an angular hash.
5. Walk only the visible portion of the hull.
6. Add new triangles along that portion.
7. Legalize the affected edges.

You can see this directly in Delaunator's
[implementation](https://github.com/mapbox/delaunator/blob/main/index.js).

Because points are inserted radially outward, Delaunator does not perform
general point location inside an increasingly large triangulation. Every new
point is handled relative to the current hull. Sorting costs `O(n log n)`, while
most of the construction after sorting is close to linear/amortized local work.

That is inherently better for batch triangulation than arbitrary-order interior
insertion.

## Edge legalization is also much tighter

Delaunator's legalization loop:

- Uses a fixed reusable edge stack.
- Receives an exact half-edge.
- Locates all four vertices using index arithmetic.
- Flips an edge with a few array assignments.
- Updates three half-edge links directly.

This crate currently:

- Creates a new `Vec` stack for each `lawson_swapping` call.
- Calls Lawson swapping three times after an interior split.
- Checks for stale triangles and neighbor relationships.
- Searches for shared vertices.
- Allocates a temporary `Vec` during every diagonal swap.
- Recomputes vertex positions and neighbor positions.
- Updates neighbor back-references by scanning.

See [`lawson_swapping`](../src/topology.rs#L220) and
[`swap_diagonal`](../src/topology.rs#L83).

There is also a predicate trade-off: this crate uses robust adaptive `incircle`;
Delaunator uses robust orientation but a compact ordinary floating-point
in-circle calculation in its hot legalization loop. That gives Delaunator
another speed advantage, although it is unlikely to explain most of the 18× gap
for uniform points.

## Likely importance

Pending profiling, the likely order is:

1. **Arbitrary-order point location:** dominant difference.
2. **More expensive edge-flip/topology bookkeeping:** substantial.
3. **Half-sized, directly navigable half-edge arrays:** meaningful.
4. **Repeated small allocations and lack of preallocation:** smaller but
   fixable.
5. **Predicate differences and finalization:** secondary for uniform data.

The best optimization path would be:

1. Use `triangle_neighbors[current][opposite_vertex]` directly during point
   location; the current neighbor search appears unnecessary.
2. Remove the allocation from `swap_diagonal`.
3. Reuse one legalization stack and invoke legalization once per insertion.
4. Reserve point and triangle capacities up front.
5. Spatially sort batch inputs while preserving original point indices.
6. Eventually represent adjacency as exact `u32` half-edge indices.
7. For maximum batch performance, implement a dedicated sweep-hull path while
   retaining the existing incremental path for updates and constraints.

Steps 1–4 should improve the current implementation without changing its public
topology model. Steps 5–7 address the fundamental algorithmic gap.
