# pycdt-rs

[![CI](https://github.com/Nicolaus93/pycdt-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/Nicolaus93/pycdt-rs/actions/workflows/ci.yml)
[![CodSpeed](https://img.shields.io/endpoint?url=https://codspeed.io/badge.json)](https://app.codspeed.io/Nicolaus93/pycdt-rs?utm_source=badge)

A Rust implementation of 2D Delaunay and constrained Delaunay triangulation.

The repository is now split into two projects:

- the **core Rust crate** at the repository root
- the **Python bindings project** in [`python/`](python)

This lets other Rust projects depend on the triangulation library without pulling in `pyo3`, `numpy`, or other Python-side tooling.

## Repository layout

- `src/` — core Rust implementation
- `tests/` — Rust tests
- `benches/` — Rust benchmarks, tracked with CodSpeed
- `python/` — standalone Python extension project built with `maturin`
- `python/examples/` — Python examples

## Rust usage

Build or test the core crate:

```bash
cargo build --release
cargo test
```

Lint/format helpers are available through `just`:

```bash
just fmt
just fmt-check
just clippy
just check
```

## Benchmarks

Benchmarks live in [`benches/`](benches) and use
[divan](https://github.com/nvzqz/divan) through the
[CodSpeed compatibility layer](https://codspeed.io/docs/benchmarks/rust/divan).
They are run on every pull request and the results are reported by CodSpeed.

Run them locally with the [CodSpeed CLI](https://codspeed.io/docs/cli):

```bash
cargo codspeed build
codspeed run --mode simulation -- cargo codspeed run
```

Two suites are available:

- `triangulation` — end-to-end triangulation, incremental updates, constrained
  edge insertion and hole removal
- `geometry` — the robust predicates and topology helpers on the hot path

Use it from another Rust project as a normal dependency:

```toml
[dependencies]
pycdt-rs = { path = "../pycdt-rs" }
```

Example:

```rust
use pycdt_rs::build::triangulate;

fn main() {
    let points = vec![
        [0.0, 0.0],
        [1.0, 0.0],
        [1.0, 1.0],
        [0.0, 1.0],
        [0.5, 0.5],
    ];

    let t = triangulate(&points);
    println!("{} points, {} triangles", t.num_points(), t.num_triangles());
}
```

## Python usage

The Python bindings live in [`python/`](python).

Build and install them with:

```bash
cd python
maturin develop --release
```

Or from the repository root:

```bash
maturin develop --release --manifest-path python/Cargo.toml
```

Then import:

```python
import pycdt_rs
```

## Features

- Delaunay triangulation via incremental point insertion
- constrained edge insertion
- incremental updates
- hole removal by polygon edges
- robust predicates via the `robust` crate
- Rust API for native use
- Python bindings as a separate project

## Python API

The Python extension exposes:

- `triangulate(points)`
- `update_triangulation(triangulation, new_points)`
- `add_constraints(triangulation, edges)`
- `remove_holes_by_edges(triangulation, edges)`
- `remove_super_triangle(triangulation)`
- `build_polygons_from_edges(triangulation, edges)`

`triangulate()` returns a `PyTriangulation` object with:

- `points`
- `triangle_vertices`
- `triangle_neighbors`
- `constrained_edges`
- `num_points`
- `num_triangles`

## Notes

- constraint indices refer to rows in the input point array
- Python neighbor arrays use `-1` for missing neighbors
- `update_triangulation()` should be used before removing the super triangle

## Limitations

The optimized inserter intentionally supports non-intersecting, non-overlapping
2D constraints over an existing point set. The following features require
changes to the data model and public API; they should not be approximated with
epsilon offsets or untracked topology edits.

### Intersecting constraints

A proper crossing must be constructed as a new stable point at the exact
intersection of the two supporting lines. Both physical constrained edges must
be removed and replaced by four edges incident to that point. The affected
triangles then need a local cavity retriangulation, reciprocal-neighbor repair,
incidence-map updates, and constrained Delaunay restoration outside the new
constraint edges.

One logical constraint can consequently become a chain of physical edges. The
triangulation must retain logical constraint IDs and an ordered chain per ID so
splitting an old constraint updates every owner instead of merely changing a
set of vertex pairs. Tests need multiple crossings, endpoint-on-edge
T-junctions, crossing order independence, exact rationally representable and
ill-conditioned intersections, neighbor reciprocity, chain ownership, and
Delaunay legality of every unconstrained edge.

### Overlapping constraints

Collinear overlap cannot be represented faithfully by the current
`HashSet<(usize, usize)>`: one physical edge may belong to several logical
constraints, and a partial overlap introduces split points at both overlap
boundaries. Supporting it requires edge-to-constraint ownership sets (or
reference counts), ordered chains for each logical constraint, exact projection
ordering, and atomic splitting of all owners. Removing or querying one logical
constraint must not discard an edge still owned by another.

Tests should cover equal reversed segments, containment, partial overlap on
either end, several coincident constraints, overlap combined with an existing
vertex, ownership removal, and insertion-order independence.

### Stable point identities and API changes

Intersection construction adds points after triangulation, so callers need a
stable point/vertex identity that survives internal storage changes. A future
API should return a constraint ID and expose its physical vertex chain, define
whether generated vertices are visible in `points`, and report unsupported or
invalid input with a structured error rather than a boolean. Serialization and
Python bindings must preserve generated-point and logical-constraint IDs.

An atomic insertion transaction is also needed: validate and construct all
split points and chains first, then commit topology and ownership together, or
restore the prior triangulation on any failure.

### Symbolic perturbation

Exact zero predicates currently select documented degenerate fallbacks.
Full symbolic perturbation needs a consistent simulation-of-simplicity ordering
based on stable point identities, applied uniformly to orientation, incircle,
fan selection, intersection ordering, and cavity retriangulation. Applying a
local tie-break in only one predicate can produce contradictory topology.

Tests must permute input and insertion order over collinear/cocircular fixtures,
verify deterministic physical topology and chains, and exercise generated
intersection points whose algebraic coordinates compare equal.

### Predicate caching

General orientation/incircle caching is deferred pending profiling. The current
one-orientation flip classifier carries only the combinatorial side labels
needed by a single insertion; it is not a global predicate cache. A future cache
must define keys for stable point identities (including generated points),
invalidation rules, memory bounds, and measured hit rates before adding its
complexity.

## Pre-commit

This repo includes a `pre-commit` config that runs Rust formatting and clippy checks:

```bash
pip install pre-commit
pre-commit install
pre-commit run --all-files
```

The hooks run:

- `just fmt-check`
- `just clippy`

## License

MIT License.
