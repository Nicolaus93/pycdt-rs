# pycdt-rs

[![CI](https://github.com/Nicolaus93/pycdt-rs/actions/workflows/tests.yml/badge.svg)](https://github.com/Nicolaus93/pycdt-rs/actions/workflows/tests.yml)

A Rust implementation of 2D Delaunay and constrained Delaunay triangulation.

The repository is now split into two projects:

- the **core Rust crate** at the repository root
- the **Python bindings project** in [`python/`](python)

This lets other Rust projects depend on the triangulation library without pulling in `pyo3`, `numpy`, or other Python-side tooling.

## Repository layout

- `src/` — core Rust implementation
- `tests/` — Rust tests
- `python/` — standalone Python extension project built with `maturin`
- `python/examples/` — Python examples

## Rust usage

Build or test the core crate:

```bash
cargo build --release
cargo test
```

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
- hole removal by seed points or polygon edges
- robust predicates via the `robust` crate
- Rust API for native use
- Python bindings as a separate project

## Python API

The Python extension exposes:

- `triangulate(points)`
- `update_triangulation(triangulation, new_points)`
- `add_constraints(triangulation, edges)`
- `remove_holes(triangulation, holes)`
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

## License

MIT License.
