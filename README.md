# Constrained Delaunay Triangulation

[![CI](https://github.com/Nicolaus93/pycdt-rs/actions/workflows/tests.yml/badge.svg)](https://github.com/Nicolaus93/pycdt-rs/actions/workflows/tests.yml)

A Rust implementation of 2D Delaunay and constrained Delaunay triangulation, with Python bindings for experimentation, visualization, and educational use.

## Overview

This project implements incremental Delaunay triangulation with robust geometric predicates for numerical stability. It aims to stay reasonably clear and approachable while still providing practical constrained-triangulation functionality.

The repository contains:

- a native Rust library
- a Python extension module exposed as `pycdt_rs`
- tests and property tests
- Python examples in `examples/`

## Features

- **Delaunay triangulation** via incremental point insertion
- **Constrained edge insertion** for fixed edges and boundaries
- **Incremental updates** with `update_triangulation`
- **Hole removal** by seed points or polygon edges
- **Robust predicates** via the `robust` crate for:
  - `orient2d`
  - `incircle`
- **Python bindings** built with `pyo3` and `maturin`
- **Rust API** for direct native use

## Status

- ✅ **Delaunay triangulation**: implemented and tested
- ✅ **Constrained triangulation**: implemented and tested
- ✅ **Incremental updates**: implemented and tested
- ✅ **Hole removal**: implemented and tested

## Installation

### Python

Build and install the extension locally with `maturin`:

```bash
maturin develop --release
```

Or build a wheel:

```bash
maturin build --release
```

Then import:

```python
import pycdt_rs
```

### Rust

Build or test the crate directly:

```bash
cargo build --release
cargo test
```

## Quick Start

### Python

```python
import numpy as np
import pycdt_rs

points = np.random.rand(50, 2) * 100.0
tri = pycdt_rs.triangulate(points)

print(tri.num_points)
print(tri.num_triangles)
print(tri.triangle_vertices)
```

### Rust

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

See `examples/` for more complete usage.

## Python API

Functions exposed by `pycdt_rs`:

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

Notes:

- constraint indices refer to rows in the input point array
- Python neighbor arrays use `-1` for missing neighbors
- `update_triangulation()` should be used before removing the super triangle

## Algorithm Details

The implementation uses incremental insertion:

1. Create an initial super triangle
2. Insert points one by one
3. Locate the containing triangle or edge
4. Split affected triangles
5. Restore the Delaunay property with local edge flipping
6. Remove the super triangle when finalizing the mesh

For constrained triangulation, fixed edges are inserted into the triangulation and preserved during subsequent operations.

## Robust Predicates

Computational geometry is sensitive to floating-point error, especially for nearly collinear points and near-degenerate configurations. This project uses robust predicates to improve numerical reliability:

- **orient2d**: orientation of three points
- **incircle**: whether a point lies inside the circumcircle of a triangle

These predicates are central to triangle location, edge legality checks, and topological updates.

## Examples

The `examples/` directory includes:

- `examples/basic_triangulation.py`
- `examples/constrained_triangulation.py`
- `examples/points_on_edge.py`
- `examples/circle.py`
- `examples/plane_face.py`
- `examples/plane_face_2.py`
- `examples/cyl_face.py`
- `examples/toroidal_face.py`

Run an example after installing the module locally:

```bash
python3 examples/basic_triangulation.py
```

## References

- Sloan, S.W. (1987). *A fast algorithm for constructing Delaunay triangulations in the plane*
- Shewchuk, J.R. (1997). *Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates*
- de Berg, M., et al. (2008). *Computational Geometry: Algorithms and Applications*

## License

MIT License.

## Acknowledgments

- Jonathan Shewchuk for robust geometric predicates
- S.W. Sloan for the incremental triangulation algorithm
