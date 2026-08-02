//! Benchmarks for the end-to-end triangulation pipeline: building a Delaunay
//! triangulation, updating it incrementally, inserting constraints and removing
//! holes.

mod common;

use common::{circle_points, grid_points, polygon_with_holes, uniform_points, uniform_points_in};
use divan::Bencher;
use pycdt_rs::build::{remove_holes_by_edges, triangulate, update_triangulation};
use pycdt_rs::constrained::add_constraints;

fn main() {
    divan::main();
}

/// Delaunay triangulation of uniformly distributed points, the most common
/// workload.
#[divan::bench(args = [1_000, 10_000])]
fn triangulate_uniform(bencher: Bencher, n: usize) {
    let points = uniform_points(n, 0x5EED);
    bencher.bench(|| divan::black_box(triangulate(divan::black_box(&points))));
}

/// Regular grid: many collinear points and cocircular quadruples, exercising the
/// on-edge insertion path and the robust predicates.
#[divan::bench(args = [1_000, 10_000])]
fn triangulate_grid(bencher: Bencher, n: usize) {
    let points = grid_points(n);
    bencher.bench(|| divan::black_box(triangulate(divan::black_box(&points))));
}

/// Points on a circle: fully degenerate for the in-circle predicate.
#[divan::bench(args = [1_000, 10_000])]
fn triangulate_circle(bencher: Bencher, n: usize) {
    let points = circle_points(n, 100.0);
    bencher.bench(|| divan::black_box(triangulate(divan::black_box(&points))));
}

/// Incremental insertion of new points into an already built triangulation.
///
/// The super triangle is dropped by `triangulate`, so the extra points are
/// drawn from a smaller centred box to guarantee they land inside the convex
/// hull of the base point set.
#[divan::bench(args = [100, 1_000])]
fn update_triangulation_incremental(bencher: Bencher, added: usize) {
    let base = uniform_points(256, 0x5EED);
    let extra = uniform_points_in(added, 0xBEEF, 250.0, 750.0);
    let triangulation = triangulate(&base);

    bencher
        .with_inputs(|| triangulation.clone())
        .bench_values(|mut t| {
            update_triangulation(&mut t, divan::black_box(&extra));
            divan::black_box(t.num_triangles())
        });
}

/// Constrained edge insertion: each constraint forces edge flips along the
/// segment it crosses.
#[divan::bench(args = [100, 1_000])]
fn add_constraints_crossing(bencher: Bencher, n: usize) {
    // A diameter of the circle crosses the triangulation interior, forcing
    // edge flips without introducing intersections between constraints.
    let points = circle_points(n, 100.0);
    let triangulation = triangulate(&points);
    let constraints = [(0, n / 2)];

    bencher
        .with_inputs(|| triangulation.clone())
        .bench_values(|mut t| {
            divan::black_box(add_constraints(&mut t, divan::black_box(&constraints)))
        });
}

/// Hole removal: insert the boundary constraints, then discard every triangle
/// outside the outer polygon or inside one of four holes.
#[divan::bench(args = [64, 128])]
fn remove_holes(bencher: Bencher, n: usize) {
    const HOLES_PER_ROW: usize = 2;
    const BOUNDARY_POINTS: usize = 4 + 4 * HOLES_PER_ROW * HOLES_PER_ROW;
    let (points, edges) = polygon_with_holes(HOLES_PER_ROW, n - BOUNDARY_POINTS);
    debug_assert_eq!(points.len(), n);
    let triangulation = triangulate(&points);

    bencher
        .with_inputs(|| triangulation.clone())
        .bench_values(|mut t| {
            remove_holes_by_edges(&mut t, divan::black_box(&edges));
            divan::black_box(t.num_triangles())
        });
}
