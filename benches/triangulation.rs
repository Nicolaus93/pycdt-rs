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
#[divan::bench(args = [64, 256, 512])]
fn triangulate_uniform(bencher: Bencher, n: usize) {
    let points = uniform_points(n, 0x5EED);
    bencher.bench(|| divan::black_box(triangulate(divan::black_box(&points))));
}

/// Regular grid: many collinear points and cocircular quadruples, exercising the
/// on-edge insertion path and the robust predicates.
#[divan::bench(args = [8, 16, 22])]
fn triangulate_grid(bencher: Bencher, side: usize) {
    let points = grid_points(side);
    bencher.bench(|| divan::black_box(triangulate(divan::black_box(&points))));
}

/// Points on a circle: fully degenerate for the in-circle predicate.
#[divan::bench(args = [64, 256])]
fn triangulate_circle(bencher: Bencher, n: usize) {
    let points = circle_points(n, 100.0);
    bencher.bench(|| divan::black_box(triangulate(divan::black_box(&points))));
}

/// Incremental insertion of new points into an already built triangulation.
///
/// The super triangle is dropped by `triangulate`, so the extra points are
/// drawn from a smaller centred box to guarantee they land inside the convex
/// hull of the base point set.
#[divan::bench(args = [16, 64])]
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
#[divan::bench(args = [1, 8, 32])]
fn add_constraints_crossing(bencher: Bencher, count: usize) {
    // Points on a circle, constraints are long chords crossing the interior so
    // every insertion has to remove and re-triangulate intersecting edges.
    let n = 128;
    let points = circle_points(n, 100.0);
    let triangulation = triangulate(&points);
    let constraints: Vec<(usize, usize)> = (0..count).map(|i| (i, (i + n / 2) % n)).collect();

    bencher
        .with_inputs(|| triangulation.clone())
        .bench_values(|mut t| {
            divan::black_box(add_constraints(&mut t, divan::black_box(&constraints)))
        });
}

/// Hole removal: insert the boundary constraints, then discard every triangle
/// outside the outer polygon or inside a hole. `side` is the number of holes per
/// row of the `side x side` hole grid.
#[divan::bench(args = [1, 2, 4])]
fn remove_holes(bencher: Bencher, side: usize) {
    let (points, edges) = polygon_with_holes(side, 32 * side);
    let triangulation = triangulate(&points);

    bencher
        .with_inputs(|| triangulation.clone())
        .bench_values(|mut t| {
            remove_holes_by_edges(&mut t, divan::black_box(&edges));
            divan::black_box(t.num_triangles())
        });
}
