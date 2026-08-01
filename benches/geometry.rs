//! Micro-benchmarks for the geometric predicates and topological helpers that
//! sit on the hot path of the triangulation algorithms.

mod common;

use common::{uniform_points, Rng};
use divan::Bencher;
use pycdt_rs::build::{
    build_polygons_from_edges, find_containing_triangle, polygon_area, triangulate,
};
use pycdt_rs::constrained::segments_intersect;
use pycdt_rs::geometry::{incircle, is_point_inside_polygon, orient2d, point_in_triangle};
use pycdt_rs::types::Point;

fn main() {
    divan::main();
}

const BATCH: usize = 1024;

fn quadruples(n: usize, seed: u64) -> Vec<[Point; 4]> {
    let mut rng = Rng::new(seed);
    (0..n)
        .map(|_| {
            let mut point = || [rng.next_f64() * 100.0, rng.next_f64() * 100.0];
            [point(), point(), point(), point()]
        })
        .collect()
}

#[divan::bench]
fn orient2d_batch(bencher: Bencher) {
    let data = quadruples(BATCH, 0x0121);
    bencher.bench(|| {
        let mut acc = 0.0;
        for [a, b, c, _] in divan::black_box(&data) {
            acc += orient2d(a, b, c);
        }
        divan::black_box(acc)
    });
}

#[divan::bench]
fn incircle_batch(bencher: Bencher) {
    let data = quadruples(BATCH, 0x0122);
    bencher.bench(|| {
        let mut acc = 0.0;
        for [a, b, c, d] in divan::black_box(&data) {
            acc += incircle(a, b, c, d);
        }
        divan::black_box(acc)
    });
}

#[divan::bench]
fn point_in_triangle_batch(bencher: Bencher) {
    let data = quadruples(BATCH, 0x0123);
    bencher.bench(|| {
        let mut hits = 0usize;
        for [p, a, b, c] in divan::black_box(&data) {
            hits += point_in_triangle(p, a, b, c) as usize;
        }
        divan::black_box(hits)
    });
}

#[divan::bench]
fn segments_intersect_batch(bencher: Bencher) {
    let data = quadruples(BATCH, 0x0124);
    bencher.bench(|| {
        let mut hits = 0usize;
        for [a, b, c, d] in divan::black_box(&data) {
            hits += usize::from(segments_intersect(a, b, c, d));
        }
        divan::black_box(hits)
    });
}

/// Point-in-polygon against a convex polygon with many vertices.
#[divan::bench(args = [16, 256])]
fn point_inside_polygon(bencher: Bencher, vertices: usize) {
    let polygon = common::circle_points(vertices, 50.0);
    let queries = uniform_points(256, 0x0125);
    bencher.bench(|| {
        let mut hits = 0usize;
        for query in divan::black_box(&queries) {
            hits += usize::from(is_point_inside_polygon(query, divan::black_box(&polygon)));
        }
        divan::black_box(hits)
    });
}

/// Linear scan used to locate the triangle containing a point during insertion.
#[divan::bench]
fn find_containing_triangle_scan(bencher: Bencher) {
    let points = uniform_points(256, 0x5EED);
    let triangulation = triangulate(&points);
    let queries = uniform_points(64, 0x0126);
    bencher.bench(|| {
        let mut found = 0usize;
        for query in divan::black_box(&queries) {
            found += usize::from(matches!(
                find_containing_triangle(&triangulation, query),
                pycdt_rs::types::PointLocation::Interior(_)
            ));
        }
        divan::black_box(found)
    });
}

/// Rebuild closed polygons out of an unordered edge soup.
#[divan::bench(args = [4, 32])]
fn build_polygons(bencher: Bencher, loops: usize) {
    let mut edges = Vec::with_capacity(loops * 8);
    for l in 0..loops {
        let base = l * 8;
        for i in 0..8 {
            edges.push((base + i, base + (i + 1) % 8));
        }
    }
    bencher.bench(|| divan::black_box(build_polygons_from_edges(divan::black_box(&edges))));
}

#[divan::bench]
fn polygon_area_large(bencher: Bencher) {
    let points = common::circle_points(1024, 50.0);
    let polygon: Vec<usize> = (0..points.len()).collect();
    bencher.bench(|| divan::black_box(polygon_area(divan::black_box(&points), &polygon)));
}
