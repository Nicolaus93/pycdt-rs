use pycdt_rs::build::triangulate;
use pycdt_rs::types::NO_NEIGHBOR;
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::time::Instant;
fn convex_hull_size(points: &[[f64; 2]]) -> usize {
    let mut sorted = points.to_vec();
    sorted.sort_unstable_by(|a, b| a[0].total_cmp(&b[0]).then(a[1].total_cmp(&b[1])));
    sorted.dedup();
    let turn = |a: [f64; 2], b: [f64; 2], c: [f64; 2]| {
        robust::orient2d(
            robust::Coord { x: a[0], y: a[1] },
            robust::Coord { x: b[0], y: b[1] },
            robust::Coord { x: c[0], y: c[1] },
        )
    };
    let build_half = |iter: Box<dyn Iterator<Item = [f64; 2]>>| {
        let mut hull = Vec::new();
        for point in iter {
            while hull.len() >= 2 && turn(hull[hull.len() - 2], hull[hull.len() - 1], point) <= 0.0
            {
                hull.pop();
            }
            hull.push(point);
        }
        hull.len()
    };
    build_half(Box::new(sorted.iter().copied()))
        + build_half(Box::new(sorted.iter().rev().copied()))
        - 2
}

fn main() {
    let mut args = std::env::args().skip(1);
    let point_count = args
        .next()
        .map(|arg| arg.parse().expect("point count must be an integer"))
        .unwrap_or(100_000);
    let measured_runs = args
        .next()
        .map(|arg| arg.parse().expect("run count must be an integer"))
        .unwrap_or(9);

    let mut rng = StdRng::seed_from_u64(0x5EED);
    let points: Vec<[f64; 2]> = (0..point_count)
        .map(|_| [rng.gen::<f64>() * 1000.0, rng.gen::<f64>() * 1000.0])
        .collect();

    // Warm code and data paths before collecting wall-clock samples.
    for _ in 0..2 {
        std::hint::black_box(triangulate(std::hint::black_box(&points)));
    }

    let mut samples_ms = Vec::with_capacity(measured_runs);
    let mut triangles = 0;
    let mut boundary_edges = 0;
    for _ in 0..measured_runs {
        let started = Instant::now();
        let triangulation = triangulate(std::hint::black_box(&points));
        samples_ms.push(started.elapsed().as_secs_f64() * 1000.0);
        triangles = triangulation.num_triangles();
        boundary_edges = triangulation
            .triangle_neighbors
            .iter()
            .flatten()
            .filter(|&&neighbor| neighbor == NO_NEIGHBOR)
            .count();
        std::hint::black_box(triangulation);
    }

    samples_ms.sort_by(f64::total_cmp);
    let median = samples_ms[samples_ms.len() / 2];
    println!("points={point_count}");
    println!("runs={measured_runs}");
    println!("median_ms={median:.3}");
    println!("min_ms={:.3}", samples_ms[0]);
    println!("max_ms={:.3}", samples_ms[samples_ms.len() - 1]);
    println!("triangles={triangles}");
    println!("boundary_edges={boundary_edges}");
    println!("expected_hull_edges={}", convex_hull_size(&points));
    println!(
        "expected_triangles={}",
        2 * point_count - 2 - convex_hull_size(&points)
    );
    println!("samples_ms={samples_ms:?}");
}
