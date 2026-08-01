//! Shared deterministic point generators used by the benchmark targets.
//!
//! The generators avoid any external RNG dependency so the benchmark inputs are
//! byte-for-byte reproducible across machines and CI runs.

// Each benchmark target includes this module and only uses part of it.
#![allow(dead_code)]

use pycdt_rs::types::Point;

/// Small xorshift64* generator: deterministic and dependency free.
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Self(seed | 1)
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// Uniform value in `[0, 1)`.
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// `n` uniformly distributed points inside a 1000x1000 box.
pub fn uniform_points(n: usize, seed: u64) -> Vec<Point> {
    uniform_points_in(n, seed, 0.0, 1000.0)
}

/// `n` uniformly distributed points inside the `[min, max]^2` box.
pub fn uniform_points_in(n: usize, seed: u64, min: f64, max: f64) -> Vec<Point> {
    let mut rng = Rng::new(seed);
    let span = max - min;
    (0..n)
        .map(|_| [min + rng.next_f64() * span, min + rng.next_f64() * span])
        .collect()
}

/// A `side x side` regular grid. Highly degenerate input: many collinear points
/// and cocircular quadruples, which stresses the robust predicates and the
/// on-edge insertion path.
pub fn grid_points(side: usize) -> Vec<Point> {
    let mut points = Vec::with_capacity(side * side);
    for i in 0..side {
        for j in 0..side {
            points.push([i as f64, j as f64]);
        }
    }
    points
}

/// `n` points evenly spaced on a circle: every point is cocircular with every
/// other one, the worst case for the in-circle test.
pub fn circle_points(n: usize, radius: f64) -> Vec<Point> {
    (0..n)
        .map(|i| {
            let angle = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
            [radius * angle.cos(), radius * angle.sin()]
        })
        .collect()
}

/// A square outer boundary with a `side x side` arrangement of square holes
/// punched in it, plus `filler` uniformly distributed interior points.
///
/// Returns the points and the constrained edges describing the boundary and the
/// holes, ready to be passed to `remove_holes_by_edges`.
pub fn polygon_with_holes(side: usize, filler: usize) -> (Vec<Point>, Vec<(usize, usize)>) {
    let extent = 10.0 * side as f64;

    let mut points: Vec<Point> = vec![[0.0, 0.0], [extent, 0.0], [extent, extent], [0.0, extent]];
    let mut edges: Vec<(usize, usize)> = vec![(0, 1), (1, 2), (2, 3), (3, 0)];

    for row in 0..side {
        for col in 0..side {
            let x0 = 10.0 * col as f64 + 3.0;
            let y0 = 10.0 * row as f64 + 3.0;
            let base = points.len();
            points.push([x0, y0]);
            points.push([x0 + 4.0, y0]);
            points.push([x0 + 4.0, y0 + 4.0]);
            points.push([x0, y0 + 4.0]);
            edges.push((base, base + 1));
            edges.push((base + 1, base + 2));
            edges.push((base + 2, base + 3));
            edges.push((base + 3, base));
        }
    }

    let mut rng = Rng::new(0xC0D5_5EED);
    for _ in 0..filler {
        points.push([rng.next_f64() * extent, rng.next_f64() * extent]);
    }

    (points, edges)
}
