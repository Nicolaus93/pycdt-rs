use criterion::{black_box, criterion_group, criterion_main, BatchSize, Criterion};
use pycdt_rs::{build::triangulate, constrained::add_constraints, Triangulation};
use rand::{rngs::StdRng, Rng, SeedableRng};

const RANDOM_SEED: u64 = 0x05ee_dc0d_ed15_ca11;
// The 100,000-point case is temporarily disabled because its setup makes the
// simulated CI benchmark prohibitively slow.
const RANDOM_POINT_COUNTS: [usize; 2] = [1_000, 10_000];

fn random_points(count: usize, seed: u64) -> Vec<[f64; 2]> {
    assert!(count >= 2);
    let mut generator = StdRng::seed_from_u64(seed);
    let mut points = Vec::with_capacity(count);
    // All random points lie strictly above and between the fixed endpoints, so
    // (0,1) is a guaranteed physical hull edge for every fixture size.
    points.push([0.0, 0.0]);
    points.push([1.0, 0.0]);
    points.extend((2..count).map(|_| {
        [
            generator.gen_range(0.001..0.999),
            generator.gen_range(0.001..0.999),
        ]
    }));
    points
}

fn grid(width: usize, height: usize, perturb: bool) -> Vec<[f64; 2]> {
    (0..height)
        .flat_map(|y| {
            (0..width).map(move |x| {
                let offset = if perturb && x > 0 && x + 1 < width && y > 0 && y + 1 < height {
                    (((x * 17 + y * 31) % 11) as f64 - 5.0) * 0.013
                } else {
                    0.0
                };
                [x as f64, y as f64 + offset]
            })
        })
        .collect()
}

fn insertion_benchmark(
    group: &mut criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>,
    name: &str,
    base: &Triangulation,
    constraints: &[(usize, usize)],
) {
    group.bench_function(name, |bencher| {
        bencher.iter_batched(
            || base.clone(),
            |mut triangulation| {
                let inserted = add_constraints(&mut triangulation, black_box(constraints));
                assert!(inserted, "benchmark constraints must be valid");
                black_box(triangulation);
            },
            // Cloning is fixture setup and deliberately excluded from the timed region.
            BatchSize::SmallInput,
        )
    });
}

fn constraint_insertion(c: &mut Criterion) {
    let mut group = c.benchmark_group("constraint_insertion");

    // The perturbed rows avoid collinear intermediate vertices, producing a
    // long corridor whose crossed edges must be flipped and restored.
    let width = 24;
    let height = 24;
    let perturbed = triangulate(&grid(width, height, true));
    let long_constraint = [(width * 5, width * 18 + width - 1)];
    insertion_benchmark(
        &mut group,
        "long_corridor_24x24",
        &perturbed,
        &long_constraint,
    );

    // Several disjoint corridors share one add_constraints call, covering the
    // batch-lifetime incidence map, queues, and generation-marker storage.
    let batched_constraints = [
        (width * 3, width * 7 + width - 1),
        (width * 9, width * 12 + width - 1),
        (width * 14, width * 17 + width - 1),
        (width * 19, width * 21 + width - 1),
    ];
    insertion_benchmark(
        &mut group,
        "four_disjoint_corridors_24x24",
        &perturbed,
        &batched_constraints,
    );

    // A regular-grid diagonal is split into its physical constrained subedges
    // at every existing vertex on the segment.
    let regular = triangulate(&grid(width, height, false));
    let through_vertices = [(0, width * height - 1)];
    insertion_benchmark(
        &mut group,
        "through_existing_vertices_24x24",
        &regular,
        &through_vertices,
    );

    group.finish();
}

fn random_constraint_insertion(c: &mut Criterion) {
    let mut group = c.benchmark_group("constraint_insertion_random_points");

    for point_count in RANDOM_POINT_COUNTS {
        let points = random_points(point_count, RANDOM_SEED);
        let base = std::sync::OnceLock::new();
        let constraints = [(0, 1)];

        group.bench_function(point_count.to_string(), |bencher| {
            bencher.iter_batched(
                || base.get_or_init(|| triangulate(&points)).clone(),
                |mut triangulation| {
                    assert!(add_constraints(&mut triangulation, black_box(&constraints)));
                    black_box(triangulation);
                },
                BatchSize::SmallInput,
            )
        });
    }

    group.finish();
}

fn triangulation_construction(c: &mut Criterion) {
    let mut group = c.benchmark_group("triangulation_random_points");

    for point_count in RANDOM_POINT_COUNTS {
        // Point generation is deterministic fixture setup. Only construction
        // of the complete triangulation is inside the measured region.
        let points = random_points(point_count, RANDOM_SEED);
        group.bench_function(point_count.to_string(), |bencher| {
            bencher.iter(|| black_box(triangulate(black_box(&points))))
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    constraint_insertion,
    random_constraint_insertion,
    triangulation_construction
);
criterion_main!(benches);
