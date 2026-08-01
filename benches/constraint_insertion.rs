use criterion::{black_box, criterion_group, criterion_main, BatchSize, Criterion};
use pycdt_rs::{build::triangulate, constrained::add_constraints, Triangulation};

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

criterion_group!(benches, constraint_insertion);
criterion_main!(benches);
