"""Constrained Delaunay triangulation examples using pycdt_rs.

Demonstrates: triangulate, add_constraints, remove_super_triangle.
"""

import numpy as np
import pycdt_rs


def grid_diagonal_constraint():
    """Insert a long diagonal constraint across a grid."""
    print("=" * 60)
    print("GRID DIAGONAL CONSTRAINT")
    print("=" * 60)

    # Create a 6x6 grid
    n = 6
    x = np.linspace(0, 5, n)
    y = np.linspace(0, 5, n)
    xx, yy = np.meshgrid(x, y)
    points = np.column_stack([xx.ravel(), yy.ravel()])

    tri = pycdt_rs.triangulate(points)
    print(f"Points: {tri.num_points}, Triangles: {tri.num_triangles}")

    # Diagonal constraint from top-left to bottom-right
    p_idx, q_idx = 30, 5
    print(f"Adding constraint {p_idx} -> {q_idx}")
    print(f"  From {points[p_idx]} to {points[q_idx]}")

    pycdt_rs.add_constraints(tri, [(p_idx, q_idx)])
    print(f"After constraint: {tri.num_triangles} triangles")
    print(f"Constrained edges: {tri.constrained_edges}")


def polygon_boundary_constraints():
    """Add multiple constraints forming an L-shaped boundary."""
    print("\n" + "=" * 60)
    print("POLYGON BOUNDARY CONSTRAINTS")
    print("=" * 60)

    boundary = np.array(
        [
            [0.0, 0.0],  # 0
            [3.0, 0.0],  # 1
            [3.0, 2.0],  # 2
            [1.0, 2.0],  # 3
            [1.0, 3.0],  # 4
            [0.0, 3.0],  # 5
            [2.0, 0.0],  # 6
            [1.0, 0.0],  # 7
        ]
    )

    interior = np.array(
        [
            [0.5, 0.5],
            [2.5, 0.5],
            [2.0, 1.5],
            [0.5, 2.5],
        ]
    )

    points = np.vstack([boundary, interior])
    tri = pycdt_rs.triangulate(points)
    print(f"Points: {tri.num_points}, Triangles: {tri.num_triangles}")

    constraints = [
        (5, 0),
        (1, 2),
        (5, 3),
        (0, 1),
        (2, 3),
        (0, 7),
        (7, 6),
        (1, 2),
    ]
    print(f"Adding {len(constraints)} constraints...")
    pycdt_rs.add_constraints(tri, constraints)

    print(f"After constraints: {tri.num_triangles} triangles")
    print(f"Constrained edges: {tri.constrained_edges}")


def build_polygons_example():
    """Demonstrate build_polygons_from_edges utility."""
    print("\n" + "=" * 60)
    print("BUILD POLYGONS FROM EDGES")
    print("=" * 60)

    points = np.array(
        [
            [0.0, 0.0],
            [3.0, 0.0],
            [3.0, 2.0],
            [1.0, 2.0],
            [1.0, 3.0],
            [0.0, 3.0],
        ],
        dtype=np.float64,
    )

    tri = pycdt_rs.triangulate(points)

    # Define edges that form a closed polygon
    edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]

    polygons = pycdt_rs.build_polygons_from_edges(tri, edges)
    print(f"Found {len(polygons)} polygon(s):")
    for i, poly in enumerate(polygons):
        print(f"  Polygon {i}: vertices {poly}")


if __name__ == "__main__":
    grid_diagonal_constraint()
    polygon_boundary_constraints()
    build_polygons_example()
