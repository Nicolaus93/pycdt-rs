"""Basic triangulation examples using pycdt_rs.

Demonstrates: triangulate, remove_super_triangle, update_triangulation.
"""

import numpy as np
import pycdt_rs


def circle_example():
    """Triangulate points on a circle."""
    arr = np.array(
        [
            [24.311, -7.358],
            [23.574, -9.456],
            [22.657, -11.481],
            [21.566, -13.419],
            [20.31, -15.253],
            [18.899, -16.971],
            [17.342, -18.558],
            [15.653, -20.004],
            [13.844, -21.296],
            [11.929, -22.425],
            [6.53, -24.546],
            [0.791, -25.388],
            [-4.989, -24.905],
            [-10.509, -23.124],
            [-15.48, -20.137],
            [-19.645, -16.1],
            [-22.786, -11.224],
            [-24.738, -5.762],
            [-25.4, 0.0],
            [-23.868, 8.687],
            [-19.458, 16.327],
            [-12.7, 21.997],
            [-4.411, 25.014],
            [4.411, 25.014],
            [12.7, 21.997],
            [19.458, 16.327],
            [23.868, 8.687],
            [25.4, 0.0],
            [25.386, -0.829],
            [25.346, -1.658],
            [25.278, -2.485],
            [25.184, -3.309],
            [25.062, -4.129],
            [24.914, -4.945],
            [24.739, -5.756],
            [24.538, -6.561],
        ]
    )

    tri = pycdt_rs.triangulate(arr)
    print(f"Circle: {tri.num_points} points, {tri.num_triangles} triangles")
    print(f"  points shape: {tri.points.shape}")
    print(f"  triangle_vertices shape: {tri.triangle_vertices.shape}")

    # Remove the super triangle to get only "real" triangles
    pycdt_rs.remove_super_triangle(tri)
    print(f"  After remove_super_triangle: {tri.num_triangles} triangles")


def concave_polygon_example():
    """Triangulate an L-shaped concave polygon."""
    points = np.array(
        [
            [0.0, 0.0],
            [4.0, 0.0],
            [4.0, 4.0],
            [3.0, 4.0],
            [3.0, 1.0],
            [1.0, 1.0],
            [1.0, 4.0],
            [0.0, 4.0],
        ],
        dtype=np.float64,
    )

    tri = pycdt_rs.triangulate(points)
    print(f"\nConcave polygon: {tri.num_points} points, {tri.num_triangles} triangles")

    pycdt_rs.remove_super_triangle(tri)
    print(f"  After remove_super_triangle: {tri.num_triangles} triangles")
    print(f"  Triangle vertices:\n{tri.triangle_vertices}")


def update_example():
    """Demonstrate incremental point insertion.

    Note: update_triangulation must be called BEFORE remove_super_triangle,
    since new points need the super triangle for insertion.
    """
    initial = np.array(
        [[0.0, 0.0], [10.0, 0.0], [5.0, 10.0], [10.0, 10.0], [0.0, 10.0]],
        dtype=np.float64,
    )
    tri = pycdt_rs.triangulate(initial)
    print(f"\nUpdate: initial {tri.num_points} points, {tri.num_triangles} triangles")

    new_pts = np.array([[5.0, 3.0], [3.0, 7.0], [7.0, 7.0]], dtype=np.float64)
    pycdt_rs.update_triangulation(tri, new_pts)
    print(f"  After adding 3 points: {tri.num_points} points, {tri.num_triangles} triangles")

    pycdt_rs.remove_super_triangle(tri)
    print(f"  After remove_super_triangle: {tri.num_triangles} triangles")


if __name__ == "__main__":
    circle_example()
    concave_polygon_example()
    update_example()
