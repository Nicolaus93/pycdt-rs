import numpy as np
import matplotlib.pyplot as plt
import pycdt_rs


def plot_triangulation(tri, title="Triangulation"):
    points = tri.points
    triangles = tri.triangle_vertices

    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    ax.set_title(title)
    ax.set_aspect("equal")
    ax.triplot(points[:, 0], points[:, 1], triangles, linewidth=0.5, color="black")
    ax.plot(points[:, 0], points[:, 1], "o", markersize=3, color="red")
    plt.tight_layout()
    plt.savefig("circle.png", dpi=150)
    plt.show()


if __name__ == "__main__":
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
    pycdt_rs.remove_super_triangle(tri)
    plot_triangulation(tri, title="Circle")
