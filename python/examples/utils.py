"""Shared plotting and utility functions for pycdt_rs examples."""

import pycdt_rs


def plot_triangulation(tri, title="Triangulation", filename=None, figsize=(8, 8), markersize=3):
    points = tri.points
    triangles = tri.triangle_vertices

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax.set_title(title)
    ax.set_aspect("equal")
    ax.triplot(points[:, 0], points[:, 1], triangles, linewidth=0.5, color="black")
    ax.plot(points[:, 0], points[:, 1], "o", markersize=markersize, color="red")
    plt.tight_layout()
    if filename:
        plt.savefig(filename, dpi=150)
    if plt.get_backend().lower() != "agg":
        plt.show()


def remove_holes_by_edges(tri, constrained_edges):
    pycdt_rs.remove_holes_by_edges(
        tri,
        [(int(edge[0]), int(edge[1])) for edge in constrained_edges],
    )
