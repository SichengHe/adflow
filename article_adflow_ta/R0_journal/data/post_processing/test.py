#!/usr/bin/env python
"""Prototype script: keep field as cell data, interpolate to a dense structured grid, and plot with Matplotlib."""
from __future__ import annotations

from pathlib import Path
from typing import Tuple

import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv

# ---------------------------------------------------------------------------
# User controls
# ---------------------------------------------------------------------------
SURFACE_CGNS_PATH = Path(
    r"C:\Users\Rohit\Desktop\article_resolvent_wing\data\resolvent_analysis_cylinder\output_cyl\cyl_000_surf.cgns"
)
BLOCK_PATH = ("BaseSurfaceSol", "SymmetryBCZone3", "Internal")  # nested block hierarchy
FIELD_NAME = "Velocity"
FIELD_COMPONENT = 0  # None -> magnitude
GRID_RESOLUTION = 400  # number of samples per axis
X_LIMITS = (-3.0, 3.0)
Y_LIMITS = (-3.0, 3.0)
COLORMAP = "RdBu_r"
LEVELS = 300
OUTPUT_FILE = Path("figures") / "test_cell_interpolated.pdf"
SHOW_FIG = True
CYLINDER_RADIUS = 0.5
CYLINDER_EDGE_COLOR = "black"
CYLINDER_EDGE_WIDTH = 2.0
# ---------------------------------------------------------------------------


def _get_nested_block(root: pv.MultiBlock, path: Tuple[str, ...]) -> pv.DataSet:
    block = root
    for key in path:
        if not isinstance(block, pv.MultiBlock):
            raise ValueError(f"Expected MultiBlock while traversing '{key}', got {type(block)}")
        if key not in block.keys():
            raise KeyError(f"Key '{key}' not found in block: available keys = {list(block.keys())}")
        block = block[key]
    if isinstance(block, pv.MultiBlock):
        # Use the first non-null leaf
        for child in block:
            if child is not None:
                block = child
                break
    if not isinstance(block, pv.DataSet):
        raise TypeError(f"Final block is not a dataset: {type(block)}")
    return block


def interpolate_cell_field(block: pv.DataSet) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Sample the cell-centered field onto a structured grid."""
    surface = block.extract_surface().copy()
    if FIELD_NAME not in surface.cell_data:
        raise KeyError(f"Field '{FIELD_NAME}' not found in cell_data.")
    values = np.asarray(surface.cell_data[FIELD_NAME])
    if values.ndim == 2:
        if FIELD_COMPONENT is None:
            values = np.linalg.norm(values, axis=1)
        else:
            values = values[:, FIELD_COMPONENT]

    x = np.linspace(X_LIMITS[0], X_LIMITS[1], GRID_RESOLUTION)
    y = np.linspace(Y_LIMITS[0], Y_LIMITS[1], GRID_RESOLUTION)
    xv, yv = np.meshgrid(x, y, indexing="xy")
    zv = np.zeros_like(xv)
    plane = pv.StructuredGrid(xv, yv, zv)

    tmp_field = FIELD_NAME if values.ndim == 1 else f"{FIELD_NAME}_scalar"
    if tmp_field in surface.cell_data:
        surface.cell_data.remove(tmp_field)
    surface.cell_data[tmp_field] = values
    sampled = plane.sample(surface)
    if tmp_field not in sampled.point_data:
        raise KeyError(f"Sampled grid missing '{tmp_field}' in point_data.")
    dense = np.asarray(sampled.point_data[tmp_field]).reshape(xv.shape)
    return xv, yv, dense


def plot_field(xv: np.ndarray, yv: np.ndarray, values: np.ndarray) -> None:
    fig, ax = plt.subplots(figsize=(6, 5), dpi=300)
    contour = ax.contourf(xv, yv, values, levels=LEVELS, cmap=COLORMAP)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"{FIELD_NAME}[{FIELD_COMPONENT if FIELD_COMPONENT is not None else '|.|'}]")
    circle = plt.Circle((0.0, 0.0), CYLINDER_RADIUS, color=CYLINDER_EDGE_COLOR, fill=False, linewidth=CYLINDER_EDGE_WIDTH)
    ax.add_patch(circle)
    fig.colorbar(contour, ax=ax)
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_FILE, bbox_inches="tight")
    if SHOW_FIG:
        plt.show()
    else:
        plt.close(fig)


def main() -> None:
    if not SURFACE_CGNS_PATH.exists():
        raise FileNotFoundError(SURFACE_CGNS_PATH)
    dataset = pv.read(SURFACE_CGNS_PATH)
    block = _get_nested_block(dataset, BLOCK_PATH)
    xv, yv, dense_values = interpolate_cell_field(block)
    plot_field(xv, yv, dense_values)


if __name__ == "__main__":
    main()
