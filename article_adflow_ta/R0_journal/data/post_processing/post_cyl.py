#!/usr/bin/env python
"""
Post-processing utilities for the cylinder resolvent case.

Step 1: Inspect the multiblock CGNS surface file (`out_surf.cgns`) and print
        a concise tree of the available datasets and fields.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Tuple

import numpy as np
import re
import matplotlib.pyplot as plt
from matplotlib import tri as mtri
from matplotlib.patches import Circle

try:
    import niceplots  # type: ignore
except ImportError:  # pragma: no cover
    niceplots = None

try:
    import pyvista as pv
except ImportError:  # pragma: no cover - runtime dependency
    pv = None


# Edit this path if you want a different default CGNS file without passing CLI args.
DEFAULT_SURFACE_PATH = Path(
    r"/home/rohit/Desktop/resolvent_adflow/article_resolvent_wing/data/resolvent_analysis_cylinder/output_cyl/out_surf_forcing.cgns"
)
if not DEFAULT_SURFACE_PATH.exists():
    # Fallback to relative path (keeps script working if repo is moved).
    DEFAULT_SURFACE_PATH = (
        Path(__file__)
        .resolve()
        .parents[1]
        / "resolvent_analysis_cylinder"
        / "output_cyl"
        / "cyl_000_surf.cgns"
    )

# ---------------------------------------------------------------------------
# Plot/visualization controls (edit these constants instead of hunting below)
# ---------------------------------------------------------------------------
SURFACE_CGNS_PATH = DEFAULT_SURFACE_PATH
BLOCK_FILTER = "SymmetryBCZone3"
OUTPUT_DIR = Path("figures") / "cylinder_post"
USE_TEX = True  # toggle LaTeX rendering globally
SHOW_DATASET_SUMMARY = True
LIST_FIELDS_ON_LOAD = True
PERFORM_PLOT = True
SHOW_PLOT_WINDOWS = True
FIELD_TO_PLOT = "DeltaVelYSurf"
FIELD_COMPONENT_INDEX: int | None = 1  # None -> magnitude for vector/tensor
OUTPUT_FILENAME_OVERRIDE = None
MANUAL_COLORBAR_LIMITS = None  # overrides COLORBAR_LIMITS if not None
USE_NICEPLOTS_STYLE = True
FONT_FAMILY = "CMU Serif"
AXES_LABELSIZE = 20
TICK_LABELSIZE = 16
LEGEND_LABELSIZE = 18
TITLE_FONT_SIZE = 22
AXIS_LABELS = (r"$x$", r"$y$")
Y_LABEL_ROTATION = 0
Y_LABEL_PAD = None  # e.g., 30
FIGURE_SIZE = (10, 4)
FIGURE_DPI = 300
PLOT_LEVELS = 50
PLOT_COLORMAP = "RdBu_r"
PLOT_COLORMAP_EXTEND = "both"
FORCE_EQUAL_ASPECT = True
X_LIMITS = (-10,10)#(-1,7)  # e.g., (0.0, 1.0) or None to auto
Y_LIMITS = (-2,2)

COLORBAR_ORIENTATION = "vertical"
COLORBAR_FRACTION = 0.04
COLORBAR_PAD = 0.05
COLORBAR_LABELPAD = 10
COLORBAR_LABEL_ROTATION = 0
COLORBAR_TICKSIZE = 16
COLORBAR_ANCHOR = None  # e.g., (0.5, -0.1) or None to let Matplotlib decide
COLORBAR_LABEL_TEMPLATE = r"$\delta \mathbf{{v}}$"
COLORBAR_LIMITS = (-0.05 , 0.05) #None  # e.g., (-0.1, 0.1)
COLORBAR_NUM_TICKS = 5
COLORBAR_TICK_FORMAT = "{:.4f}"

CYLINDER_OUTLINE_ENABLED = True
CYLINDER_CENTER = (0.0, 0.0)
CYLINDER_RADIUS = 0.5  # Diameter = 1
CYLINDER_EDGE_COLOR = "black"
CYLINDER_EDGE_WIDTH = 2.0

SURFACE_SUBDIVIDE_LEVELS = 1  # >0 subdivides triangles to smooth contours

TITLE_TEMPLATE = ""
OUTPUT_FILENAME_TEMPLATE = "{field_label_safe}_forcing.pdf"
# ---------------------------------------------------------------------------


def _collect_blocks(dataset: pv.DataSet, prefix: str = "") -> Iterable[Tuple[str, pv.DataSet]]:
    """Yield (path, block) pairs while drilling through a MultiBlock hierarchy."""
    if isinstance(dataset, pv.MultiBlock):
        keys = dataset.keys()
        for idx, child in enumerate(dataset):
            name = keys[idx] or f"block_{idx}"
            child_prefix = f"{prefix}/{name}" if prefix else name
            yield from _collect_blocks(child, child_prefix)
    else:
        yield prefix or dataset.__class__.__name__, dataset


def _summarize_arrays(data_dict: pv.DataSetAttributes) -> str:
    """Return a short description of arrays inside a vtkDataSetAttributes object."""
    names = list(data_dict.keys())
    if not names:
        return "none"
    bits = []
    for name in names:
        array = np.asarray(data_dict[name])
        shape = "x".join(str(dim) for dim in array.shape)
        bits.append(f"{name} (shape={shape}, dtype={array.dtype})")
    return "; ".join(bits)


def describe_surface_file(dataset: pv.DataSet, path: Path) -> None:
    """Print dataset information."""
    print("=" * 80)
    print("Cylinder Resolvent Surface Inspector")
    print("=" * 80)
    print(f"CGNS file: {path}")
    print()
    blocks = list(_collect_blocks(dataset))
    print(f"Detected {len(blocks)} drawable block(s):")
    for idx, (name, block) in enumerate(blocks):
        if block is None:
            print(f"  [{idx:02d}] {name}: <empty block>")
            continue
        n_points = getattr(block, "n_points", 0)
        n_cells = getattr(block, "n_cells", 0)
        print(f"  [{idx:02d}] {name}: {block.__class__.__name__} "
              f"(points={n_points}, cells={n_cells})")
        # Provide a quick look at the array names to guide later plotting.
        if n_points:
            point_info = _summarize_arrays(block.point_data)
            print(f"       point_data: {point_info}")
        if n_cells:
            cell_info = _summarize_arrays(block.cell_data)
            print(f"       cell_data: {cell_info}")


def _apply_plot_style(use_tex: bool) -> None:
    """Configure Matplotlib to mimic visualize_aircraft_surface_data.py aesthetics."""
    if USE_NICEPLOTS_STYLE and niceplots is not None:
        plt.style.use(niceplots.get_style())
    plt.rcParams.update(
        {
            "axes.unicode_minus": False,
            "font.family": FONT_FAMILY,
            "text.usetex": use_tex,
            "mathtext.fontset": "cm",
            "axes.labelsize": AXES_LABELSIZE,
            "xtick.labelsize": TICK_LABELSIZE,
            "ytick.labelsize": TICK_LABELSIZE,
            "legend.fontsize": LEGEND_LABELSIZE,
        }
    )
    if use_tex:
        plt.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"


def _select_block(dataset: pv.DataSet, keyword: str) -> Tuple[str, pv.DataSet]:
    """Return the first block whose path contains `keyword`."""
    blocks = list(_collect_blocks(dataset))
    for name, block in blocks:
        if keyword.lower() in name.lower():
            return name, block
    raise ValueError(f"Could not find block containing '{keyword}'.")


def _prepare_surface(block: pv.DataSet) -> pv.PolyData:
    """Convert a structured block to a triangulated surface."""
    surface = block
    if not isinstance(surface, pv.PolyData):
        surface = block.extract_surface()
    surface = surface.triangulate()
    if SURFACE_SUBDIVIDE_LEVELS and SURFACE_SUBDIVIDE_LEVELS > 0:
        surface = surface.subdivide(SURFACE_SUBDIVIDE_LEVELS, subfilter="loop")
    return surface


def _triangulation_from_surface(surface: pv.PolyData) -> Tuple[mtri.Triangulation, np.ndarray]:
    """Return Matplotlib triangulation and point-data dictionary from a PyVista surface."""
    faces = surface.faces.reshape(-1, 4)[:, 1:]
    points = surface.points
    tri = mtri.Triangulation(points[:, 0], points[:, 1], triangles=faces)
    return tri, points


def _ensure_point_field(surface: pv.PolyData, field_name: str) -> Tuple[pv.PolyData, np.ndarray]:
    """Return a surface (possibly promoted) with the requested field on point_data."""
    if field_name in surface.point_data:
        return surface, np.asarray(surface.point_data[field_name])
    if field_name in surface.cell_data:
        promoted = surface.cell_data_to_point_data(pass_cell_data=True)
        if field_name not in promoted.point_data:
            raise KeyError(f"Field '{field_name}' not found after promotion to point data.")
        return promoted, np.asarray(promoted.point_data[field_name])
    raise KeyError(f"Field '{field_name}' not found in point_data or cell_data.")


def _prepare_field_values(values: np.ndarray, field_name: str, component: int | None) -> Tuple[np.ndarray, str]:
    """Select a component or magnitude for plotting and return display label."""
    arr = np.asarray(values)
    if arr.ndim == 1:
        return arr, field_name
    if arr.ndim == 2:
        n_comp = arr.shape[1]
        if component is None:
            magnitude = np.linalg.norm(arr, axis=1)
            label = f"|{field_name}|"
            return magnitude, label
        if component < 0 or component >= n_comp:
            raise ValueError(f"Component {component} out of bounds for field '{field_name}' (size {n_comp}).")
        return arr[:, component], f"{field_name}[{component}]"
    raise ValueError(f"Unsupported field shape for '{field_name}': {arr.shape}")


def _print_available_fields(block: pv.DataSet) -> None:
    """List plottable fields for the selected block."""
    def format_fields(data: pv.DataSetAttributes) -> str:
        if not data.keys():
            return "  (none)"
        lines = []
        for name in data.keys():
            array = np.asarray(data[name])
            lines.append(f"  - {name}: shape={array.shape}, dtype={array.dtype}")
        return "\n".join(lines)

    print("Available point_data fields:")
    print(format_fields(block.point_data))
    print("Available cell_data fields:")
    print(format_fields(block.cell_data))


def _slugify(text: str) -> str:
    """Generate a filesystem-friendly string."""
    safe = re.sub(r"[^0-9A-Za-z._-]+", "_", text).strip("_")
    return safe or "field"


def _compute_value_limits(values: np.ndarray, manual_limits):
    if manual_limits is not None:
        vmin, vmax = manual_limits
    else:
        vmin = float(np.nanmin(values))
        vmax = float(np.nanmax(values))
    if vmin == vmax:
        vmax = vmin + 1e-12
    return vmin, vmax


def _set_colorbar_ticks(cbar, vmin: float, vmax: float, num_ticks: int | None) -> None:
    if num_ticks is None or num_ticks < 2:
        ticks = np.asarray(cbar.get_ticks(), dtype=float)
        if ticks.size < 2:
            ticks = np.linspace(vmin, vmax, 5)
    else:
        ticks = np.linspace(vmin, vmax, num_ticks)
    ticks[0] = vmin
    ticks[-1] = vmax
    cbar.set_ticks(ticks)


def _plot_component(
    tri: mtri.Triangulation,
    values: np.ndarray,
    title: str,
    out_path: Path,
    colorbar_label: str,
    value_limits,
    num_ticks: int,
    use_tex: bool,
    show: bool = True,
) -> plt.Figure:
    """Create and save a tricontour plot for a single velocity component."""
    _apply_plot_style(use_tex=use_tex)
    fig, ax = plt.subplots(figsize=FIGURE_SIZE, dpi=FIGURE_DPI)
    vmin, vmax = _compute_value_limits(values, value_limits)
    if isinstance(PLOT_LEVELS, int):
        levels = np.linspace(vmin, vmax, PLOT_LEVELS)
    else:
        levels = PLOT_LEVELS
    contour = ax.tricontourf(
        tri,
        values,
        levels=levels,
        cmap=PLOT_COLORMAP,
        extend=PLOT_COLORMAP_EXTEND,
        vmin=vmin,
        vmax=vmax,
    )
    if FORCE_EQUAL_ASPECT:
        ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(AXIS_LABELS[0])
    if Y_LABEL_PAD is not None:
        ax.set_ylabel(AXIS_LABELS[1], rotation=Y_LABEL_ROTATION, labelpad=Y_LABEL_PAD)
    else:
        ax.set_ylabel(AXIS_LABELS[1], rotation=Y_LABEL_ROTATION)
    ax.set_title(title, fontsize=TITLE_FONT_SIZE)
    if X_LIMITS is not None:
        ax.set_xlim(*X_LIMITS)
    if Y_LIMITS is not None:
        ax.set_ylim(*Y_LIMITS)
    if CYLINDER_OUTLINE_ENABLED:
        circle = Circle(
            CYLINDER_CENTER,
            CYLINDER_RADIUS,
            fill=False,
            color=CYLINDER_EDGE_COLOR,
            linewidth=CYLINDER_EDGE_WIDTH,
        )
        ax.add_patch(circle)

    cbar_kwargs = dict(
        orientation=COLORBAR_ORIENTATION,
        fraction=COLORBAR_FRACTION,
        pad=COLORBAR_PAD,
    )
    if COLORBAR_ANCHOR is not None:
        cbar_kwargs["anchor"] = COLORBAR_ANCHOR
    cbar = fig.colorbar(contour, ax=ax, **cbar_kwargs)
    cbar.set_label(colorbar_label, rotation=COLORBAR_LABEL_ROTATION, labelpad=COLORBAR_LABELPAD)
    cbar.ax.tick_params(labelsize=COLORBAR_TICKSIZE)
    _set_colorbar_ticks(cbar, vmin, vmax, num_ticks)
    tick_labels = [COLORBAR_TICK_FORMAT.format(tick) for tick in cbar.get_ticks()]
    cbar.set_ticklabels(tick_labels)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    if not show:
        plt.close(fig)
    return fig


def plot_field(
    block_name: str,
    block: pv.DataSet,
    field_name: str,
    component: int | None,
    output_dir: Path,
    output_name: str | None,
    use_tex: bool,
    show: bool,
    manual_limits,
) -> None:
    """Extract block data and plot the requested field."""
    surface = _prepare_surface(block)
    surface, field_values = _ensure_point_field(surface, field_name)
    tri, _points = _triangulation_from_surface(surface)

    values, field_label = _prepare_field_values(field_values, field_name, component)

    title = TITLE_TEMPLATE.format(field_label=field_label, block=block_name)
    filename = output_name or OUTPUT_FILENAME_TEMPLATE.format(
        field_label_safe=_slugify(field_label)
    )
    component_label = component if component is not None else "mag"
    colorbar_label = COLORBAR_LABEL_TEMPLATE.format(
        field_label=field_label,
        field=field_name,
        component=component_label,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    _plot_component(
        tri,
        values,
        title=title,
        out_path=output_dir / filename,
        colorbar_label=colorbar_label,
        value_limits=manual_limits,
        num_ticks=COLORBAR_NUM_TICKS,
        use_tex=use_tex,
        show=show,
    )
    if show:
        plt.show()


def main() -> None:
    if pv is None:
        raise ImportError(
            "pyvista is required for CGNS inspection. "
            "Install it with `pip install pyvista` in your environment."
        )
    if not SURFACE_CGNS_PATH.exists():
        raise FileNotFoundError(f"Surface CGNS file not found: {SURFACE_CGNS_PATH}")

    dataset = pv.read(SURFACE_CGNS_PATH)
    if SHOW_DATASET_SUMMARY:
        describe_surface_file(dataset, SURFACE_CGNS_PATH)
    block_name, block = _select_block(dataset, BLOCK_FILTER)
    print(f"Selected block: {block_name}")

    if LIST_FIELDS_ON_LOAD or not PERFORM_PLOT or FIELD_TO_PLOT is None:
        _print_available_fields(block)
    if not PERFORM_PLOT or FIELD_TO_PLOT is None:
        return

    plot_field(
        block_name=block_name,
        block=block,
        field_name=FIELD_TO_PLOT,
        component=FIELD_COMPONENT_INDEX,
        output_dir=OUTPUT_DIR,
        output_name=OUTPUT_FILENAME_OVERRIDE,
        use_tex=USE_TEX,
        show=SHOW_PLOT_WINDOWS,
        manual_limits=MANUAL_COLORBAR_LIMITS if MANUAL_COLORBAR_LIMITS is not None else COLORBAR_LIMITS,
    )


if __name__ == "__main__":
    main()
