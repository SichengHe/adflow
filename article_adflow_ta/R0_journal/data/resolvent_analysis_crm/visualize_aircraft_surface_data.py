"""
CGNS surface viewer with:
  - menu of plottable arrays (indexed 0..N-1)
  - selection by --var <index>
  - Matplotlib 3D plotting
  - controllable Matplotlib camera: elev/azim/roll
  - controllable axis limits (auto_equal / auto_tight / manual)

Matplotlib 3D default camera (commonly seen by users):
  elev=30, azim=-60, roll=0
"""

from __future__ import annotations

import argparse
from collections import OrderedDict
from io import BytesIO
from pathlib import Path
from typing import Dict, List, Tuple
import warnings

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.transforms import Bbox
import numpy as np
import pyvista as pv
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

try:
    from PIL import Image
except ImportError:  # pragma: no cover - optional dependency
    Image = None

try:
    import niceplots
except ImportError:  # pragma: no cover - optional dependency
    niceplots = None


# -----------------------
# Defaults (IDE-friendly)
# -----------------------
CGNS_FILE = Path("out_surf_base_flow_plane_body.cgns")

# Variable selection
VAR_INDEX = 1          # 0, 1, 2 ... (overridden by --var)
PREFER = "cell"        # "cell" or "point" (controls listing order too. 
                        #  See console output for the plottable quantities.)

# Typography / style (inspired by plot_obj_func1_niceplots.py)
USE_NICEPLOTS_STYLE = True
FONT_FAMILY = "CMU Serif"
USE_TEX = True
LABEL_FONT_SIZE = 20
TICK_FONT_SIZE = 20
COLORBAR_LABEL_FONT_SIZE = 20

# Matplotlib camera controls (defaults most users see)
MPL_ELEV_DEFAULT = 90 # 30   # degrees
MPL_AZIM_DEFAULT = -90 #-60  # degrees
MPL_ROLL_DEFAULT = 0    # degrees
MPL_SHOW_AXES = False    # Toggle Matplotlib axes/pane rendering
MPL_SHOW_GRID = False    # Toggle Matplotlib background grid

# Colorbar controls
COLORBAR_MIN = 0#None
COLORBAR_MAX = 1# None
COLORBAR_ORIENTATION = "horizontal"  # "vertical" | "horizontal"
COLORBAR_SHRINK = 0.4
COLORBAR_ASPECT = 30.0
COLORBAR_PAD = 0.02
COLORBAR_ANCHOR = (0.5, -55.5)
COLORBAR_LABEL = r"$C_p$" #  r"$\delta \boldsymbol{\rho}$" #r"$\delta \mathbf{u}$" #None  # Defaults to array name when None. Use $xyz$ for latex
COLORBAR_LABELPAD = 12.0
COLORBAR_LABEL_ROTATION = 0
COLORBAR_FLUSH_END_TICKS = True
COLORBAR_NUM_TICKS = 5
COLORBAR_TICK_FORMAT = "%.5f"

# Final figure cropping (performed post-render via Pillow)
# Needs a little playing around but gets the job done
FIGURE_CROP_ENABLED = True
FIGURE_CROP_TRIM_LEFT = 0.2   # Fraction trimmed from left edge (0..1)
FIGURE_CROP_TRIM_RIGHT = 0.2  # Fraction trimmed from right edge (0..1)
FIGURE_CROP_TRIM_BOTTOM = 0.15
FIGURE_CROP_TRIM_TOP = 0.1
FIGURE_CROP_SAVE_PATH: Path | None = Path("cp_top_view.pdf") #None  # e.g., Path("cropped_visualization.png")
FIGURE_CROP_FORMAT = "pdf"  # "png", "jpeg", or "pdf"

# Axis limits / aspect controls
#   - "auto_equal": equal aspect cube limits (preserves geometry shape)
#   - "auto_tight": tight mins/maxs (not equalized)
#   - "manual": use X_LIM/Y_LIM/Z_LIM below
AX_LIMITS_MODE = "manual"  # "auto_equal" | "auto_tight" | "manual"

# Manual limits (used only when AX_LIMITS_MODE == "manual")
X_LIM = (5, 65.0)
Y_LIM = (0.0, 30.0)
Z_LIM = (0,1)

# If you are confused where to start for camera positioning, try these combos:
# https://matplotlib.org/stable/_images/api-toolkits-mplot3d-view_planes_3d_00_00.2x.png


_STYLE_APPLIED = False


def apply_global_plot_style() -> None:
    """Apply niceplots + typography controls once."""
    global _STYLE_APPLIED
    if _STYLE_APPLIED:
        return

    if USE_NICEPLOTS_STYLE and niceplots is not None:
        plt.style.use(niceplots.get_style())
    elif USE_NICEPLOTS_STYLE and niceplots is None:
        warnings.warn("niceplots is not available; using Matplotlib defaults.")

    plt.rcParams.update(
        {
            "axes.unicode_minus": False,
            "font.family": FONT_FAMILY,
            "font.size": LABEL_FONT_SIZE,
            "axes.labelsize": LABEL_FONT_SIZE,
            "xtick.labelsize": TICK_FONT_SIZE,
            "ytick.labelsize": TICK_FONT_SIZE,
        }
    )
    plt.rcParams["text.usetex"] = USE_TEX
    if USE_TEX:
        plt.rcParams["mathtext.fontset"] = "cm"
        plt.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"

    _STYLE_APPLIED = True


def _ensure_colorbar_end_ticks(
    cbar, vmin: float, vmax: float, num_ticks: int | None
) -> None:
    if not COLORBAR_FLUSH_END_TICKS:
        return
    ticks = np.asarray(cbar.get_ticks(), dtype=float)
    if num_ticks is not None and num_ticks >= 2:
        ticks = np.linspace(vmin, vmax, num_ticks)
    elif ticks.size < 2:
        fallback = num_ticks if num_ticks is not None and num_ticks >= 2 else 5
        ticks = np.linspace(vmin, vmax, fallback)
    ticks[0] = vmin
    ticks[-1] = vmax
    cbar.set_ticks(ticks)


def _trim_image(image) -> "Image.Image":
    """Trim edges from a PIL image using fractional settings."""
    width, height = image.size

    left_trim = float(np.clip(FIGURE_CROP_TRIM_LEFT, 0.0, 0.49))
    right_trim = float(np.clip(FIGURE_CROP_TRIM_RIGHT, 0.0, 0.49))
    top_trim = float(np.clip(FIGURE_CROP_TRIM_TOP, 0.0, 0.49))
    bottom_trim = float(np.clip(FIGURE_CROP_TRIM_BOTTOM, 0.0, 0.49))

    left_px = int(round(left_trim * width))
    right_px = int(round(width - right_trim * width))
    top_px = int(round(top_trim * height))
    bottom_px = int(round(height - bottom_trim * height))

    if right_px <= left_px:
        right_px = left_px + 1
    if bottom_px <= top_px:
        bottom_px = top_px + 1

    return image.crop((left_px, top_px, right_px, bottom_px))


def _finalize_figure(fig: plt.Figure) -> None:
    """Show the Matplotlib figure, optionally cropping via Pillow."""
    if not FIGURE_CROP_ENABLED:
        plt.show()
        return

    fmt = FIGURE_CROP_FORMAT.lower()
    if fmt not in {"png", "jpeg", "pdf"}:
        warnings.warn(
            f"Unsupported FIGURE_CROP_FORMAT '{FIGURE_CROP_FORMAT}'. "
            "Falling back to PNG."
        )
        fmt = "png"

    if fmt == "pdf":
        _save_cropped_pdf(fig)
        return

    if Image is None:
        warnings.warn(
            "FIGURE_CROP_ENABLED is True, but Pillow is not installed. "
            "Disable cropping or install pillow to enable post-processing."
        )
        plt.show()
        return

    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=fig.dpi, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    with Image.open(buf) as im:
        cropped = _trim_image(im)
        if FIGURE_CROP_SAVE_PATH:
            out_path = FIGURE_CROP_SAVE_PATH.expanduser()
            out_path.parent.mkdir(parents=True, exist_ok=True)
            save_image = cropped
            if fmt == "jpeg" and save_image.mode in {"RGBA", "P"}:
                save_image = save_image.convert("RGB")
            save_kwargs = {}
            save_kwargs["format"] = fmt.upper()
            save_image.save(out_path, **save_kwargs)
            print(f"Cropped figure saved to {out_path} ({save_kwargs['format']})")
        else:
            plt.figure(figsize=fig.get_size_inches())
            plt.imshow(cropped)
            plt.axis("off")
            plt.show()


def _save_cropped_pdf(fig: plt.Figure) -> None:
    """Save a vector PDF cropped using bbox operations."""
    if FIGURE_CROP_SAVE_PATH is None:
        warnings.warn(
            "FIGURE_CROP_SAVE_PATH must be set when exporting cropped PDFs. "
            "Showing figure without cropping."
        )
        plt.show()
        return

    bbox = fig.bbox_inches
    width = bbox.width
    height = bbox.height

    trim_left = float(np.clip(FIGURE_CROP_TRIM_LEFT, 0.0, 0.49))
    trim_right = float(np.clip(FIGURE_CROP_TRIM_RIGHT, 0.0, 0.49))
    trim_bottom = float(np.clip(FIGURE_CROP_TRIM_BOTTOM, 0.0, 0.49))
    trim_top = float(np.clip(FIGURE_CROP_TRIM_TOP, 0.0, 0.49))

    left = bbox.x0 + width * trim_left
    right = bbox.x1 - width * trim_right
    bottom = bbox.y0 + height * trim_bottom
    top = bbox.y1 - height * trim_top

    if right <= left:
        right = left + 1e-3
    if top <= bottom:
        top = bottom + 1e-3

    target_bbox = Bbox.from_extents(left, bottom, right, top)

    out_path = FIGURE_CROP_SAVE_PATH.expanduser()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format="pdf", bbox_inches=target_bbox)
    plt.close(fig)
    print(f"Cropped PDF saved to {out_path}")


def _collect_datasets(ds) -> List[Tuple[str, pv.DataSet]]:
    out: List[Tuple[str, pv.DataSet]] = []

    def rec(obj, prefix: str):
        if obj is None:
            return
        if isinstance(obj, pv.MultiBlock):
            keys = obj.keys()
            for i, child in enumerate(obj):
                name = keys[i] or f"Block_{i}"
                rec(child, f"{prefix}/{name}" if prefix else name)
            return
        out.append((prefix or obj.__class__.__name__, obj))

    rec(ds, "")
    return out


def _prepare_surface(mesh: pv.DataSet) -> pv.PolyData:
    if isinstance(mesh, pv.PolyData):
        surf = mesh
    else:
        try:
            surf = mesh.extract_surface()
        except Exception:
            surf = mesh.extract_geometry()
    return surf.triangulate()


def _summarize_array(arr: np.ndarray) -> str:
    a = np.asarray(arr)
    comps = 1 if a.ndim == 1 else a.shape[1]
    return f"dtype={a.dtype}, comps={comps}"


def build_plottables(
    entries: List[Tuple[str, pv.DataSet]], prefer: str = "cell"
) -> List[Tuple[str, str]]:
    """
    Return a stable, indexed list of plottables: [(location, name), ...]
      location is "cell" or "point"
    """
    point_names = OrderedDict()
    cell_names = OrderedDict()

    for _, block in entries:
        if block is None or getattr(block, "n_points", 0) == 0:
            continue
        surf = _prepare_surface(block)
        if surf.n_points == 0 or surf.n_cells == 0:
            continue

        for name in surf.point_data.keys():
            if name not in point_names:
                try:
                    point_names[name] = _summarize_array(surf.point_data[name])
                except Exception:
                    point_names[name] = "unreadable"

        for name in surf.cell_data.keys():
            if name not in cell_names:
                try:
                    cell_names[name] = _summarize_array(surf.cell_data[name])
                except Exception:
                    cell_names[name] = "unreadable"

    plottables: List[Tuple[str, str]] = []
    if prefer == "cell":
        for n in cell_names.keys():
            plottables.append(("cell", n))
        for n in point_names.keys():
            plottables.append(("point", n))
    else:
        for n in point_names.keys():
            plottables.append(("point", n))
        for n in cell_names.keys():
            plottables.append(("cell", n))

    return plottables


def print_plottables(
    plottables: List[Tuple[str, str]], entries: List[Tuple[str, pv.DataSet]]
) -> None:
    summaries: Dict[Tuple[str, str], str] = {(loc, name): "" for loc, name in plottables}

    for _, block in entries:
        if block is None or getattr(block, "n_points", 0) == 0:
            continue
        surf = _prepare_surface(block)
        if surf.n_points == 0 or surf.n_cells == 0:
            continue

        for (loc, name) in list(summaries.keys()):
            if summaries[(loc, name)]:
                continue
            data = surf.point_data if loc == "point" else surf.cell_data
            if name in data:
                try:
                    summaries[(loc, name)] = _summarize_array(data[name])
                except Exception:
                    summaries[(loc, name)] = "unreadable"

    print("\nPlottable quantities (select with --var INDEX):")
    for i, (loc, name) in enumerate(plottables):
        meta = summaries[(loc, name)]
        meta_str = f" [{meta}]" if meta else ""
        print(f"  {i:>3d}: ({loc}) {name}{meta_str}")
    print("")


def merge_surfaces(entries: List[Tuple[str, pv.DataSet]]) -> pv.PolyData:
    """Merge all drawable surfaces into one PolyData."""
    merged: pv.PolyData | None = None

    for _, block in entries:
        if block is None or getattr(block, "n_points", 0) == 0:
            continue
        surf = _prepare_surface(block)
        if surf.n_points == 0 or surf.n_cells == 0:
            continue
        merged = surf if merged is None else merged.merge(surf, merge_points=False)

    if merged is None:
        raise RuntimeError("No drawable surfaces found in CGNS.")
    return merged


def apply_axes_limits(ax, pts: np.ndarray) -> None:
    """
    Apply axis limits based on AX_LIMITS_MODE and header settings.
    """
    mins = pts.min(axis=0)
    maxs = pts.max(axis=0)

    mode = AX_LIMITS_MODE.lower().strip()

    if mode == "manual":
        ax.set_xlim(*X_LIM)
        ax.set_ylim(*Y_LIM)
        ax.set_zlim(*Z_LIM)
    
        # Use proportional aspect from the limits (keeps geometry "true")
        xr = X_LIM[1] - X_LIM[0]
        yr = Y_LIM[1] - Y_LIM[0]
        zr = Z_LIM[1] - Z_LIM[0]
    
        # Avoid zero ranges
        xr = xr if xr != 0 else 1.0
        yr = yr if yr != 0 else 1.0
        zr = zr if zr != 0 else 1.0
    
        ax.set_box_aspect((xr, yr, zr))
        return


    if mode == "auto_tight":
        ax.set_xlim(mins[0], maxs[0])
        ax.set_ylim(mins[1], maxs[1])
        ax.set_zlim(mins[2], maxs[2])
        ranges = (maxs - mins)
        # Use proportional box aspect for data extents
        ax.set_box_aspect((ranges[0], ranges[1], ranges[2]))
        return

    # default: auto_equal
    ranges = maxs - mins
    cx, cy, cz = (mins + maxs) / 2.0
    R = ranges.max() / 2.0
    ax.set_xlim(cx - R, cx + R)
    ax.set_ylim(cy - R, cy + R)
    ax.set_zlim(cz - R, cz + R)
    ax.set_box_aspect((1, 1, 1))


def plot_with_matplotlib(
    surface: pv.PolyData,
    loc: str,
    name: str,
    elev: float,
    azim: float,
    roll: float,
    cbar_min: float | None,
    cbar_max: float | None,
    cbar_label: str | None,
    cbar_labelpad: float,
    cbar_label_rotation: float,
    cbar_shrink: float,
    cbar_aspect: float,
    cbar_orientation: str,
    cbar_anchor: tuple[float, float],
    cbar_pad: float,
    cbar_label_fontsize: float,
    cbar_num_ticks: int | None,
    cbar_tick_format: str,
) -> None:
    """
    Matplotlib 3D rendering with scalar coloring.
    For cell-data, we convert to point-data first.
    """
    apply_global_plot_style()

    surf = surface

    if loc == "cell":
        # Convert cell -> point so each vertex has a scalar (better for mpl trisurf)
        surf = surf.cell_data_to_point_data(pass_cell_data=True)
        loc = "point"

    if name not in surf.point_data:
        raise KeyError(
            f"Selected array '{name}' not found in point_data after conversion."
        )

    tri = surf.triangulate()
    tri = tri.subdivide(1, subfilter="linear")
    faces = tri.faces.reshape(-1, 4)[:, 1:]  # (ntri, 3)
    pts = tri.points
    scal = np.asarray(tri.point_data[name]).ravel()

    fig = plt.figure(figsize=(11, 7), dpi=400)
    ax = fig.add_subplot(111, projection="3d")
    ax.tick_params(labelsize=TICK_FONT_SIZE)
    
    

    # Apply axis limits / aspect first (so geometry doesn't look skewed)
    apply_axes_limits(ax, pts)

    # Geometry
    tpc = ax.plot_trisurf(
        pts[:, 0],
        pts[:, 1],
        pts[:, 2],
        triangles=faces,
        cmap="inferno",
        linewidth=0,
        antialiased=True,
        shade=True,
        #edgecolor="k",
    )

    # Color by scalar (Matplotlib trisurf uses per-face colors)
    face_vals = scal[faces].mean(axis=1)
    tpc.set_array(face_vals)

    auto_vmin = float(face_vals.min())
    auto_vmax = float(face_vals.max())
    vmin = auto_vmin if cbar_min is None else float(cbar_min)
    vmax = auto_vmax if cbar_max is None else float(cbar_max)
    if vmin == vmax:
        vmax = vmin + 1e-12
    tpc.set_clim(vmin, vmax)

    cbar = fig.colorbar(
        tpc,
        ax=ax,
        shrink=cbar_shrink,
        pad=cbar_pad,
        aspect=cbar_aspect,
        orientation=cbar_orientation,
        anchor=cbar_anchor,
    )
    label_text = name if not cbar_label else cbar_label
    if cbar_orientation == "vertical":
        cbar.ax.set_ylabel(
            label_text,
            rotation=cbar_label_rotation,
            labelpad=cbar_labelpad,
            fontsize=cbar_label_fontsize,
        )
    else:
        cbar.ax.set_xlabel(
            label_text,
            rotation=cbar_label_rotation,
            labelpad=cbar_labelpad,
            fontsize=cbar_label_fontsize,
        )
    cbar.ax.tick_params(labelsize=TICK_FONT_SIZE)
    formatter = FormatStrFormatter(cbar_tick_format)
    if cbar_orientation == "vertical":
        cbar.ax.yaxis.set_major_formatter(formatter)
    else:
        cbar.ax.xaxis.set_major_formatter(formatter)

    _ensure_colorbar_end_ticks(cbar, vmin, vmax, cbar_num_ticks)

    # Labels
    #ax.set_title(f"{name} ({'cell' if loc=='point' else loc} data shown as point-colored)")
    if MPL_SHOW_AXES:
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.grid(MPL_SHOW_GRID)
    else:
        # Remove pane shading/axes while keeping the title/colorbar
        ax.set_axis_off()

    # Camera (native roll supported in your setup)
    ax.view_init(elev=elev, azim=azim, roll=roll)

    _finalize_figure(fig)



def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("cgns", nargs="?", default=str(CGNS_FILE), help="Path to CGNS")
    p.add_argument("--list", action="store_true", help="Print plottable array list and exit")
    p.add_argument("--var", type=int, default=VAR_INDEX, help="Index into plottable list (0..N-1)")
    p.add_argument("--prefer", choices=("cell", "point"), default=PREFER, help="Ordering preference for listing")

    # Optional camera overrides (otherwise header defaults apply)
    p.add_argument("--elev", type=float, default=MPL_ELEV_DEFAULT, help="Matplotlib elev (deg)")
    p.add_argument("--azim", type=float, default=MPL_AZIM_DEFAULT, help="Matplotlib azim (deg)")
    p.add_argument("--roll", type=float, default=MPL_ROLL_DEFAULT, help="Matplotlib roll (deg)")

    # Colorbar overrides
    p.add_argument("--cbar-min", type=float, default=COLORBAR_MIN, help="Colorbar lower bound / clim min")
    p.add_argument("--cbar-max", type=float, default=COLORBAR_MAX, help="Colorbar upper bound / clim max")
    p.add_argument("--cbar-label", default=COLORBAR_LABEL, help="Override colorbar label text")
    p.add_argument("--cbar-labelpad", type=float, default=COLORBAR_LABELPAD, help="Colorbar label padding")
    p.add_argument(
        "--cbar-label-rotation",
        type=float,
        default=COLORBAR_LABEL_ROTATION,
        help="Colorbar label rotation (deg)",
    )
    p.add_argument(
        "--cbar-label-fontsize",
        type=float,
        default=COLORBAR_LABEL_FONT_SIZE,
        help="Colorbar label font size",
    )
    p.add_argument("--cbar-shrink", type=float, default=COLORBAR_SHRINK, help="Colorbar shrink factor")
    p.add_argument("--cbar-aspect", type=float, default=COLORBAR_ASPECT, help="Colorbar aspect ratio")
    p.add_argument("--cbar-pad", type=float, default=COLORBAR_PAD, help="Padding between axes and colorbar")
    p.add_argument(
        "--cbar-orientation",
        choices=("vertical", "horizontal"),
        default=COLORBAR_ORIENTATION,
        help="Orientation of colorbar",
    )
    p.add_argument(
        "--cbar-tick-format",
        default=COLORBAR_TICK_FORMAT,
        help="Format string used for colorbar tick labels (e.g. %.3f)",
    )
    p.add_argument(
        "--cbar-num-ticks",
        type=int,
        default=COLORBAR_NUM_TICKS,
        help="Total number of ticks on the colorbar (>=2)",
    )
    p.add_argument(
        "--cbar-anchor",
        type=float,
        nargs=2,
        metavar=("ANCHOR_X", "ANCHOR_Y"),
        default=list(COLORBAR_ANCHOR),
        help="Anchor position for colorbar when shrink != 1",
    )

    return p.parse_args()


def main() -> None:
    args = parse_args()
    file_path = Path(args.cgns).expanduser().resolve()

    ds = pv.read(file_path)
    entries = _collect_datasets(ds)

    plottables = build_plottables(entries, prefer=args.prefer)
    if not plottables:
        raise RuntimeError("No plottable arrays found (point_data/cell_data empty).")

    print_plottables(plottables, entries)
    if args.list:
        return

    if args.var < 0 or args.var >= len(plottables):
        raise ValueError(f"--var must be in [0, {len(plottables)-1}]")

    loc, name = plottables[args.var]
    print(f"Selected: index={args.var} -> ({loc}) {name}")

    surface = merge_surfaces(entries)
    anchor = tuple(args.cbar_anchor) if args.cbar_anchor else COLORBAR_ANCHOR
    plot_with_matplotlib(
        surface=surface,
        loc=loc,
        name=name,
        elev=args.elev,
        azim=args.azim,
        roll=args.roll,
        cbar_min=args.cbar_min,
        cbar_max=args.cbar_max,
        cbar_label=args.cbar_label,
        cbar_labelpad=args.cbar_labelpad,
        cbar_label_rotation=args.cbar_label_rotation,
        cbar_shrink=args.cbar_shrink,
        cbar_aspect=args.cbar_aspect,
        cbar_orientation=args.cbar_orientation,
        cbar_tick_format=args.cbar_tick_format,
        cbar_anchor=anchor,
        cbar_pad=args.cbar_pad,
        cbar_label_fontsize=args.cbar_label_fontsize,
        cbar_num_ticks=args.cbar_num_ticks,
    )


if __name__ == "__main__":
    main()
