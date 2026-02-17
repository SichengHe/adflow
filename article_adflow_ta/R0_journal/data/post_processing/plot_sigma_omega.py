#!/usr/bin/env python
"""Plot σ₁(ω) from a resolvent sweep .npz file using consistent styling."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

try:
    import niceplots  # type: ignore
except ImportError:  # pragma: no cover
    niceplots = None


# ---------------------------------------------------------------------------
# User controls
# ---------------------------------------------------------------------------
NPZ_PATH = Path(
    r"/home/rohit/Desktop/resolvent_adflow/article_resolvent_wing/data/resolvent_analysis_airfoil_naca64a010/sweep_omega_0.005_to_0.06.npz"
)
OUTPUT_PATH = Path("figures") / "sigma_vs_omega.pdf"
USE_TEX = True
USE_NICEPLOTS_STYLE = True
FONT_FAMILY = "CMU Serif"
LABEL_FONT_SIZE = 20
TICK_FONT_SIZE = 18
LINE_WIDTH = 3.0
MARKER_STYLE = None
MARKER_SIZE = 5
X_LABEL = r"$\omega^*$"
Y_LABEL = r"$\sigma_1$"
X_LABEL_PAD = 10
Y_LABEL_PAD = 20
X_LIMITS = None  # e.g., (0.0, 0.06)
Y_LIMITS = None
SHOW_GRID = False
GRID_ALPHA = 0.3
FIGURE_SIZE = (7, 4.5)
FIGURE_DPI = 300
SHOW_FIGURE = True
# ---------------------------------------------------------------------------


def _apply_style() -> None:
    if USE_NICEPLOTS_STYLE and niceplots is not None:
        plt.style.use(niceplots.get_style())
    plt.rcParams.update(
        {
            "axes.unicode_minus": False,
            "font.family": FONT_FAMILY,
            "text.usetex": USE_TEX,
            "axes.labelsize": LABEL_FONT_SIZE,
            "xtick.labelsize": TICK_FONT_SIZE,
            "ytick.labelsize": TICK_FONT_SIZE,
        }
    )
    if USE_TEX:
        plt.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"


def _load_sigma_data(path: Path) -> tuple[np.ndarray, np.ndarray]:
    if not path.exists():
        raise FileNotFoundError(path)
    data = np.load(path)
    omega = np.asarray(data["omega"])
    sigma = np.asarray(data["sigma"])
    if sigma.ndim == 1:
        sigma_mode1 = sigma
    else:
        sigma_mode1 = sigma[:, 0]
    return omega, sigma_mode1


def main() -> None:
    omega, sigma1 = _load_sigma_data(NPZ_PATH)
    _apply_style()

    fig, ax = plt.subplots(figsize=FIGURE_SIZE, dpi=FIGURE_DPI)
    ax.plot(
        omega,
        sigma1,
        linewidth=LINE_WIDTH,
        marker=MARKER_STYLE,
        markersize=MARKER_SIZE,
    )
    ax.set_xlabel(X_LABEL, labelpad=X_LABEL_PAD)
    ax.set_ylabel(Y_LABEL, rotation=0, labelpad=Y_LABEL_PAD)
    if X_LIMITS is not None:
        ax.set_xlim(*X_LIMITS)
    if Y_LIMITS is not None:
        ax.set_ylim(*Y_LIMITS)
    if SHOW_GRID:
        ax.grid(True, alpha=GRID_ALPHA)

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_PATH, bbox_inches="tight")
    if SHOW_FIGURE:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()
