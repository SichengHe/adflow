#!/usr/bin/env python
"""Plot resolvent frequency sweep results from a .npz file."""
import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description="Plot resolvent sweep data")
    parser.add_argument("npz", help="Path to sweep .npz file")
    parser.add_argument("--out", default=None, help="Output image file (png)")
    args = parser.parse_args()

    npz_path = Path(args.npz)
    data = np.load(npz_path)
    omega = data["omega"]
    sigma1 = data["sigma1"]

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(omega, sigma1, marker='o', linewidth=1.5)
    ax.set_xlabel(r"$\omega$")
    ax.set_ylabel(r"$\sigma_1(\omega)$")
    ax.set_title("Resolvent Sweep")
    ax.grid(True, alpha=0.3)

    out = args.out
    if out is None:
        out = npz_path.with_suffix('.png')
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print(f"Saved plot: {out}")


if __name__ == "__main__":
    main()
