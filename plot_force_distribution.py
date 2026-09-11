#!/usr/bin/env python3
"""Plot a spatial force distribution produced by 2Dmain.cpp.

The historical simulator writes rows of:
    row  column  mean_force  uncertainty

This helper uses only NumPy and Matplotlib and avoids the removed
matplotlib.mlab.griddata API used by the original plotting script.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_force_grid(path: Path):
    data = np.loadtxt(path, ndmin=2)
    if data.shape[1] < 3:
        raise ValueError(
            f"{path} must contain at least 3 columns: row, column, force"
        )

    row = data[:, 0].astype(int)
    col = data[:, 1].astype(int)
    force = data[:, 2]

    rows = np.unique(row)
    cols = np.unique(col)
    expected = len(rows) * len(cols)
    if data.shape[0] != expected:
        raise ValueError(
            "Input is not a complete rectangular grid: "
            f"found {data.shape[0]} rows but expected {expected}."
        )

    grid = np.full((len(rows), len(cols)), np.nan, dtype=float)
    r_index = {value: i for i, value in enumerate(rows)}
    c_index = {value: i for i, value in enumerate(cols)}

    for r, c, f in zip(row, col, force):
        grid[r_index[r], c_index[c]] = f

    if np.isnan(grid).any():
        raise ValueError("Input contains duplicate/missing grid coordinates.")

    return rows, cols, grid


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("Force_Dist.txt"),
        help="Path to Force_Dist.txt (default: %(default)s)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional image output path, e.g. force_distribution.png",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Do not open an interactive plotting window.",
    )
    args = parser.parse_args()

    rows, cols, force = load_force_grid(args.input)
    X, Y = np.meshgrid(cols, rows)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(X, Y, force, alpha=0.75)

    ax.set_xlabel("Column")
    ax.set_ylabel("Row")
    ax.set_zlabel("Mean force (pN)")
    ax.set_title("Actin filament force distribution")
    fig.tight_layout()

    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.output, dpi=200, bbox_inches="tight")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
