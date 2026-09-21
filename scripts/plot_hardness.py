#!/usr/bin/env python3
"""
plot_hardness.py -- render output/hardness_map.csv and
output/case_depth_profile.csv from laser_hardening_sim as PNGs.

Usage:
    python3 plot_hardness.py [output_dir]   (default: ./output)

Requires: numpy, matplotlib (not needed to build or run the C++ core).
"""
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def load_hardness_map(path):
    xs, ys, hv = [], [], []
    with open(path) as f:
        for row in csv.DictReader(f):
            xs.append(float(row["x_mm"]))
            ys.append(float(row["y_mm"]))
            hv.append(float(row["hardness_HV"]))
    xs, ys, hv = map(np.array, (xs, ys, hv))
    nx = len(np.unique(xs))
    ny = len(np.unique(ys))
    return xs.reshape(ny, nx), ys.reshape(ny, nx), hv.reshape(ny, nx)


def main():
    out_dir = Path(sys.argv[1] if len(sys.argv) > 1 else "./output")

    X, Y, HV = load_hardness_map(out_dir / "hardness_map.csv")
    fig, ax = plt.subplots(figsize=(9, 3))
    mesh = ax.pcolormesh(X, Y, HV, shading="auto", cmap="inferno")
    ax.invert_yaxis()
    ax.set_xlabel("x (scan direction) [mm]")
    ax.set_ylabel("depth [mm]")
    ax.set_title("Estimated hardness after laser scan [HV]")
    fig.colorbar(mesh, ax=ax, label="HV")
    fig.tight_layout()
    fig.savefig(out_dir / "hardness_map.png", dpi=150)
    print(f"wrote {out_dir / 'hardness_map.png'}")

    xs, depth = [], []
    with open(out_dir / "case_depth_profile.csv") as f:
        for row in csv.DictReader(f):
            xs.append(float(row["x_mm"]))
            depth.append(float(row["case_depth_mm"]))

    fig2, ax2 = plt.subplots(figsize=(9, 3))
    ax2.plot(xs, depth)
    ax2.set_xlabel("x (scan direction) [mm]")
    ax2.set_ylabel("case depth [mm]")
    ax2.set_title("Hardened case depth along the scan")
    fig2.tight_layout()
    fig2.savefig(out_dir / "case_depth_profile.png", dpi=150)
    print(f"wrote {out_dir / 'case_depth_profile.png'}")


if __name__ == "__main__":
    main()
