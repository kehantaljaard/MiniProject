#!/usr/bin/env python3
"""Animate miniproject marginals as a GIF.

Top panel (optional): binary detections from the dataset.
Bottom panel: marginal probabilities, with green `+` at the MAP estimate and
blue `x` at the ground-truth wumpus position (if supplied).
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter


def _numbered_files(directory: Path, pattern: re.Pattern):
    found = []
    for p in directory.iterdir():
        m = pattern.match(p.name)
        if m:
            found.append((int(m.group(1)), p))
    found.sort()
    return [p for _, p in found]


def load_marginals(out_dir: Path):
    files = _numbered_files(out_dir, re.compile(r"marginal_(\d+)\.txt$"))
    if not files:
        sys.exit(f"No marginal_*.txt files found in {out_dir}")
    return [np.loadtxt(p) for p in files]


def load_observations(dataset_dir: Path):
    files = _numbered_files(dataset_dir, re.compile(r"data_file(\d+)\.txt$"))
    if not files:
        sys.exit(f"No data_file*.txt files found in {dataset_dir}")
    return [np.loadtxt(p) for p in files]


def load_trajectory(path: Path):
    if not path.exists():
        return None
    return np.loadtxt(path, dtype=int).reshape(-1, 2)  # (x, y) per row


def setup_grid_axis(ax, W, H, title_text):
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks(np.arange(W) + 0.5)
    ax.set_yticks(np.arange(H) + 0.5)
    ax.set_xticklabels(np.arange(W))
    ax.set_yticklabels(np.arange(H))
    ax.tick_params(length=0)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    return ax.set_title(title_text)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "dir",
        help="Directory containing marginal_*.txt and trajectory.txt; visualisation.gif is written here",
    )
    ap.add_argument(
        "--observations",
        help="Dataset directory containing data_file*.txt to render on the top panel",
    )
    ap.add_argument(
        "--ground-truth",
        help="Path to wumpus_trajectory.txt with the true (x, y) per time step",
    )
    ap.add_argument("--fps", type=float, default=2.0)
    ap.add_argument("--cmap", default="Reds", help="Marginal panel colormap")
    args = ap.parse_args()

    out_dir = Path(args.dir).expanduser().resolve()
    frames = load_marginals(out_dir)
    traj = load_trajectory(out_dir / "trajectory.txt")
    obs = load_observations(Path(args.observations).expanduser()) if args.observations else None
    truth = load_trajectory(Path(args.ground_truth).expanduser()) if args.ground_truth else None

    H, W = frames[0].shape
    vmax = max(f.max() for f in frames)
    edges_x = np.arange(W + 1)
    edges_y = np.arange(H + 1)

    if obs is not None:
        fig, (ax_obs, ax_marg) = plt.subplots(2, 1, figsize=(6, 10))
    else:
        fig, ax_marg = plt.subplots()
        ax_obs = None

    obs_mesh = None
    obs_title = None
    if ax_obs is not None:
        obs_mesh = ax_obs.pcolormesh(
            edges_x, edges_y, obs[0],
            cmap="Greys", vmin=0.0, vmax=1.0,
            edgecolors="0.6", linewidth=0.5, shading="flat",
        )
        obs_title = setup_grid_axis(ax_obs, W, H, "observations t=0")

    marg_mesh = ax_marg.pcolormesh(
        edges_x, edges_y, frames[0],
        cmap=args.cmap, vmin=0.0, vmax=vmax,
        edgecolors="0.6", linewidth=0.5, shading="flat",
    )
    marg_title = setup_grid_axis(ax_marg, W, H, "marginals t=0")

    (map_marker,) = ax_marg.plot(
        [], [], marker="+", color="green", markersize=12, markeredgewidth=2,
        linestyle="None", label="MAP",
    )
    (truth_marker,) = ax_marg.plot(
        [], [], marker="x", color="blue", markersize=10, markeredgewidth=2,
        linestyle="None", label="truth",
    )
    if truth is not None:
        ax_marg.legend(loc="upper right", fontsize=8)

    def argmax_xy(i):
        if traj is not None and i < len(traj):
            return int(traj[i, 0]), int(traj[i, 1])
        r, c = np.unravel_index(np.argmax(frames[i]), frames[i].shape)
        return int(c), int(r)

    def update(i):
        marg_mesh.set_array(frames[i].ravel())
        marg_title.set_text(f"marginals t={i}")
        x, y = argmax_xy(i)
        map_marker.set_data([x + 0.5], [y + 0.5])

        if truth is not None and i < len(truth):
            tx, ty = int(truth[i, 0]), int(truth[i, 1])
            truth_marker.set_data([tx + 0.5], [ty + 0.5])

        if obs_mesh is not None and i < len(obs):
            obs_mesh.set_array(obs[i].ravel())
            obs_title.set_text(f"observations t={i}")

        return marg_mesh, map_marker, truth_marker

    anim = FuncAnimation(fig, update, frames=len(frames), interval=1000 / args.fps, blit=False)

    gif_path = out_dir / "visualisation.gif"
    anim.save(gif_path, writer=PillowWriter(fps=args.fps))
    print(f"Wrote {gif_path}")

    images_dir = out_dir / "images"
    images_dir.mkdir(exist_ok=True)
    pad = max(3, len(str(len(frames) - 1)))
    for i in range(len(frames)):
        update(i)
        fig.savefig(images_dir / f"frame_{i:0{pad}d}.png")
    print(f"Wrote {len(frames)} frames to {images_dir}")


if __name__ == "__main__":
    main()
