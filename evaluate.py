#!/usr/bin/env python3
"""
Evaluate inferred wumpus trajectories against ground truth.

Implements the three quantitative metrics described in the report:
  - Per-step accuracy
  - Mean l2 cell distance
  - Top-k accuracy (k = 1, 3, 5)

Outputs:
  * output/dataset_<N>/marginal_<ttt>.txt  -- H x W probability grid per step
  * Ground truths/dataset_<N>/wumpus_trajectory.txt  -- one "x y" per step
    where x = col, y = row (matching the report's coordinate convention).
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = ROOT / "output"
GT_DIR = ROOT / "Ground truths"
TOP_KS = (1, 3, 5)
MARGINAL_RE = re.compile(r"marginal_(\d+)\.txt$")


def load_marginals(dataset_dir: Path) -> np.ndarray:
    """Return array of shape (T+1, H, W) of posterior marginals."""
    files = sorted(
        (p for p in dataset_dir.iterdir() if MARGINAL_RE.search(p.name)),
        key=lambda p: int(MARGINAL_RE.search(p.name).group(1)),
    )
    if not files:
        raise FileNotFoundError(f"No marginal_*.txt files in {dataset_dir}")
    grids = [np.loadtxt(f) for f in files]
    return np.stack(grids, axis=0)


def load_ground_truth(gt_file: Path) -> np.ndarray:
    """Return array of shape (T+1, 2) with (row, col) per step.

    Each line of the ground-truth file is "x y" with x = col, y = row.
    """
    xy = np.loadtxt(gt_file, dtype=int)
    if xy.ndim == 1:
        xy = xy[None, :]
    cols, rows = xy[:, 0], xy[:, 1]
    return np.stack([rows, cols], axis=1)


def evaluate(marginals: np.ndarray, gt_rc: np.ndarray) -> dict:
    T, H, W = marginals.shape
    if gt_rc.shape[0] != T:
        raise ValueError(
            f"Mismatch: {T} marginal files vs {gt_rc.shape[0]} GT steps"
        )

    flat = marginals.reshape(T, H * W)
    gt_flat = gt_rc[:, 0] * W + gt_rc[:, 1]

    argmax_flat = flat.argmax(axis=1)
    argmax_rc = np.stack([argmax_flat // W, argmax_flat % W], axis=1)

    per_step_acc = float((argmax_flat == gt_flat).mean())
    dists = np.linalg.norm(argmax_rc - gt_rc, axis=1)
    mean_dist = float(dists.mean())

    # Top-k: rank cells by descending probability, check GT rank.
    order = np.argsort(-flat, axis=1)
    gt_rank = (order == gt_flat[:, None]).argmax(axis=1)
    topk = {k: float((gt_rank < k).mean()) for k in TOP_KS}

    return {
        "T": T,
        "H": H,
        "W": W,
        "per_step_acc": per_step_acc,
        "mean_l2_dist": mean_dist,
        "top_k": topk,
    }


def print_row(name: str, r: dict) -> None:
    top = r["top_k"]
    print(
        f"  {name:<12s} "
        f"acc={r['per_step_acc']*100:5.1f}%  "
        f"mean d={r['mean_l2_dist']:5.2f}  "
        f"top-1={top[1]*100:5.1f}%  "
        f"top-3={top[3]*100:5.1f}%  "
        f"top-5={top[5]*100:5.1f}%  "
        f"({r['H']}x{r['W']}, T+1={r['T']})"
    )


def main() -> None:
    datasets = sorted(
        d.name for d in OUTPUT_DIR.iterdir()
        if d.is_dir() and d.name.startswith("dataset")
    )
    if not datasets:
        raise SystemExit(f"No dataset_* dirs in {OUTPUT_DIR}")

    print("Wumpus tracking evaluation\n")
    results = {}
    for name in datasets:
        out_dir = OUTPUT_DIR / name
        gt_file = GT_DIR / name / "wumpus_trajectory.txt"
        if not gt_file.exists():
            print(f"  {name}: skipped (no ground truth at {gt_file})")
            continue
        marginals = load_marginals(out_dir)
        gt_rc = load_ground_truth(gt_file)
        results[name] = evaluate(marginals, gt_rc)
        print_row(name, results[name])

    print("\nLaTeX table row format (acc | mean dist | top-3 | top-5):")
    for name, r in results.items():
        top = r["top_k"]
        idx = name.replace("dataset", "")
        print(
            f"  {idx} & {r['per_step_acc']*100:.1f}\\% "
            f"& {r['mean_l2_dist']:.2f} "
            f"& {top[3]*100:.1f}\\% "
            f"& {top[5]*100:.1f}\\% \\\\"
        )


if __name__ == "__main__":
    main()
