# code/manhattan_genomewide.py  (기존 파일에 아래처럼 --suptitle 옵션만 추가된 전체본으로 덮어쓰기)
#!/usr/bin/env python3
import argparse
import math
from typing import List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd


def read_assoc_linear(path: str, test: str = "ADD") -> pd.DataFrame:
    df = pd.read_table(path, sep=r"\s+", dtype=str)
    need = {"CHR", "BP", "P"}
    if not need.issubset(df.columns):
        raise SystemExit(f"[ERR] missing columns {sorted(need - set(df.columns))}: {path}")

    if "TEST" in df.columns and test:
        df = df[df["TEST"] == test].copy()

    df["CHR"] = pd.to_numeric(df["CHR"], errors="coerce")
    df["BP"] = pd.to_numeric(df["BP"], errors="coerce")
    df["P"] = pd.to_numeric(df["P"], errors="coerce")
    df = df.dropna(subset=["CHR", "BP", "P"])
    df = df[(df["CHR"] >= 1) & (df["CHR"] <= 22)]
    df = df[df["BP"] > 0]
    df.loc[df["P"] <= 0, "P"] = 1e-300
    df["mlogp"] = -np.log10(df["P"].astype(float).values)
    df["CHR"] = df["CHR"].astype(int)
    df["BP"] = df["BP"].astype(int)
    return df[["CHR", "BP", "P", "mlogp"]]


def thin_points(df: pd.DataFrame, max_points: int, seed: int) -> pd.DataFrame:
    if max_points <= 0 or len(df) <= max_points:
        return df
    rng = np.random.default_rng(seed)
    idx = rng.choice(df.index.values, size=max_points, replace=False)
    return df.loc[idx].copy()


def manhattan_coords(df: pd.DataFrame) -> Tuple[pd.DataFrame, np.ndarray, List[int]]:
    chr_max = df.groupby("CHR")["BP"].max().reindex(range(1, 23)).fillna(0).astype(int).values
    chr_offsets = np.concatenate([[0], np.cumsum(chr_max[:-1])])
    d = df.copy()
    d["x"] = d["BP"].values + chr_offsets[d["CHR"].values - 1]
    centers = chr_offsets + (chr_max / 2.0)
    chrs = list(range(1, 23))
    return d, centers, chrs


def plot_one(ax, df: pd.DataFrame, title: str, gw_p: float,
             max_points: int, seed: int, xtick_step: int = 1):
    df = thin_points(df, max_points=max_points, seed=seed)
    d, centers, chrs = manhattan_coords(df)

    for chr_i in range(1, 23):
        sub = d[d["CHR"] == chr_i]
        if len(sub) == 0:
            continue
        color = "0.2" if (chr_i % 2 == 1) else "0.7"
        ax.scatter(sub["x"].values, sub["mlogp"].values, s=1.2, c=color, linewidths=0)

    ax.axhline(-math.log10(gw_p), linewidth=1.0)

    show = np.array(chrs)[::xtick_step]
    ax.set_xticks(centers[::xtick_step])
    ax.set_xticklabels([str(c) for c in show], fontsize=9)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel(r"$-\log_{10}(P)$")
    ax.set_title(title)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assoc", required=True, help="assoc.linear file(s), comma-separated for multi-panel")
    ap.add_argument("--out-png", required=True)
    ap.add_argument("--title", default=None, help="single title")
    ap.add_argument("--titles", default=None, help="comma-separated titles for multi-panel")
    ap.add_argument("--suptitle", default=None, help="overall title for multi-panel")
    ap.add_argument("--gw-threshold", type=float, default=5e-8)
    ap.add_argument("--test", default="ADD")
    ap.add_argument("--max-points", type=int, default=0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--xtick-step", type=int, default=1)
    args = ap.parse_args()

    assoc_list = [x.strip() for x in args.assoc.split(",") if x.strip()]
    n = len(assoc_list)

    if n == 1:
        df = read_assoc_linear(assoc_list[0], test=args.test)
        title = args.title if args.title else "Genome-wide eQTL scan"
        fig = plt.figure(figsize=(12, 4.6), dpi=150)
        ax = fig.add_subplot(111)
        plot_one(ax, df, title, args.gw_threshold, args.max_points, args.seed, xtick_step=1)
        fig.tight_layout()
        fig.savefig(args.out_png, bbox_inches="tight")
        plt.close(fig)
        return

    titles = None
    if args.titles:
        titles = [x.strip() for x in args.titles.split(",")]
        if len(titles) != n:
            raise SystemExit("[ERR] --titles count mismatch with --assoc")

    fig_w = 8.2 * n
    fig_h = 4.6
    fig = plt.figure(figsize=(fig_w, fig_h), dpi=150)
    axes = [fig.add_subplot(1, n, i + 1) for i in range(n)]

    for i in range(n):
        df = read_assoc_linear(assoc_list[i], test=args.test)
        t = titles[i] if titles else f"Panel {i+1}"
        plot_one(axes[i], df, t, args.gw_threshold, args.max_points, args.seed + i, xtick_step=args.xtick_step)
        if i != 0:
            axes[i].set_ylabel("")

    if args.suptitle:
        fig.suptitle(args.suptitle, y=0.98)
        fig.tight_layout(rect=[0, 0, 1, 0.94])
    else:
        fig.tight_layout()

    fig.savefig(args.out_png, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
