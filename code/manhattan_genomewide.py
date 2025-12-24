#!/usr/bin/env python3
import argparse
import math
import random
from typing import List, Optional, Tuple, Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# GRCh37/hg19 chr1-22 lengths (bp)
CHR_LEN_HG19 = {
    1: 249250621, 2: 243199373, 3: 198022430, 4: 191154276, 5: 180915260,
    6: 171115067, 7: 159138663, 8: 146364022, 9: 141213431, 10: 135534747,
    11: 135006516, 12: 133851895, 13: 115169878, 14: 107349540, 15: 102531392,
    16:  90354753, 17:  81195210, 18:  78077248, 19:  59128983, 20:  63025520,
    21:  48129895, 22:  51304566
}

def _chr_offsets() -> Tuple[List[int], Dict[int, int], List[int]]:
    chrs = list(range(1, 23))
    offsets = {}
    cum = 0
    mids = []
    for c in chrs:
        offsets[c] = cum
        mids.append(cum + CHR_LEN_HG19[c] // 2)
        cum += CHR_LEN_HG19[c]
    return chrs, offsets, mids


def _read_plink_assoc_linear_downsample(
    assoc_path: str,
    test: str = "ADD",
    max_points: int = 1_200_000,
    keep_p: float = 1e-5,
    seed: int = 1,
    highlight_snps: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Reads PLINK *.assoc.linear in chunks and downsamples for plotting.
    Keeps all P <= keep_p and all highlight_snps; reservoir-samples the rest up to max_points.
    """
    rng = random.Random(seed)
    highlight = set([s.strip() for s in (highlight_snps or []) if s.strip()])

    keep_rows = []
    reservoir = []
    n_other_seen = 0

    it = pd.read_table(assoc_path, delim_whitespace=True, dtype=str, chunksize=800_000)
    for chunk in it:
        for c in ["CHR", "SNP", "BP", "P"]:
            if c not in chunk.columns:
                raise SystemExit(f"[ERR] missing column {c} in {assoc_path}")

        if "TEST" in chunk.columns and test is not None:
            chunk = chunk[chunk["TEST"] == test]

        chunk["CHR"] = pd.to_numeric(chunk["CHR"], errors="coerce")
        chunk["BP"]  = pd.to_numeric(chunk["BP"], errors="coerce")
        chunk["P"]   = pd.to_numeric(chunk["P"], errors="coerce")
        chunk = chunk.dropna(subset=["CHR", "BP", "P", "SNP"])
        chunk = chunk[(chunk["CHR"] >= 1) & (chunk["CHR"] <= 22)]
        chunk = chunk[chunk["BP"] > 0]
        chunk.loc[chunk["P"] <= 0, "P"] = 1e-300

        is_keep = (chunk["P"] <= keep_p) | (chunk["SNP"].isin(highlight))
        keep_part = chunk.loc[is_keep, ["CHR", "SNP", "BP", "P"]]
        other_part = chunk.loc[~is_keep, ["CHR", "SNP", "BP", "P"]]

        if len(keep_part) > 0:
            keep_rows.append(keep_part)

        for row in other_part.itertuples(index=False):
            n_other_seen += 1
            if len(reservoir) < max_points:
                reservoir.append((int(row.CHR), str(row.SNP), int(row.BP), float(row.P)))
            else:
                j = rng.randrange(n_other_seen)
                if j < max_points:
                    reservoir[j] = (int(row.CHR), str(row.SNP), int(row.BP), float(row.P))

    df_keep = pd.concat(keep_rows, ignore_index=True) if keep_rows else pd.DataFrame(columns=["CHR","SNP","BP","P"])
    df_res  = pd.DataFrame(reservoir, columns=["CHR","SNP","BP","P"]) if reservoir else pd.DataFrame(columns=["CHR","SNP","BP","P"])
    df = pd.concat([df_keep, df_res], ignore_index=True)

    df["CHR"] = df["CHR"].astype(int)
    df["BP"] = df["BP"].astype(int)
    df["P"] = df["P"].astype(float)
    df["mlogp"] = -np.log10(df["P"].values)
    return df


def _plot_one(ax, df: pd.DataFrame, title: str, gw_p: float, point_size: float,
              highlight_snps: Optional[List[str]] = None):
    chrs, offsets, mids = _chr_offsets()
    df = df.copy()
    df["x"] = df.apply(lambda r: offsets[int(r["CHR"])] + int(r["BP"]), axis=1)

    # Plot chromosome-by-chromosome to force fixed grayscale (no colormap)
    for c in chrs:
        d = df[df["CHR"] == c]
        if len(d) == 0:
            continue
        col = "0.25" if (c % 2 == 1) else "0.65"   # dark gray vs light gray
        ax.scatter(d["x"].values, d["mlogp"].values, s=point_size, c=col, alpha=0.7, linewidths=0)

    # optional highlight (OFF by default)
    hl = set([s.strip() for s in (highlight_snps or []) if s.strip()])
    if hl:
        d2 = df[df["SNP"].isin(hl)]
        if len(d2) > 0:
            ax.scatter(d2["x"].values, d2["mlogp"].values, s=max(point_size * 8, 18), c="red", alpha=0.9, linewidths=0)

    ax.axhline(-math.log10(gw_p), linewidth=1.0)
    ax.set_title(title)
    ax.set_ylabel(r"$-\log_{10}(P)$")
    ax.set_xticks(mids)
    ax.set_xticklabels([str(c) for c in chrs])
    ax.set_xlabel("Chromosome")
    ax.margins(x=0.01)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assoc", required=True,
                    help="One file or comma-separated list of PLINK .assoc.linear files (multi-panel if >1).")
    ap.add_argument("--out-png", required=True)
    ap.add_argument("--title", default=None, help="Single panel title.")
    ap.add_argument("--titles", default=None, help="Comma-separated panel titles (multi-panel).")
    ap.add_argument("--gw-threshold", type=float, default=5e-8)
    ap.add_argument("--max-points", type=int, default=1_200_000)
    ap.add_argument("--keep-p", type=float, default=1e-5)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--highlight-snps", default="", help="Comma-separated rsIDs to highlight (optional). Default none.")
    ap.add_argument("--point-size", type=float, default=1.0)
    args = ap.parse_args()

    assoc_list = [x.strip() for x in args.assoc.split(",") if x.strip()]
    hl = [x.strip() for x in args.highlight_snps.split(",") if x.strip()]

    if len(assoc_list) == 1:
        df = _read_plink_assoc_linear_downsample(
            assoc_list[0],
            max_points=args.max_points,
            keep_p=args.keep_p,
            seed=args.seed,
            highlight_snps=hl
        )
        fig = plt.figure(figsize=(16, 5), dpi=150)
        ax = fig.add_subplot(111)
        _plot_one(ax, df, args.title or "Genome-wide eQTL", args.gw_threshold, args.point_size, hl)
        fig.tight_layout()
        fig.savefig(args.out_png, bbox_inches="tight")
        plt.close(fig)
        return

    titles = None
    if args.titles:
        titles = [x.strip() for x in args.titles.split(",")]
        if len(titles) != len(assoc_list):
            raise SystemExit("[ERR] --titles count must match --assoc count")

    n = len(assoc_list)
    fig = plt.figure(figsize=(18, 5), dpi=150)
    axes = [fig.add_subplot(1, n, i+1) for i in range(n)]

    for i, (ax, path) in enumerate(zip(axes, assoc_list)):
        df = _read_plink_assoc_linear_downsample(
            path,
            max_points=args.max_points,
            keep_p=args.keep_p,
            seed=args.seed,
            highlight_snps=hl
        )
        t = titles[i] if titles else f"Panel {i+1}"
        _plot_one(ax, df, t, args.gw_threshold, args.point_size, hl)
        if i != 0:
            ax.set_ylabel("")

    fig.tight_layout()
    fig.savefig(args.out_png, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
