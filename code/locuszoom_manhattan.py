#!/usr/bin/env python3
import argparse
import math
from typing import Dict, List

import matplotlib
matplotlib.use("Agg")  # headless (no Qt/wayland)
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd


def read_assoc_linear(path: str, test: str = "ADD") -> pd.DataFrame:
    df = pd.read_table(path, sep=r"\s+", dtype=str)
    need = {"CHR", "SNP", "BP", "P"}
    miss = need - set(df.columns)
    if miss:
        raise SystemExit(f"[ERR] assoc missing columns {sorted(miss)}: {path}")

    if "TEST" in df.columns and test:
        df = df[df["TEST"] == test].copy()

    df["CHR"] = pd.to_numeric(df["CHR"], errors="coerce")
    df["BP"] = pd.to_numeric(df["BP"], errors="coerce")
    df["P"] = pd.to_numeric(df["P"], errors="coerce")
    df = df.dropna(subset=["CHR", "BP", "P", "SNP"])
    df = df[(df["CHR"] >= 1) & (df["CHR"] <= 22)]
    df = df[df["BP"] > 0]
    df.loc[df["P"] <= 0, "P"] = 1e-300
    df["mlogp"] = -np.log10(df["P"].astype(float).values)
    df["BP"] = df["BP"].astype(int)
    df["SNP"] = df["SNP"].astype(str)
    return df[["CHR", "SNP", "BP", "P", "mlogp"]]


def read_ld_r2(path: str, lead: str) -> Dict[str, float]:
    if not path:
        return {}
    df = pd.read_table(path, sep=r"\s+", dtype=str)
    if not {"SNP_A", "SNP_B", "R2"} <= set(df.columns):
        raise SystemExit(f"[ERR] LD file missing SNP_A/SNP_B/R2: {path}")

    df = df[df["SNP_A"] == lead].copy()
    if len(df) == 0:
        out = {lead: 1.0}
        return out

    df["R2"] = pd.to_numeric(df["R2"], errors="coerce").fillna(0.0).clip(0, 1)
    out = dict(zip(df["SNP_B"].astype(str), df["R2"].astype(float)))
    out[lead] = 1.0
    return out


def plot_panel(ax, assoc: pd.DataFrame, r2map: Dict[str, float], lead: str,
               title: str, gw_p: float, point_size: float = 18.0):
    d = assoc.copy()
    d["r2"] = d["SNP"].map(r2map).fillna(0.0).astype(float).clip(0, 1)

    x = d["BP"].values / 1e6
    y = d["mlogp"].values

    sc = ax.scatter(x, y, c=d["r2"].values, s=point_size, linewidths=0)
    ax.axhline(-math.log10(gw_p), linewidth=1.0)
    ax.set_title(title)
    ax.set_xlabel("Position (Mb)")
    ax.set_ylabel(r"$-\log_{10}(P)$")

    # lead: red dot only (no label)
    ld = d[d["SNP"] == lead]
    if len(ld) > 0:
        ax.scatter((ld["BP"].values / 1e6), ld["mlogp"].values,
                   s=point_size * 2.5, c="red", linewidths=0, zorder=5)
    return sc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assoc", required=True, help="assoc.linear file(s), comma-separated for multi-panel")
    ap.add_argument("--ld", required=True, help="plink .ld file(s), comma-separated (same count as --assoc)")
    ap.add_argument("--lead", required=True, help="lead SNP rsID(s), comma-separated (same count as --assoc)")
    ap.add_argument("--out-png", required=True)

    ap.add_argument("--title", default=None, help="single-panel title")
    ap.add_argument("--titles", default=None, help="comma-separated titles for multi-panel")

    ap.add_argument("--gw-threshold", type=float, default=5e-8)
    ap.add_argument("--test", default="ADD")
    ap.add_argument("--point-size", type=float, default=18.0)
    args = ap.parse_args()

    assoc_list = [x.strip() for x in args.assoc.split(",") if x.strip()]
    ld_list = [x.strip() for x in args.ld.split(",") if x.strip()]
    lead_list = [x.strip() for x in args.lead.split(",") if x.strip()]
    if not (len(assoc_list) == len(ld_list) == len(lead_list)):
        raise SystemExit("[ERR] --assoc/--ld/--lead count mismatch")

    n = len(assoc_list)

    if n == 1:
        assoc = read_assoc_linear(assoc_list[0], test=args.test)
        r2map = read_ld_r2(ld_list[0], lead_list[0])
        t = args.title if args.title else "Locus (±500 kb)"

        fig = plt.figure(figsize=(10, 7), dpi=150)
        ax = fig.add_subplot(111)
        sc = plot_panel(ax, assoc, r2map, lead_list[0], t, args.gw_threshold, args.point_size)
        cb = fig.colorbar(sc, ax=ax)
        cb.set_label(r"LD $r^2$ to lead")
        fig.tight_layout()
        fig.savefig(args.out_png, bbox_inches="tight")
        plt.close(fig)
        return

    titles = None
    if args.titles:
        titles = [x.strip() for x in args.titles.split(",")]
        if len(titles) != n:
            raise SystemExit("[ERR] --titles count mismatch with --assoc")

    fig = plt.figure(figsize=(5.8 * n, 4.2), dpi=150)
    axes = [fig.add_subplot(1, n, i + 1) for i in range(n)]

    mappable = None
    for i in range(n):
        assoc = read_assoc_linear(assoc_list[i], test=args.test)
        r2map = read_ld_r2(ld_list[i], lead_list[i])
        t = titles[i] if titles else f"Panel {i+1}"
        sc = plot_panel(axes[i], assoc, r2map, lead_list[i], t, args.gw_threshold, args.point_size)
        mappable = sc
        if i != 0:
            axes[i].set_ylabel("")

    cb = fig.colorbar(mappable, ax=axes, fraction=0.025, pad=0.02)
    cb.set_label(r"LD $r^2$ to lead")
    fig.tight_layout()
    fig.savefig(args.out_png, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
