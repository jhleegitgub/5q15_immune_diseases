#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd

def logsumexp(a: np.ndarray) -> float:
    m = np.nanmax(a)
    if not np.isfinite(m):
        return float("nan")
    return float(m + np.log(np.nansum(np.exp(a - m))))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, help="PLINK .assoc.linear")
    ap.add_argument("--out", dest="out", required=True, help="output pip table tsv")
    ap.add_argument("--prefix", dest="prefix", required=True, help="prefix for credible set")
    ap.add_argument("--prior-sd", dest="prior_sd", type=float, required=True, help="prior SD in phenotype units")
    ap.add_argument("--credible", dest="credible", type=float, default=0.95, help="credible set threshold")
    args = ap.parse_args()

    df = pd.read_csv(args.inp, sep=r"\s+|\t", engine="python", dtype=str)
    # normalize columns
    for col in ["CHR","SNP","BP","A1","TEST","NMISS","BETA","STAT","P"]:
        if col not in df.columns:
            raise SystemExit(f"[ERR] missing column {col} in {args.inp}")

    # keep ADD only if present
    if "TEST" in df.columns:
        df = df[df["TEST"] == "ADD"].copy()

    # numeric
    for c in ["CHR","BP","NMISS","BETA","STAT","P"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # drop rows where plink prints NA for conditioned SNPs etc.
    df = df.dropna(subset=["SNP","CHR","BP","A1","BETA","STAT","P"])
    df = df[df["STAT"] != 0].copy()

    # SE: PLINK linear has no SE; use SE = |BETA/STAT|
    df["SE"] = (df["BETA"].abs() / df["STAT"].abs()).replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=["SE"])
    df = df[df["SE"] > 0].copy()

    # Z
    df["Z"] = df["BETA"] / df["SE"]

    W = float(args.prior_sd) ** 2
    V = df["SE"].values ** 2
    r = W / V

    z = df["Z"].values
    # Wakefield log(ABF): -0.5*log(1+r) + (z^2*r)/(2*(1+r))
    logabf = -0.5 * np.log1p(r) + (z*z*r) / (2.0*(1.0 + r))
    df["logABF"] = logabf

    lse = logsumexp(logabf)
    pip = np.exp(logabf - lse)
    df["PIP"] = pip

    # sort by PIP desc
    df = df.sort_values("PIP", ascending=False).reset_index(drop=True)
    df["CUM_PIP"] = df["PIP"].cumsum()

    # ABF (optional) - clip to avoid inf printing mess, but PIP uses logABF anyway
    df["ABF"] = np.exp(np.clip(df["logABF"].values, -700, 700))

    out_cols = ["SNP","CHR","BP","A1","BETA","SE","STAT","P","logABF","ABF","PIP","CUM_PIP"]
    df[out_cols].to_csv(args.out, sep="\t", index=False)

    cred = df[df["CUM_PIP"] <= args.credible].copy()
    # ensure the first row that crosses threshold is included
    if len(df) > 0 and (len(cred) == 0 or cred["CUM_PIP"].max() < args.credible):
        k = int(np.searchsorted(df["CUM_PIP"].values, args.credible, side="left"))
        cred = df.iloc[:k+1].copy()

    cred_path = f"{args.prefix}_credible95.tsv"
    cred[out_cols].to_csv(cred_path, sep="\t", index=False)

    print(f"[OK] wrote {args.out}")
    print(f"[OK] wrote {cred_path} (n={cred.shape[0]})")

if __name__ == "__main__":
    main()
