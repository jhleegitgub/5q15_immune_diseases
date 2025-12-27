#!/usr/bin/env python3
import argparse
import os
import re
import glob
import pandas as pd

def bh_fdr(pvals):
    """Benjamini-Hochberg FDR; returns q-values aligned to input order."""
    p = pd.Series(pvals, dtype="float64")
    m = p.notna().sum()
    q = pd.Series([float("nan")] * len(p), index=p.index, dtype="float64")
    if m == 0:
        return q

    p_non = p[p.notna()]
    order = p_non.sort_values().index
    ranks = pd.Series(range(1, len(order) + 1), index=order, dtype="float64")
    q_raw = (p_non.loc[order] * m / ranks)

    # enforce monotonicity
    q_mon = q_raw.copy()
    for i in range(len(q_mon) - 2, -1, -1):
        q_mon.iloc[i] = min(q_mon.iloc[i], q_mon.iloc[i + 1])

    q.loc[order] = q_mon.clip(upper=1.0).values
    return q

def parse_chr_from_path(path: str):
    m = re.search(r"chr(\d+)", os.path.basename(path))
    return int(m.group(1)) if m else None

def read_assoc_linear(path: str):
    if os.path.getsize(path) == 0:
        return None
    try:
        df = pd.read_csv(path, sep=r"\s+", dtype=str, engine="python")
    except Exception:
        return None
    if df is None or df.empty:
        return None
    return df

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, help="directory with *.assoc.linear")
    ap.add_argument("--pattern", default="*.assoc.linear", help="glob pattern")
    ap.add_argument("--p-raw", default="5e-6")
    ap.add_argument("--q-fdr", default="0.1")
    ap.add_argument("--out-all", required=True)
    ap.add_argument("--out-sig", required=True)
    ap.add_argument("--out-top", required=True)
    args = ap.parse_args()

    p_raw = float(args.p_raw)
    q_fdr = float(args.q_fdr)

    files = sorted(glob.glob(os.path.join(args.dir, args.pattern)))
    if not files:
        raise SystemExit(f"[ERR] no assoc files in: {args.dir}")

    rows = []
    for f in files:
        df = read_assoc_linear(f)
        if df is None:
            continue

        # required columns
        # PLINK --linear: CHR SNP BP A1 TEST NMISS BETA STAT P (SE may not exist)
        for c in ["CHR", "SNP", "BP", "A1", "TEST", "NMISS", "BETA", "STAT", "P"]:
            if c not in df.columns:
                # tolerate missing STAT/BETA in some cases? but usually exist
                pass

        # PHENO column appears when --all-pheno used
        if "PHENO" not in df.columns:
            df["PHENO"] = "PHENO"

        # filter ADD
        if "TEST" in df.columns:
            df = df[df["TEST"].astype(str) == "ADD"]
        if df.empty:
            continue

        # numeric P
        df["P"] = pd.to_numeric(df["P"], errors="coerce")
        df = df[df["P"].notna()]
        if df.empty:
            continue

        pheno_chr = parse_chr_from_path(f)
        df["pheno_chr"] = pheno_chr
        df["pheno_gene"] = df["PHENO"].astype(str)

        # normalize SNP_CHR
        if "CHR" in df.columns:
            df["SNP_CHR"] = pd.to_numeric(df["CHR"], errors="coerce").astype("Int64")
        else:
            df["SNP_CHR"] = pd.NA

        # cis/trans label (pheno_chr 없으면 NA)
        def cis_trans(row):
            if pd.isna(row.get("pheno_chr")) or pd.isna(row.get("SNP_CHR")):
                return "NA"
            return "cis" if int(row["pheno_chr"]) == int(row["SNP_CHR"]) else "trans"

        df["cis_trans"] = df.apply(cis_trans, axis=1)
        df["source_file"] = os.path.basename(f)

        keep_cols = ["pheno_gene","pheno_chr","cis_trans","SNP_CHR","SNP","BP","A1","BETA","STAT","P","NMISS","source_file"]
        for c in keep_cols:
            if c not in df.columns:
                df[c] = pd.NA
        rows.append(df[keep_cols])

    if not rows:
        raise SystemExit("[ERR] no valid rows parsed from assoc files")

    out = pd.concat(rows, ignore_index=True)

    out["q_bh"] = bh_fdr(out["P"].values)

    # write all
    out.to_csv(args.out_all, sep="\t", index=False)

    # significant
    sig = out[(out["P"] <= p_raw) & (out["q_bh"] <= q_fdr)].copy()
    sig.to_csv(args.out_sig, sep="\t", index=False)

    # top1 per phenotype among significant
    if not sig.empty:
        sig_sorted = sig.sort_values(["pheno_gene","P"], ascending=[True, True])
        top1 = sig_sorted.groupby("pheno_gene", as_index=False).head(1)
    else:
        top1 = sig.head(0)

    top1.to_csv(args.out_top, sep="\t", index=False)
    print(f"[OK] wrote {args.out_all}")
    print(f"[OK] wrote {args.out_sig}")
    print(f"[OK] wrote {args.out_top}")

if __name__ == "__main__":
    main()
