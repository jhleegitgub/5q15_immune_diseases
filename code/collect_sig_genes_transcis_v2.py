#!/usr/bin/env python3
# code/collect_sig_genes_transcis_v2.py
# - reads many PLINK *.assoc.linear (produced by --all-pheno) under --dir
# - keeps TEST=ADD
# - outputs columns like your old format:
#   pheno_gene pheno_chr SNP_CHR SNP BP A1 BETA STAT P NMISS q_bh
# - BH-FDR is computed across ALL collected rows (after optional trans-only filters)
import argparse, glob, os, re, sys
import pandas as pd

PATTERNS = [
    re.compile(r"_chr(\d+)\.([^.]+)\.assoc\.linear$"),
    re.compile(r"chr(\d+)\.([^.]+)\.assoc\.linear$"),
]

def parse_meta(path: str):
    base = os.path.basename(path)
    for pat in PATTERNS:
        m = pat.search(base)
        if m:
            return m.group(1), m.group(2)
    return None, None

def load_gene_tss_from_gtf(gtf: str):
    import gzip
    op = gzip.open if gtf.endswith(".gz") else open
    out = {}
    with op(gtf, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("#"):
                continue
            toks = line.rstrip("\n").split("\t")
            if len(toks) < 9:
                continue
            if toks[2] != "gene":
                continue
            chrom = toks[0].replace("chr", "")
            start = int(toks[3]); end = int(toks[4])
            strand = toks[6]
            attrs = toks[8]
            gname = None
            for seg in attrs.split(";"):
                seg = seg.strip()
                if seg.startswith("gene_name"):
                    q = seg.split('"')
                    if len(q) >= 2:
                        gname = q[1].upper()
                        break
            if not gname:
                continue
            tss = start if strand == "+" else end
            out[gname] = (chrom, tss)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, help="folder containing many *.assoc.linear (can be nested)")
    ap.add_argument("--p-raw", type=float, default=5e-6)
    ap.add_argument("--q-fdr", type=float, default=0.10)
    ap.add_argument("--test", default="ADD", help="which TEST to keep (default ADD)")
    ap.add_argument("--trans-only", action="store_true", help="keep only trans (pheno_chr != SNP_CHR)")
    ap.add_argument("--gtf", default=None, help="(optional) gtf for same-chr distance filter")
    ap.add_argument("--samechr-min-mb", type=float, default=None,
                    help="treat same-chr as 'trans-like' if |BP-TSS| >= this many Mb (requires --gtf)")
    ap.add_argument("--out-all", required=True)
    ap.add_argument("--out-sig", required=True)
    args = ap.parse_args()

    paths = sorted(glob.glob(os.path.join(args.dir, "**", "*.assoc.linear"), recursive=True))
    if not paths:
        print(f"[ERR] no assoc.linear found under: {args.dir}", file=sys.stderr)
        sys.exit(1)

    frames = []
    for p in paths:
        pheno_chr, pheno_gene = parse_meta(p)
        if not pheno_chr or not pheno_gene:
            # skip files that don't match expected naming
            continue
        try:
            df = pd.read_table(p, sep=r"\s+", dtype=str, engine="python")
        except Exception:
            continue
        if df.empty:
            continue
        if "TEST" in df.columns:
            df = df[df["TEST"] == args.test].copy()
        if df.empty:
            continue

        # standardize columns
        if "CHR" not in df.columns or "SNP" not in df.columns or "P" not in df.columns:
            continue
        df = df.rename(columns={"CHR": "SNP_CHR"})
        df["pheno_gene"] = pheno_gene
        df["pheno_chr"] = str(pheno_chr)

        keep_cols = ["pheno_gene", "pheno_chr", "SNP_CHR", "SNP", "BP", "A1", "BETA", "STAT", "P", "NMISS"]
        for c in keep_cols:
            if c not in df.columns:
                df[c] = pd.NA
        df = df[keep_cols].copy()
        frames.append(df)

    if not frames:
        print("[ERR] no usable ADD rows parsed (check filename patterns and PLINK outputs)", file=sys.stderr)
        sys.exit(1)

    all_df = pd.concat(frames, ignore_index=True)

    # numeric casting for filtering/FDR
    all_df["P"] = pd.to_numeric(all_df["P"], errors="coerce")
    all_df = all_df.dropna(subset=["P"]).copy()

    # trans-only filter
    if args.trans_only:
        all_df["pheno_chr"] = all_df["pheno_chr"].astype(str)
        all_df["SNP_CHR"] = all_df["SNP_CHR"].astype(str)
        trans_mask = (all_df["pheno_chr"] != all_df["SNP_CHR"])

        if args.samechr_min_mb is not None and args.gtf:
            gene_tss = load_gene_tss_from_gtf(args.gtf)
            min_bp = int(args.samechr_min_mb * 1_000_000)

            def far_samechr(row):
                g = str(row["pheno_gene"]).upper()
                if g not in gene_tss:
                    return False
                gchr, tss = gene_tss[g]
                if str(gchr) != str(row["SNP_CHR"]):
                    return False
                try:
                    bp = int(float(row["BP"]))
                except Exception:
                    return False
                return abs(bp - int(tss)) >= min_bp

            far_mask = all_df.apply(far_samechr, axis=1)
            all_df = all_df[trans_mask | far_mask].copy()
        else:
            all_df = all_df[trans_mask].copy()

    # BH-FDR
    try:
        from statsmodels.stats.multitest import multipletests
    except Exception:
        print("[ERR] statsmodels needed: pip install statsmodels", file=sys.stderr)
        sys.exit(1)

    all_df["q_bh"] = multipletests(all_df["P"].values, method="fdr_bh")[1]

    # output ALL
    os.makedirs(os.path.dirname(args.out_all) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(args.out_sig) or ".", exist_ok=True)

    out_cols = ["pheno_gene", "pheno_chr", "SNP_CHR", "SNP", "BP", "A1", "BETA", "STAT", "P", "NMISS", "q_bh"]
    # keep order + stable formatting
    all_df[out_cols].to_csv(args.out_all, sep="\t", index=False, float_format="%.12g")

    # SIG subset
    sig = all_df[(all_df["P"] <= args.p_raw) & (all_df["q_bh"] <= args.q_fdr)].copy()
    sig = sig.sort_values(["q_bh", "P", "pheno_gene", "SNP"], kind="mergesort")
    sig[out_cols].to_csv(args.out_sig, sep="\t", index=False, float_format="%.12g")

    print(f"[OK] all  -> {args.out_all} (n={len(all_df)})")
    print(f"[OK] sig  -> {args.out_sig} (n={len(sig)})")

if __name__ == "__main__":
    main()
