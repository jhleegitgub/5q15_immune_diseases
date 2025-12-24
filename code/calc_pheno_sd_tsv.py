#!/usr/bin/env python3
import sys
import pandas as pd

def main():
    if len(sys.argv) < 3:
        print("usage: calc_pheno_sd_tsv.py <pheno.tsv> GENE1 [GENE2 ...]", file=sys.stderr)
        sys.exit(1)

    pheno = sys.argv[1]
    genes = sys.argv[2:]

    df = pd.read_csv(pheno, sep=r"\s+|\t", engine="python")
    # expect FID IID + gene columns
    out = []
    for g in genes:
        if g not in df.columns:
            raise SystemExit(f"[ERR] gene column not found: {g}")
        x = pd.to_numeric(df[g], errors="coerce").dropna()
        out.append((g, float(x.std(ddof=1))))
    print("gene\tsd")
    for g, sd in out:
        print(f"{g}\t{sd}")

if __name__ == "__main__":
    main()
