#!/usr/bin/env python3
import argparse
import gzip
import pandas as pd
from typing import Dict, Optional

def read_ld_r2(ld_gz: str, lead: str) -> Dict[str, float]:
    """
    Read plink --r2 gz output and return r2 to lead for each partner SNP.
    Expected columns include SNP_A SNP_B R2 (and others).
    """
    with gzip.open(ld_gz, "rt") as f:
        df = pd.read_csv(f, sep=r"\s+", engine="python")
    # normalize column names
    cols = {c.lower(): c for c in df.columns}
    snp_a = cols.get("snp_a") or "SNP_A"
    snp_b = cols.get("snp_b") or "SNP_B"
    r2c = cols.get("r2") or "R2"

    r2map: Dict[str, float] = {lead: 1.0}
    for _, r in df.iterrows():
        a = str(r[snp_a])
        b = str(r[snp_b])
        r2 = float(r[r2c])
        if a == lead and b != lead:
            r2map[b] = max(r2map.get(b, 0.0), r2)
        elif b == lead and a != lead:
            r2map[a] = max(r2map.get(a, 0.0), r2)
    return r2map

def consequence_score(cons: str, impact: str) -> float:
    """
    Simple severity scoring (논문용: '가중치 기반 우선순위' 정도로만 사용).
    """
    c = (cons or "").lower()
    i = (impact or "").upper()

    # hard hits
    if any(x in c for x in ["stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant"]):
        return 5.0
    # splice region-ish
    if "splice" in c:
        return 3.4
    # coding
    if "missense_variant" in c:
        return 1.2
    # UTR / regulatory-ish
    if "3_prime_utr_variant" in c or "5_prime_utr_variant" in c:
        return 1.2

    # impact fallback
    if i == "HIGH":
        return 5.0
    if i == "MODERATE":
        return 1.2
    if i == "LOW":
        return 0.8
    return 0.4

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ld", required=True, help="plink .ld.gz (r2 output)")
    ap.add_argument("--vep", required=True, help="VEP TSV from quick_vep_grch37_v2.py")
    ap.add_argument("--lead", required=True, help="lead SNP rsID")
    ap.add_argument("--pip", required=False, help="*_pip.tsv from finemap step (SNP,PIP columns)")
    ap.add_argument("--out", required=True, help="ranked output TSV")
    ap.add_argument("--top", type=int, default=20, help="top N to export")
    ap.add_argument("--top-out", required=False, help="top N TSV path")
    args = ap.parse_args()

    r2map = read_ld_r2(args.ld, args.lead)

    vep = pd.read_csv(args.vep, sep="\t", dtype=str).fillna("")
    if "SNP" not in vep.columns:
        raise SystemExit("[ERR] VEP file must contain column 'SNP'")

    # pip merge (optional)
    if args.pip:
        pip = pd.read_csv(args.pip, sep="\t")
        # tolerate different header cases
        if "SNP" not in pip.columns:
            # try first column
            pip.columns = ["SNP"] + list(pip.columns[1:])
        if "PIP" not in pip.columns:
            # try second col as PIP
            if len(pip.columns) >= 2:
                pip = pip.rename(columns={pip.columns[1]: "PIP"})
            else:
                pip["PIP"] = 0.0
        pip["SNP"] = pip["SNP"].astype(str)
        pip["PIP"] = pd.to_numeric(pip["PIP"], errors="coerce").fillna(0.0)
        vep = vep.merge(pip[["SNP", "PIP"]], on="SNP", how="left")
        vep["PIP"] = vep["PIP"].fillna(0.0)
    else:
        vep["PIP"] = 0.0

    # r2 to lead
    vep["R2_to_lead"] = vep["SNP"].map(lambda s: float(r2map.get(str(s), float("nan"))))
    # score
    vep["score"] = vep.apply(lambda r: consequence_score(r.get("most_severe_consequence",""), r.get("impact","")), axis=1)
    vep["is_lead"] = (vep["SNP"] == args.lead).astype(int)

    # sort: score -> PIP -> R2 -> lead first (ties)
    vep_sorted = vep.sort_values(
        by=["score", "PIP", "R2_to_lead", "is_lead"],
        ascending=[False, False, False, False],
        na_position="last"
    )

    # write full ranked
    vep_sorted.to_csv(args.out, sep="\t", index=False)

    # write top N
    if args.top_out:
        vep_sorted.head(args.top).to_csv(args.top_out, sep="\t", index=False)

if __name__ == "__main__":
    main()
