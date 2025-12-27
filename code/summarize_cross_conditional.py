#!/usr/bin/env python3
import argparse
import os
import re
import pandas as pd

def read_signals(path: str) -> pd.DataFrame:
    # tolerate tabs/whitespace + lead/snp/SNP column names
    df = pd.read_csv(path, sep=r"\s+|\t", engine="python")
    cols = {c.lower(): c for c in df.columns}
    if "lead" not in cols:
        if "snp" in cols:
            df = df.rename(columns={cols["snp"]: "lead"})
        elif "SNP" in df.columns:
            df = df.rename(columns={"SNP": "lead"})
        else:
            # fallback: 3rd column is lead
            if df.shape[1] < 3:
                raise ValueError(f"signals file needs >=3 cols: {path}")
            df = df.rename(columns={df.columns[2]: "lead"})
    if "gene" not in cols and "gene" not in df.columns:
        df = df.rename(columns={df.columns[0]: "gene"})
    if "signal_id" not in cols and "signal_id" not in df.columns:
        df = df.rename(columns={df.columns[1]: "signal_id"})
    return df[["gene", "signal_id", "lead"]].copy()

def read_runs(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    for c in ["run_id", "outcome", "assoc"]:
        if c not in df.columns:
            raise ValueError(f"runs.tsv missing column '{c}': {path}")
    return df

def read_assoc_linear(path: str) -> pd.DataFrame:
    df = pd.read_table(path, sep=r"\s+")
    if "TEST" in df.columns:
        df = df[df["TEST"] == "ADD"].copy()
    return df

def extract_beta_p(assoc_path: str, snp: str):
    if not os.path.exists(assoc_path):
        return (None, None)
    df = read_assoc_linear(assoc_path)
    if "SNP" not in df.columns:
        return (None, None)
    hit = df[df["SNP"] == snp]
    if hit.empty:
        return (None, None)
    beta = pd.to_numeric(hit.iloc[0].get("BETA", None), errors="coerce")
    p = pd.to_numeric(hit.iloc[0].get("P", None), errors="coerce")
    return (None if pd.isna(beta) else float(beta),
            None if pd.isna(p) else float(p))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--signals", required=True)
    ap.add_argument("--runs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--erap1-cond-mode", choices=["sig1","all"], default="sig1",
                    help="label ERAP1 SNP-conditioning as sig1 or sig1+sig2+sig3")
    args = ap.parse_args()

    signals = read_signals(args.signals)
    runs = read_runs(args.runs)

    # lead SNPs for each outcome gene = signal1 lead
    def lead_of(g):
        x = signals[(signals["gene"]==g) & (signals["signal_id"]=="signal1")]
        if x.empty:
            raise ValueError(f"missing signal1 for {g} in {args.signals}")
        return str(x.iloc[0]["lead"])

    lead_ERAP2 = lead_of("ERAP2")
    lead_ERAP1 = lead_of("ERAP1")
    lead_LNPEP = lead_of("LNPEP")

    # if ERAP1_COND_MODE=all, collect extra leads for label
    erap1_extra = ""
    if args.erap1_cond_mode == "all":
        s2 = signals[(signals["gene"]=="ERAP1") & (signals["signal_id"]=="signal2")]
        s3 = signals[(signals["gene"]=="ERAP1") & (signals["signal_id"]=="signal3")]
        if not s2.empty and not s3.empty:
            erap1_extra = f"+{s2.iloc[0]['lead']}+{s3.iloc[0]['lead']}"

    # conditioning label mapping
    def cond_label(run_id: str):
        if run_id.startswith("cond_snp_"):
            snp = run_id.replace("cond_snp_","",1)
            if snp == lead_ERAP2:
                return f"ERAP2 top ({lead_ERAP2})"
            if snp == lead_ERAP1:
                if args.erap1_cond_mode == "all":
                    return f"ERAP1 signals ({lead_ERAP1}{erap1_extra})"
                return f"ERAP1 top ({lead_ERAP1})"
            if snp == lead_LNPEP:
                return f"LNPEP top ({lead_LNPEP})"
            return f"SNP ({snp})"
        if run_id == "cov_expr_ERAP2": return "ERAP2 expression"
        if run_id == "cov_expr_ERAP1": return "ERAP1 expression"
        if run_id == "cov_expr_LNPEP": return "LNPEP expression"
        if run_id == "cov_expr_self":  return "self expression (covariate)"
        if run_id == "baseline":       return "baseline"
        return run_id

    out_rows = []
    for outcome in ["ERAP2","ERAP1","LNPEP"]:
        lead = {"ERAP2":lead_ERAP2, "ERAP1":lead_ERAP1, "LNPEP":lead_LNPEP}[outcome]

        base = runs[(runs["outcome"]==outcome) & (runs["run_id"]=="baseline")]
        if base.empty:
            raise ValueError(f"missing baseline for {outcome} in {args.runs}")
        base_assoc = base.iloc[0]["assoc"]
        b_beta, b_p = extract_beta_p(base_assoc, lead)

        # compare against the 6 cross-condition runs (excluding baseline/self)
        for rid in ["cond_snp_"+lead_ERAP2, "cond_snp_"+lead_ERAP1, "cond_snp_"+lead_LNPEP,
                    "cov_expr_ERAP2", "cov_expr_ERAP1", "cov_expr_LNPEP"]:
            rr = runs[(runs["outcome"]==outcome) & (runs["run_id"]==rid)]
            if rr.empty:
                continue
            a_beta, a_p = extract_beta_p(rr.iloc[0]["assoc"], lead)

            # percent change (guard div0)
            def pct(new, old):
                if old is None or new is None or old == 0:
                    return None
                return 100.0*(new-old)/old

            def dneglogp(newp, oldp):
                import math
                if oldp is None or newp is None or oldp<=0 or newp<=0:
                    return None
                return 100.0*((-math.log10(newp)) - (-math.log10(oldp))) / (-math.log10(oldp))

            out_rows.append({
                "Outcome_gene": outcome,
                "Conditioning / covariate": cond_label(rid),
                "lead_snp_tested": lead,
                "baseline_beta": b_beta,
                "baseline_p": b_p,
                "cond_beta": a_beta,
                "cond_p": a_p,
                "Δβ (%)": pct(a_beta, b_beta),
                "Δ−log10(P) (%)": dneglogp(a_p, b_p),
            })

    out = pd.DataFrame(out_rows)
    out.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    main()

