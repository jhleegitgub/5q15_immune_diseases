#!/usr/bin/env python3
"""Summarise cross-conditional attenuation at lead SNP(s).

This script is intentionally forgiving about input column names because
people often change runs.tsv formats during pipeline iterations.

Inputs
------
signals_summary.tsv : must contain at least gene + signal_id + lead SNP
runs.tsv : must contain at least outcome + run_id + assoc path

Accepted column aliases
-----------------------
runs.tsv assoc path : assoc_path OR assoc OR path

Output
------
TSV with per (outcome, signal_id, lead_snp, run_id) delta% of beta and -log10(p).
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


def _read_table(path: str) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    return pd.read_table(path)


def _pick_col(df: pd.DataFrame, candidates: List[str], label: str, path: str) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    raise ValueError(f"{path} missing column for {label}. Tried: {candidates}. Found: {list(df.columns)}")


def _read_signals(path: str) -> pd.DataFrame:
    df = _read_table(path)
    gene_col = _pick_col(df, ["gene", "pheno_gene", "outcome_gene"], "gene", path)
    sig_col = _pick_col(df, ["signal_id", "signal", "sig", "signal_idx"], "signal_id", path)
    lead_col = _pick_col(df, ["lead", "lead_snp", "lead_rsid", "rsid"], "lead", path)
    out = df[[gene_col, sig_col, lead_col]].copy()
    out.columns = ["gene", "signal_id", "lead"]
    out["gene"] = out["gene"].astype(str)
    out["signal_id"] = out["signal_id"].astype(str)
    out["lead"] = out["lead"].astype(str)
    return out


def _read_runs(path: str) -> pd.DataFrame:
    df = _read_table(path)
    outcome_col = _pick_col(df, ["outcome", "pheno", "pheno_gene"], "outcome", path)
    runid_col = _pick_col(df, ["run_id", "run", "label"], "run_id", path)
    assoc_col = _pick_col(df, ["assoc_path", "assoc", "path"], "assoc path", path)

    # Optional columns (nice to have)
    cond_type = df["cond_type"].astype(str) if "cond_type" in df.columns else "NA"
    cond_label = df["cond_label"].astype(str) if "cond_label" in df.columns else "NA"

    out = pd.DataFrame(
        {
            "outcome": df[outcome_col].astype(str),
            "run_id": df[runid_col].astype(str),
            "cond_type": cond_type,
            "cond_label": cond_label,
            "assoc_path": df[assoc_col].astype(str),
        }
    )
    return out


def _read_plink_assoc(path: str) -> pd.DataFrame:
    df = pd.read_table(path, sep=r"\s+", engine="python")
    if "TEST" not in df.columns or "SNP" not in df.columns:
        raise ValueError(f"Unexpected assoc file format: {path}")
    df = df[df["TEST"] == "ADD"].copy()
    for c in ["P", "BETA", "BP", "STAT"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def _extract_at_snp(assoc: pd.DataFrame, snp: str) -> Dict[str, float]:
    row = assoc.loc[assoc["SNP"] == snp]
    if row.empty:
        return {"P": float("nan"), "BETA": float("nan"), "BP": float("nan"), "STAT": float("nan")}
    row = row.iloc[0]
    return {
        "P": float(row["P"]) if "P" in row else float("nan"),
        "BETA": float(row["BETA"]) if "BETA" in row else float("nan"),
        "BP": float(row["BP"]) if "BP" in row else float("nan"),
        "STAT": float(row["STAT"]) if "STAT" in row else float("nan"),
    }


def _pct_change(new: float, base: float) -> float:
    if pd.isna(new) or pd.isna(base) or base == 0:
        return float("nan")
    return 100.0 * (new - base) / base


def _safe_logp(p: float) -> float:
    if pd.isna(p) or p <= 0:
        return float("nan")
    return float(-np.log10(p))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--signals", required=True, help="signals_summary.tsv")
    ap.add_argument("--runs", required=True, help="runs.tsv")
    ap.add_argument("--out", required=True, help="output tsv")
    ap.add_argument(
        "--baseline_run_id",
        default="baseline",
        help="run_id that represents baseline per outcome (default: baseline)",
    )
    args = ap.parse_args()

    sig = _read_signals(args.signals)
    runs = _read_runs(args.runs)

    # Preload assoc files to avoid repeated IO
    assoc_cache: Dict[str, pd.DataFrame] = {}
    for p in sorted(set(runs["assoc_path"])):
        if not Path(p).exists():
            raise FileNotFoundError(f"Missing assoc file in runs.tsv: {p}")
        assoc_cache[p] = _read_plink_assoc(p)

    out_rows = []
    for outcome, sub in runs.groupby("outcome", sort=False):
        base_rows = sub[sub["run_id"] == args.baseline_run_id]
        if base_rows.empty:
            raise ValueError(
                f"No baseline run for outcome={outcome}. Expected run_id='{args.baseline_run_id}'."
            )
        base_path = base_rows.iloc[0]["assoc_path"]
        base_assoc = assoc_cache[base_path]

        # lead SNP(s) for this outcome: include all signals where gene==outcome
        leads = sig[sig["gene"] == outcome]
        if leads.empty:
            continue

        for _, lead_row in leads.iterrows():
            signal_id = lead_row["signal_id"]
            lead_snp = lead_row["lead"]
            base_val = _extract_at_snp(base_assoc, lead_snp)
            base_beta = base_val["BETA"]
            base_logp = _safe_logp(base_val["P"])

            for _, r in sub.iterrows():
                assoc = assoc_cache[r["assoc_path"]]
                val = _extract_at_snp(assoc, lead_snp)
                logp = _safe_logp(val["P"])

                out_rows.append(
                    {
                        "outcome": outcome,
                        "signal_id": signal_id,
                        "lead_snp": lead_snp,
                        "run_id": r["run_id"],
                        "cond_type": r["cond_type"],
                        "cond_label": r["cond_label"],
                        "beta_base": base_beta,
                        "beta": val["BETA"],
                        "logp_base": base_logp,
                        "logp": logp,
                        "delta_beta_pct": _pct_change(val["BETA"], base_beta),
                        "delta_logp_pct": _pct_change(logp, base_logp),
                    }
                )

    out_df = pd.DataFrame(out_rows)
    out_df.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
