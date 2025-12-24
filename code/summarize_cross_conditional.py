#!/usr/bin/env python3
import sys, csv, math
from collections import defaultdict

def read_tsv(path):
    with open(path, "r", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        rows = list(r)
        return rows

def read_assoc_linear(path):
    # PLINK .assoc.linear (whitespace-delimited)
    rows = []
    with open(path, "r") as f:
        header = f.readline().split()
        idx = {k: i for i, k in enumerate(header)}
        for line in f:
            sp = line.split()
            if not sp:
                continue
            d = {k: sp[idx[k]] for k in idx.keys() if idx[k] < len(sp)}
            rows.append(d)
    return rows

def safe_float(x):
    try:
        if x is None:
            return None
        if x in ("NA", "nan", "NaN", ".", ""):
            return None
        return float(x)
    except Exception:
        return None

def neglog10p(p):
    if p is None or p <= 0:
        return None
    return -math.log10(p)

def main():
    if len(sys.argv) != 4:
        print("usage: summarize_cross_conditional.py <signals_summary.tsv> <runs.tsv> <out.tsv>", file=sys.stderr)
        sys.exit(2)

    signals_path, runs_path, out_path = sys.argv[1:]

    signals = read_tsv(signals_path)
    runs = read_tsv(runs_path)

    # Expect columns in signals_summary.tsv: gene, signal_id, lead (at least)
    gene_leads = defaultdict(list)
    for s in signals:
        g = s.get("gene")
        sid = s.get("signal_id")
        lead = s.get("lead")
        if g and sid and lead:
            gene_leads[g].append((sid, lead))

    # cache assoc parses
    assoc_cache = {}
    for r in runs:
        ap = r["assoc_path"]
        if ap not in assoc_cache:
            assoc_cache[ap] = read_assoc_linear(ap)

    # baseline path per outcome
    baseline_path = {}
    for r in runs:
        if r["run_id"] == "baseline":
            baseline_path[r["outcome"]] = r["assoc_path"]

    def index_add(rows):
        m = {}
        for d in rows:
            if d.get("TEST") != "ADD":
                continue
            snp = d.get("SNP")
            if snp:
                m[snp] = d
        return m

    out_fields = [
        "outcome","signal_id","lead_snp",
        "run_id","cond_type","cond_label",
        "beta_base","p_base","neglog10p_base",
        "beta_cond","p_cond","neglog10p_cond",
        "delta_beta_pct","delta_neglog10p_pct",
        "assoc_path"
    ]

    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=out_fields)
        w.writeheader()

        for r in runs:
            outcome = r["outcome"]
            if outcome not in baseline_path:
                continue

            base_rows = assoc_cache[baseline_path[outcome]]
            cond_rows = assoc_cache[r["assoc_path"]]
            base_m = index_add(base_rows)
            cond_m = index_add(cond_rows)

            for sig_id, lead in gene_leads.get(outcome, []):
                b = base_m.get(lead)
                c = cond_m.get(lead)

                beta_b = safe_float(b.get("BETA") if b else None)
                p_b = safe_float(b.get("P") if b else None)
                nl_b = neglog10p(p_b)

                beta_c = safe_float(c.get("BETA") if c else None)
                p_c = safe_float(c.get("P") if c else None)
                nl_c = neglog10p(p_c)

                db = None
                dn = None
                if beta_b is not None and beta_c is not None and beta_b != 0:
                    db = 100.0 * (beta_c - beta_b) / beta_b
                if nl_b is not None and nl_c is not None and nl_b != 0:
                    dn = 100.0 * (nl_c - nl_b) / nl_b

                w.writerow({
                    "outcome": outcome,
                    "signal_id": sig_id,
                    "lead_snp": lead,
                    "run_id": r["run_id"],
                    "cond_type": r["cond_type"],
                    "cond_label": r["cond_label"],
                    "beta_base": beta_b,
                    "p_base": p_b,
                    "neglog10p_base": nl_b,
                    "beta_cond": beta_c,
                    "p_cond": p_c,
                    "neglog10p_cond": nl_c,
                    "delta_beta_pct": db,
                    "delta_neglog10p_pct": dn,
                    "assoc_path": r["assoc_path"]
                })

if __name__ == "__main__":
    main()
