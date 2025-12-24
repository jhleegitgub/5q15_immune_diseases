#!/usr/bin/env python3
import sys
import time
import json
import requests
from typing import List, Dict, Any, Optional

VEP_URL = "https://grch37.rest.ensembl.org/vep/human/id"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}

# Ensembl REST는 과도 호출하면 429가 뜸 -> 배치 + backoff
BATCH = 150
SLEEP_SEC = 0.25
MAX_RETRY = 6

def read_rsids(path: str) -> List[str]:
    rsids = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            rsids.append(s)
    return sorted(set(rsids))

def pick_gene_symbol(tc: Dict[str, Any]) -> str:
    # prefer "gene_symbol" if present, else "gene_id"
    return tc.get("gene_symbol") or tc.get("gene_id") or ""

def flatten_one(entry: Dict[str, Any]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    out["SNP"] = entry.get("id", "")
    out["most_severe_consequence"] = entry.get("most_severe_consequence", "")

    # transcript_consequences 우선 1개(impact/gene) 뽑기
    tcs = entry.get("transcript_consequences") or []
    impact = ""
    gene = ""
    biotype = ""
    consequence_terms = ""

    if tcs:
        # impact 우선 HIGH>MODERATE>LOW>MODIFIER 순으로 강한 것
        rank = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}
        tcs_sorted = sorted(
            tcs,
            key=lambda x: (rank.get(x.get("impact", ""), 0),
                           1 if x.get("canonical") else 0),
            reverse=True
        )
        tc = tcs_sorted[0]
        impact = tc.get("impact", "") or ""
        gene = pick_gene_symbol(tc)
        biotype = tc.get("biotype", "") or ""
        consequence_terms = ",".join(tc.get("consequence_terms") or []) or ""

        out["impact"] = impact
        out["gene"] = gene
        out["biotype"] = biotype
        out["consequence_terms"] = consequence_terms

        # optional fields (may be absent)
        out["amino_acids"] = tc.get("amino_acids", "") or ""
        out["protein_start"] = tc.get("protein_start", "") or ""
        out["protein_end"] = tc.get("protein_end", "") or ""
        out["sift_prediction"] = tc.get("sift_prediction", "") or ""
        out["polyphen_prediction"] = tc.get("polyphen_prediction", "") or ""
    else:
        out["impact"] = ""
        out["gene"] = ""
        out["biotype"] = ""
        out["consequence_terms"] = ""
        out["amino_acids"] = ""
        out["protein_start"] = ""
        out["protein_end"] = ""
        out["sift_prediction"] = ""
        out["polyphen_prediction"] = ""

    # keep a placeholder column for downstream scripts
    out["CADD_PHRED_max"] = ""  # REST 기본 응답엔 보통 없음(빈칸 유지)
    return out

def post_ids(ids: List[str]) -> List[Dict[str, Any]]:
    payload = {"ids": ids}
    backoff = 1.0
    for attempt in range(MAX_RETRY):
        r = requests.post(VEP_URL, headers=HEADERS, data=json.dumps(payload), timeout=60)
        if r.status_code == 200:
            return r.json()
        if r.status_code in (429, 500, 502, 503, 504):
            time.sleep(backoff)
            backoff *= 1.8
            continue
        raise RuntimeError(f"VEP request failed: status={r.status_code} body={r.text[:300]}")
    raise RuntimeError("VEP request failed after retries (rate limit / server error)")

def write_tsv(rows: List[Dict[str, Any]], out_path: str) -> None:
    cols = [
        "SNP",
        "gene",
        "impact",
        "most_severe_consequence",
        "consequence_terms",
        "biotype",
        "amino_acids",
        "protein_start",
        "protein_end",
        "sift_prediction",
        "polyphen_prediction",
        "CADD_PHRED_max",
    ]
    with open(out_path, "w") as w:
        w.write("\t".join(cols) + "\n")
        for row in rows:
            w.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <rsids.txt> <out.tsv>", file=sys.stderr)
        sys.exit(2)

    rsids_path = sys.argv[1]
    out_tsv = sys.argv[2]

    rsids = read_rsids(rsids_path)
    rows: List[Dict[str, Any]] = []

    for i in range(0, len(rsids), BATCH):
        batch = rsids[i:i+BATCH]
        data = post_ids(batch)
        for entry in data:
            rows.append(flatten_one(entry))
        time.sleep(SLEEP_SEC)

    # stable order
    rows.sort(key=lambda x: x.get("SNP",""))
    write_tsv(rows, out_tsv)

if __name__ == "__main__":
    main()
