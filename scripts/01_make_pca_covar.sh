#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_config.sh"

OUTDIR="${OUTDIR:-$RESULTDIR/01_pca}"
mkdir -p "$OUTDIR"

[[ -x "$PLINK" ]] || { echo "[ERR] PLINK not executable: $PLINK" >&2; exit 1; }
[[ -s "${BFILE}.bed" ]] || { echo "[ERR] BFILE not found: $BFILE(.bed/.bim/.fam)" >&2; exit 1; }

$PLINK --bfile "$BFILE" --indep-pairwise 50 5 0.2 --out "$OUTDIR/indepSNP"
$PLINK --bfile "$BFILE" --extract "$OUTDIR/indepSNP.prune.in" --pca 10 --out "$OUTDIR/pca10"

awk 'BEGIN{
  OFS="\t";
  print "FID","IID","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10"
}
NR==1 && ($1=="FID" || $1=="#FID"){next}
{
  print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12
}' "$OUTDIR/pca10.eigenvec" > "$OUTDIR/covar_pca10.tsv"

echo "[OK] wrote $OUTDIR/covar_pca10.tsv"
