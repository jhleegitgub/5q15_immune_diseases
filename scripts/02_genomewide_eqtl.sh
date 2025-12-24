#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CODEDIR="$ROOT/code"
RESULTDIR="$ROOT/result/02_eqtl_genomewide"
FIGDIR="$ROOT/fig"
LOGDIR="$RESULTDIR/logs"
mkdir -p "$RESULTDIR" "$FIGDIR" "$LOGDIR"

need() { command -v "$1" >/dev/null 2>&1 || { echo "[ERR] missing: $1" >&2; exit 1; }; }
need awk; need sed; need grep

# ---- inputs ----
BFILE="${1:-${BFILE:-$ROOT/../GenotypeData/GW.E-GEUV-3.EUR.MAF005.HWE1e-06}}"

PHENO="${2:-${PHENO:-}}"
if [[ -z "${PHENO}" ]]; then
  shopt -s nullglob
  CANDS=(
    "$ROOT/PhenotypeData/chr5_GD462.signalGeneQuantRPKM_plink.txt"
    "$ROOT/PhenotypeData/chr5_GD462.signalGeneQuantRPKM_plink_input.txt"
    "$ROOT/PhenotypeData/"*signalGeneQuantRPKM*plink*.txt
    "$ROOT/PhenotypeData/"*signalGeneQuantRPKM*plink*.tsv
    "$ROOT/PhenotypeData/"*pheno*.txt
    "$ROOT/PhenotypeData/"*pheno*.tsv

    "$ROOT/../PhenotypeData/chr5_GD462.signalGeneQuantRPKM_plink.txt"
    "$ROOT/../PhenotypeData/chr5_GD462.signalGeneQuantRPKM_plink_input.txt"
    "$ROOT/../PhenotypeData/"*signalGeneQuantRPKM*plink*.txt
    "$ROOT/../PhenotypeData/"*signalGeneQuantRPKM*plink*.tsv
    "$ROOT/../PhenotypeData/"*pheno*.txt
    "$ROOT/../PhenotypeData/"*pheno*.tsv
  )
  shopt -u nullglob
  for f in "${CANDS[@]}"; do
    if [[ -f "$f" ]]; then PHENO="$f"; break; fi
  done
fi

COVAR_DEFAULT_A="$ROOT/result/01_pca/covar_pca10.tsv"
COVAR_DEFAULT_B="$ROOT/covar_pca10.tsv"
if [[ -f "$COVAR_DEFAULT_A" ]]; then
  COVAR_DEFAULT="$COVAR_DEFAULT_A"
else
  COVAR_DEFAULT="$COVAR_DEFAULT_B"
fi
COVAR="${3:-${COVAR:-$COVAR_DEFAULT}}"

PLINK="${4:-${PLINK:-$HOME/Software/Plink/plink}}"
COVAR_NAMES="${5:-${COVAR_NAMES:-C1 C2 C3 C4 C5 C6 C7 C8 C9 C10}}"

PY="${PY:-python3}"
MANHATTAN_PY="$CODEDIR/manhattan_genomewide.py"

[[ -x "$PLINK" ]] || { echo "[ERR] PLINK not executable: $PLINK" >&2; exit 1; }
[[ -f "${BFILE}.bim" ]] || { echo "[ERR] BFILE not found: ${BFILE}.bim" >&2; exit 1; }
[[ -n "$PHENO" && -f "$PHENO" ]] || { echo "[ERR] PHENO not found: $PHENO" >&2; exit 1; }
[[ -f "$COVAR" ]] || { echo "[ERR] COVAR not found: $COVAR" >&2; exit 1; }
[[ -f "$MANHATTAN_PY" ]] || { echo "[ERR] missing code: $MANHATTAN_PY" >&2; exit 1; }

read -r -a COVAR_ARR <<< "$COVAR_NAMES"

run_plink_eqtl() {
  local gene="$1"
  local outprefix="$2"
  local logfile="$3"

  "$PLINK" \
    --allow-no-sex \
    --bfile "$BFILE" \
    --pheno "$PHENO" --pheno-name "$gene" \
    --covar "$COVAR" --covar-name "${COVAR_ARR[@]}" \
    --linear hide-covar \
    --out "$outprefix" \
    >"$logfile" 2>&1
}

plot_manhattan() {
  local gene="$1"
  local assoc="$2"
  local outpng="$3"
  "$PY" "$MANHATTAN_PY" \
    --assoc "$assoc" \
    --out-png "$outpng" \
    --title "Genome-wide eQTL for ${gene} expression" \
    --gw-threshold 5e-8 \
    --max-points 1200000 \
    --keep-p 1e-5 \
    --seed 1 \
    --point-size 1.0
}

for gene in ERAP2 ERAP1 LNPEP; do
  outprefix="$RESULTDIR/${gene}_genomewide"
  assoc="${outprefix}.assoc.linear"
  log="$LOGDIR/${gene}_genomewide.plink.log"

  if [[ ! -s "$assoc" ]]; then
    echo "[RUN] genome-wide eQTL: $gene"
    run_plink_eqtl "$gene" "$outprefix" "$log"
  else
    echo "[SKIP] exists: $assoc"
  fi

  [[ -s "$assoc" ]] || { echo "[ERR] missing assoc: $assoc" >&2; exit 1; }

  outpng="$FIGDIR/Fig_genomewide_${gene}_manhattan.png"
  echo "[RUN] plot: $outpng"
  plot_manhattan "$gene" "$assoc" "$outpng"
done

# combined 3-panel
COMBO="$FIGDIR/Fig_genomewide_3genes_manhattan.png"
echo "[RUN] plot (3-panel): $COMBO"
"$PY" "$MANHATTAN_PY" \
  --assoc "$RESULTDIR/ERAP2_genomewide.assoc.linear,$RESULTDIR/ERAP1_genomewide.assoc.linear,$RESULTDIR/LNPEP_genomewide.assoc.linear" \
  --out-png "$COMBO" \
  --titles "ERAP2,ERAP1,LNPEP" \
  --gw-threshold 5e-8 \
  --max-points 1200000 \
  --keep-p 1e-5 \
  --seed 1 \
  --point-size 1.0

echo "[OK] done"
echo " - PHENO:  $PHENO"
echo " - result: $RESULTDIR"
echo " - figs:   $FIGDIR"
