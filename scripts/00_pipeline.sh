#!/usr/bin/env bash
# scripts/00_pipeline.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_config.sh"

log(){ echo "[$(date '+%F %T')] $*"; }

# outputs
OUT_PCA="$RESULTDIR/01_pca"
OUT_GW="$RESULTDIR/02_eqtl_genomewide"
OUT_SIG="$RESULTDIR/03_signal_check_5q15"
OUT_FM="$RESULTDIR/04_finemap_pip"
OUT_FA="$RESULTDIR/05_functional_annot"
OUT_CC="$RESULTDIR/06_cross_conditional"

mkdir -p "$OUT_PCA" "$OUT_GW" "$OUT_SIG" "$OUT_FM" "$OUT_FA" "$OUT_CC"
mkdir -p "$FIGDIR" "$TABLEDIR"

# finemap priors
PRIOR_MULT="${PRIOR_MULT:-0.15}"
PRIOR_MULT_LIST="${PRIOR_MULT_LIST:-0.05 0.10 0.15 0.20 0.30}"
WIN="${WIN:-500000}"
WINKB="${WINKB:-500}"

log "01) PCA covariates"
OUTDIR="$OUT_PCA" bash "$SCRIPTDIR/01_make_pca_covar.sh"

log "02) Genome-wide eQTL"
OUTDIR="$OUT_GW" bash "$SCRIPTDIR/02_genomewide_eqtl.sh"

log "03) 5q15 signal check"
OUTDIR="$OUT_SIG" bash "$SCRIPTDIR/03_signal_check_5q15.sh"

log "04) Fine-mapping (ABF->PIP)"
OUTDIR="$OUT_FM" bash "$SCRIPTDIR/04_finemap_pip_run.sh" \
  "$BFILE" \
  "$PHENO5" \
  "$OUT_PCA/covar_pca10.tsv" \
  "$PLINK" \
  "$COVAR_NAMES" \
  "$WIN" \
  "$ERAP2_LEAD" \
  "$LNPEP_LEAD" \
  "$ERAP1_S1" \
  "$ERAP1_S2" \
  "$ERAP1_S3" \
  "$PRIOR_MULT" \
  "$PRIOR_MULT_LIST"

log "05) Functional annotation"
OUTDIR="$OUT_FA" \
FINEMAP_DIR="$OUT_FM" \
SIGNALS_TSV="$OUT_SIG/signals_summary.tsv" \
bash "$SCRIPTDIR/05_functional_annot.sh"

log "06) Cross-conditional"
OUTDIR="$OUT_CC" \
SIGNALS_TSV="$OUT_SIG/signals_summary.tsv" \
COVAR="$OUT_PCA/covar_pca10.tsv" \
bash "$SCRIPTDIR/06_cross_conditional.sh"

log "07) 27-panel figure"
CCDIR="$OUT_CC" \
SIGNALS_TSV="$OUT_SIG/signals_summary.tsv" \
bash "$SCRIPTDIR/07_make_crossconditional_figs.sh"

log "DONE"
log "result: $RESULTDIR"
log "fig   : $FIGDIR"
log "table : $TABLEDIR"
