#!/usr/bin/env bash
set -euo pipefail

# this file lives in scripts/
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

SCRIPTDIR="$ROOT/scripts"
CODEDIR="$ROOT/code"

RESULTDIR="$ROOT/result"
FIGDIR="$ROOT/fig"
TABLEDIR="$ROOT/table"

mkdir -p "$RESULTDIR" "$FIGDIR" "$TABLEDIR"

# tools / inputs (env override OK)
PLINK="${PLINK:-$HOME/Software/Plink/plink}"

# adjust these if your project layout differs
BFILE="${BFILE:-$ROOT/../GenotypeData/GW.E-GEUV-3.EUR.MAF005.HWE1e-06}"
PHENO5="${PHENO5:-$ROOT/../PhenotypeData/chr5_GD462.signalGeneQuantRPKM_plink.txt}"

# default covar output from step01
COVAR="${COVAR:-$RESULTDIR/01_pca/covar_pca10.tsv}"
COVAR_NAMES="${COVAR_NAMES:-C1 C2 C3 C4 C5 C6 C7 C8 C9 C10}"

# window size (bp for fine-map/signals; kb for cross-conditional/locus)
WIN="${WIN:-500000}"
WINKB="${WINKB:-500}"

# leads (override if needed)
ERAP2_LEAD="${ERAP2_LEAD:-rs2910686}"
LNPEP_LEAD="${LNPEP_LEAD:-rs248215}"
ERAP1_S1="${ERAP1_S1:-rs30379}"
ERAP1_S2="${ERAP1_S2:-rs27039}"
ERAP1_S3="${ERAP1_S3:-rs1065407}"
