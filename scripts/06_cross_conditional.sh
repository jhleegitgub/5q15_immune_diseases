#!/usr/bin/env bash
# scripts/06_cross_conditional.sh
# Cross-conditional cis-eQTL runs for ERAP2/ERAP1/LNPEP at 5q15 (±WIN)
# Supports: ERAP1_COND_MODE=sig1 (default) or all (sig1+sig2+sig3)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# optional project config
if [[ -f "$SCRIPT_DIR/_config.sh" ]]; then
  # shellcheck disable=SC1091
  source "$SCRIPT_DIR/_config.sh"
fi

# ---- defaults (env override OK) ----
RESULTDIR="${RESULTDIR:-$ROOT/result}"
CODEDIR="${CODEDIR:-$ROOT/code}"

BFILE="${BFILE:-$ROOT/GenotypeData/GW.E-GEUV-3.EUR.MAF005.HWE1e-06}"
PHENO="${PHENO:-${PHENO5:-$ROOT/PhenotypeData/chr5_GD462.signalGeneQuantRPKM_plink.txt}}"
COVAR="${COVAR:-$RESULTDIR/01_pca/covar_pca10.tsv}"
COVAR_NAMES="${COVAR_NAMES:-C1 C2 C3 C4 C5 C6 C7 C8 C9 C10}"

PLINK="${PLINK:-$HOME/Software/Plink/plink}"
WIN="${WIN:-500000}"   # bp

SIGNALS_RAW="${SIGNALS_RAW:-$RESULTDIR/03_signal_check_5q15/signals_summary.tsv}"
ERAP1_COND_MODE="${ERAP1_COND_MODE:-sig1}"   # sig1 | all

MAKE_COVAR_EXPR_PY="${MAKE_COVAR_EXPR_PY:-$CODEDIR/make_covar_plus_expr.py}"

OUTDIR="${OUTDIR:-$RESULTDIR/06_cross_conditional}"
mkdir -p "$OUTDIR"

die(){ echo "[ERR] $*" 1>&2; exit 1; }
need(){ [[ -e "$1" ]] || die "missing: $1"; }

[[ -x "$PLINK" ]] || die "PLINK not executable: $PLINK"
need "${BFILE}.bed"; need "${BFILE}.bim"; need "${BFILE}.fam"
need "$PHENO"; need "$COVAR"
need "$SIGNALS_RAW"
need "$MAKE_COVAR_EXPR_PY"

# ---- 0) normalize signals_summary.tsv -> signals_min.tsv (gene/signal_id/lead only) ----
SIGNALS_MIN="$OUTDIR/signals_min.tsv"
awk -v OFS="\t" '
  NR==1{
    # find columns
    for(i=1;i<=NF;i++){
      if($i=="gene") g=i
      if($i=="signal_id") s=i
      if($i=="lead" || $i=="snp" || $i=="SNP") l=i
    }
    if(!g) g=1; if(!s) s=2; if(!l) l=3;
    print "gene","signal_id","lead";
    next
  }
  {print $g,$s,$l}
' "$SIGNALS_RAW" > "$SIGNALS_MIN"

get_lead(){
  local gene="$1" sid="$2"
  awk -F'\t' -v g="$gene" -v s="$sid" '
    NR==1{next}
    $1==g && $2==s {print $3; exit}
  ' "$SIGNALS_MIN"
}

ERAP2_LEAD="$(get_lead ERAP2 signal1 || true)"
LNPEP_LEAD="$(get_lead LNPEP signal1 || true)"
ERAP1_S1="$(get_lead ERAP1 signal1 || true)"
ERAP1_S2="$(get_lead ERAP1 signal2 || true)"
ERAP1_S3="$(get_lead ERAP1 signal3 || true)"

[[ -n "${ERAP2_LEAD:-}" && -n "${LNPEP_LEAD:-}" && -n "${ERAP1_S1:-}" ]] \
  || die "failed to read leads from $SIGNALS_RAW (need ERAP2/LNPEP signal1 + ERAP1 signal1). Check $SIGNALS_MIN"

if [[ "$ERAP1_COND_MODE" == "all" ]]; then
  [[ -n "${ERAP1_S2:-}" && -n "${ERAP1_S3:-}" ]] \
    || die "ERAP1_COND_MODE=all needs ERAP1 signal2+signal3 in signals_summary.tsv"
fi

# ---- 1) region: centered on ERAP2 lead (same x-range for all panels) ----
get_chr_bp(){
  local snp="$1"
  awk -v s="$snp" '$2==s{print $1, $4; exit}' "${BFILE}.bim"
}

read -r CHR CENTER_BP < <(get_chr_bp "$ERAP2_LEAD" || true)
[[ -n "${CHR:-}" && -n "${CENTER_BP:-}" ]] || die "SNP not found in BIM: $ERAP2_LEAD"

FROM=$((CENTER_BP - WIN)); ((FROM<1)) && FROM=1
TO=$((CENTER_BP + WIN))

# ---- 2) covar+expr files (created once) ----
COVAR_E2="$OUTDIR/covar_plus_expr_ERAP2.tsv"
COVAR_E1="$OUTDIR/covar_plus_expr_ERAP1.tsv"
COVAR_LN="$OUTDIR/covar_plus_expr_LNPEP.tsv"

python3 "$MAKE_COVAR_EXPR_PY" --covar "$COVAR" --pheno "$PHENO" --gene ERAP2 --out "$COVAR_E2"
python3 "$MAKE_COVAR_EXPR_PY" --covar "$COVAR" --pheno "$PHENO" --gene ERAP1 --out "$COVAR_E1"
python3 "$MAKE_COVAR_EXPR_PY" --covar "$COVAR" --pheno "$PHENO" --gene LNPEP --out "$COVAR_LN"

# ---- 3) ERAP1 condition list (if mode=all) ----
ERAP1_CONDLIST="$OUTDIR/cond_ERAP1_sig123.list"
if [[ "$ERAP1_COND_MODE" == "all" ]]; then
  printf "%s\n%s\n%s\n" "$ERAP1_S1" "$ERAP1_S2" "$ERAP1_S3" > "$ERAP1_CONDLIST"
fi

run_plink(){
  local outcome="$1"
  local run_id="$2"
  shift 2

  local od="$OUTDIR/$outcome"
  mkdir -p "$od"
  local pref="$od/$run_id"
  local assoc="${pref}.assoc.linear"

  if [[ -s "$assoc" ]]; then
    echo "[SKIP] $outcome $run_id"
    return 0
  fi

  "$PLINK" \
    --bfile "$BFILE" \
    --chr "$CHR" --from-bp "$FROM" --to-bp "$TO" \
    --pheno "$PHENO" --pheno-name "$outcome" \
    "$@" \
    --linear hide-covar --allow-no-sex \
    --out "$pref" >/dev/null

  [[ -s "$assoc" ]] || die "PLINK failed to produce: $assoc"
  echo "[OK] $assoc"
}

# ---- 4) build runs.tsv (for R + summary) ----
RUNS_TSV="$OUTDIR/runs.tsv"
echo -e "run_id\toutcome\tcond_gene\tcov_type\tassoc" > "$RUNS_TSV"

emit_run(){
  local run_id="$1" outcome="$2" cond_gene="$3" cov_type="$4"
  local assoc="$OUTDIR/$outcome/${run_id}.assoc.linear"
  echo -e "${run_id}\t${outcome}\t${cond_gene}\t${cov_type}\t${assoc}" >> "$RUNS_TSV"
}

# outcomes
OUTCOMES=(ERAP2 ERAP1 LNPEP)

for outcome in "${OUTCOMES[@]}"; do
  # baseline
  run_plink "$outcome" "baseline" --covar "$COVAR" --covar-name $COVAR_NAMES
  emit_run "baseline" "$outcome" "NA" "baseline"

  # cond on ERAP2 SNP
  run_plink "$outcome" "cond_snp_${ERAP2_LEAD}" --covar "$COVAR" --covar-name $COVAR_NAMES --condition "$ERAP2_LEAD"
  emit_run "cond_snp_${ERAP2_LEAD}" "$outcome" "ERAP2" "cond_snp"

  # cond on ERAP1 SNP (sig1 or sig123)
  if [[ "$ERAP1_COND_MODE" == "all" ]]; then
    run_plink "$outcome" "cond_snp_${ERAP1_S1}" --covar "$COVAR" --covar-name $COVAR_NAMES --condition-list "$ERAP1_CONDLIST"
  else
    run_plink "$outcome" "cond_snp_${ERAP1_S1}" --covar "$COVAR" --covar-name $COVAR_NAMES --condition "$ERAP1_S1"
  fi
  emit_run "cond_snp_${ERAP1_S1}" "$outcome" "ERAP1" "cond_snp"

  # cond on LNPEP SNP
  run_plink "$outcome" "cond_snp_${LNPEP_LEAD}" --covar "$COVAR" --covar-name $COVAR_NAMES --condition "$LNPEP_LEAD"
  emit_run "cond_snp_${LNPEP_LEAD}" "$outcome" "LNPEP" "cond_snp"

  # covariate: ERAP2 expression
  run_plink "$outcome" "cov_expr_ERAP2" --covar "$COVAR_E2" --covar-name $COVAR_NAMES expr_ERAP2
  emit_run "cov_expr_ERAP2" "$outcome" "ERAP2" "cov_expr"

  # covariate: ERAP1 expression
  run_plink "$outcome" "cov_expr_ERAP1" --covar "$COVAR_E1" --covar-name $COVAR_NAMES expr_ERAP1
  emit_run "cov_expr_ERAP1" "$outcome" "ERAP1" "cov_expr"

  # covariate: LNPEP expression
  run_plink "$outcome" "cov_expr_LNPEP" --covar "$COVAR_LN" --covar-name $COVAR_NAMES expr_LNPEP
  emit_run "cov_expr_LNPEP" "$outcome" "LNPEP" "cov_expr"

  # covariate: self expression (optional panel)
  case "$outcome" in
    ERAP2) run_plink "$outcome" "cov_expr_self" --covar "$COVAR_E2" --covar-name $COVAR_NAMES expr_ERAP2 ;;
    ERAP1) run_plink "$outcome" "cov_expr_self" --covar "$COVAR_E1" --covar-name $COVAR_NAMES expr_ERAP1 ;;
    LNPEP) run_plink "$outcome" "cov_expr_self" --covar "$COVAR_LN" --covar-name $COVAR_NAMES expr_LNPEP ;;
  esac
  emit_run "cov_expr_self" "$outcome" "$outcome" "cov_expr_self"
done

# ---- 5) LD to ref SNPs (for coloring by column) ----
WINKB=$((WIN/1000))

make_ld(){
  local tag="$1" ref="$2"
  local outpref="$OUTDIR/ld_to_${tag}"
  local outgz="${outpref}.ld.gz"
  if [[ -s "$outgz" ]]; then
    echo "[SKIP] LD $tag"
    return 0
  fi
  "$PLINK" \
    --bfile "$BFILE" \
    --chr "$CHR" --from-bp "$FROM" --to-bp "$TO" \
    --ld-snp "$ref" \
    --r2 gz --ld-window 99999 --ld-window-kb "$WINKB" --ld-window-r2 0 \
    --out "$outpref" >/dev/null
  [[ -s "$outgz" ]] || die "failed LD: $outgz"
  echo "[OK] $outgz"
}

make_ld "ERAP2" "$ERAP2_LEAD"
make_ld "ERAP1" "$ERAP1_S1"
make_ld "LNPEP" "$LNPEP_LEAD"

echo "[OK] runs: $RUNS_TSV"
echo "[OK] signals_min: $SIGNALS_MIN"
echo "[OK] LD: $OUTDIR/ld_to_{ERAP2,ERAP1,LNPEP}.ld.gz"
