#!/usr/bin/env bash
# scripts/06_cross_conditional.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_config.sh"

SIGNALS_TSV="${SIGNALS_TSV:-$RESULTDIR/03_signal_check_5q15/signals_summary.tsv}"
CENTER_SNP="${CENTER_SNP:-$ERAP2_LEAD}"
OUTDIR="${OUTDIR:-$RESULTDIR/06_cross_conditional}"

MAKE_COVAR_PY="${MAKE_COVAR_PY:-$CODEDIR/make_covar_plus_expr.py}"
SUM_CC_PY="${SUM_CC_PY:-$CODEDIR/summarize_cross_conditional.py}"

mkdir -p "$OUTDIR"/{assoc,tmp,covar}

die(){ echo "[ERR] $*" 1>&2; exit 1; }
need(){ [[ -e "$1" ]] || die "missing: $1"; }

[[ -x "$PLINK" ]] || die "PLINK not executable: $PLINK"
need "${BFILE}.bim"
need "$PHENO5"
need "$COVAR"
need "$SIGNALS_TSV"
need "$MAKE_COVAR_PY"
need "$SUM_CC_PY"

# ---------------------------------------
# region window centered at CENTER_SNP
# ---------------------------------------
get_chr_bp() { awk -v s="$1" '$2==s {print $1, $4; exit}' "${BFILE}.bim"; }

read -r CHR0 BP0 < <(get_chr_bp "$CENTER_SNP")
[[ -n "${CHR0:-}" && -n "${BP0:-}" ]] || die "CENTER_SNP not found in BIM: $CENTER_SNP"

FROM=$((BP0 - WINKB*1000)); TO=$((BP0 + WINKB*1000)); ((FROM<1)) && FROM=1

# ---------------------------------------
# load signals
# ---------------------------------------
mapfile -t ERAP2_LEADS < <(awk -F'\t' 'NR>1 && $1=="ERAP2"{print $3}' "$SIGNALS_TSV" | sed '/^$/d' | sort -u)
mapfile -t ERAP1_LEADS < <(awk -F'\t' 'NR>1 && $1=="ERAP1"{print $3}' "$SIGNALS_TSV" | sed '/^$/d' | sort -u)
mapfile -t LNPEP_LEADS < <(awk -F'\t' 'NR>1 && $1=="LNPEP"{print $3}' "$SIGNALS_TSV" | sed '/^$/d' | sort -u)

RUNS_TSV="$OUTDIR/runs.tsv"
echo -e "outcome\trun_id\tcond_type\tcond_label\tassoc_path" > "$RUNS_TSV"

# ---------------------------------------
# helper: run plink in window
# IMPORTANT: no --covar-name; use all columns in covar file
# ---------------------------------------
plink_run() {
  local outcome="$1"
  local run_id="$2"
  local cond_type="$3"
  local cond_label="$4"
  local covar_file="$5"
  shift 5
  local -a cond_snps=("$@")

  local outpref="$OUTDIR/assoc/${outcome}.${run_id}"
  local assoc="${outpref}.assoc.linear"

  local -a cond_args=()
  if [[ ${#cond_snps[@]} -gt 0 ]]; then
    local cf="$OUTDIR/tmp/cond_${outcome}_${run_id}.txt"
    printf "%s\n" "${cond_snps[@]}" > "$cf"
    cond_args=(--condition-list "$cf")
  fi

  "$PLINK" \
    --bfile "$BFILE" \
    --chr "$CHR0" --from-bp "$FROM" --to-bp "$TO" \
    --pheno "$PHENO5" --pheno-name "$outcome" \
    --covar "$covar_file" \
    --linear hide-covar --allow-no-sex \
    "${cond_args[@]}" \
    --out "$outpref" >/dev/null

  echo -e "${outcome}\t${run_id}\t${cond_type}\t${cond_label}\t${assoc}" >> "$RUNS_TSV"
}

# covar + expression column
make_expr_covar() {
  local expr_gene="$1"
  local out="$OUTDIR/covar/covar_plus_expr_${expr_gene}.tsv"
  python3 "$MAKE_COVAR_PY" "$COVAR" "$PHENO5" "$expr_gene" "$out"
  echo "$out"
}

# ---------------------------------------
# baseline
# ---------------------------------------
for g in ERAP2 ERAP1 LNPEP; do
  plink_run "$g" "baseline" "none" "baseline" "$COVAR"
done

# ---------------------------------------
# SNP-conditional (signal SNPs)
# ---------------------------------------
for snp in "${ERAP2_LEADS[@]}"; do plink_run "ERAP2" "cond_snp_${snp}" "snp" "cond:${snp}" "$COVAR" "$snp"; done
for snp in "${ERAP1_LEADS[@]}"; do plink_run "ERAP2" "cond_snp_${snp}" "snp" "cond:${snp}" "$COVAR" "$snp"; done
for snp in "${LNPEP_LEADS[@]}"; do plink_run "ERAP2" "cond_snp_${snp}" "snp" "cond:${snp}" "$COVAR" "$snp"; done

for snp in "${ERAP1_LEADS[@]}"; do plink_run "ERAP1" "cond_snp_${snp}" "snp" "cond:${snp}" "$COVAR" "$snp"; done
for snp in "${ERAP2_LEADS[@]}"; do plink_run "ERAP1" "cond_snp_${snp}" "snp" "cond:${snp}" "$COVAR" "$snp"; done
for snp in "${LNPEP_LEADS[@]}"; do plink_run "ERAP1" "cond_snp_${snp}" "snp" "cond:${snp}" "$COVAR" "$snp"; done
if [[ ${#ERAP1_LEADS[@]} -ge 2 ]]; then
  plink_run "ERAP1" "cond_ERAP1_allSignals" "snp" "cond:ERAP1_allSignals" "$COVAR" "${ERAP1_LEADS[@]}"
fi

for snp in "${LNPEP_LEADS[@]}"; do plink_run "LNPEP" "cond_snp_${snp}" "snp" "cond:${snp}" "$COVAR" "$snp"; done
for snp in "${ERAP2_LEADS[@]}"; do plink_run "LNPEP" "cond_snp_${snp}" "snp" "cond:${snp}" "$COVAR" "$snp"; done
for snp in "${ERAP1_LEADS[@]}"; do plink_run "LNPEP" "cond_snp_${snp}" "snp" "cond:${snp}" "$COVAR" "$snp"; done

# ---------------------------------------
# Expression-covariate runs (27-panel의 expr 조건)
# ---------------------------------------
for expr in ERAP1 LNPEP; do
  covf="$(make_expr_covar "$expr")"
  plink_run "ERAP2" "cov_expr_${expr}" "expr" "covar:expr_${expr}" "$covf"
done
for expr in ERAP2 LNPEP; do
  covf="$(make_expr_covar "$expr")"
  plink_run "ERAP1" "cov_expr_${expr}" "expr" "covar:expr_${expr}" "$covf"
done
for expr in ERAP2 ERAP1; do
  covf="$(make_expr_covar "$expr")"
  plink_run "LNPEP" "cov_expr_${expr}" "expr" "covar:expr_${expr}" "$covf"
done

# ---------------------------------------
# summarize attenuation table
# ---------------------------------------
python3 "$SUM_CC_PY" "$SIGNALS_TSV" "$RUNS_TSV" "$OUTDIR/attenuation_summary.tsv"

echo "[OK] $RUNS_TSV"
echo "[OK] $OUTDIR/attenuation_summary.tsv"
