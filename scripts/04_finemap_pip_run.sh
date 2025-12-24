#!/usr/bin/env bash
# scripts/04_finemap_pip_run.sh
# - runs PLINK 1Mb window association (once)
# - computes phenotype SD
# - runs fine-mapping (ABF->PIP) with prior_sd = SD * prior_mult
# - writes MAIN prior outputs to OUTDIR root (for downstream steps)
# - writes sensitivity outputs to OUTDIR/sensitivity/mult_<...>/
set -euo pipefail

# ----------------------------
# args (compatible with pipeline)
# ----------------------------
BFILE="${1:-${BFILE:-../GenotypeData/GW.E-GEUV-3.EUR.MAF005.HWE1e-06}}"
PHENO="${2:-${PHENO:-${PHENO5:-../PhenotypeData/chr5_GD462.signalGeneQuantRPKM_plink.txt}}}"
COVAR="${3:-${COVAR:-./covar_pca10.tsv}}"
PLINK="${4:-${PLINK:-$HOME/Software/Plink/plink}}"
COVAR_NAMES="${5:-${COVAR_NAMES:-C1 C2 C3 C4 C5 C6 C7 C8 C9 C10}}"
WIN="${6:-${WIN:-500000}}"

ERAP2_LEAD="${7:-${ERAP2_LEAD:-rs2910686}}"
LNPEP_LEAD="${8:-${LNPEP_LEAD:-rs248215}}"
ERAP1_S1="${9:-${ERAP1_S1:-rs30379}}"
ERAP1_S2="${10:-${ERAP1_S2:-rs27039}}"
ERAP1_S3="${11:-${ERAP1_S3:-rs1065407}}"

PRIOR_MULT="${12:-${PRIOR_MULT:-0.15}}"
PRIOR_MULT_LIST_RAW="${13:-${PRIOR_MULT_LIST:-0.05 0.10 0.15 0.20 0.30}}"
CREDIBLE="${CREDIBLE:-0.95}"

# ----------------------------
# project paths (script-relative)
# ----------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CODEDIR="${CODEDIR:-$ROOT/code}"

OUTDIR="${OUTDIR:-$ROOT/result/04_finemap_pip}"
mkdir -p "$OUTDIR" "$OUTDIR/sensitivity"

CALC_SD_PY="${CALC_SD_PY:-$CODEDIR/calc_pheno_sd_tsv.py}"
FINEMAP_PY="${FINEMAP_PY:-$CODEDIR/finemap_pip.py}"

# ----------------------------
# checks
# ----------------------------
die(){ echo "[ERR] $*" 1>&2; exit 1; }
need(){ [[ -e "$1" ]] || die "missing: $1"; }

[[ -x "$PLINK" ]] || die "PLINK not executable: $PLINK"
need "${BFILE}.bed"; need "${BFILE}.bim"; need "${BFILE}.fam"
need "$PHENO"; need "$COVAR"
need "$CALC_SD_PY"; need "$FINEMAP_PY"

# ----------------------------
# helpers
# ----------------------------
py_mul() { python3 -c "print(float('$1')*float('$2'))"; }

get_chr_bp() {
  local snp="$1"
  awk -v s="$snp" '$2==s {print $1, $4; exit}' "${BFILE}.bim"
}

plink_window_baseline() {
  local gene="$1" center_snp="$2" outprefix="$3"
  local chr bp
  read -r chr bp < <(get_chr_bp "$center_snp" || true)
  [[ -n "${chr:-}" && -n "${bp:-}" ]] || die "SNP not in BIM: $center_snp"
  local from=$((bp - WIN)); local to=$((bp + WIN)); ((from<0)) && from=0

  "$PLINK" --bfile "$BFILE" \
    --chr "$chr" --from-bp "$from" --to-bp "$to" \
    --pheno "$PHENO" --pheno-name "$gene" \
    --covar "$COVAR" --covar-name $COVAR_NAMES \
    --linear hide-covar --allow-no-sex \
    --out "$outprefix" >/dev/null
}

plink_window_condlist() {
  local gene="$1" center_snp="$2" cond_list="$3" outprefix="$4"
  local chr bp
  read -r chr bp < <(get_chr_bp "$center_snp" || true)
  [[ -n "${chr:-}" && -n "${bp:-}" ]] || die "SNP not in BIM: $center_snp"
  local from=$((bp - WIN)); local to=$((bp + WIN)); ((from<0)) && from=0

  "$PLINK" --bfile "$BFILE" \
    --chr "$chr" --from-bp "$from" --to-bp "$to" \
    --pheno "$PHENO" --pheno-name "$gene" \
    --covar "$COVAR" --covar-name $COVAR_NAMES \
    --condition-list "$cond_list" \
    --linear hide-covar --allow-no-sex \
    --out "$outprefix" >/dev/null
}

pip_stats_line() {
  # args: mult label lead_snp prior_sd pip_tsv cred_tsv
  local mult="$1" label="$2" lead="$3" prior_sd="$4" pip_tsv="$5" cred_tsv="$6"
  local top_snp top_pip lead_pip cs_n
  top_snp="$(awk -F'\t' 'NR==2{print $1; exit}' "$pip_tsv" 2>/dev/null || true)"
  top_pip="$(awk -F'\t' 'NR==2{print $11; exit}' "$pip_tsv" 2>/dev/null || true)"
  lead_pip="$(awk -F'\t' -v s="$lead" 'NR==1{next} $1==s{print $11; found=1; exit} END{if(!found) print "NA"}' "$pip_tsv")"
  cs_n="$(awk 'NR>1{n++} END{print (n+0)}' "$cred_tsv")"
  echo -e "${mult}\t${label}\t${lead}\t${prior_sd}\t${lead_pip}\t${top_snp}\t${top_pip}\t${cs_n}"
}

run_finemap_once() {
  # args: assoc prior_sd outprefix out_pip
  local assoc="$1" prior_sd="$2" prefix="$3" outpip="$4"
  python3 "$FINEMAP_PY" \
    --in "$assoc" \
    --out "$outpip" \
    --prefix "$prefix" \
    --prior-sd "$prior_sd" \
    --credible "$CREDIBLE" >/dev/null
}

# ----------------------------
# 0) phenotype SD
# ----------------------------
SD_TSV="$OUTDIR/pheno_sd.tsv"
python3 "$CALC_SD_PY" "$PHENO" ERAP2 ERAP1 LNPEP > "$SD_TSV"

SD_ERAP2="$(awk -F'\t' '$1=="ERAP2"{print $2; exit}' "$SD_TSV")"
SD_ERAP1="$(awk -F'\t' '$1=="ERAP1"{print $2; exit}' "$SD_TSV")"
SD_LNPEP="$(awk -F'\t' '$1=="LNPEP"{print $2; exit}' "$SD_TSV")"
[[ -n "${SD_ERAP2:-}" && -n "${SD_ERAP1:-}" && -n "${SD_LNPEP:-}" ]] || die "failed to compute phenotype SDs"

# ----------------------------
# 1) PLINK association inputs (generated once)
# ----------------------------
ERAP2_PREF="$OUTDIR/ERAP2_base_pm${WIN}"
ERAP2_ASSOC="${ERAP2_PREF}.assoc.linear"
[[ -s "$ERAP2_ASSOC" ]] || { echo "[RUN] PLINK window baseline ERAP2 ($ERAP2_LEAD)"; plink_window_baseline ERAP2 "$ERAP2_LEAD" "$ERAP2_PREF"; }

LNPEP_PREF="$OUTDIR/LNPEP_base_pm${WIN}"
LNPEP_ASSOC="${LNPEP_PREF}.assoc.linear"
[[ -s "$LNPEP_ASSOC" ]] || { echo "[RUN] PLINK window baseline LNPEP ($LNPEP_LEAD)"; plink_window_baseline LNPEP "$LNPEP_LEAD" "$LNPEP_PREF"; }

# ERAP1: isolate each signal by conditioning on the other two
L_S1="$OUTDIR/cond_ERAP1_sig1.txt"; printf "%s\n%s\n" "$ERAP1_S2" "$ERAP1_S3" > "$L_S1"
L_S2="$OUTDIR/cond_ERAP1_sig2.txt"; printf "%s\n%s\n" "$ERAP1_S1" "$ERAP1_S3" > "$L_S2"
L_S3="$OUTDIR/cond_ERAP1_sig3.txt"; printf "%s\n%s\n" "$ERAP1_S1" "$ERAP1_S2" > "$L_S3"

ERAP1_S1_PREF="$OUTDIR/ERAP1_sig1_isolated_pm${WIN}"
ERAP1_S1_ASSOC="${ERAP1_S1_PREF}.assoc.linear"
[[ -s "$ERAP1_S1_ASSOC" ]] || { echo "[RUN] PLINK ERAP1 sig1 isolated"; plink_window_condlist ERAP1 "$ERAP1_S1" "$L_S1" "$ERAP1_S1_PREF"; }

ERAP1_S2_PREF="$OUTDIR/ERAP1_sig2_isolated_pm${WIN}"
ERAP1_S2_ASSOC="${ERAP1_S2_PREF}.assoc.linear"
[[ -s "$ERAP1_S2_ASSOC" ]] || { echo "[RUN] PLINK ERAP1 sig2 isolated"; plink_window_condlist ERAP1 "$ERAP1_S2" "$L_S2" "$ERAP1_S2_PREF"; }

ERAP1_S3_PREF="$OUTDIR/ERAP1_sig3_isolated_pm${WIN}"
ERAP1_S3_ASSOC="${ERAP1_S3_PREF}.assoc.linear"
[[ -s "$ERAP1_S3_ASSOC" ]] || { echo "[RUN] PLINK ERAP1 sig3 isolated"; plink_window_condlist ERAP1 "$ERAP1_S3" "$L_S3" "$ERAP1_S3_PREF"; }

# ----------------------------
# 2) MAIN prior (writes to OUTDIR root for downstream)
# ----------------------------
PRIOR_ERAP2_MAIN="$(py_mul "$SD_ERAP2" "$PRIOR_MULT")"
PRIOR_ERAP1_MAIN="$(py_mul "$SD_ERAP1" "$PRIOR_MULT")"
PRIOR_LNPEP_MAIN="$(py_mul "$SD_LNPEP" "$PRIOR_MULT")"

run_finemap_once "$ERAP2_ASSOC"    "$PRIOR_ERAP2_MAIN" "$OUTDIR/ERAP2"      "$OUTDIR/ERAP2_pip.tsv"
run_finemap_once "$LNPEP_ASSOC"    "$PRIOR_LNPEP_MAIN" "$OUTDIR/LNPEP"      "$OUTDIR/LNPEP_pip.tsv"
run_finemap_once "$ERAP1_S1_ASSOC" "$PRIOR_ERAP1_MAIN" "$OUTDIR/ERAP1_sig1" "$OUTDIR/ERAP1_sig1_pip.tsv"
run_finemap_once "$ERAP1_S2_ASSOC" "$PRIOR_ERAP1_MAIN" "$OUTDIR/ERAP1_sig2" "$OUTDIR/ERAP1_sig2_pip.tsv"
run_finemap_once "$ERAP1_S3_ASSOC" "$PRIOR_ERAP1_MAIN" "$OUTDIR/ERAP1_sig3" "$OUTDIR/ERAP1_sig3_pip.tsv"

echo "[OK] finemap main outputs in $OUTDIR"

# ----------------------------
# 3) sensitivity loop (writes to OUTDIR/sensitivity)
# ----------------------------
# parse list (single string with spaces)
IFS=' ' read -r -a MULTS <<< "$PRIOR_MULT_LIST_RAW"

SENS_SUM="$OUTDIR/sensitivity/prior_sensitivity_summary.tsv"
echo -e "prior_mult\tlabel\tlead_snp\tprior_sd\tlead_pip\ttop_snp\ttop_pip\tcredible_n" > "$SENS_SUM"

for m in "${MULTS[@]}"; do
  [[ -n "${m:-}" ]] || continue
  tag="$(echo "$m" | tr '.' 'p')"              # 0.15 -> 0p15
  sub="$OUTDIR/sensitivity/mult_${tag}"
  mkdir -p "$sub"

  prior_e2="$(py_mul "$SD_ERAP2" "$m")"
  prior_e1="$(py_mul "$SD_ERAP1" "$m")"
  prior_ln="$(py_mul "$SD_LNPEP" "$m")"

  # outputs
  p1="$sub/ERAP2_pip.tsv"
  p2="$sub/LNPEP_pip.tsv"
  p3="$sub/ERAP1_sig1_pip.tsv"
  p4="$sub/ERAP1_sig2_pip.tsv"
  p5="$sub/ERAP1_sig3_pip.tsv"

  run_finemap_once "$ERAP2_ASSOC"    "$prior_e2" "$sub/ERAP2"      "$p1"
  run_finemap_once "$LNPEP_ASSOC"    "$prior_ln" "$sub/LNPEP"      "$p2"
  run_finemap_once "$ERAP1_S1_ASSOC" "$prior_e1" "$sub/ERAP1_sig1" "$p3"
  run_finemap_once "$ERAP1_S2_ASSOC" "$prior_e1" "$sub/ERAP1_sig2" "$p4"
  run_finemap_once "$ERAP1_S3_ASSOC" "$prior_e1" "$sub/ERAP1_sig3" "$p5"

  # collect summary lines
  echo "$(pip_stats_line "$m" "ERAP2"      "$ERAP2_LEAD" "$prior_e2" "$p1" "$sub/ERAP2_credible95.tsv")" >> "$SENS_SUM"
  echo "$(pip_stats_line "$m" "LNPEP"      "$LNPEP_LEAD" "$prior_ln" "$p2" "$sub/LNPEP_credible95.tsv")" >> "$SENS_SUM"
  echo "$(pip_stats_line "$m" "ERAP1_sig1" "$ERAP1_S1"   "$prior_e1" "$p3" "$sub/ERAP1_sig1_credible95.tsv")" >> "$SENS_SUM"
  echo "$(pip_stats_line "$m" "ERAP1_sig2" "$ERAP1_S2"   "$prior_e1" "$p4" "$sub/ERAP1_sig2_credible95.tsv")" >> "$SENS_SUM"
  echo "$(pip_stats_line "$m" "ERAP1_sig3" "$ERAP1_S3"   "$prior_e1" "$p5" "$sub/ERAP1_sig3_credible95.tsv")" >> "$SENS_SUM"
done

echo "[OK] sensitivity summary: $SENS_SUM"
echo "[OK] sensitivity outputs: $OUTDIR/sensitivity/mult_*"
