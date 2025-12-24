#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ===== paths =====
PLINK="${PLINK:-$HOME/Software/Plink/plink}"
PYTHON="${PYTHON:-python3}"

BFILE="${BFILE:-$ROOT/../GenotypeData/GW.E-GEUV-3.EUR.MAF005.HWE1e-06}"
PHENO="${PHENO:-$ROOT/../PhenotypeData/chr5_GD462.signalGeneQuantRPKM_plink.txt}"
COVAR="${COVAR:-$ROOT/result/01_pca/covar_pca10.tsv}"
COVAR_NAMES="${COVAR_NAMES:-C1 C2 C3 C4 C5 C6 C7 C8 C9 C10}"

WINDOW_BP="${WINDOW_BP:-500000}"

OUTDIR="$ROOT/result/03_signal_check_5q15"
FIGDIR="$ROOT/fig"
CODE_LOCUS="$ROOT/code/locuszoom_manhattan.py"

mkdir -p "$OUTDIR" "$FIGDIR"
export MPLBACKEND=Agg

log(){ printf "%s\n" "$*" >&2; }
die(){ log "[ERR] $*"; exit 1; }
req(){ [[ -f "$1" ]] || die "not found: $1"; }

# ===== checks =====
[[ -x "$PLINK" ]] || die "PLINK not executable: $PLINK"
req "${BFILE}.bed"; req "${BFILE}.bim"; req "${BFILE}.fam"
req "$PHENO"; req "$COVAR"; req "$CODE_LOCUS"

get_bp() {
  local rsid="$1"
  local bp
  bp="$(awk -v snp="$rsid" '$2==snp {print $4; exit}' "${BFILE}.bim" || true)"
  [[ -n "$bp" ]] || die "SNP not in BIM: $rsid"
  echo "$bp"
}

region_bounds() {
  local bp="$1"
  local from=$((bp - WINDOW_BP))
  local to=$((bp + WINDOW_BP))
  (( from < 1 )) && from=1
  echo "$from $to"
}

# stdout: ONLY "assoc<TAB>ld"
plink_window() {
  local gene="$1"
  local tag="$2"
  local lead="$3"
  local from_bp="$4"
  local to_bp="$5"
  local cond_mode="${6:-none}"   # none | snp | list
  local cond_arg="${7:-}"        # rsid | file

  local out_prefix="$OUTDIR/${tag}"
  local assoc="${out_prefix}.assoc.linear"
  local ld_prefix="${out_prefix}.ld_to_${lead}"
  local ld_file="${ld_prefix}.ld"

  if [[ ! -f "$assoc" ]]; then
    log "[RUN] PLINK assoc: $tag"
    if [[ "$cond_mode" == "none" ]]; then
      "$PLINK" --bfile "$BFILE" \
        --chr 5 --from-bp "$from_bp" --to-bp "$to_bp" \
        --pheno "$PHENO" --pheno-name "$gene" \
        --covar "$COVAR" --covar-name $COVAR_NAMES \
        --linear hide-covar --allow-no-sex \
        --out "$out_prefix" >/dev/null
    elif [[ "$cond_mode" == "snp" ]]; then
      "$PLINK" --bfile "$BFILE" \
        --chr 5 --from-bp "$from_bp" --to-bp "$to_bp" \
        --pheno "$PHENO" --pheno-name "$gene" \
        --covar "$COVAR" --covar-name $COVAR_NAMES \
        --condition "$cond_arg" \
        --linear hide-covar --allow-no-sex \
        --out "$out_prefix" >/dev/null
    else
      "$PLINK" --bfile "$BFILE" \
        --chr 5 --from-bp "$from_bp" --to-bp "$to_bp" \
        --pheno "$PHENO" --pheno-name "$gene" \
        --covar "$COVAR" --covar-name $COVAR_NAMES \
        --condition-list "$cond_arg" \
        --linear hide-covar --allow-no-sex \
        --out "$out_prefix" >/dev/null
    fi
  else
    log "[SKIP] exists: $assoc"
  fi

  if [[ ! -f "$ld_file" ]]; then
    log "[RUN] PLINK LD r2 to lead: $tag (lead=$lead)"
    "$PLINK" --bfile "$BFILE" \
      --chr 5 --from-bp "$from_bp" --to-bp "$to_bp" \
      --r2 --ld-snp "$lead" \
      --ld-window 999999 --ld-window-kb 2000 --ld-window-r2 0 \
      --out "$ld_prefix" >/dev/null
  else
    log "[SKIP] exists: $ld_file"
  fi

  [[ -f "$assoc" ]] || die "assoc not produced: $assoc"
  [[ -f "$ld_file" ]] || die "ld not produced: $ld_file"

  printf "%s\t%s\n" "$assoc" "$ld_file"
}

# ===== lead SNPs =====
ERAP2_LEAD="rs2910686"
LNPEP_LEAD="rs248215"
ERAP1_SIG1="rs30379"
ERAP1_SIG2="rs27039"
ERAP1_SIG3="rs1065407"

# ===== regions =====
bp_er2="$(get_bp "$ERAP2_LEAD")"; read -r er2_from er2_to <<<"$(region_bounds "$bp_er2")"
bp_ln="$(get_bp "$LNPEP_LEAD")"; read -r ln_from ln_to <<<"$(region_bounds "$bp_ln")"
bp_s1="$(get_bp "$ERAP1_SIG1")"; read -r s1_from s1_to <<<"$(region_bounds "$bp_s1")"
bp_s2="$(get_bp "$ERAP1_SIG2")"; read -r s2_from s2_to <<<"$(region_bounds "$bp_s2")"
bp_s3="$(get_bp "$ERAP1_SIG3")"; read -r s3_from s3_to <<<"$(region_bounds "$bp_s3")"

COND_S1S2="$OUTDIR/ERAP1_sig3_condition_list.txt"
printf "%s\n%s\n" "$ERAP1_SIG1" "$ERAP1_SIG2" > "$COND_S1S2"

# ===== run windows (capture ONLY paths) =====
IFS=$'\t' read -r ER2_ASSOC ER2_LD < <(plink_window "ERAP2" "ERAP2_sig1_pm${WINDOW_BP}" "$ERAP2_LEAD" "$er2_from" "$er2_to" "none")
IFS=$'\t' read -r LN_ASSOC  LN_LD  < <(plink_window "LNPEP" "LNPEP_sig1_pm${WINDOW_BP}" "$LNPEP_LEAD" "$ln_from" "$ln_to" "none")
IFS=$'\t' read -r S1_ASSOC  S1_LD  < <(plink_window "ERAP1" "ERAP1_sig1_pm${WINDOW_BP}" "$ERAP1_SIG1" "$s1_from" "$s1_to" "none")
IFS=$'\t' read -r S2_ASSOC  S2_LD  < <(plink_window "ERAP1" "ERAP1_sig2_pm${WINDOW_BP}" "$ERAP1_SIG2" "$s2_from" "$s2_to" "snp"  "$ERAP1_SIG1")
IFS=$'\t' read -r S3_ASSOC  S3_LD  < <(plink_window "ERAP1" "ERAP1_sig3_pm${WINDOW_BP}" "$ERAP1_SIG3" "$s3_from" "$s3_to" "list" "$COND_S1S2")

# ===== summary =====
SUMMARY="$OUTDIR/signals_summary.tsv"
{
  echo -e "gene\tsignal\tlead\tbp\tfrom_bp\tto_bp\tassoc\tld"
  echo -e "ERAP2\tsig1\t$ERAP2_LEAD\t$bp_er2\t$er2_from\t$er2_to\t$ER2_ASSOC\t$ER2_LD"
  echo -e "LNPEP\tsig1\t$LNPEP_LEAD\t$bp_ln\t$ln_from\t$ln_to\t$LN_ASSOC\t$LN_LD"
  echo -e "ERAP1\tsig1\t$ERAP1_SIG1\t$bp_s1\t$s1_from\t$s1_to\t$S1_ASSOC\t$S1_LD"
  echo -e "ERAP1\tsig2\t$ERAP1_SIG2\t$bp_s2\t$s2_from\t$s2_to\t$S2_ASSOC\t$S2_LD"
  echo -e "ERAP1\tsig3\t$ERAP1_SIG3\t$bp_s3\t$s3_from\t$s3_to\t$S3_ASSOC\t$S3_LD"
} > "$SUMMARY"
log "[OK] wrote $SUMMARY"

# ===== plots =====
# single panels
"$PYTHON" "$CODE_LOCUS" \
  --assoc "$ER2_ASSOC" --ld "$ER2_LD" --lead "$ERAP2_LEAD" \
  --out-png "$FIGDIR/Fig_locus_ERAP2_sig1_pm${WINDOW_BP}.png" \
  --title "Regional cis-eQTL association at 5q15 for ERAP2 (±500 kb)" \
  --gw-threshold 5e-8 --test ADD

"$PYTHON" "$CODE_LOCUS" \
  --assoc "$LN_ASSOC" --ld "$LN_LD" --lead "$LNPEP_LEAD" \
  --out-png "$FIGDIR/Fig_locus_LNPEP_sig1_pm${WINDOW_BP}.png" \
  --title "Regional cis-eQTL association at 5q15 for LNPEP (±500 kb)" \
  --gw-threshold 5e-8 --test ADD

"$PYTHON" "$CODE_LOCUS" \
  --assoc "$S1_ASSOC" --ld "$S1_LD" --lead "$ERAP1_SIG1" \
  --out-png "$FIGDIR/Fig_locus_ERAP1_sig1_pm${WINDOW_BP}.png" \
  --title "Regional cis-eQTL association at 5q15 for ERAP1 signal1 (±500 kb)" \
  --gw-threshold 5e-8 --test ADD

"$PYTHON" "$CODE_LOCUS" \
  --assoc "$S2_ASSOC" --ld "$S2_LD" --lead "$ERAP1_SIG2" \
  --out-png "$FIGDIR/Fig_locus_ERAP1_sig2_pm${WINDOW_BP}.png" \
  --title "Regional cis-eQTL association at 5q15 for ERAP1 signal2 (±500 kb)" \
  --gw-threshold 5e-8 --test ADD

"$PYTHON" "$CODE_LOCUS" \
  --assoc "$S3_ASSOC" --ld "$S3_LD" --lead "$ERAP1_SIG3" \
  --out-png "$FIGDIR/Fig_locus_ERAP1_sig3_pm${WINDOW_BP}.png" \
  --title "Regional cis-eQTL association at 5q15 for ERAP1 signal3 (±500 kb)" \
  --gw-threshold 5e-8 --test ADD

# 3-panel: ERAP2 + ERAP1 sig1 + LNPEP
"$PYTHON" "$CODE_LOCUS" \
  --assoc "$ER2_ASSOC,$S1_ASSOC,$LN_ASSOC" \
  --ld    "$ER2_LD,$S1_LD,$LN_LD" \
  --lead  "$ERAP2_LEAD,$ERAP1_SIG1,$LNPEP_LEAD" \
  --titles "ERAP2,ERAP1 signal1,LNPEP" \
  --suptitle "Comparative cis-eQTL profiles at 5q15 for ERAP2, ERAP1 signal1 and LNPEP (±500 kb)" \
  --out-png "$FIGDIR/Fig_locus_3panel_ERAP2sig1_ERAP1sig1_LNPEPsig1_pm${WINDOW_BP}.png" \
  --gw-threshold 5e-8 --test ADD

# 3-panel: ERAP1 sig1/sig2/sig3
"$PYTHON" "$CODE_LOCUS" \
  --assoc "$S1_ASSOC,$S2_ASSOC,$S3_ASSOC" \
  --ld    "$S1_LD,$S2_LD,$S3_LD" \
  --lead  "$ERAP1_SIG1,$ERAP1_SIG2,$ERAP1_SIG3" \
  --titles "signal1,signal2 (cond S1),signal3 (cond S1+S2)" \
  --suptitle "Independent ERAP1 cis-eQTL signals at 5q15 revealed by stepwise conditional analysis (±500 kb)" \
  --out-png "$FIGDIR/Fig_locus_3panel_ERAP1_sig123_pm${WINDOW_BP}.png" \
  --gw-threshold 5e-8 --test ADD

log "[OK] 03 done."
