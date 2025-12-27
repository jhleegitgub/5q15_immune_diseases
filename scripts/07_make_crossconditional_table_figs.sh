#!/usr/bin/env bash
# scripts/07_make_crossconditional_table_figs.sh
# - 5q15 window cross-conditional associations + attenuation tables + 3x3/2x2 plots
set -euo pipefail

ts(){ date +"[%H:%M:%S]"; }
log(){ echo "$(ts) $*"; }
die(){ echo "[ERR] $*" 1>&2; exit 1; }

# ----------------------------
# inputs (env override OK)
# ----------------------------
BFILE="${BFILE:-../GenotypeData/GW.E-GEUV-3.EUR.MAF005.HWE1e-06}"
PHENO5="${PHENO5:-../PhenotypeData/chr5_GD462.signalGeneQuantRPKM_plink.txt}"
COVAR="${COVAR:-./result/01_pca/covar_pca10.tsv}"
PLINK="${PLINK:-$HOME/Software/Plink/plink}"
COVAR_NAMES="${COVAR_NAMES:-C1 C2 C3 C4 C5 C6 C7 C8 C9 C10}"

# 5q15 window definition (center = ERAP1 sig3 by default)
CENTER_SNP="${CENTER_SNP:-rs1065407}"
WIN="${WIN:-500000}"

# lead SNPs
ERAP2_LEAD="${ERAP2_LEAD:-rs2910686}"
LNPEP_LEAD="${LNPEP_LEAD:-rs248215}"

ERAP1_S1="${ERAP1_S1:-rs30379}"
ERAP1_S2="${ERAP1_S2:-rs27039}"
ERAP1_S3="${ERAP1_S3:-rs1065407}"

# project paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CODEDIR="${CODEDIR:-$ROOT/code}"
RESULTDIR="${RESULTDIR:-$ROOT/result}"
FIGDIR="${FIGDIR:-$ROOT/fig}"
TABLEDIR="${TABLEDIR:-$ROOT/table}"

OUTDIR="${OUTDIR:-$RESULTDIR/07_crossconditional}"
TMPDIR="$OUTDIR/tmp"
ASSOCDIR="$OUTDIR/assoc"
LOGDIR="$OUTDIR/logs"
mkdir -p "$TMPDIR" "$ASSOCDIR" "$LOGDIR" "$FIGDIR" "$TABLEDIR"

# code deps
MAKE_COVAR_PLUS_EXPR_PY="${MAKE_COVAR_PLUS_EXPR_PY:-$CODEDIR/make_covar_plus_expr.py}"
PLOT_R="${PLOT_R:-$CODEDIR/locus_grid_multi.R}"
SUMMARISE_PY="${SUMMARISE_PY:-$CODEDIR/summarize_cross_conditional_v2.py}"

# ----------------------------
# checks
# ----------------------------
[[ -x "$PLINK" ]] || die "PLINK not executable: $PLINK"
[[ -f "${BFILE}.bim" ]] || die "missing: ${BFILE}.bim"
[[ -f "$PHENO5" ]] || die "missing: $PHENO5"
[[ -f "$COVAR" ]] || die "missing: $COVAR"
[[ -f "$MAKE_COVAR_PLUS_EXPR_PY" ]] || die "missing: $MAKE_COVAR_PLUS_EXPR_PY"
[[ -f "$PLOT_R" ]] || die "missing: $PLOT_R"
[[ -f "$SUMMARISE_PY" ]] || die "missing: $SUMMARISE_PY"

# ----------------------------
# helpers
# ----------------------------
get_chr_bp() {
  local snp="$1"
  awk -v s="$snp" '$2==s{print $1"\t"$4; exit}' "${BFILE}.bim"
}

mk_window() {
  local snp="$1"
  local chr bp
  read -r chr bp < <(get_chr_bp "$snp" || true)
  [[ -n "${chr:-}" && -n "${bp:-}" ]] || die "SNP not found in BIM: $snp"
  local from=$((bp - WIN)); ((from<1)) && from=1
  local to=$((bp + WIN))
  echo -e "$chr\t$from\t$to"
}

get_best_snp() {
  # best SNP in an assoc.linear (ADD rows, min P)
  local assoc="$1"
  awk '
    BEGIN{best=""; bp=0; p=1}
    NR==1{next}
    $5=="ADD"{
      if($9+0 < p){ p=$9+0; best=$2; bp=$3 }
    }
    END{
      if(best==""){ print "NA\tNA" } else { print best"\t"bp }
    }
  ' "$assoc"
}

plink_window_baseline() {
  local gene="$1" center="$2" outprefix="$3"
  local chr from to
  read -r chr from to < <(mk_window "$center")
  "$PLINK" --bfile "$BFILE" \
    --chr "$chr" --from-bp "$from" --to-bp "$to" \
    --pheno "$PHENO5" --pheno-name "$gene" \
    --covar "$COVAR" --covar-name $COVAR_NAMES \
    --linear hide-covar --allow-no-sex \
    --out "$outprefix" >"$LOGDIR/${gene}_baseline.log" 2>&1
}

plink_window_cond_snp() {
  local gene="$1" center="$2" cond_snp="$3" outprefix="$4"
  local chr from to
  read -r chr from to < <(mk_window "$center")
  "$PLINK" --bfile "$BFILE" \
    --chr "$chr" --from-bp "$from" --to-bp "$to" \
    --pheno "$PHENO5" --pheno-name "$gene" \
    --covar "$COVAR" --covar-name $COVAR_NAMES \
    --condition "$cond_snp" \
    --linear hide-covar --allow-no-sex \
    --out "$outprefix" >"$LOGDIR/${gene}_cond_${cond_snp}.log" 2>&1
}

plink_window_cond_list() {
  local gene="$1" center="$2" cond_list="$3" tag="$4" outprefix="$5"
  local chr from to
  read -r chr from to < <(mk_window "$center")
  "$PLINK" --bfile "$BFILE" \
    --chr "$chr" --from-bp "$from" --to-bp "$to" \
    --pheno "$PHENO5" --pheno-name "$gene" \
    --covar "$COVAR" --covar-name $COVAR_NAMES \
    --condition-list "$cond_list" \
    --linear hide-covar --allow-no-sex \
    --out "$outprefix" >"$LOGDIR/${gene}_cond_${tag}.log" 2>&1
}

plink_window_cond_expr() {
  # IMPORTANT: expr_col must match make_covar_plus_expr.py output = expr_<GENE>
  local gene="$1" center="$2" covar_tsv="$3" expr_col="$4" outprefix="$5"
  local chr from to
  read -r chr from to < <(mk_window "$center")
  "$PLINK" --bfile "$BFILE" \
    --chr "$chr" --from-bp "$from" --to-bp "$to" \
    --pheno "$PHENO5" --pheno-name "$gene" \
    --covar "$covar_tsv" --covar-name $COVAR_NAMES "$expr_col" \
    --linear hide-covar --allow-no-sex \
    --out "$outprefix" >"$LOGDIR/${gene}_condexpr_${expr_col}.log" 2>&1
}

# ----------------------------
# 0) snplist + LD matrix for plotting
# ----------------------------
log "center SNP: $CENTER_SNP"
read -r CHR FROM TO < <(mk_window "$CENTER_SNP")

SNPLIST="$TMPDIR/snplist_${CENTER_SNP}_pm${WIN}.snplist"
LDGZ="$TMPDIR/ld_${CENTER_SNP}_pm${WIN}.ld.gz"

log "make snplist: $SNPLIST"
awk -v chr="$CHR" -v from="$FROM" -v to="$TO" '$1==chr && $4>=from && $4<=to{print $2}' "${BFILE}.bim" > "$SNPLIST"
[[ -s "$SNPLIST" ]] || die "SNPLIST empty: $SNPLIST"

if [[ -s "$LDGZ" ]]; then
  log "LD matrix exists: $LDGZ"
else
  log "make LD matrix: $LDGZ"
  PREF_LD="$TMPDIR/_plink_ld_${CENTER_SNP}_pm${WIN}"
  "$PLINK" --bfile "$BFILE" \
    --chr "$CHR" --from-bp "$FROM" --to-bp "$TO" \
    --extract "$SNPLIST" --keep-allele-order \
    --r square gz --out "$PREF_LD" >"$LOGDIR/ld_matrix.log" 2>&1
  mv -f "${PREF_LD}.ld.gz" "$LDGZ"
fi

# ----------------------------
# 1) covar + expression (expr_<GENE> 컬럼 생성)
# ----------------------------
COVAR_PLUS_ERAP2="$TMPDIR/covar_plus_ERAP2.tsv"
COVAR_PLUS_ERAP1="$TMPDIR/covar_plus_ERAP1.tsv"
COVAR_PLUS_LNPEP="$TMPDIR/covar_plus_LNPEP.tsv"
COVAR_PLUS_CSF2="$TMPDIR/covar_plus_CSF2.tsv"

log "covar+expr: ERAP2 -> $COVAR_PLUS_ERAP2"
python3 "$MAKE_COVAR_PLUS_EXPR_PY" "$COVAR" "$PHENO5" ERAP2 "$COVAR_PLUS_ERAP2"
log "covar+expr: ERAP1 -> $COVAR_PLUS_ERAP1"
python3 "$MAKE_COVAR_PLUS_EXPR_PY" "$COVAR" "$PHENO5" ERAP1 "$COVAR_PLUS_ERAP1"
log "covar+expr: LNPEP -> $COVAR_PLUS_LNPEP"
python3 "$MAKE_COVAR_PLUS_EXPR_PY" "$COVAR" "$PHENO5" LNPEP "$COVAR_PLUS_LNPEP"
log "covar+expr: CSF2 -> $COVAR_PLUS_CSF2"
python3 "$MAKE_COVAR_PLUS_EXPR_PY" "$COVAR" "$PHENO5" CSF2 "$COVAR_PLUS_CSF2"

# expr column names (make_covar_plus_expr.py output)
EXPR_ERAP2="expr_ERAP2"
EXPR_ERAP1="expr_ERAP1"
EXPR_LNPEP="expr_LNPEP"
EXPR_CSF2="expr_CSF2"

# ----------------------------
# 2) condition lists
# ----------------------------
COND_ERAP1_SIG123="$TMPDIR/cond_ERAP1_sig123.txt"
printf "%s\n%s\n%s\n" "$ERAP1_S1" "$ERAP1_S2" "$ERAP1_S3" > "$COND_ERAP1_SIG123"

# ----------------------------
# 3) baseline (5q15 window)
# ----------------------------
for g in ERAP2 ERAP1 LNPEP CSF2; do
  if [[ ! -s "$ASSOCDIR/${g}_base.assoc.linear" ]]; then
    log "baseline: $g"
    plink_window_baseline "$g" "$CENTER_SNP" "$ASSOCDIR/${g}_base"
  else
    log "baseline exists: $g"
  fi
done

# CSF2 lead inside window (for ERAP1-CSF2 pair)
read -r CSF2_LEAD_WIN CSF2_LEAD_BP < <(get_best_snp "$ASSOCDIR/CSF2_base.assoc.linear")
[[ "$CSF2_LEAD_WIN" != "NA" ]] || die "could not pick CSF2 lead in window"

log "CSF2 local lead in window: $CSF2_LEAD_WIN"

# ----------------------------
# 4) SNP-conditionals
# ----------------------------
# outcome ERAP2
log "ERAP2 | cond ERAP1(sig1+2+3)"
plink_window_cond_list ERAP2 "$CENTER_SNP" "$COND_ERAP1_SIG123" "ERAP1sig123" "$ASSOCDIR/ERAP2_cond_ERAP1sig123"
log "ERAP2 | cond ERAP2($ERAP2_LEAD)"
plink_window_cond_snp  ERAP2 "$CENTER_SNP" "$ERAP2_LEAD" "$ASSOCDIR/ERAP2_cond_${ERAP2_LEAD}"
log "ERAP2 | cond LNPEP($LNPEP_LEAD)"
plink_window_cond_snp  ERAP2 "$CENTER_SNP" "$LNPEP_LEAD" "$ASSOCDIR/ERAP2_cond_${LNPEP_LEAD}"

# outcome ERAP1
log "ERAP1 | cond ERAP1(sig1+2+3)"
plink_window_cond_list ERAP1 "$CENTER_SNP" "$COND_ERAP1_SIG123" "ERAP1sig123" "$ASSOCDIR/ERAP1_cond_ERAP1sig123"
log "ERAP1 | cond ERAP2($ERAP2_LEAD)"
plink_window_cond_snp  ERAP1 "$CENTER_SNP" "$ERAP2_LEAD" "$ASSOCDIR/ERAP1_cond_${ERAP2_LEAD}"
log "ERAP1 | cond LNPEP($LNPEP_LEAD)"
plink_window_cond_snp  ERAP1 "$CENTER_SNP" "$LNPEP_LEAD" "$ASSOCDIR/ERAP1_cond_${LNPEP_LEAD}"
log "ERAP1 | cond CSF2($CSF2_LEAD_WIN)"
plink_window_cond_snp  ERAP1 "$CENTER_SNP" "$CSF2_LEAD_WIN" "$ASSOCDIR/ERAP1_cond_${CSF2_LEAD_WIN}"

# outcome LNPEP
log "LNPEP | cond ERAP1(sig1+2+3)"
plink_window_cond_list LNPEP "$CENTER_SNP" "$COND_ERAP1_SIG123" "ERAP1sig123" "$ASSOCDIR/LNPEP_cond_ERAP1sig123"
log "LNPEP | cond ERAP2($ERAP2_LEAD)"
plink_window_cond_snp  LNPEP "$CENTER_SNP" "$ERAP2_LEAD" "$ASSOCDIR/LNPEP_cond_${ERAP2_LEAD}"
log "LNPEP | cond LNPEP($LNPEP_LEAD)"
plink_window_cond_snp  LNPEP "$CENTER_SNP" "$LNPEP_LEAD" "$ASSOCDIR/LNPEP_cond_${LNPEP_LEAD}"

# outcome CSF2
log "CSF2 | cond ERAP1(sig1+2+3)"
plink_window_cond_list CSF2 "$CENTER_SNP" "$COND_ERAP1_SIG123" "ERAP1sig123" "$ASSOCDIR/CSF2_cond_ERAP1sig123"
log "CSF2 | cond CSF2($CSF2_LEAD_WIN)"
plink_window_cond_snp  CSF2 "$CENTER_SNP" "$CSF2_LEAD_WIN" "$ASSOCDIR/CSF2_cond_${CSF2_LEAD_WIN}"

# ----------------------------
# 5) expression-conditionals (FIXED: expr_<GENE>)
# ----------------------------
log "ERAP2 | cond on ERAP1 expression"
plink_window_cond_expr ERAP2 "$CENTER_SNP" "$COVAR_PLUS_ERAP1" "$EXPR_ERAP1" "$ASSOCDIR/ERAP2_cond_on_ERAP1expr"
log "ERAP2 | cond on LNPEP expression"
plink_window_cond_expr ERAP2 "$CENTER_SNP" "$COVAR_PLUS_LNPEP" "$EXPR_LNPEP" "$ASSOCDIR/ERAP2_cond_on_LNPEPexpr"

log "ERAP1 | cond on ERAP2 expression"
plink_window_cond_expr ERAP1 "$CENTER_SNP" "$COVAR_PLUS_ERAP2" "$EXPR_ERAP2" "$ASSOCDIR/ERAP1_cond_on_ERAP2expr"
log "ERAP1 | cond on LNPEP expression"
plink_window_cond_expr ERAP1 "$CENTER_SNP" "$COVAR_PLUS_LNPEP" "$EXPR_LNPEP" "$ASSOCDIR/ERAP1_cond_on_LNPEPexpr"
log "ERAP1 | cond on CSF2 expression"
plink_window_cond_expr ERAP1 "$CENTER_SNP" "$COVAR_PLUS_CSF2" "$EXPR_CSF2" "$ASSOCDIR/ERAP1_cond_on_CSF2expr"

log "LNPEP | cond on ERAP2 expression"
plink_window_cond_expr LNPEP "$CENTER_SNP" "$COVAR_PLUS_ERAP2" "$EXPR_ERAP2" "$ASSOCDIR/LNPEP_cond_on_ERAP2expr"
log "LNPEP | cond on ERAP1 expression"
plink_window_cond_expr LNPEP "$CENTER_SNP" "$COVAR_PLUS_ERAP1" "$EXPR_ERAP1" "$ASSOCDIR/LNPEP_cond_on_ERAP1expr"

log "CSF2 | cond on ERAP1 expression"
plink_window_cond_expr CSF2 "$CENTER_SNP" "$COVAR_PLUS_ERAP1" "$EXPR_ERAP1" "$ASSOCDIR/CSF2_cond_on_ERAP1expr"

# ----------------------------
# 6) attenuation tables + plots (3x3 + 2x2 x3)
# ----------------------------
SIG_MAIN="$TMPDIR/signals_main.tsv"
{
  echo -e "gene\tsignal_id\tlead"
  echo -e "ERAP2\tsignal1\t$ERAP2_LEAD"
  echo -e "ERAP1\tsignal1\t$ERAP1_S1"
  echo -e "LNPEP\tsignal1\t$LNPEP_LEAD"
} > "$SIG_MAIN"

RUNS_MAIN="$TMPDIR/runs_main.tsv"
{
  echo -e "outcome\trun_id\tcond_type\tcond_label\tbaseline_path\tassoc_path"
  for g in ERAP2 ERAP1 LNPEP; do
    echo -e "${g}\tbaseline\tbaseline\tbaseline\t$ASSOCDIR/${g}_base.assoc.linear\t$ASSOCDIR/${g}_base.assoc.linear"
  done

  echo -e "ERAP2\tcond_ERAP1sig123\tsnp\tERAP1_sig123\t$ASSOCDIR/ERAP2_base.assoc.linear\t$ASSOCDIR/ERAP2_cond_ERAP1sig123.assoc.linear"
  echo -e "ERAP2\tcond_${ERAP2_LEAD}\tsnp\tERAP2\t$ASSOCDIR/ERAP2_base.assoc.linear\t$ASSOCDIR/ERAP2_cond_${ERAP2_LEAD}.assoc.linear"
  echo -e "ERAP2\tcond_${LNPEP_LEAD}\tsnp\tLNPEP\t$ASSOCDIR/ERAP2_base.assoc.linear\t$ASSOCDIR/ERAP2_cond_${LNPEP_LEAD}.assoc.linear"
  echo -e "ERAP2\tcondexpr_ERAP1\texpr\tERAP1_expr\t$ASSOCDIR/ERAP2_base.assoc.linear\t$ASSOCDIR/ERAP2_cond_on_ERAP1expr.assoc.linear"
  echo -e "ERAP2\tcondexpr_LNPEP\texpr\tLNPEP_expr\t$ASSOCDIR/ERAP2_base.assoc.linear\t$ASSOCDIR/ERAP2_cond_on_LNPEPexpr.assoc.linear"

  echo -e "ERAP1\tcond_${ERAP2_LEAD}\tsnp\tERAP2\t$ASSOCDIR/ERAP1_base.assoc.linear\t$ASSOCDIR/ERAP1_cond_${ERAP2_LEAD}.assoc.linear"
  echo -e "ERAP1\tcond_${LNPEP_LEAD}\tsnp\tLNPEP\t$ASSOCDIR/ERAP1_base.assoc.linear\t$ASSOCDIR/ERAP1_cond_${LNPEP_LEAD}.assoc.linear"
  echo -e "ERAP1\tcondexpr_ERAP2\texpr\tERAP2_expr\t$ASSOCDIR/ERAP1_base.assoc.linear\t$ASSOCDIR/ERAP1_cond_on_ERAP2expr.assoc.linear"
  echo -e "ERAP1\tcondexpr_LNPEP\texpr\tLNPEP_expr\t$ASSOCDIR/ERAP1_base.assoc.linear\t$ASSOCDIR/ERAP1_cond_on_LNPEPexpr.assoc.linear"

  echo -e "LNPEP\tcond_ERAP1sig123\tsnp\tERAP1_sig123\t$ASSOCDIR/LNPEP_base.assoc.linear\t$ASSOCDIR/LNPEP_cond_ERAP1sig123.assoc.linear"
  echo -e "LNPEP\tcond_${ERAP2_LEAD}\tsnp\tERAP2\t$ASSOCDIR/LNPEP_base.assoc.linear\t$ASSOCDIR/LNPEP_cond_${ERAP2_LEAD}.assoc.linear"
  echo -e "LNPEP\tcond_${LNPEP_LEAD}\tsnp\tLNPEP\t$ASSOCDIR/LNPEP_base.assoc.linear\t$ASSOCDIR/LNPEP_cond_${LNPEP_LEAD}.assoc.linear"
  echo -e "LNPEP\tcondexpr_ERAP2\texpr\tERAP2_expr\t$ASSOCDIR/LNPEP_base.assoc.linear\t$ASSOCDIR/LNPEP_cond_on_ERAP2expr.assoc.linear"
  echo -e "LNPEP\tcondexpr_ERAP1\texpr\tERAP1_expr\t$ASSOCDIR/LNPEP_base.assoc.linear\t$ASSOCDIR/LNPEP_cond_on_ERAP1expr.assoc.linear"
} > "$RUNS_MAIN"

OUT_MAIN="$TABLEDIR/crossconditional_ERAP1_ERAP2_LNPEP_attenuation.tsv"
log "summarise 3x3 => $OUT_MAIN"
python3 "$SUMMARISE_PY" --signals "$SIG_MAIN" --runs "$RUNS_MAIN" --out "$OUT_MAIN"

log "plot 3x3: ERAP2/ERAP1/LNPEP"
( cd "$ASSOCDIR" && Rscript "$PLOT_R" "$LDGZ" "$SNPLIST" "$FIGDIR/Fig_crossconditional_3x3_ERAP1_ERAP2_LNPEP.png" \
    --genes ERAP2,ERAP1,LNPEP \
    --ref-snps ERAP2=$ERAP2_LEAD,ERAP1=$ERAP1_S1,LNPEP=$LNPEP_LEAD \
    --cond-suffix ERAP1=ERAP1sig123 \
    --cond-labels ERAP1=sig1+2+3 \
    --width 18 --height 26 --dpi 300 )

log "plot 2x2: ERAP1/ERAP2"
( cd "$ASSOCDIR" && Rscript "$PLOT_R" "$LDGZ" "$SNPLIST" "$FIGDIR/Fig_crossconditional_2x2_ERAP1_ERAP2.png" \
    --genes ERAP1,ERAP2 \
    --ref-snps ERAP1=$ERAP1_S1,ERAP2=$ERAP2_LEAD \
    --cond-suffix ERAP1=ERAP1sig123 \
    --cond-labels ERAP1=sig1+2+3 \
    --width 12 --height 18 --dpi 300 )

log "plot 2x2: ERAP2/LNPEP"
( cd "$ASSOCDIR" && Rscript "$PLOT_R" "$LDGZ" "$SNPLIST" "$FIGDIR/Fig_crossconditional_2x2_ERAP2_LNPEP.png" \
    --genes ERAP2,LNPEP \
    --ref-snps ERAP2=$ERAP2_LEAD,LNPEP=$LNPEP_LEAD \
    --width 12 --height 18 --dpi 300 )

# ERAP1-CSF2 table + plot
SIG_CSF2="$TMPDIR/signals_ERAP1_CSF2.tsv"
{
  echo -e "gene\tsignal_id\tlead"
  echo -e "ERAP1\tsignal1\t$ERAP1_S1"
  echo -e "CSF2\tsignal1\t$CSF2_LEAD_WIN"
} > "$SIG_CSF2"

RUNS_CSF2="$TMPDIR/runs_ERAP1_CSF2.tsv"
{
  echo -e "outcome\trun_id\tcond_type\tcond_label\tbaseline_path\tassoc_path"
  echo -e "ERAP1\tbaseline\tbaseline\tbaseline\t$ASSOCDIR/ERAP1_base.assoc.linear\t$ASSOCDIR/ERAP1_base.assoc.linear"
  echo -e "CSF2\tbaseline\tbaseline\tbaseline\t$ASSOCDIR/CSF2_base.assoc.linear\t$ASSOCDIR/CSF2_base.assoc.linear"

  echo -e "ERAP1\tcond_${CSF2_LEAD_WIN}\tsnp\tCSF2\t$ASSOCDIR/ERAP1_base.assoc.linear\t$ASSOCDIR/ERAP1_cond_${CSF2_LEAD_WIN}.assoc.linear"
  echo -e "ERAP1\tcondexpr_CSF2\texpr\tCSF2_expr\t$ASSOCDIR/ERAP1_base.assoc.linear\t$ASSOCDIR/ERAP1_cond_on_CSF2expr.assoc.linear"

  echo -e "CSF2\tcond_ERAP1sig123\tsnp\tERAP1_sig123\t$ASSOCDIR/CSF2_base.assoc.linear\t$ASSOCDIR/CSF2_cond_ERAP1sig123.assoc.linear"
  echo -e "CSF2\tcondexpr_ERAP1\texpr\tERAP1_expr\t$ASSOCDIR/CSF2_base.assoc.linear\t$ASSOCDIR/CSF2_cond_on_ERAP1expr.assoc.linear"
} > "$RUNS_CSF2"

OUT_CSF2="$TABLEDIR/crossconditional_ERAP1_CSF2_attenuation.tsv"
log "summarise ERAP1-CSF2 => $OUT_CSF2"
python3 "$SUMMARISE_PY" --signals "$SIG_CSF2" --runs "$RUNS_CSF2" --out "$OUT_CSF2"

log "plot 2x2: ERAP1/CSF2"
( cd "$ASSOCDIR" && Rscript "$PLOT_R" "$LDGZ" "$SNPLIST" "$FIGDIR/Fig_crossconditional_2x2_ERAP1_CSF2.png" \
    --genes ERAP1,CSF2 \
    --ref-snps ERAP1=$ERAP1_S1,CSF2=$CSF2_LEAD_WIN \
    --cond-suffix ERAP1=ERAP1sig123 \
    --cond-labels ERAP1=sig1+2+3 \
    --width 12 --height 18 --dpi 300 )

log "done"
log "figs: $FIGDIR/Fig_crossconditional_*"
log "tables: $TABLEDIR/crossconditional_*"
