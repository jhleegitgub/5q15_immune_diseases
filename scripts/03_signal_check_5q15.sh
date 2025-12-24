#!/usr/bin/env bash
set -euo pipefail

export MPLBACKEND=Agg

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CODEDIR="$ROOT/code"
RESULTDIR="$ROOT/result/03_signal_check_5q15"
FIGDIR="$ROOT/fig"
LOGDIR="$RESULTDIR/logs"
mkdir -p "$RESULTDIR" "$FIGDIR" "$LOGDIR"

PY="${PY:-python3}"
PLOT_PY="$CODEDIR/locuszoom_manhattan.py"

BFILE="${BFILE:-$ROOT/../GenotypeData/GW.E-GEUV-3.EUR.MAF005.HWE1e-06}"
PLINK="${PLINK:-$HOME/Software/Plink/plink}"

PHENO="${PHENO:-}"
if [[ -z "${PHENO}" ]]; then
  for f in \
    "$ROOT/PhenotypeData/chr5_GD462.signalGeneQuantRPKM_plink.txt" \
    "$ROOT/../PhenotypeData/chr5_GD462.signalGeneQuantRPKM_plink.txt" \
    "$ROOT/PhenotypeData/chr5_GD462.signalGeneQuantRPKM_plink_input.txt" \
    "$ROOT/../PhenotypeData/chr5_GD462.signalGeneQuantRPKM_plink_input.txt"
  do
    [[ -f "$f" ]] && PHENO="$f" && break
  done
fi

COVAR="${COVAR:-}"
if [[ -z "${COVAR}" ]]; then
  [[ -f "$ROOT/result/01_pca/covar_pca10.tsv" ]] && COVAR="$ROOT/result/01_pca/covar_pca10.tsv"
  [[ -z "${COVAR}" && -f "$ROOT/covar_pca10.tsv" ]] && COVAR="$ROOT/covar_pca10.tsv"
fi
COVAR_NAMES="${COVAR_NAMES:-C1 C2 C3 C4 C5 C6 C7 C8 C9 C10}"

WIN="${WIN:-500000}"
GW_P="${GW_P:-5e-8}"

ERAP2_LEAD="${ERAP2_LEAD:-rs2910686}"
LNPEP_LEAD="${LNPEP_LEAD:-rs248215}"
ERAP1_S1="${ERAP1_S1:-rs30379}"
ERAP1_S2="${ERAP1_S2:-rs27039}"
ERAP1_S3="${ERAP1_S3:-rs1065407}"

[[ -x "$PLINK" ]] || { echo "[ERR] PLINK not executable: $PLINK" >&2; exit 1; }
[[ -f "${BFILE}.bim" ]] || { echo "[ERR] BFILE not found: ${BFILE}.bim" >&2; exit 1; }
[[ -n "${PHENO}" && -f "${PHENO}" ]] || { echo "[ERR] PHENO not found: $PHENO" >&2; exit 1; }
[[ -f "${COVAR}" ]] || { echo "[ERR] COVAR not found: $COVAR" >&2; exit 1; }
[[ -f "${PLOT_PY}" ]] || { echo "[ERR] plotter missing: $PLOT_PY" >&2; exit 1; }

get_chr_bp() {
  local snp="$1"
  awk -v s="$snp" '$2==s {print $1"\t"$4; exit}' "${BFILE}.bim"
}

calc_window() {
  local snp="$1"
  local chr bp
  read -r chr bp < <(get_chr_bp "$snp")
  [[ -n "${chr:-}" && -n "${bp:-}" ]] || { echo "[ERR] SNP not in BIM: $snp" >&2; exit 1; }
  local from=$((bp - WIN)); ((from < 1)) && from=1
  local to=$((bp + WIN))
  echo -e "${chr}\t${from}\t${to}"
}

run_assoc_window() {
  local gene="$1" center="$2" outprefix="$3" condlist="${4:-}"
  local chr from to
  read -r chr from to < <(calc_window "$center")
  local logfile="$LOGDIR/$(basename "$outprefix").plink.log"

  if [[ -n "${condlist}" ]]; then
    "$PLINK" --bfile "$BFILE" \
      --chr "$chr" --from-bp "$from" --to-bp "$to" \
      --pheno "$PHENO" --pheno-name "$gene" \
      --covar "$COVAR" --covar-name $COVAR_NAMES \
      --condition-list "$condlist" \
      --linear hide-covar --allow-no-sex \
      --out "$outprefix" \
      >"$logfile" 2>&1
  else
    "$PLINK" --bfile "$BFILE" \
      --chr "$chr" --from-bp "$from" --to-bp "$to" \
      --pheno "$PHENO" --pheno-name "$gene" \
      --covar "$COVAR" --covar-name $COVAR_NAMES \
      --linear hide-covar --allow-no-sex \
      --out "$outprefix" \
      >"$logfile" 2>&1
  fi
}

run_ld() {
  local lead="$1" center="$2" outprefix="$3"
  local chr from to
  read -r chr from to < <(calc_window "$center")
  local kb=$(( (WIN + 999) / 1000 + 10 ))

  "$PLINK" --bfile "$BFILE" \
    --chr "$chr" --from-bp "$from" --to-bp "$to" \
    --r2 \
    --ld-snp "$lead" \
    --ld-window 999999 \
    --ld-window-kb "$kb" \
    --ld-window-r2 0 \
    --out "$outprefix" \
    >"$LOGDIR/$(basename "$outprefix").ld.log" 2>&1
}

plot_single() {
  local assoc="$1" ld="$2" lead="$3" outpng="$4" title="$5"
  "$PY" "$PLOT_PY" \
    --assoc "$assoc" \
    --ld "$ld" \
    --lead "$lead" \
    --out-png "$outpng" \
    --title "$title" \
    --gw-threshold "$GW_P" \
    --point-size 18
}

echo "[RUN] ERAP2 sig1 baseline"
run_assoc_window ERAP2 "$ERAP2_LEAD" "$RESULTDIR/ERAP2_sig1"
run_ld "$ERAP2_LEAD" "$ERAP2_LEAD" "$RESULTDIR/ERAP2_sig1_ld"

echo "[RUN] LNPEP sig1 baseline"
run_assoc_window LNPEP "$LNPEP_LEAD" "$RESULTDIR/LNPEP_sig1"
run_ld "$LNPEP_LEAD" "$LNPEP_LEAD" "$RESULTDIR/LNPEP_sig1_ld"

echo "[RUN] ERAP1 sig1 baseline"
run_assoc_window ERAP1 "$ERAP1_S1" "$RESULTDIR/ERAP1_sig1"
run_ld "$ERAP1_S1" "$ERAP1_S1" "$RESULTDIR/ERAP1_sig1_ld"

echo "[RUN] ERAP1 sig2 (cond on sig1)"
COND_S1="$RESULTDIR/cond_ERAP1_S1.txt"
printf "%s\n" "$ERAP1_S1" > "$COND_S1"
run_assoc_window ERAP1 "$ERAP1_S1" "$RESULTDIR/ERAP1_sig2_condS1" "$COND_S1"
run_ld "$ERAP1_S2" "$ERAP1_S1" "$RESULTDIR/ERAP1_sig2_ld"

echo "[RUN] ERAP1 sig3 (cond on sig1+sig2)"
COND_S1S2="$RESULTDIR/cond_ERAP1_S1S2.txt"
printf "%s\n%s\n" "$ERAP1_S1" "$ERAP1_S2" > "$COND_S1S2"
run_assoc_window ERAP1 "$ERAP1_S1" "$RESULTDIR/ERAP1_sig3_condS1S2" "$COND_S1S2"
run_ld "$ERAP1_S3" "$ERAP1_S1" "$RESULTDIR/ERAP1_sig3_ld"

echo "[RUN] locus plots (single)"
plot_single "$RESULTDIR/ERAP2_sig1.assoc.linear" "$RESULTDIR/ERAP2_sig1_ld.ld" "$ERAP2_LEAD" \
  "$FIGDIR/Fig_locus_ERAP2_signal1_pm500kb.png" "ERAP2 locus (±500 kb)"

plot_single "$RESULTDIR/LNPEP_sig1.assoc.linear" "$RESULTDIR/LNPEP_sig1_ld.ld" "$LNPEP_LEAD" \
  "$FIGDIR/Fig_locus_LNPEP_signal1_pm500kb.png" "LNPEP locus (±500 kb)"

plot_single "$RESULTDIR/ERAP1_sig1.assoc.linear" "$RESULTDIR/ERAP1_sig1_ld.ld" "$ERAP1_S1" \
  "$FIGDIR/Fig_locus_ERAP1_signal1_pm500kb.png" "ERAP1 signal1 (±500 kb)"

plot_single "$RESULTDIR/ERAP1_sig2_condS1.assoc.linear" "$RESULTDIR/ERAP1_sig2_ld.ld" "$ERAP1_S2" \
  "$FIGDIR/Fig_locus_ERAP1_signal2_pm500kb.png" "ERAP1 signal2 (±500 kb; conditioned on signal1)"

plot_single "$RESULTDIR/ERAP1_sig3_condS1S2.assoc.linear" "$RESULTDIR/ERAP1_sig3_ld.ld" "$ERAP1_S3" \
  "$FIGDIR/Fig_locus_ERAP1_signal3_pm500kb.png" "ERAP1 signal3 (±500 kb; conditioned on signal1+2)"

echo "[RUN] combined (ERAP2 sig1 + ERAP1 sig1 + LNPEP sig1)"
"$PY" "$PLOT_PY" \
  --assoc "$RESULTDIR/ERAP2_sig1.assoc.linear,$RESULTDIR/ERAP1_sig1.assoc.linear,$RESULTDIR/LNPEP_sig1.assoc.linear" \
  --ld    "$RESULTDIR/ERAP2_sig1_ld.ld,$RESULTDIR/ERAP1_sig1_ld.ld,$RESULTDIR/LNPEP_sig1_ld.ld" \
  --lead  "$ERAP2_LEAD,$ERAP1_S1,$LNPEP_LEAD" \
  --titles "ERAP2,ERAP1 signal1,LNPEP" \
  --out-png "$FIGDIR/Fig_locus_3panel_ERAP2sig1_ERAP1sig1_LNPEPsig1_pm500kb.png" \
  --gw-threshold "$GW_P" \
  --point-size 16

echo "[RUN] combined (ERAP1 sig1+sig2+sig3)"
"$PY" "$PLOT_PY" \
  --assoc "$RESULTDIR/ERAP1_sig1.assoc.linear,$RESULTDIR/ERAP1_sig2_condS1.assoc.linear,$RESULTDIR/ERAP1_sig3_condS1S2.assoc.linear" \
  --ld    "$RESULTDIR/ERAP1_sig1_ld.ld,$RESULTDIR/ERAP1_sig2_ld.ld,$RESULTDIR/ERAP1_sig3_ld.ld" \
  --lead  "$ERAP1_S1,$ERAP1_S2,$ERAP1_S3" \
  --titles "signal1,signal2 (cond S1),signal3 (cond S1+S2)" \
  --out-png "$FIGDIR/Fig_locus_3panel_ERAP1_sig123_pm500kb.png" \
  --gw-threshold "$GW_P" \
  --point-size 16

echo "[OK] done"
echo " - result: $RESULTDIR"
echo " - figs:   $FIGDIR"
