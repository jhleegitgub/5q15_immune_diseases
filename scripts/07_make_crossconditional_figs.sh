#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_config.sh"

SIGNALS_TSV="${SIGNALS_TSV:-$RESULTDIR/03_signal_check_5q15/signals_summary.tsv}"
CCDIR="${CCDIR:-$RESULTDIR/06_cross_conditional}"

PLOT_R="${PLOT_R:-$CODEDIR/plot_crossconditional_27panel.R}"

OUTLD="${OUTLD:-$CCDIR/ld}"
mkdir -p "$OUTLD" "$FIGDIR" "$TABLEDIR"

[[ -x "$PLINK" ]] || { echo "[ERR] PLINK not executable: $PLINK" >&2; exit 1; }
[[ -f "${BFILE}.bim" ]] || { echo "[ERR] missing BIM: ${BFILE}.bim" >&2; exit 1; }
[[ -f "$SIGNALS_TSV" ]] || { echo "[ERR] missing SIGNALS_TSV: $SIGNALS_TSV" >&2; exit 1; }
[[ -f "$CCDIR/runs.tsv" ]] || { echo "[ERR] missing 06 output: $CCDIR/runs.tsv" >&2; exit 1; }
[[ -f "$CCDIR/attenuation_summary.tsv" ]] || { echo "[ERR] missing 06 output: $CCDIR/attenuation_summary.tsv" >&2; exit 1; }
[[ -f "$PLOT_R" ]] || { echo "[ERR] missing: $PLOT_R" >&2; exit 1; }

ERAP2_LEAD1="$(awk -F'\t' 'NR>1 && $1=="ERAP2" && $2=="signal1"{print $3; exit}' "$SIGNALS_TSV")"
ERAP1_LEAD1="$(awk -F'\t' 'NR>1 && $1=="ERAP1" && $2=="signal1"{print $3; exit}' "$SIGNALS_TSV")"
LNPEP_LEAD1="$(awk -F'\t' 'NR>1 && $1=="LNPEP" && $2=="signal1"{print $3; exit}' "$SIGNALS_TSV")"

get_chr_bp() { awk -v s="$1" '$2==s{print $1"\t"$4; found=1; exit} END{if(!found) exit 2}' "${BFILE}.bim"; }
read -r CHR0 BP0 < <(get_chr_bp "$ERAP2_LEAD1" | awk '{print $1, $2}')
FROM=$((BP0 - WINKB*1000)); TO=$((BP0 + WINKB*1000)); [[ $FROM -lt 1 ]] && FROM=1

make_ld() {
  local ref="$1"
  local out="$OUTLD/${ref}_pm${WINKB}kb"
  "$PLINK" \
    --bfile "$BFILE" \
    --chr "$CHR0" --from-bp "$FROM" --to-bp "$TO" \
    --ld-snp "$ref" \
    --r2 gz --ld-window 99999 --ld-window-kb "$WINKB" --ld-window-r2 0 \
    --out "$out" >/dev/null
  echo "${out}.ld.gz"
}

LD_ERAP2="$(make_ld "$ERAP2_LEAD1")"
LD_ERAP1="$(make_ld "$ERAP1_LEAD1")"
LD_LNPEP="$(make_ld "$LNPEP_LEAD1")"

OUTPDF="$FIGDIR/FigS4_crossconditional_27panel_pm${WINKB}kb.pdf"
OUTPNG="$FIGDIR/FigS4_crossconditional_27panel_pm${WINKB}kb.png"
OUTTSV="$TABLEDIR/TableS2_cross_attenuation_signal1_pm${WINKB}kb.tsv"

Rscript "$PLOT_R" \
  "$CCDIR/runs.tsv" \
  "$SIGNALS_TSV" \
  "$LD_ERAP2" "$LD_ERAP1" "$LD_LNPEP" \
  "$WINKB" \
  "$OUTPDF" "$OUTPNG" "$OUTTSV"

echo "[OK] $OUTPDF"
echo "[OK] $OUTPNG"
echo "[OK] $OUTTSV"
