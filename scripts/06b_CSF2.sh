#!/usr/bin/env bash
# scripts/06b_CSF2.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_config.sh"

export MPLBACKEND=Agg
export QT_QPA_PLATFORM=offscreen

# --------------------------
# inputs / params
# --------------------------
GENE="${GENE:-CSF2}"

CENTER_SNP="${CENTER_SNP:-rs1065407}"   # ERAP1 signal3
WIN_BP="${WIN_BP:-500000}"              # ±500kb

PHENO="${PHENO:-${PHENO5:-$PHENODIR/chr5_GD462.signalGeneQuantRPKM_plink.txt}}"
PLINK="${PLINK:-$HOME/Software/Plink/plink}"

GW_PY="${GW_PY:-$CODEDIR/manhattan_genomewide.py}"
LOCUS_PY="${LOCUS_PY:-$CODEDIR/locuszoom_manhattan.py}"

OUTDIR="${OUTDIR:-$RESULTDIR/06b_CSF2}"
mkdir -p "$OUTDIR" "$FIGDIR"

die(){ echo "[ERR] $*" 1>&2; exit 1; }
need(){ [[ -e "$1" ]] || die "missing: $1"; }

[[ -x "$PLINK" ]] || die "PLINK not executable: $PLINK"
need "${BFILE}.bed"; need "${BFILE}.bim"; need "${BFILE}.fam"
need "$PHENO"; need "$COVAR"
need "$GW_PY"; need "$LOCUS_PY"

get_chr_bp(){
  local snp="$1"
  awk -v s="$snp" '$2==s{print $1"\t"$4; exit}' "${BFILE}.bim"
}

find_local_lead(){
  local assoc="$1"
  awk '
    NR==1{
      for(i=1;i<=NF;i++){
        if($i=="SNP") snp=i;
        else if($i=="P") p=i;
        else if($i=="TEST") t=i;
      }
      next
    }
    ($t=="ADD" && $p!="NA"){
      pv=$p+0
      if(pv>0 && (min==0 || pv<min)){min=pv; lead=$snp}
    }
    END{ if(lead=="") print "NA"; else print lead }
  ' "$assoc"
}

# --------------------------
# 1) genome-wide CSF2 eQTL + manhattan
# --------------------------
GW_PREF="$OUTDIR/${GENE}_genomewide"
GW_ASSOC="${GW_PREF}.assoc.linear"

if [[ ! -s "$GW_ASSOC" ]]; then
  echo "[RUN] PLINK genome-wide: $GENE"
  "$PLINK" --bfile "$BFILE" \
    --pheno "$PHENO" --pheno-name "$GENE" \
    --covar "$COVAR" --covar-name $COVAR_NAMES \
    --linear hide-covar --allow-no-sex \
    --out "$GW_PREF" >/dev/null
else
  echo "[SKIP] exists: $GW_ASSOC"
fi

GW_FIG="$FIGDIR/Fig_genomewide_${GENE}_manhattan.png"
echo "[RUN] plot: $GW_FIG"
python3 "$GW_PY" \
  --assoc "$GW_ASSOC" \
  --out-png "$GW_FIG" \
  --title "Genome-wide eQTL scan for ${GENE} expression" \
  --gw-threshold 5e-8 \
  --test ADD \
  --max-points 2000000 \
  --xtick-step 1 >/dev/null

# --------------------------
# 2) 5q15 window around ERAP1 signal3: CSF2 baseline
# --------------------------
read -r CHR BP < <(get_chr_bp "$CENTER_SNP" || true)
[[ -n "${CHR:-}" && -n "${BP:-}" ]] || die "CENTER_SNP not in BIM: $CENTER_SNP"

FROM=$((BP - WIN_BP)); TO=$((BP + WIN_BP))
((FROM<1)) && FROM=1

REG_PREF="$OUTDIR/${GENE}_at_5q15_center_${CENTER_SNP}_pm${WIN_BP}"
REG_ASSOC="${REG_PREF}.assoc.linear"

echo "[RUN] PLINK window baseline: $GENE (center=$CENTER_SNP ±${WIN_BP}bp)"
"$PLINK" --bfile "$BFILE" \
  --chr "$CHR" --from-bp "$FROM" --to-bp "$TO" \
  --pheno "$PHENO" --pheno-name "$GENE" \
  --covar "$COVAR" --covar-name $COVAR_NAMES \
  --linear hide-covar --allow-no-sex \
  --out "$REG_PREF" >/dev/null

LEAD="$(find_local_lead "$REG_ASSOC")"
[[ "$LEAD" != "NA" ]] || die "failed to find local lead in: $REG_ASSOC"
echo "[OK] local lead in window: $LEAD"

# --------------------------
# 2.5) make LD file for locuszoom (required)
# --------------------------
WINKB=$((WIN_BP/1000))
LD_PREF="$OUTDIR/${GENE}_pm${WIN_BP}.ld_to_${LEAD}"
LD_FILE="${LD_PREF}.ld"

echo "[RUN] PLINK LD (lead=$LEAD, window=${WINKB}kb)"
"$PLINK" --bfile "$BFILE" \
  --chr "$CHR" --from-bp "$FROM" --to-bp "$TO" \
  --ld-snp "$LEAD" \
  --r2 --ld-window 99999 --ld-window-kb "$WINKB" --ld-window-r2 0 \
  --out "$LD_PREF" >/dev/null
need "$LD_FILE"

# plot baseline locus (needs --ld)
REG_FIG="$FIGDIR/Fig_locus_${GENE}_center_${CENTER_SNP}_pm500kb.png"
python3 "$LOCUS_PY" \
  --assoc "$REG_ASSOC" \
  --ld "$LD_FILE" \
  --lead "$LEAD" \
  --out-png "$REG_FIG" \
  --title "Regional cis-eQTL association at 5q15 for ${GENE} (±500 kb)" \
  --test ADD \
  --gw-threshold 5e-8 \
  --point-size 10 >/dev/null

# --------------------------
# 3) conditional on local lead, re-run in same window
# --------------------------
COND_LIST="$OUTDIR/${GENE}_local_lead.txt"
printf "%s\n" "$LEAD" > "$COND_LIST"

REGC_PREF="$OUTDIR/${GENE}_at_5q15_center_${CENTER_SNP}_pm${WIN_BP}_cond_${LEAD}"
REGC_ASSOC="${REGC_PREF}.assoc.linear"

echo "[RUN] PLINK window conditional: $GENE | cond($LEAD)"
"$PLINK" --bfile "$BFILE" \
  --chr "$CHR" --from-bp "$FROM" --to-bp "$TO" \
  --pheno "$PHENO" --pheno-name "$GENE" \
  --covar "$COVAR" --covar-name $COVAR_NAMES \
  --condition-list "$COND_LIST" \
  --linear hide-covar --allow-no-sex \
  --out "$REGC_PREF" >/dev/null

REGC_FIG="$FIGDIR/Fig_locus_${GENE}_center_${CENTER_SNP}_pm500kb_cond_${LEAD}.png"
python3 "$LOCUS_PY" \
  --assoc "$REGC_ASSOC" \
  --ld "$LD_FILE" \
  --lead "$LEAD" \
  --out-png "$REGC_FIG" \
  --title "Regional cis-eQTL association at 5q15 for ${GENE} conditioned on ${LEAD} (±500 kb)" \
  --test ADD \
  --gw-threshold 5e-8 \
  --point-size 10 >/dev/null

echo "[OK] outputs:"
echo "  - $GW_ASSOC"
echo "  - $REG_ASSOC"
echo "  - $REGC_ASSOC"
echo "  - $LD_FILE"
echo "  - $GW_FIG"
echo "  - $REG_FIG"
echo "  - $REGC_FIG"
