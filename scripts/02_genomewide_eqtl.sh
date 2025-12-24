# scripts/02_genomewide_eqtl.sh
#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

PLINK="${PLINK:-$HOME/Software/Plink/plink}"
PYTHON="${PYTHON:-python3}"

BFILE="${BFILE:-$ROOT/../GenotypeData/GW.E-GEUV-3.EUR.MAF005.HWE1e-06}"
PHENO="${PHENO:-$ROOT/../PhenotypeData/chr5_GD462.signalGeneQuantRPKM_plink.txt}"
COVAR="${COVAR:-$ROOT/result/01_pca/covar_pca10.tsv}"
COVAR_NAMES="${COVAR_NAMES:-C1 C2 C3 C4 C5 C6 C7 C8 C9 C10}"

OUTDIR="$ROOT/result/02_eqtl_genomewide"
FIGDIR="$ROOT/fig"
mkdir -p "$OUTDIR" "$FIGDIR"

CODE_MANH="$ROOT/code/manhattan_genomewide.py"

req() { [[ -f "$1" ]] || { echo "[ERR] not found: $1" >&2; exit 1; }; }
req "${BFILE}.bed"; req "${BFILE}.bim"; req "${BFILE}.fam"
req "$PHENO"
req "$COVAR"
req "$CODE_MANH"

run_plink_gene() {
  local gene="$1"
  local out_prefix="$OUTDIR/${gene}_genomewide"
  local assoc="${out_prefix}.assoc.linear"

  if [[ -f "$assoc" ]]; then
    echo "[SKIP] exists: $assoc"
  else
    echo "[RUN] PLINK genome-wide: $gene"
    "$PLINK" \
      --bfile "$BFILE" \
      --pheno "$PHENO" --pheno-name "$gene" \
      --covar "$COVAR" --covar-name $COVAR_NAMES \
      --linear hide-covar \
      --allow-no-sex \
      --out "$out_prefix"
  fi

  local out_png="$FIGDIR/Fig_genomewide_${gene}_manhattan.png"
  echo "[RUN] plot: $out_png"
  "$PYTHON" "$CODE_MANH" \
    --assoc "$assoc" \
    --out-png "$out_png" \
    --title "Genome-wide eQTL scan for ${gene} expression" \
    --gw-threshold 5e-8 \
    --test ADD
}

run_plink_gene "ERAP2"
run_plink_gene "ERAP1"
run_plink_gene "LNPEP"

# 3-panel combined
COMBO_PNG="$FIGDIR/Fig_genomewide_3genes_manhattan.png"
echo "[RUN] combined plot: $COMBO_PNG"
"$PYTHON" "$CODE_MANH" \
  --assoc "$OUTDIR/ERAP2_genomewide.assoc.linear,$OUTDIR/ERAP1_genomewide.assoc.linear,$OUTDIR/LNPEP_genomewide.assoc.linear" \
  --out-png "$COMBO_PNG" \
  --titles "ERAP2,ERAP1,LNPEP" \
  --suptitle "Genome-wide eQTL landscapes of ERAP2, ERAP1 and LNPEP expression" \
  --gw-threshold 5e-8 \
  --test ADD \
  --xtick-step 1

echo "[OK] 02 done."
