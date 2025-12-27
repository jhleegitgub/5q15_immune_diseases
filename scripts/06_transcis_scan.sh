#!/usr/bin/env bash
# scripts/06_transcis_scan.sh
# Candidate SNPs (from top20 proxy + top20 credible) -> genome-wide cis/trans scan (all phenos per chr)
# Output: result/06_transcis_scan/<TARGET>/<TARGET>_{all,sig}_genome.tsv
set -euo pipefail

# ----------------------------
# project paths (script-relative)
# ----------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# ----------------------------
# configurable inputs (env override OK)
# ----------------------------
BFILE="${BFILE:-$ROOT/../GenotypeData/GW.E-GEUV-3.EUR.MAF005.HWE1e-06}"
PHENO_DIR="${PHENO_DIR:-$ROOT/../PhenotypeData}"
PHENO_PATTERN="${PHENO_PATTERN:-chr%d_GD462.signalGeneQuantRPKM_plink.txt}"   # %d=chr
COVAR="${COVAR:-$ROOT/result/01_pca/covar_pca10.tsv}"
COVAR_NAMES="${COVAR_NAMES:-C1 C2 C3 C4 C5 C6 C7 C8 C9 C10}"
PLINK="${PLINK:-$HOME/Software/Plink/plink}"

TABLEDIR="${TABLEDIR:-$ROOT/table}"
CAND_DIR="${CAND_DIR:-$TABLEDIR/functional_candidates}"   # from step05 outputs

R2TH="${R2TH:-0.8}"
TOPN="${TOPN:-20}"

P_RAW="${P_RAW:-5e-6}"
Q_FDR="${Q_FDR:-0.1}"

# ERAP1 candidate set mode:
#   sig1  : ERAP1_sig1 proxy+credible
#   all   : (ERAP1_sig1+sig2+sig3) proxy+credible (union)
ERAP1_MODE="${ERAP1_MODE:-sig1}"

OUTROOT="${OUTROOT:-$ROOT/result/06_transcis_scan}"
CODEDIR="${CODEDIR:-$ROOT/code}"
COLLECT_PY="${COLLECT_PY:-$CODEDIR/collect_sig_genes_transcis_v2.py}"

# which targets to run (space-separated)
TARGETS_RAW="${TARGETS_RAW:-ERAP2 ERAP1 LNPEP}"

# ----------------------------
die(){ echo "[ERR] $*" 1>&2; exit 1; }
need(){ [[ -e "$1" ]] || die "missing: $1"; }

[[ -x "$PLINK" ]] || die "PLINK not executable: $PLINK"
need "${BFILE}.bed"; need "${BFILE}.bim"; need "${BFILE}.fam"
need "$COVAR"
need "$COLLECT_PY"
[[ -d "$PHENO_DIR" ]] || die "PHENO_DIR not found: $PHENO_DIR"
[[ -d "$CAND_DIR" ]] || die "CAND_DIR not found: $CAND_DIR"

mkdir -p "$OUTROOT"

# ----------------------------
# helpers
# ----------------------------
extract_rsids_from_top_tsv() {
  # args: in_top.tsv out_rsids.txt
  local in="$1" out="$2"
  need "$in"
  awk -F'\t' '
    NR==1{
      c=0
      for(i=1;i<=NF;i++){
        if($i=="SNP" || $i=="rsid" || $i=="RSID" || $i=="ID"){c=i}
      }
      if(c==0){c=1}
      next
    }
    {print $c}
  ' "$in" | sed '/^$/d' | sort -u > "$out"
}

first_existing() {
  # print first existing file among args
  for f in "$@"; do
    [[ -f "$f" ]] && { echo "$f"; return 0; }
  done
  return 1
}

build_candidates() {
  # args: target out_rsids
  local target="$1" out_rsids="$2"
  mkdir -p "$(dirname "$out_rsids")"

  local proxy cred
  local tag="$(echo "$R2TH" | tr '.' 'p')"   # 0.8 -> 0p8

  if [[ "$target" == "ERAP2" ]]; then
    proxy="$(first_existing \
      "$CAND_DIR/ERAP2_proxy_r2${R2TH}.top${TOPN}.tsv" \
      "$CAND_DIR/ERAP2_proxy_r2${tag}.top${TOPN}.tsv" \
      "$CAND_DIR/ERAP2_proxy_r2${R2TH}.top20.tsv" \
      "$CAND_DIR/ERAP2_proxy_r2${tag}.top20.tsv" \
    )" || die "ERAP2 proxy top table not found in $CAND_DIR"

    cred="$(first_existing \
      "$CAND_DIR/ERAP2_credible.top${TOPN}.tsv" \
      "$CAND_DIR/ERAP2_sig1_credible.top${TOPN}.tsv" \
      "$CAND_DIR/ERAP2_credible.top20.tsv" \
      "$CAND_DIR/ERAP2_sig1_credible.top20.tsv" \
    )" || die "ERAP2 credible top table not found in $CAND_DIR"

    tmp1="$(mktemp)"; tmp2="$(mktemp)"
    extract_rsids_from_top_tsv "$proxy" "$tmp1"
    extract_rsids_from_top_tsv "$cred"  "$tmp2"
    cat "$tmp1" "$tmp2" | sort -u > "$out_rsids"
    rm -f "$tmp1" "$tmp2"
    return 0
  fi

  if [[ "$target" == "LNPEP" ]]; then
    proxy="$(first_existing \
      "$CAND_DIR/LNPEP_proxy_r2${R2TH}.top${TOPN}.tsv" \
      "$CAND_DIR/LNPEP_proxy_r2${tag}.top${TOPN}.tsv" \
      "$CAND_DIR/LNPEP_proxy_r2${R2TH}.top20.tsv" \
      "$CAND_DIR/LNPEP_proxy_r2${tag}.top20.tsv" \
    )" || die "LNPEP proxy top table not found in $CAND_DIR"

    cred="$(first_existing \
      "$CAND_DIR/LNPEP_credible.top${TOPN}.tsv" \
      "$CAND_DIR/LNPEP_sig1_credible.top${TOPN}.tsv" \
      "$CAND_DIR/LNPEP_credible.top20.tsv" \
      "$CAND_DIR/LNPEP_sig1_credible.top20.tsv" \
    )" || die "LNPEP credible top table not found in $CAND_DIR"

    tmp1="$(mktemp)"; tmp2="$(mktemp)"
    extract_rsids_from_top_tsv "$proxy" "$tmp1"
    extract_rsids_from_top_tsv "$cred"  "$tmp2"
    cat "$tmp1" "$tmp2" | sort -u > "$out_rsids"
    rm -f "$tmp1" "$tmp2"
    return 0
  fi

  if [[ "$target" == "ERAP1" ]]; then
    local labels=()
    if [[ "$ERAP1_MODE" == "all" ]]; then
      labels=("ERAP1_sig1" "ERAP1_sig2" "ERAP1_sig3")
    else
      labels=("ERAP1_sig1")
    fi

    tmpU="$(mktemp)"
    : > "$tmpU"
    for lab in "${labels[@]}"; do
      proxy="$(first_existing \
        "$CAND_DIR/${lab}_proxy_r2${R2TH}.top${TOPN}.tsv" \
        "$CAND_DIR/${lab}_proxy_r2${tag}.top${TOPN}.tsv" \
        "$CAND_DIR/${lab}_proxy_r2${R2TH}.top20.tsv" \
        "$CAND_DIR/${lab}_proxy_r2${tag}.top20.tsv" \
      )" || die "$lab proxy top table not found in $CAND_DIR"

      cred="$(first_existing \
        "$CAND_DIR/${lab}_credible.top${TOPN}.tsv" \
        "$CAND_DIR/${lab}_sig1_credible.top${TOPN}.tsv" \
        "$CAND_DIR/${lab}_credible.top20.tsv" \
        "$CAND_DIR/${lab}_sig1_credible.top20.tsv" \
      )" || die "$lab credible top table not found in $CAND_DIR"

      t1="$(mktemp)"; t2="$(mktemp)"
      extract_rsids_from_top_tsv "$proxy" "$t1"
      extract_rsids_from_top_tsv "$cred"  "$t2"
      cat "$t1" "$t2" >> "$tmpU"
      rm -f "$t1" "$t2"
    done
    sort -u "$tmpU" > "$out_rsids"
    rm -f "$tmpU"
    return 0
  fi

  die "Unknown target: $target"
}

run_plink_chr() {
  # args: target chr rsids outprefix log
  local target="$1" chr="$2" rsids="$3" outprefix="$4" log="$5"

  local phenofile="$PHENO_DIR/$(printf "$PHENO_PATTERN" "$chr")"
  [[ -f "$phenofile" ]] || die "PHENO not found: $phenofile"

  "$PLINK" \
    --bfile "$BFILE" \
    --pheno "$phenofile" \
    --all-pheno \
    --extract "$rsids" \
    --covar "$COVAR" --covar-name $COVAR_NAMES \
    --linear hide-covar \
    --allow-no-sex \
    --out "$outprefix" \
    >"$log" 2>&1
}

# ----------------------------
# main
# ----------------------------
IFS=' ' read -r -a TARGETS <<< "$TARGETS_RAW"

for target in "${TARGETS[@]}"; do
  [[ -n "${target:-}" ]] || continue
  echo "[RUN] trans/cis scan target=$target"

  TDIR="$OUTROOT/$target"
  ADIR="$TDIR/assoc_by_chr"
  LDIR="$TDIR/logs"
  RDIR="$OUTROOT/rsids"
  mkdir -p "$ADIR" "$LDIR" "$RDIR"

  RSIDS="$RDIR/${target}_union_proxy${R2TH}_cred_top${TOPN}.rsids.txt"
  build_candidates "$target" "$RSIDS"
  echo "[OK] rsids: $RSIDS (n=$(wc -l < "$RSIDS" | tr -d ' '))"

  # plink per chr
  for chr in $(seq 1 22); do
    pref="$ADIR/${target}_cand_chr${chr}"
    # any assoc already exists?
    if ls "$ADIR/${target}_cand_chr${chr}".*.assoc.linear >/dev/null 2>&1; then
      echo "[SKIP] exists: $pref.*.assoc.linear"
      continue
    fi
    echo "  [RUN] chr$chr"
    run_plink_chr "$target" "$chr" "$RSIDS" "$pref" "$LDIR/${target}_chr${chr}.plink.log"
  done

  # collect + FDR
  OUT_ALL="$TDIR/${target}_all_genome.tsv"
  OUT_SIG="$TDIR/${target}_sig_genome.tsv"

  python3 "$COLLECT_PY" \
    --dir "$ADIR" \
    --p-raw "$P_RAW" \
    --q-fdr "$Q_FDR" \
    --out-all "$OUT_ALL" \
    --out-sig "$OUT_SIG"

  echo "[OK] $target => $OUT_SIG"
done

echo "[OK] done: $OUTROOT"
