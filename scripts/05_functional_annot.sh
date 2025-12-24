#!/usr/bin/env bash
# scripts/05_functional_annot.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_config.sh"

# --------------------------
# Inputs (env override OK)
# --------------------------
FINEMAP_DIR="${FINEMAP_DIR:-$RESULTDIR/04_finemap_pip}"
SIGNALS_TSV="${SIGNALS_TSV:-$RESULTDIR/03_signal_check_5q15/signals_summary.tsv}"

OUTDIR="${OUTDIR:-$RESULTDIR/05_functional_annot}"
WINKB="${WINKB:-500}"          # +/- 500kb
R2TH="${R2TH:-0.8}"            # proxy threshold
TOPN="${TOPN:-20}"             # top candidates to export

VEP_PY="${VEP_PY:-$CODEDIR/quick_vep_grch37_v2.py}"
RANK_PY="${RANK_PY:-$CODEDIR/rank_candidates_from_vep37_v2.py}"

mkdir -p "$OUTDIR"/{credible,proxy,ld,rsids,vep,rank,logs}
mkdir -p "$TABLEDIR"/functional_candidates

# --------------------------
# helpers
# --------------------------
die(){ echo "[ERR] $*" 1>&2; exit 1; }

need_cmd(){
  command -v "$1" >/dev/null 2>&1 || die "missing command: $1"
}

zcat_auto(){
  # usage: zcat_auto file.gz
  if command -v zcat >/dev/null 2>&1; then
    zcat "$1"
  else
    gzip -cd "$1"
  fi
}

get_chr_bp(){
  # prints: chr bp
  local snp="$1"
  awk -v s="$snp" '$2==s{print $1"\t"$4; exit}' "${BFILE}.bim"
}

mk_window(){
  # prints: chr from to
  local snp="$1"
  local chr bp
  read -r chr bp < <(get_chr_bp "$snp")
  [[ -n "${chr:-}" && -n "${bp:-}" ]] || die "SNP not found in .bim: $snp"
  local from=$((bp - WINKB*1000))
  local to=$((bp + WINKB*1000))
  if [[ $from -lt 1 ]]; then from=1; fi
  echo -e "${chr}\t${from}\t${to}"
}

extract_credible_rsids(){
  # supports TSV with header containing SNP column, or first column
  local credible_tsv="$1"
  local out_rsids="$2"
  awk -F'\t' '
    NR==1{
      for(i=1;i<=NF;i++){
        if($i=="SNP"){snpcol=i}
      }
      if(snpcol==0){snpcol=1}
      next
    }
    {print $snpcol}
  ' "$credible_tsv" | sed '/^$/d' | sort -u > "$out_rsids"
}

make_proxy_rsids(){
  local label="$1"
  local lead="$2"

  local chr from to
  read -r chr from to < <(mk_window "$lead")

  local pref="$OUTDIR/ld/${label}_proxy_r2${R2TH}"
  "$PLINK" \
    --bfile "$BFILE" \
    --chr "$chr" --from-bp "$from" --to-bp "$to" \
    --ld-snp "$lead" \
    --r2 gz --ld-window 99999 --ld-window-kb "$WINKB" --ld-window-r2 "$R2TH" \
    --out "$pref" >/dev/null

  local rsids="$OUTDIR/rsids/${label}_proxy_r2${R2TH}.rsids.txt"
  {
    echo "$lead"
    zcat_auto "${pref}.ld.gz" | awk 'NR>1{print $3"\n"$6}'
  } | sed '/^$/d' | sort -u > "$rsids"

  echo "$rsids"
}

make_ld_window(){
  # LD (lead vs all SNPs in window): output .ld.gz
  local label="$1"
  local lead="$2"
  local tag="$3"

  local chr from to
  read -r chr from to < <(mk_window "$lead")

  local pref="$OUTDIR/ld/${label}_${tag}"
  "$PLINK" \
    --bfile "$BFILE" \
    --chr "$chr" --from-bp "$from" --to-bp "$to" \
    --ld-snp "$lead" \
    --r2 gz --ld-window 99999 --ld-window-kb "$WINKB" --ld-window-r2 0 \
    --out "$pref" >/dev/null

  echo "${pref}.ld.gz"
}

find_finemap_paths(){
  # returns: credible_path \t pip_path
  local label="$1"
  local cred1="$FINEMAP_DIR/${label}_credible95.tsv"
  local pip1="$FINEMAP_DIR/${label}_pip.tsv"

  local cred="$cred1"
  local pip="$pip1"

  if [[ ! -f "$cred" ]]; then
    cred="$(ls -1 "$FINEMAP_DIR"/"${label}"*credible*95*.tsv 2>/dev/null | head -1 || true)"
  fi
  if [[ ! -f "$pip" ]]; then
    pip="$(ls -1 "$FINEMAP_DIR"/"${label}"*pip*.tsv 2>/dev/null | head -1 || true)"
  fi

  [[ -f "$cred" ]] || die "missing credible95 for $label in $FINEMAP_DIR"
  [[ -f "$pip"  ]] || die "missing pip for $label in $FINEMAP_DIR"

  echo -e "${cred}\t${pip}"
}

run_one_set(){
  local label="$1" lead="$2" rsids="$3" ld_gz="$4" pip_tsv="$5" out_sub="$6"

  local vep_tsv="$OUTDIR/vep/${label}_${out_sub}.vep.tsv"
  local ranked_tsv="$OUTDIR/rank/${label}_${out_sub}.ranked.tsv"
  local top_tsv="$OUTDIR/rank/${label}_${out_sub}.top${TOPN}.tsv"

  python3 "$VEP_PY" "$rsids" "$vep_tsv" >"$OUTDIR/logs/${label}_${out_sub}.vep.log" 2>&1

  python3 "$RANK_PY" \
    --ld "$ld_gz" \
    --vep "$vep_tsv" \
    --lead "$lead" \
    --pip "$pip_tsv" \
    --out "$ranked_tsv" \
    --top "$TOPN" \
    --top-out "$top_tsv" \
    >"$OUTDIR/logs/${label}_${out_sub}.rank.log" 2>&1

  # copy top list to table dir (stable name)
  cp -f "$top_tsv" "$TABLEDIR/functional_candidates/${label}_${out_sub}.top${TOPN}.tsv"
}

# --------------------------
# main
# --------------------------
need_cmd awk
need_cmd sed
need_cmd sort
need_cmd gzip
[[ -x "$PLINK" ]] || die "PLINK not executable: $PLINK"
[[ -f "${BFILE}.bim" ]] || die "BFILE not found: ${BFILE}.bim"
[[ -f "$SIGNALS_TSV" ]] || die "signals_summary.tsv not found: $SIGNALS_TSV"
[[ -f "$VEP_PY" ]] || die "missing: $VEP_PY"
[[ -f "$RANK_PY" ]] || die "missing: $RANK_PY"

MANIFEST="$OUTDIR/functional_candidates_manifest.tsv"
echo -e "gene\tsignal_id\tlabel\tlead_snp\tset_type\trsids_path\tpip_path\tcredible_path\tld_gz\ttop_table\toutdir" > "$MANIFEST"

echo "[RUN] functional annotation => $OUTDIR"
echo "[INFO] FINEMAP_DIR=$FINEMAP_DIR"
echo "[INFO] SIGNALS_TSV=$SIGNALS_TSV"

# expected columns: gene  signal_id  snp
while IFS=$'\t' read -r gene sid lead; do
  [[ -n "${gene:-}" && -n "${sid:-}" && -n "${lead:-}" ]] || continue

  # finemap prefix naming
  label="$gene"
  if [[ "$gene" == "ERAP1" ]]; then
    case "$sid" in
      signal1) label="ERAP1_sig1" ;;
      signal2) label="ERAP1_sig2" ;;
      signal3) label="ERAP1_sig3" ;;
      *) label="ERAP1_${sid}" ;;
    esac
  fi

  read -r credible pip < <(find_finemap_paths "$label" | awk -F'\t' '{print $1, $2}')

  echo "[RUN] $label lead=$lead"

  # ---- A) credible-set 기반 ----
  rs_cred="$OUTDIR/rsids/${label}_credible.rsids.txt"
  extract_credible_rsids "$credible" "$rs_cred"
  { echo "$lead"; cat "$rs_cred"; } | sed '/^$/d' | sort -u > "${rs_cred}.tmp" && mv "${rs_cred}.tmp" "$rs_cred"

  ld_cred="$(make_ld_window "$label" "$lead" "credible")"
  run_one_set "$label" "$lead" "$rs_cred" "$ld_cred" "$pip" "credible"

  topA="$TABLEDIR/functional_candidates/${label}_credible.top${TOPN}.tsv"
  echo -e "${gene}\t${sid}\t${label}\t${lead}\tcredible\t${rs_cred}\t${pip}\t${credible}\t${ld_cred}\t${topA}\t${OUTDIR}" >> "$MANIFEST"

  # ---- B) proxy(r2>=R2TH) 기반 ----
  rs_proxy="$(make_proxy_rsids "$label" "$lead")"
  ld_proxy="$(make_ld_window "$label" "$lead" "proxy_r2${R2TH}")"
  run_one_set "$label" "$lead" "$rs_proxy" "$ld_proxy" "$pip" "proxy_r2${R2TH}"

  topB="$TABLEDIR/functional_candidates/${label}_proxy_r2${R2TH}.top${TOPN}.tsv"
  echo -e "${gene}\t${sid}\t${label}\t${lead}\tproxy_r2${R2TH}\t${rs_proxy}\t${pip}\t${credible}\t${ld_proxy}\t${topB}\t${OUTDIR}" >> "$MANIFEST"

done < <(tail -n +2 "$SIGNALS_TSV")

# copy manifest to table dir too
cp -f "$MANIFEST" "$TABLEDIR/functional_candidates_manifest.tsv"

echo "[OK] done."
echo "[OK] manifest: $MANIFEST"
echo "[OK] top tables: $TABLEDIR/functional_candidates/"
