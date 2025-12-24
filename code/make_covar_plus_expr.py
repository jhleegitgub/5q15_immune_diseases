#!/usr/bin/env python3
import sys, csv

def detect_delim(path: str):
    with open(path, "r", newline="") as f:
        line = f.readline()
    return "\t" if "\t" in line else None  # None -> whitespace split

def read_table(path: str):
    delim = detect_delim(path)
    rows = []
    with open(path, "r", newline="") as f:
        if delim:
            r = csv.reader(f, delimiter=delim)
            for row in r:
                if row:
                    rows.append(row)
        else:
            for line in f:
                line = line.strip()
                if line:
                    rows.append(line.split())
    return rows

def main():
    if len(sys.argv) != 5:
        print("usage: make_covar_plus_expr.py <base_covar.tsv> <pheno.txt> <expr_gene> <out.tsv>", file=sys.stderr)
        sys.exit(2)

    base, pheno, gene, out = sys.argv[1:]

    base_rows = read_table(base)
    if len(base_rows) < 2:
        raise SystemExit("base covar too small")

    header = base_rows[0]
    has_header = (len(header) >= 3 and not header[2].replace('.', '', 1).replace('-', '', 1).isdigit())

    if has_header:
        base_header = header
        data_rows = base_rows[1:]
    else:
        base_header = ["FID", "IID"] + [f"C{i}" for i in range(1, len(base_rows[0]) - 1)]
        data_rows = base_rows

    base_map = {}
    for row in data_rows:
        if len(row) < 3:
            continue
        key = (row[0], row[1])
        base_map[key] = row

    ph_rows = read_table(pheno)
    ph_header = ph_rows[0]
    try:
        idx = ph_header.index(gene)
    except ValueError:
        raise SystemExit(f"gene column not found in pheno header: {gene}")

    expr_map = {}
    for row in ph_rows[1:]:
        if len(row) <= idx:
            continue
        key = (row[0], row[1])
        val = row[idx]
        if val in ("NA", "", "nan", "NaN", "."):
            continue
        expr_map[key] = val

    merged = []
    for key, brow in base_map.items():
        if key in expr_map:
            merged.append(brow + [expr_map[key]])

    with open(out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(base_header + [f"expr_{gene}"])
        for row in merged:
            w.writerow(row)

if __name__ == "__main__":
    main()
