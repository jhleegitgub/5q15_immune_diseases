# manhattan_genomewide.py (수정 포인트만 포함된 "전체 덮어쓰기" 버전)
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def read_assoc_linear(path: str, test: str = "ADD") -> pd.DataFrame:
    df = pd.read_table(path, sep=r"\s+", dtype=str)
    # PLINK --linear output: CHR SNP BP A1 TEST NMISS BETA STAT P
    df = df[df["TEST"] == test].copy()
    df["CHR"] = pd.to_numeric(df["CHR"], errors="coerce")
    df["BP"] = pd.to_numeric(df["BP"], errors="coerce")
    df["P"] = pd.to_numeric(df["P"], errors="coerce")
    df = df.dropna(subset=["CHR", "BP", "P"])
    df = df[(df["CHR"] >= 1) & (df["CHR"] <= 22)]
    return df

def manhattan_ax(ax, df: pd.DataFrame, gw_threshold: float, max_points: int, seed: int, xtick_step: int):
    df = df.sort_values(["CHR", "BP"]).copy()
    df["LOGP"] = -np.log10(np.clip(df["P"].values, 1e-300, 1.0))

    # downsample (랜덤이 아니라 균일 샘플링 + 상위 신호 보존)
    if max_points and len(df) > max_points:
        rng = np.random.default_rng(seed)
        top = df.nlargest(min(20000, len(df)), "LOGP")
        rest = df.drop(top.index)
        k = max_points - len(top)
        if k > 0 and len(rest) > k:
            rest = rest.sample(n=k, random_state=seed)
        df = pd.concat([top, rest], ignore_index=True).sort_values(["CHR", "BP"])

    # cumulative x
    chr_max = df.groupby("CHR")["BP"].max().to_dict()
    chr_order = list(range(1, 23))
    offsets = {}
    cum = 0
    for c in chr_order:
        offsets[c] = cum
        cum += chr_max.get(c, 0) + 1_000_000  # 1Mb gap

    df["X"] = df.apply(lambda r: r["BP"] + offsets[int(r["CHR"])], axis=1)

    # tick at chr centers
    ticks = []
    labels = []
    for c in chr_order:
        if c not in chr_max:
            continue
        start = offsets[c]
        end = offsets[c] + chr_max[c]
        ticks.append((start + end) / 2)
        labels.append(str(c))

    # alternate color by chr (단색/회색톤)
    colors = []
    for c in df["CHR"].astype(int).values:
        colors.append("0.2" if (c % 2 == 1) else "0.65")

    ax.scatter(df["X"].values, df["LOGP"].values, s=1.0, c=colors, linewidths=0, rasterized=True)

    ax.axhline(-np.log10(gw_threshold), linewidth=1.0)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel(r"$-\log_{10}(P)$")

    # xtick-step 적용
    if xtick_step and xtick_step > 1:
        ticks2 = ticks[::xtick_step]
        labels2 = labels[::xtick_step]
    else:
        ticks2, labels2 = ticks, labels

    ax.set_xticks(ticks2)
    ax.set_xticklabels(labels2)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assoc", required=True, help="assoc.linear path, or comma-separated list")
    ap.add_argument("--out-png", required=True)
    ap.add_argument("--title", default=None)
    ap.add_argument("--titles", default=None, help="comma-separated titles for multi-panel")
    ap.add_argument("--suptitle", default=None, help="overall title for multi-panel")
    ap.add_argument("--gw-threshold", type=float, default=5e-8)
    ap.add_argument("--test", default="ADD")
    ap.add_argument("--max-points", type=int, default=800000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--xtick-step", type=int, default=1)
    args = ap.parse_args()

    assoc_list = [x.strip() for x in args.assoc.split(",") if x.strip()]
    titles = None
    if args.titles:
        titles = [x.strip() for x in args.titles.split(",")]
        if len(titles) != len(assoc_list):
            raise SystemExit("[ERR] --titles count mismatch with --assoc")

    dfs = [read_assoc_linear(p, test=args.test) for p in assoc_list]
    n = len(dfs)

    if n == 1:
        fig_w, fig_h = 12, 4.8
        fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h), dpi=150)
        manhattan_ax(ax, dfs[0], args.gw_threshold, args.max_points, args.seed, args.xtick_step)
        if args.title:
            ax.set_title(args.title)
        # x축 라벨 크기
        ax.tick_params(axis="x", labelsize=9)
        fig.tight_layout()
    else:
        # ✅ 여기 핵심: 패널당 가로를 키움 (9.5~10 권장)
        panel_w = 10.0
        fig_w, fig_h = panel_w * n, 4.6
        fig, axes = plt.subplots(1, n, figsize=(fig_w, fig_h), dpi=150)
        if n == 1:
            axes = [axes]
        for i, ax in enumerate(axes):
            manhattan_ax(ax, dfs[i], args.gw_threshold, args.max_points, args.seed, args.xtick_step)
            ax.tick_params(axis="x", labelsize=8)  # ✅ 글자 조금 줄여서 21/22 겹침 방지
            ax.set_title(titles[i] if titles else f"panel{i+1}")
            if i != 0:
                ax.set_ylabel("")
        if args.suptitle:
            fig.suptitle(args.suptitle, y=1.02)
        fig.subplots_adjust(wspace=0.18)
        fig.tight_layout()

    fig.savefig(args.out_png, bbox_inches="tight")
    plt.close(fig)

if __name__ == "__main__":
    main()
