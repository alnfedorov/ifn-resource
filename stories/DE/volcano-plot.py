from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from matplotlib.lines import Line2D

import ld

plt.rcParams['svg.fonttype'] = 'none'

PALETTE = {
    "Significant up": "#a02428", "Significant down": "#393b92", "Not significant": "#D3D3D3"
}
PALETTE_ORDER = ["Not significant", "Significant up", "Significant down"]
XLIM = (-12, 12)
XTICKS = [-12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12]
YLIM = (0, 300)
YTICKS = [0, 100, 200, 300]


def volcano_plot(df: pd.DataFrame, title: str, saveto: Path):
    X = df['log2FoldChange']

    zeropval = df['pvalue'] == 0
    Y = np.repeat(YLIM[1] + 1, len(df))
    Y[~zeropval] = -np.log10(df.loc[~zeropval, 'pvalue'])

    fig, ax = plt.subplots(figsize=(8, 7))

    # Clip values and use triangle markers for clipped values
    marker = np.asarray(['o' for _ in range(len(df))])
    marker[X < XLIM[0]] = '<'
    marker[X > XLIM[1]] = '>'
    marker[Y < YLIM[0]] = 'v'
    marker[(Y > YLIM[1]) | zeropval] = '^'
    X, Y = np.clip(X, *XLIM), np.clip(Y, *YLIM)

    # Plot each marker and category separately. First plot the 'Not significant' category to ensure it is at the bottom.
    for mk in "o<>^v":
        for cat in PALETTE_ORDER:
            mask = (marker == mk) & (df['category'] == cat)
            ax.scatter(
                X[mask], Y[mask], color=PALETTE[cat], marker=mk, s=35, lw=0.25, edgecolor='black', rasterized=True
            )

    # Axes
    ax.set_xlabel("$log_2$(fold~change)", fontsize=12)
    ax.set_ylabel("$-log_{10}(p~value)$", fontsize=12)

    ax.spines[['right', 'top']].set_visible(False)

    ax.set_xlim(XLIM[0], XLIM[1])
    ax.set_xticks(XTICKS)

    ax.set_ylim(YLIM[0], YLIM[1])
    ax.set_yticks(YTICKS)

    # Fold-change threshold
    ax.axvline(ld.thresholds.log2fc, ls='--', color='black')
    ax.axvline(-ld.thresholds.log2fc, ls='--', color='black')

    # FDR labels
    thr = Y[df['padj'] <= ld.thresholds.padj].min()
    ax.axhline(thr, ls='--', color='black')
    ax.text(
        ax.get_xlim()[0], thr * 1.01, f'$\\bf{{FDR~≤~{ld.thresholds.padj}}}$', ha='left', va='bottom', fontsize=10
    )

    legend, counts = [], df['category'].value_counts().to_dict()
    for description, color in PALETTE.items():
        if counts.get(description, 0) > 0:
            legend.append(Line2D(
                [0], [0], marker='o', color='white', label=f"[N={counts[description]}] {description}",
                markerfacecolor=color, markersize=10
            ))
    ax.legend(handles=legend, loc='upper left', labelspacing=0.15, handletextpad=0.02, fontsize=16)

    # Upregulated / Downregulated counts
    for category, x, ha in [('Significant up', 0.95, 'right'), ('Significant down', 0.05, 'left')]:
        count = counts.get(category, 0)
        ax.text(x, 0.5, count, color=PALETTE[category], ha=ha, va='center', fontsize=16, transform=ax.transAxes)

    ax.set_title(title, fontsize=18, loc='left')

    # fig.show()
    saveto.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(saveto, bbox_inches="tight", pad_inches=0, transparent=True, dpi=400)
    plt.close(fig)


summary = pd.read_pickle(ld.DESeq2.summary)

workload = []
for partition, pname in ("protein_coding", "mRNA"), ("non_coding", "lncRNA"), ("artifacts", "artifacts"):
    partdf = summary[summary['partition'] == partition]
    for ifn in ld.IFNS:
        columns = [(ifn, 'mock', x) for x in ['category', 'log2FoldChange', 'padj', 'pvalue']]
        df = partdf[columns].copy().rename(columns={x: x[-1] for x in columns})
        if len(df) == 0:
            print(f"No data for {ifn} vs mock in partition {pname}")
            continue
        saveto = ld.plots.volcano / partition / f"{ifn}_vs_mock.svg"
        title = f"{ifn} vs mock ({pname})"
        workload.append(delayed(volcano_plot)(df, title, saveto))

Parallel(n_jobs=-1, backend='multiprocessing')(workload)
