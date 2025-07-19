from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import ld

plt.rcParams['svg.fonttype'] = 'none'

PALETTE = {
    "Significant up": "#a02428", "Significant down": "#393b92", "Not significant": "#D3D3D3",
}
ANNOTATE = {
    "protein_coding": {
        "Significant up": 300,
        "Significant down": 75,
    },
    "non_coding": {
        "Significant up": 150,
        "Significant down": 75,
    }
}
SETTINGS = {
    "protein_coding": dict(
        xlim=(-12, 12),
        xticks=[-12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12],
        ylim=(0, 300),
        yticks=[0, 50, 100, 150, 200, 250, 300],
        x_always_label=(-5, 7),
        y_always_label=50
    ),
    "non_coding": dict(
        xlim=(-12, 12),
        xticks=[-12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12],
        ylim=(0, 60),
        yticks=[0, 10, 20, 30, 40, 50, 60],
        x_always_label=(-4, 4),
        y_always_label=15
    ),
}


def volcano_plot(df: pd.DataFrame, title: str, settings: dict[str, Any], annotate: dict[str, int], saveto: Path):
    # Replace 0 p-values with a non-zero minimum
    zeropval = df['pvalue'] == 0
    df.loc[zeropval, 'pvalue'] = df.loc[df['pvalue'] != 0, 'pvalue'].min() / 2

    X, Y = df['log2FoldChange'], -np.log10(df['pvalue'])

    fig, ax = plt.subplots(figsize=(14, 8))

    # Color based on categories
    color = df['category'].apply(lambda x: PALETTE[x])

    # Clip values and use triangle markers for clipped values
    marker = np.asarray(['o' for _ in range(len(df))])
    xlim, ylim = settings['xlim'], settings['ylim']
    marker[X < xlim[0]] = '<'
    marker[X > xlim[1]] = '>'
    marker[Y < ylim[0]] = 'v'
    marker[(Y > ylim[1]) | zeropval] = '^'
    X, Y = np.clip(X, *xlim), np.clip(Y, *ylim)

    # Plot each marker separately
    for mk in "o<>^v":
        mask = (marker == mk)
        ax.scatter(X[mask], Y[mask], color=color[mask], marker=mk, s=15, lw=0.25, edgecolor='black')

    # Axes
    ax.set_xlabel("$log_2$(fold~change)", fontsize=12)
    ax.set_ylabel("$-log_{10}(p~value)$", fontsize=12)

    ax.spines[['right', 'top']].set_visible(False)

    ax.set(**{k: settings[k] for k in ['xlim', 'ylim', 'xticks', 'yticks']})

    # Fold-change threshold
    ax.axvline(ld.thresholds.log2fc, ls='--', color='black')
    ax.axvline(-ld.thresholds.log2fc, ls='--', color='black')

    # FDR labels
    thr = Y[df['padj'] <= ld.thresholds.padj].min()
    ax.axhline(thr, ls='--', color='black')
    ax.text(
        ax.get_xlim()[1] * 0.999, thr * 1.01, f'FDR ≤ {ld.thresholds.padj}', fontweight='bold',
        ha='right', va='bottom', fontsize=10
    )

    ax.set_title(title, fontsize=18, loc='left')

    # Annotate top hits or hits of interest
    df['X'], df['Y'] = X, Y
    df['rank'] = df['X'].abs().rank(method='first', ascending=False) + df['Y'].rank(method='first', ascending=False)
    for cat, ha in [('Significant up', 'left'), ('Significant down', 'right')]:
        xalways, yalways = settings['x_always_label'], settings['y_always_label']
        top_hits = [
            *df[df['category'] == cat].nsmallest(annotate[cat], 'rank').iterrows(),
            # Always label these
            *df[(df['category'] == cat) & (~df['X'].between(*xalways))].iterrows(),
            *df[(df['category'] == cat) & (df['Y'] >= yalways)].iterrows()
        ]
        annotated = set()
        for _, row in top_hits:
            if row['name'] in annotated:
                continue
            ax.text(row['X'], row['Y'], row['name'], fontsize=6, ha=ha, va='bottom', color='black')
            annotated.add(row['name'])
    # fig.show()
    saveto.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(saveto, bbox_inches="tight", pad_inches=0, transparent=True, dpi=400)
    plt.close(fig)


summary = pd.read_pickle(ld.DESeq2.summary)

workload = []
for partition, pname in ("protein_coding", "mRNA"), ("non_coding", "lncRNA"):
    partdf = summary[summary['partition'] == partition]
    for ifn in ["IFNb"]:
        columns = [(ifn, 'mock', x) for x in ['category', 'log2FoldChange', 'padj', 'pvalue']]
        df = partdf[columns + ['name']].copy().rename(columns={x: x[-1] for x in columns})
        saveto = ld.plots.volcano / partition / f"{ifn}_vs_mock [labeled].svg"
        title = f"{ifn}-β ({pname})"
        workload.append(delayed(volcano_plot)(df, title, SETTINGS[partition], ANNOTATE[partition], saveto))

Parallel(n_jobs=-1, backend='multiprocessing')(workload)
