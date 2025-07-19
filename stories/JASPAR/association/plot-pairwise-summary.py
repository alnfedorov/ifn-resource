import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

import ld
from stories.DE import IFNS

plt.rcParams['svg.fonttype'] = 'none'

PALETTE = {'Not significant': '#93a2cb', 'Significant': '#da6434'}
ORDER = IFNS
LABEL = 100

Y = 'p-value'
YLIM = (0, 80)

X = 'Median Δ(log2 fold change)'
XLIM = (-1, 1)

# Load the results of the statistical tests
df = pd.read_pickle(ld.STAT_TESTS)

df['-log10(p-value)'] = -df[Y].apply(lambda x: np.log10(x))

# Sanity checks
assert df[X].between(XLIM[0], XLIM[1]).all(), f"X values are not within the limits: {df[X].min()} - {df[X].max()}"
assert df['-log10(p-value)'].between(YLIM[0], YLIM[1]).all(), \
    f"Y values are not within the limits: {df['-log10(p-value)'].min()} - {df['-log10(p-value)'].max()}"

controls = df['control'].unique()
rows, cols = len(IFNS), len(controls)

fig, axes = plt.subplots(rows, cols, figsize=(cols * 4 + 1, rows * 3 + 1), sharex=True, sharey=True)
axes = axes.reshape(rows, cols)
for i, target in enumerate(IFNS):
    for j, control in enumerate(controls):
        ax = axes[i, j]
        subdf = df[(df['target'] == target) & (df['control'] == control)].copy()

        # Mark siginificant motifs
        subdf['q-value'] = stats.false_discovery_control(subdf[Y], method='bh')
        subdf['Category'] = 'Not significant'
        subdf.loc[subdf['q-value'] <= ld.thresholds.qvalue, 'Category'] = 'Significant'

        sns.scatterplot(
            data=subdf, x=X, y='-log10(p-value)', hue='Category',
            s=50, edgecolor='white', linewidth=0.75, palette=PALETTE, legend=False, ax=ax
        )

        ax.spines[['top', 'right']].set_visible(False)
        ax.set(
            xlabel=X, ylabel='-log10(p-value)', xlim=XLIM, ylim=YLIM, title=f"{target} vs {control}"
        )

        # FDR line
        minpval = subdf[subdf['Category'] != 'Not significant']['-log10(p-value)']
        if not minpval.empty:
            minpval = minpval.min()
            ax.axhline(minpval, color='black', linestyle='--', linewidth=1)
            ax.text(1.0, minpval, f'FDR ≤ {ld.thresholds.qvalue}', ha='right', va='bottom', color='black',
                    transform=ax.get_yaxis_transform())

        # Label top hits
        signif = subdf[subdf['Category'] != 'Not significant']
        signif = signif.sort_values(by='-log10(p-value)', ascending=False).head(LABEL)
        for motif, delta, pval, qval in signif[['motif', X, '-log10(p-value)', 'q-value']].itertuples(index=False):
            ax.text(
                delta, pval, motif, ha='left', va='center', color='black', fontsize=4, fontweight='bold',
                zorder=1000
            )

fig.suptitle("Statistical significance of motif associations", fontsize=16)
fig.savefig(ld.RESULTS / "motif-association-significance.svg", bbox_inches='tight', pad_inches=0)
# fig.show()
plt.close(fig)
