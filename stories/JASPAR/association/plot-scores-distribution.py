import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import ld
from stories.DE import IFNS

plt.rcParams['svg.fonttype'] = 'none'

summary = pd.read_pickle(ld.TXGROUP_SUMMARY)

motif = ('cluster', 'ISRE-like')
assert motif in summary.columns, f"Motif {motif} not found in the summary"

observations = []

# Select the background category - transcripts that were unresponsive to any IFN
mask = summary['tags'].apply(lambda x: "Background" in x)
subdf = summary[mask][[motif]].copy()
subdf['Category'] = 'Background'
observations.append(subdf)

# Select the significant upregulated transcripts for each IFN
for ifn in IFNS:
    mask = summary['tags'].apply(lambda x: (ifn, "mock", "Significant up") in x)
    subdf = summary[mask][[motif]].copy()
    subdf['Category'] = ifn
    observations.append(subdf)

df = pd.concat(observations, axis=0, ignore_index=True)

# Plot the distribution of motif scores for each category
grid = sns.displot(
    df, x=motif, row='Category', kind='hist', hue='Category', stat='density',
    common_bins=True, common_norm=False
)


def plot(data, **kwargs):
    ax = plt.gca()
    for percentile in [0.1, 0.25, 0.5, 0.75, 0.9]:
        value = data.quantile(percentile)
        ax.axvline(value, color='red', linestyle='--', label=f'{int(percentile * 100)}th Percentile')
        ax.text(
            value, 1.0, f'{int(percentile * 100)}th:\n{value:.2f}',
            horizontalalignment='left', verticalalignment='top', color='red',
            transform=ax.get_xaxis_transform()
        )
    ax.set_xlim(-6, 6)


grid.map(plot, motif)
grid.figure.suptitle(f"Distribution of {motif[1]}", fontsize=16)

grid.figure.savefig(ld.RESULTS / f"scores-distribution-{motif[1]}.svg", bbox_inches='tight', pad_inches=0)
# grid.figure.show()
plt.close(grid.figure)
