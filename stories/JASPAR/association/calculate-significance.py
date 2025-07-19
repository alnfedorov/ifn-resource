import pickle

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy import stats

import ld
from stories import DE


def run(ifn, control, motif, targets):
    # Mark motifs with Z-score above the threshold as True ('positive') and below average as False
    hasmotif = targets[motif] >= ld.thresholds.zscore[1]
    nomotif = targets[motif] <= ld.thresholds.zscore[0]

    if hasmotif.sum() == 0 or nomotif.sum() == 0:
        print(f"Skipping {ifn} vs {control} for motif {motif[1]}: no targets with motif")
        return None

    hasmotif = targets[hasmotif][(ifn, control, 'log2FoldChange')].values
    nomotif = targets[nomotif][(ifn, control, 'log2FoldChange')].values

    # Wilcoxon rank-sum (Mann Whitney U) test between the two groups
    pvalue = stats.mannwhitneyu(hasmotif, nomotif, alternative='two-sided').pvalue

    # Calculate mean/median delta log2 fold change
    mean_delta_log2fc = np.average(hasmotif) - np.average(nomotif)
    median_delta_log2fc = np.median(hasmotif) - np.median(nomotif)

    return {
        'motif': motif[1], 'target': ifn, 'control': control, 'p-value': pvalue,
        'Mean Δ(log2 fold change)': mean_delta_log2fc,
        'Median Δ(log2 fold change)': median_delta_log2fc,
        'With motif': len(hasmotif), 'Without motif': len(nomotif),
    }


summary = pd.read_pickle(ld.TXGROUP_SUMMARY)
motifs = summary.filter(like='cluster').columns

workload = []
for ifn in DE.IFNS:
    for control in ['mock']:
        if ifn == control:
            continue
        mask = (summary[(ifn, 'TPM')] >= ld.thresholds.min_tpm) | (summary[(control, 'TPM')] >= ld.thresholds.min_tpm)
        print(f"Calculating significance for {ifn}-vs-{control} with {mask.sum()} ({mask.mean():.1%}) targets")
        targets = summary[mask]
        for motif in motifs:
            workload.append(delayed(run)(ifn, control, motif, targets))

results = Parallel(n_jobs=-1, backend='threading', verbose=1000, pre_dispatch='all')(workload)
results = [res for res in results if res is not None]

df = pd.DataFrame(results)
ld.STAT_TESTS.parent.mkdir(parents=True, exist_ok=True)
df.to_pickle(ld.STAT_TESTS, protocol=pickle.HIGHEST_PROTOCOL)
