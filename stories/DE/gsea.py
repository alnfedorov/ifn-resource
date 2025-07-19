import gseapy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import ld
from stories import terminus

plt.rcParams['svg.fonttype'] = 'none'

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

# Prepare target gene sets
# names = gseapy.get_library_name(organism='Human')
# GO_Biological_Process_2025 / KEGG_2021_Human

libname = 'GO_Biological_Process_2025'
library = gseapy.get_library(name=libname, organism='Human')
genes_in_library = set.union(*[set(v) for v in library.values()])
print(f"Number of genes in {libname} library: {len(genes_in_library)}")

# Load target genes - up/down-regulated by all IFNs
group2gene = pd.read_csv(terminus.GROUP2GENE, sep='\t', index_col=0)['gene_names'].to_dict()
group2gene = {k: v.split('+') for k, v in group2gene.items()}

summary = pd.read_pickle(ld.DESeq2.summary)
summary['genes'] = summary.index.map(group2gene)

gene_lists = {}
for k in ['Significant up', 'Significant down']:
    tags = {(ifn, 'mock', k) for ifn in ld.IFNS}
    mask = summary['tags'].apply(lambda x: len(x & tags) == len(tags)) & summary['partition'].eq('protein_coding')
    allgenes = summary.loc[mask, 'genes'].explode().unique()
    gene_lists[k] = sorted(set(allgenes) & genes_in_library)
    print(f"Number of {k} genes in {libname} library: {len(gene_lists[k])}")

# Select all expressed genes - background
expressed = set(summary['genes'].explode().unique()) & genes_in_library
background = sorted(expressed)
print(f"Number of expressed genes in {libname} library: {len(background)}")

# Run hypergeometric test to find enriched terms
records = []
for k, genes in gene_lists.items():
    results = gseapy.enrich(genes, library, background=background, no_plot=True).results

    results['Category'] = k
    results['Category size'] = len(genes)

    results['Hits'] = results['Overlap'].apply(lambda x: int(x.split('/')[0]))
    results['Term size'] = results['Overlap'].apply(lambda x: int(x.split('/')[1]))
    results['Hits (%)'] = results['Hits'] / len(genes) * 100

    records.append(results)

df = pd.concat(records, ignore_index=True)
df = df.drop(columns=['Overlap', 'P-value', 'Odds Ratio', 'Combined Score'])
df['-log10(Adjusted P-value)'] = -np.log10(df['Adjusted P-value'])

hue_norm = (1, round(df['Hits (%)'].max()))
xlim = (0, df['-log10(Adjusted P-value)'].max())

for cat, subdf in df.groupby('Category'):
    subdf = subdf[subdf['Adjusted P-value'] <= 0.1].copy()
    subdf = subdf.sort_values('Adjusted P-value').head(50)

    size_norm = (1, int(subdf['Hits'].max()))

    fig, ax = plt.subplots(1, 1, figsize=(8, len(subdf) * 0.2 + 1))
    sns.scatterplot(
        data=subdf, x='-log10(Adjusted P-value)', y='Term', hue='Hits (%)', size='Hits', ax=ax,
        legend='brief', hue_norm=hue_norm, palette="viridis", sizes=(10, 150), size_norm=size_norm,
        linewidth=0.5, edgecolor='black'
    )
    ax.set_title(f"{cat} - {libname} (N={len(subdf)})", fontsize=14)

    ax.set_xlabel('-log10(Adjusted P-value)', fontsize=12)
    ax.set_xlim(xlim)

    ax.set_ylabel('Term', fontsize=12)

    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    ax.axvline(1, color='red', linestyle='--', label='Adjusted p-value <= 0.1')

    # Plot the color bar
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(*hue_norm))
    cbar = plt.colorbar(sm, ax=ax, orientation='vertical', ticks=np.arange(hue_norm[0], hue_norm[1] + 1, 1))

    saveto = ld.plots.gsea / f"{cat}-{libname}.svg"
    saveto.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(saveto, bbox_inches='tight', pad_inches=0)

    # fig.show()
    plt.close(fig)
