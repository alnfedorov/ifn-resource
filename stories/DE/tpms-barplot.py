import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import ld
from stories import terminus

plt.rcParams['svg.fonttype'] = 'none'
ld.plots.barplot.mkdir(parents=True, exist_ok=True)

GROUP2GENE = pd.read_csv(terminus.GROUP2GENE, sep='\t', index_col=0)['gene_names'].to_dict()

# Load and aggregate to deduced gene names
TPM = pd.read_csv(terminus.TPMS, sep='\t')
TPM['Gene'] = TPM['ID'].map(GROUP2GENE)
TPM = TPM.drop(columns='ID')
TPM = TPM.groupby('Gene', as_index=False).sum()

TPM = TPM.melt(id_vars=['Gene'], var_name='Sample', value_name='TPM')
TPM['Donor'] = TPM['Sample'].apply(lambda x: x.split('+')[1].split('-')[0])
TPM['IFN'] = TPM['Sample'].apply(lambda x: x.split('+')[1].split('-')[1])
TPM = TPM.drop(columns=['Sample'])

########################################################################################################################
# Plot main genes in the IFN signaling pathway
########################################################################################################################
for genes, ylim in [
    (['IFNAR1', 'IFNAR2'], 100), (['JAK1', 'TYK2'], 200), (['USP18'], 10), (['IRF9'], 100),
    (['STAT1', 'STAT2', 'STAT3', 'STAT4', 'STAT5A/B', 'STAT6'], 400),
]:
    tpm = TPM.loc[TPM['Gene'].isin(genes) & (TPM['IFN'] == 'mock')].copy()
    assert len(tpm) == len(genes) * 3, f"Expected {len(genes)} genes, got {len(tpm) / 3}"

    fig, ax = plt.subplots(figsize=(len(genes) * 0.25 + 0.75, 4))
    sns.barplot(data=tpm, x='Gene', y='TPM', order=genes, ax=ax, errorbar=None, fill=False, color='black',
                linewidth=1.0)
    sns.stripplot(
        data=tpm, x='Gene', y='TPM', hue='Donor', order=genes, dodge=True, palette=ld.plots.palette.donor, ax=ax,
        legend=False, s=6
    )
    for i, gene in enumerate(genes):
        mean_tpm = tpm.loc[tpm['Gene'] == gene, 'TPM'].mean()
        ax.text(i, mean_tpm, f"{mean_tpm:.1f}", ha='center', va='bottom', fontsize=10, color='black', rotation=75)

    ax.spines[['top', 'right']].set_visible(False)
    ax.set_xlabel(None)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='right', fontsize=10)
    ax.set_yticks(range(0, ylim + 1, ylim // 5))

    name = '-'.join(genes).replace("/", "_")

    # fig.show()
    fig.savefig(ld.plots.barplot / f"{name}.svg", bbox_inches='tight', pad_inches=0)
    plt.close()

########################################################################################################################
# Plot selected RNAs
########################################################################################################################
ifnorder = ['mock', 'IFNa1', 'IFNa2a', 'IFNa10', 'IFNo', 'IFNb']
for RNA, ylim in [
    ('NRIR', 10), ('IRF6', 25), ('ZNF93', 25), ('ZFP57', 25), ('NRF1', 25),
]:
    tpm = TPM.loc[TPM['IFN'].isin(ifnorder) & (TPM['Gene'] == RNA)].copy()
    assert len(tpm) == len(ifnorder) * 3, f"Expected {len(ifnorder)} genes, got {len(tpm) / 3}"

    fig, ax = plt.subplots(figsize=(len(ifnorder) * 0.25 + 2.5, 4))
    sns.barplot(data=tpm, x='IFN', y='TPM', order=ifnorder, ax=ax, errorbar=None, fill=False, color='black',
                linewidth=1.0)
    sns.stripplot(
        data=tpm, x='IFN', y='TPM', hue='Donor', order=ifnorder, dodge=True, palette=ld.plots.palette.donor, ax=ax,
        legend=False, s=6
    )
    for i, condition in enumerate(ifnorder):
        mean_tpm = tpm.loc[tpm['IFN'] == condition, 'TPM'].mean()
        ax.text(i, mean_tpm, f"{mean_tpm:.1f}", ha='center', va='bottom', fontsize=10, color='black', rotation=75)

    ax.spines[['top', 'right']].set_visible(False)
    ax.set_xlabel(None)
    ax.set_title(f"{RNA} TPMs")
    ax.set_ylim(0, ylim)

    # fig.show()
    fig.savefig(ld.plots.barplot / f"{RNA}-TPMs-barplot.svg", bbox_inches='tight', pad_inches=0)
    plt.close()
