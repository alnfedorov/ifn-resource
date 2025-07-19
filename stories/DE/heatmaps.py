from collections import defaultdict
from math import log2
from pathlib import Path

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

import ld
from assemblies import GRCh38
from stories import terminus

plt.rcParams['svg.fonttype'] = 'none'

# Load required data
GENCODE = GRCh38.gencode.load()
SUMMARY = pd.read_pickle(ld.DESeq2.summary)
TX2GROUP = pd.read_csv(terminus.TX2GROUP, sep='\t', index_col=0)['group'].to_dict()

EXPRESSION = pd.read_csv(ld.DESeq2.rld, index_col=0)
EXPRESSION = EXPRESSION.rename(columns={k: k.split('+')[1] for k in EXPRESSION.columns if k != 'Name'})
ORDER = [
    'A-mock', 'B-mock', 'C-mock',
    'A-IFNa1', 'B-IFNa1', 'C-IFNa1',
    'A-IFNa2a', 'B-IFNa2a', 'C-IFNa2a',
    'A-IFNa10', 'B-IFNa10', 'C-IFNa10',
    'A-IFNo', 'B-IFNo', 'C-IFNo',
    'A-IFNb', 'B-IFNb', 'C-IFNb',
]


def expression(targets: list[str]) -> pd.DataFrame:
    expr = EXPRESSION[EXPRESSION.index.isin(targets)].copy().set_index('Name')
    expr = expr.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
    expr['total'] = (expr[ORDER[3:]].sum(axis=1).round(2) * 100).astype(int)
    expr = expr.reset_index()
    expr = expr.sort_values(by=['total', 'Name'], ascending=False).drop(columns=['total']).set_index("Name")
    return expr


def heatmap(expr: pd.DataFrame, title: str, label: bool = False, saveto: Path | None = None):
    expr = expr[ORDER]

    col_colors = expr.columns.to_frame(name='Sample')
    col_colors['Donor'] = col_colors.index.str.split('-').str[0].map(ld.plots.palette.donor)
    col_colors['IFN'] = col_colors.index.str.split('-').str[1].map(ld.plots.palette.ifn)
    col_colors = col_colors[['IFN', 'Donor']].astype(str)

    rowh = 1e-3 if not label else 1.25e-1
    height = 2 + rowh * len(expr.index)
    width = 6 if not label else 7
    clustermap = sns.clustermap(
        expr,
        cmap='RdYlBu_r',
        figsize=(width, height),
        col_colors=col_colors,
        colors_ratio=(0, 0.15 / height),
        row_cluster=False, col_cluster=False,
        vmin=-4, vmax=4,
        cbar_kws={'label': f'z-score', 'location': 'right', 'ticks': [-4, -2, 0, 2, 4]},
        cbar_pos=(0, 0.1 / height, 0.1 / width, 1.0 / height),
        xticklabels=True,
        yticklabels=label,
        dendrogram_ratio=(0.6 / width, 0.25 / height),
        annot=False,
        linewidths=0
    )
    clustermap.ax_heatmap.set_ylabel(None)
    clustermap.figure.suptitle(
        f"{title}\nN={len(expr)}", fontsize=16, fontweight='bold', va='top', ha='center', y=1.0, x=0.5
    )

    if saveto:
        saveto.parent.mkdir(parents=True, exist_ok=True)
        clustermap.savefig(saveto, bbox_inches='tight', pad_inches=0, dpi=1200)
    # clustermap.figure.show()
    plt.close(clustermap.figure)


########################################################################################################################
# All coding DE transcripts significant in at least one IFN experiment
########################################################################################################################
columns = [(ifn, 'mock', 'category') for ifn in ld.IFNS]
mask = (
        (SUMMARY[columns] != 'Not significant').any(axis=1) &
        (SUMMARY['partition'] == 'protein_coding')
)
targets = SUMMARY[mask].index.tolist()
heatmap(expression(targets), 'All coding DETs', saveto=ld.plots.heatmaps / 'all_coding_dets.png')

########################################################################################################################
# Same but with increased threshold for fold change values
########################################################################################################################
for fc in [5, 10]:
    thr = log2(fc)
    transcripts_up, transcripts_down = set(), set()
    for ifn in ld.IFNS:
        mask = (
                (SUMMARY['partition'] == 'protein_coding') &
                (SUMMARY[(ifn, 'mock', 'category')] != 'Not significant')
        )
        up = mask & (SUMMARY[(ifn, 'mock', 'log2FoldChange')] >= thr)
        transcripts_up.update(SUMMARY[up].index.tolist())

        down = mask & (SUMMARY[(ifn, 'mock', 'log2FoldChange')] <= -thr)
        transcripts_down.update(SUMMARY[down].index.tolist())

    print(f"Fold-change >= {fc} upregulated: {len(transcripts_up)}")
    print(f"Fold-change <= -{fc} downregulated: {len(transcripts_down - transcripts_up)}")

    transcripts = sorted(transcripts_up | transcripts_down)
    heatmap(expression(transcripts), f'All DETs with fold-change >= {fc}',
            saveto=ld.plots.heatmaps / f'all_coding_dets_fc_{fc}.svg')

########################################################################################################################
# IFN-specific coding DETs
########################################################################################################################
for ifn in ld.IFNS:
    exclude_tags = {
        (other_ifn, 'mock', f'Significant {direction}')
        for direction in ["up", "down"] for other_ifn in ld.IFNS
        if other_ifn != ifn
    }
    for direction in ['up', 'down']:
        mask = (
                (SUMMARY['partition'] == 'protein_coding') &
                SUMMARY['tags'].apply(
                    lambda tags: (ifn, 'mock', f"Significant {direction}") in tags and len(tags & exclude_tags) == 0
                )
        )
        transcripts = SUMMARY[mask].index.tolist()

        label = True if len(transcripts) <= 50 else False
        heatmap(expression(transcripts), f'{ifn}-specific coding DETs {direction}', label,
                saveto=ld.plots.heatmaps / f'{ifn}_specific_coding_dets_{direction}.svg')

########################################################################################################################
# All transcripts mapping to curated list of ISGs
########################################################################################################################
columns = [(ifn, 'mock', 'category') for ifn in ld.IFNS]
mask = (
        (SUMMARY[columns] == 'Significant up').any(axis=1) &
        (SUMMARY['partition'] == 'protein_coding')
)
updets = set(SUMMARY[mask].index.tolist())

gname2gid = {g.attrs.name: g.ind for g in GENCODE.genes.values()}

expected_no_expression = {
    "UPP2", "MKX", "C1S", "CTCFL", "CX3CL1", "IL17RB", "GPX2", "MAB21L2", "ANGPTL1", "ZNF385B", "AMPH", "ABLIM3",
    "G6PC1", "LINC01554", "NOS2", "CRP", "VEGFC", "CREB3L3", "FNDC4", "SAA1", "TBX3", "GBA3", "SERPINE1",
    "SLC1A1", "MT1H", "GEM", "NRN1", "ENPP1", "AHNAK2", "IGFBP2", "MT1G", "PRAME", "MT1M", "PPIAP10"
}

targets, noupdet = [], 0
for isg in ld.CURATED_ISGS:
    gene = GENCODE.genes[gname2gid[isg]]

    transcripts = {TX2GROUP[tid] for tid in gene.transcripts if tid in TX2GROUP}
    assert len(transcripts) > 0, f"ISG {isg} ({gene.ind}) has no groups in TX2GROUP"
    transcripts = SUMMARY[SUMMARY.index.isin(transcripts)]
    if len(transcripts) == 0:
        if isg in expected_no_expression:
            continue
        print(f"ISG {isg} ({gene.ind}) has no transcripts in DESeq2 results")
        continue
    transcripts = set(transcripts.index)
    if len(transcripts & updets) == 0:
        noupdet += 1
    targets.extend(transcripts)

print(f"ISGs with no upregulated transcripts: {noupdet} out of {len(ld.CURATED_ISGS)}")

exp = expression(targets)
heatmap(exp, 'ISG-associated coding DETs', label=False, saveto=ld.plots.heatmaps / 'isg_coding_dets.svg')

exp = exp.head(50)  # Select top 50 transcripts by the overall Z-score
heatmap(exp, 'Top 50 ISG-associated coding DETs', label=True, saveto=ld.plots.heatmaps / 'top_50_isg_coding_dets.svg')

########################################################################################################################
# All non-coding DE transcripts significant in at least one IFN experiment
########################################################################################################################
columns = [(ifn, 'mock', 'category') for ifn in ld.IFNS]
mask = (
        (SUMMARY[columns] != 'Not significant').any(axis=1) &
        (SUMMARY['partition'] == 'non_coding')
)
targets = SUMMARY[mask].index.tolist()
heatmap(expression(targets), 'All non-coding DETs', saveto=ld.plots.heatmaps / 'all_noncoding_dets.svg')

########################################################################################################################
# All non-coding DE transcripts up/down regulated by all IFNs
########################################################################################################################
columns = [(ifn, 'mock', 'category') for ifn in ld.IFNS]
mask = (
        (SUMMARY[columns] != 'Not significant').all(axis=1) &
        (SUMMARY['partition'] == 'non_coding')
)
targets = SUMMARY[mask].index.tolist()
heatmap(expression(targets), 'Core non-coding DETs', label=True, saveto=ld.plots.heatmaps / 'core_noncoding_dets.svg')

########################################################################################################################
# All MIR-hosting genes
########################################################################################################################
columns = [(ifn, 'mock', 'category') for ifn in ld.IFNS]
mask = (
        (SUMMARY[columns] != 'Not significant').any(axis=1) &
        (SUMMARY['partition'] == 'non_coding') &
        (SUMMARY['name'].str.contains('MIR'))
)
mirhg = SUMMARY[mask].index.tolist()
heatmap(expression(mirhg), 'MIR-hosting genes', label=True, saveto=ld.plots.heatmaps / 'MIR-hosting-genes.svg')

########################################################################################################################
# DT-RNAs and AS-RNAs
########################################################################################################################
columns = [(ifn, 'mock', 'category') for ifn in ld.IFNS]

for postfix, scoring in [
    ('-DT', lambda x: x.loc[x['is-target'], 'score'].mean() + x.loc[~x['is-target'], 'score'].mean()),
    ('-AS', lambda x: x.loc[x['is-target'], 'score'].mean() - x.loc[~x['is-target'], 'score'].mean())
]:
    mask = (
            SUMMARY['name'].str.contains(postfix) &
            (SUMMARY[columns] != 'Not significant').any(axis=1) &
            (SUMMARY['partition'] == 'non_coding')
    )
    targets = SUMMARY[mask].index.tolist()

    names = SUMMARY.loc[mask, 'name'].str.split('|').apply(
        lambda x: [pref.split('-')[0] for pref in x]
    ).explode().unique()
    neighbors = (
            SUMMARY['name'].apply(lambda x: any(name in x for name in names)) &
            (SUMMARY[columns] != 'Not significant').any(axis=1) &
            (SUMMARY['partition'] == 'protein_coding')
    )
    RNA = SUMMARY.loc[neighbors | mask, ['name']].copy()

    groups = defaultdict(list)
    for groupid, name in RNA['name'].items():
        for x in name.split('|'):
            prefix = x.split('-')[0].split('[')[0]
            groups[prefix].append(groupid)

    allids = set.union(*[set(v) for v in groups.values()])
    assert set(targets).issubset(allids)
    allids = sorted(allids)

    group2rna = {group: k for k, v in groups.items() for group in v}
    rnas = SUMMARY.loc[allids, ['name']].copy()
    rnas['group'] = rnas.index.map(group2rna)
    rnas = rnas.reset_index().set_index('name')

    exp = expression(list(allids))
    exp['score'] = exp[ORDER[3:]].mean(axis=1)
    exp = exp.join(rnas).reset_index()

    exp['is-target'] = exp['ID'].isin(targets)

    # Drop unrelated without any targets
    with_target = exp[['is-target', 'group']].groupby('group')['is-target'].any().to_dict()
    exp = exp[exp['group'].map(with_target)]

    # Drop singleton targets
    singletons = exp[['is-target', 'group']].groupby('group').all()['is-target'].to_dict()
    notsingle = ~exp['group'].map(singletons)
    exp = exp[notsingle].copy()
    print(f"Keeping {exp['is-target'].sum()} / {len(targets)} {postfix} target RNA groups")

    scores = exp[['is-target', 'score', 'group']].groupby('group').apply(scoring).to_dict()
    exp['score'] = exp['group'].map(scores)

    exp = exp.sort_values(by=['score', 'is-target', 'Name']).drop(
        columns=['score', 'ID', 'group', 'is-target']
    ).set_index('Name')
    heatmap(exp, f'{postfix} transcripts', label=True, saveto=ld.plots.heatmaps / f'{postfix[1:]}_transcripts.svg')
