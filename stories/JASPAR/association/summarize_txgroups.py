import pandas as pd

import ld
from stories.DE import DESeq2, IFNS
from stories.JASPAR import scoring
from stories.terminus import TX2GROUP

TX2GROUP = pd.read_csv(TX2GROUP, sep='\t', index_col=0)['group'].to_dict()

# Match transcript responses to groups
responses = scoring.clusters().drop(columns=['seqid', 'roi-norm-start', 'roi-norm-end', 'is-reference'])
responses['ID'] = responses['Transcript ID'].apply(
    lambda alltids: {TX2GROUP[tid] for tid in alltids if tid in TX2GROUP}
)
responses = responses[responses['ID'].apply(len) > 0].drop(columns=['Transcript ID']).copy()
responses = responses.explode('ID')

responses = responses.groupby('ID').max()
assert responses.index.is_unique, "There are duplicate groups in the responses!"
responses = responses.rename(columns={
    'cluster_023': 'IRF6-like',
    'cluster_050': 'ZNF135/ZNF460',
    'cluster_041': 'ISRE-like',
    'cluster_025': 'GAS-like',
})
responses.columns = [('cluster', col) for col in responses.columns]

# Match groups to DESeq2 results
deseq2 = pd.read_pickle(DESeq2.summary).reset_index()
relevant_columns = [
    'ID', 'name', 'partition', 'tags',
    *[(ifn, "mock", stat) for ifn in IFNS for stat in ['log2FoldChange', 'lfcSE', 'pvalue', 'padj']],
    ('mock', 'TPM'), *[(ifn, 'TPM') for ifn in IFNS]
]
deseq2 = deseq2[relevant_columns].set_index('ID')

summary = pd.merge(responses, deseq2, how='right', left_index=True, right_index=True)
assert summary.isna().sum().sum() == 0 and len(summary) == len(deseq2)

# Save the summary
ld.TXGROUP_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
summary.to_pickle(ld.TXGROUP_SUMMARY, protocol=-1)
