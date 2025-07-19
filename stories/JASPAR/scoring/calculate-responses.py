from collections import defaultdict

import pandas as pd

import ld
import utils
from assemblies import GRCh38
from stories import cCRE

GENCODE = GRCh38.gencode.load()

# Load mapping between ROIs and transcripts
overlaps = cCRE.overlaps()
overlaps = overlaps.loc[overlaps['roi-type'] == 'PLS', ['Transcript ID', 'imputed', 'seqid', 'roi-start', 'roi-end']]

sequences = cCRE.sequences()
sequences = sequences[sequences['roi-type'] == 'PLS']

overlaps = overlaps.merge(sequences, on=['seqid', 'roi-start', 'roi-end'], how='outer')
overlaps = overlaps[['Transcript ID', 'imputed', 'seqid', 'roi-norm-start', 'roi-norm-end']]
assert overlaps.isna().sum().sum() == 0, "There are NaN values in the overlaps!"

# Match ROI coordinates with Transcript IDs and imputed status
roi2transcripts, roi2imputed = defaultdict(set), {}
for tid, imputed, seqid, start, end in overlaps.itertuples(index=False, name=None):
    roi2transcripts[seqid, start, end].add(tid)
    if (seqid, start, end) not in roi2imputed:
        roi2imputed[seqid, start, end] = imputed
    assert roi2imputed[seqid, start, end] == imputed, (seqid, start, end, imputed)
roi2transcripts = {k: sorted(v) for k, v in roi2transcripts.items()}

# Load pre-calculated scores for all promoters and match them with Transcript IDs / Imputed status
responses = pd.read_pickle(ld.response.scores)
responses['Transcript ID'] = [
    roi2transcripts[k] for k in zip(responses['seqid'], responses['roi-norm-start'], responses['roi-norm-end'])
]
assert responses['Transcript ID'].apply(len).min() > 0, "Some ROIs do not have any transcripts!"

responses['imputed'] = [
    roi2imputed[k] for k in zip(responses['seqid'], responses['roi-norm-start'], responses['roi-norm-end'])
]
assert responses['imputed'].isna().sum() == 0, "There are NaN values in the imputed scores!"

# Find a set of reference promoters - not imputed and matching with high-quality transcripts
references = {k for k, v in GENCODE.rnas.items() if utils.rnas.is_high_quality(v)}

responses['is-reference'] = (
        (~responses['imputed']) &
        responses['Transcript ID'].apply(lambda x: any(tid in references for tid in x))
)

# Z-score the response scores for each promoter using reference promoters
motifs = [col for col in responses.columns if isinstance(col, tuple) and len(col) == 2]
references = responses.loc[responses['is-reference'], motifs]
mean, std = references.mean(), references.std()

for col in motifs:
    responses[col] = (responses[col] - mean[col]) / std[col]

# Clean up the DataFrame
responses = responses[['Transcript ID', 'seqid', 'roi-norm-start', 'roi-norm-end', 'is-reference', *motifs]].copy()
responses.to_pickle(ld.response.per_motif, protocol=-1)

# Aggregate to the level of clusters and re-standardize
jaspar = pd.read_pickle(ld.jaspar.parsed_clusters)

# Drop TF names
renaming = {x: x[0] for x in motifs}
assert len(set(renaming.values())) == len(renaming)
responses = responses.rename(columns=renaming)

for cluster, ids in jaspar[['cluster', 'id']].itertuples(index=False, name=None):
    # Check that all motifs in the cluster are present in the responses
    assert all(id in responses.columns for id in ids), \
        f"Cluster {cluster} has missing motifs: {set(ids) - set(responses.columns)}"

    responses[cluster] = responses[list(ids)].sum(axis=1)

    # Re-standardize the cluster scores using the mean and std of the reference promoters
    references = responses.loc[responses['is-reference'], cluster]
    mean, std = references.mean(), references.std()

    responses[cluster] = (responses[cluster] - mean) / std
responses = responses.drop(columns=list(renaming.values()))
responses.to_pickle(ld.response.per_cluster, protocol=-1)
