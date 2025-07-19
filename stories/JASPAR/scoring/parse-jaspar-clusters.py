import pickle

import pandas as pd

import ld

jaspar = pd.read_csv(ld.jaspar.clusters, sep="\t")

parsed_ids = []
allids = set()
for ids in jaspar['id']:
    ids = [x.split('_')[3] for x in ids.split(',')]
    unique = set(ids)
    assert len(ids) == len(unique) and not allids & unique

    parsed_ids.append(unique)
    allids |= unique
jaspar['id'] = parsed_ids

jaspar['cluster'] = jaspar['cluster'].apply(lambda x: {
    'cluster_127': 'CTCF [1]',
    'cluster_128': 'CTCF [2]',
}.get(x, x))
jaspar['cluster'] = jaspar['cluster'].where(
    ~((jaspar['id'].apply(len) == 1) & (jaspar['cluster'].str.startswith('cluster_'))), jaspar['name']
)
assert jaspar['cluster'].nunique() == len(jaspar)

ld.jaspar.parsed_clusters.parent.mkdir(parents=True, exist_ok=True)
jaspar.to_pickle(ld.jaspar.parsed_clusters, protocol=pickle.HIGHEST_PROTOCOL)
