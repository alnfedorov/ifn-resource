from functools import reduce

import pandas as pd

import ld

# Match terminus names to actual RNA groups
groups = ld.GROUPS.load()
tid2group = {}  # Map Ensembl IDs to RNA groups
for group in groups.values():
    for member in group.members:
        assert member.ind not in tid2group
        tid2group[member.ind] = group

terminus2group = {}  # Map terminus group names to RNA groups
for tid, termname in ld.terminus.groups.load().items():
    group = tid2group[tid]
    if termname in terminus2group:
        assert terminus2group[termname] == group
    else:
        terminus2group[termname] = group

# Load terminus-collapsed TPMs and NumReads and remap IDs
tpms, reads = [], []
for sample in ld.terminus.output.iterdir():
    quant = sample / "quant.sf"
    assert quant.exists()
    df = pd.read_csv(quant, sep='\t') \
        .drop(columns=['Length', 'EffectiveLength'])

    # Remap terminus ids to new group IDs
    column = []
    for termname in df['Name']:
        column.append(terminus2group[termname].ind)
    df['ID'] = column
    assert df['ID'].is_unique
    df = df.drop(columns=['Name'])

    tpms.append(df[['ID', 'TPM']].copy().rename(columns={"TPM": sample.name}))
    reads.append(df[['ID', 'NumReads']].copy().rename(columns={"NumReads": sample.name}))

# Merge all samples into a single table
reads = reduce(lambda l, r: pd.merge(l, r, on=['ID'], how='outer'), reads)
assert reads.isna().sum().sum() == 0

tpms = reduce(lambda l, r: pd.merge(l, r, on=['ID'], how='outer'), tpms)
assert tpms.isna().sum().sum() == 0

# Save the tables
order = ['ID'] + sorted([col for col in tpms.columns if col != 'ID'])
tpms[order].to_csv(ld.TPMS, sep='\t', index=False)
reads[order].to_csv(ld.READS, sep='\t', index=False)
