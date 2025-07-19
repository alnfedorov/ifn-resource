from collections import defaultdict

import pandas as pd

import ld
from assemblies import GRCh38

GENCODE = GRCh38.gencode.load()
GROUPS = ld.GROUPS.load()

# Aggregate transcript groups to gene level
gene_ids, ambiguous = {}, defaultdict(set)
for group in GROUPS.values():
    gids = {rna.gene for rna in group.members}
    if len(gids) > 1:
        for gid in gids:
            ambiguous[gid] |= gids
    assert group.ind not in gene_ids
    gene_ids[group.ind] = list(gids)

# Resolve ambiguous gene IDs by taking the full union of all ties for these genes
rescued = {}
for group, gids in gene_ids.items():
    if len(gids) == 1:
        continue
    elif any(gid in rescued for gid in gids):
        assert all(gid in rescued and rescued[gid] == rescued[gids[0]] for gid in gids)
        gene_ids[group] = rescued[gids[0]]
        continue

    queue, visited = list(gids), set()
    while queue:
        gid = queue.pop()
        if gid in visited:
            continue
        visited.add(gid)

        # If this gene is ambiguous, add all its ties to the queue
        if gid in ambiguous:
            for tied_gid in ambiguous[gid]:
                if tied_gid not in visited:
                    queue.append(tied_gid)

    solution = sorted(visited)
    for gid in solution:
        rescued[gid] = solution
    gene_ids[group] = solution

group, gids = zip(*gene_ids.items())
df = pd.DataFrame({"group": group, "gene_ids": gids}).set_index("group")

df['gene_names'] = df['gene_ids'].apply(lambda gids: '+'.join(sorted(GENCODE.genes[gid].attrs.name for gid in gids)))

# Rename selected genes
df['gene_names'] = df['gene_names'].replace({'STAT5A+STAT5B': 'STAT5A/B'})
df.to_csv(ld.GROUP2GENE, sep="\t", index=True)
