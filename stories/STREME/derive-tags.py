from collections import defaultdict

import pandas as pd

import ld
from assemblies import GRCh38
from stories import terminus, DE

TX2GROUP = pd.read_csv(terminus.TX2GROUP, sep='\t').set_index('transcript_id')['group'].to_dict()
GENCODE = GRCh38.gencode.load()

# Load tags from DESeq2 and add single cell-based annotation
DESEQ2 = pd.read_pickle(DE.DESeq2.summary)
gname2ind = {gene.attrs.name: gene.ind for gene in GENCODE.genes.values()}
synonyms = {
    "CBWD1": "ZNG1A", "CBWD2": "ZNG1B", "DDX58": "RIGI", "C19ORF66": "SHFL", "H3F3B": "H3-3B",
    "MARCH1": "MARCHF1", "ODF3B": "CIMAP1B"
}

for category, genes in [
    ("Monocyte-specific", ld.single_cell.monocyte),
    ("Monocytes & Lymphocytes", ld.single_cell.monocyte_and_lymphocytes),
    ('Lymphocyte-specific', ld.single_cell.lymphocyte)
]:
    for gene in genes:
        gene = synonyms.get(gene.upper(), gene)
        if gene in {"AC116407.2", "AC004687.1"}:
            continue
        gene = GENCODE.genes[gname2ind[gene]]

        # scRNA-seq genes were upregulated in response to all type-I IFNs in a given cell type
        # In bulk RNA-seq, which is an average of all cell types, these genes are not necessarily
        # upregulated but are most likely reasonably expressed in at least one condition.
        groups = {TX2GROUP[x] for x in gene.transcripts if x in TX2GROUP}
        assert len(groups) >= 1, f"Gene {gene.attrs.name} has no groups: {groups}"

        groups = DESEQ2[DESEQ2.index.get_level_values('ID').isin(groups)]
        expressed = (groups[[(ifn, 'TPM') for ifn in DE.IFNS]] >= ld.single_cell.min_expression_tpm).any(axis=1)
        groups = groups[expressed].index.get_level_values('ID').unique()
        assert len(groups) >= 1, f"Gene {gene.attrs.name} has no groups: {groups}"

        # Add the category to the tags
        DESEQ2.loc[groups, 'tags'].apply(lambda x: x.add(('scRNA-seq', category)))

# Drop everything except the tags and records with no tags
DESEQ2 = DESEQ2[['tags']].reset_index(drop=False)
DESEQ2 = DESEQ2[DESEQ2['tags'].apply(lambda x: len(x) > 0)].copy()

# Print the number of records in each category
cnts = defaultdict(int)
for tags in DESEQ2['tags']:
    for tag in tags:
        cnts[tag] += 1
for tag, number in sorted(cnts.items(), key=lambda x: str(x[0])):
    print(f"{tag}: {number}")

# Save the tags
ld.TAGS.parent.mkdir(exist_ok=True, parents=True)
DESEQ2.to_pickle(ld.TAGS)
