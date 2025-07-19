import pandas as pd

import ld
from stories import terminus

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

ld.SUPPLEMENTARY_TABLES.mkdir(parents=True, exist_ok=True)

SUMMARY = pd.read_pickle(ld.DESeq2.summary).reset_index()
TX2GROUP = pd.read_csv(terminus.TX2GROUP, sep='\t')[['group', 'name']].drop_duplicates().set_index('group')

# Estimated expression (TPM) for each transcript group.
TPMS = pd.read_csv(terminus.TPMS, sep='\t')
TPMS = TPMS.join(TX2GROUP, on='ID', how='outer')

group2gene = pd.read_csv(terminus.GROUP2GENE, sep='\t', index_col='group')
group2gene['gene_ids'] = group2gene['gene_ids'].apply(lambda x: x[1:-1].replace("'", "").split(', '))
group2gene['gene_ids'] = group2gene['gene_ids'].apply(lambda x: "|".join(x))
group2gene['gene_names'] = group2gene['gene_names'].apply(lambda x: x.replace("+", "|"))
TPMS = TPMS.join(group2gene, on='ID', how='outer')

columns = ['ID', 'name', 'gene_ids', 'gene_names']
columns += [col for col in TPMS.columns if col not in columns]
TPMS = TPMS[columns].sort_values(by='ID')

# Reorder columns and save
columns = [
    "ID", "name", "gene_ids", "gene_names",
    "ERX10466138+A-mock", "ERX10466137+B-mock", "ERX10476454+C-mock",
    "ERX10487730+A-IFNa1", "ERX10476248+B-IFNa1", "ERX10476470+C-IFNa1",
    "ERX10475515+A-IFNa2a", "ERX10479637+B-IFNa2a", "ERX10476466+C-IFNa2a",
    "ERX10475516+A-IFNa10", "ERX10475581+B-IFNa10", "ERX10476334+C-IFNa10",
    "ERX10487731+A-IFNo", "ERX10476340+B-IFNo", "ERX10476471+C-IFNo",
    "ERX10475529+A-IFNb", "ERX10475594+B-IFNb", "ERX10476465+C-IFNb",
]
TPMS = TPMS[columns].rename(columns={
    "ID": "Group ID", "name": "Group name", "gene_ids": "Gene IDs", "gene_names": "Gene names"
})
TPMS.to_csv(ld.SUPPLEMENTARY_TABLES / "Supplementary Table 9: Transcript TPMs.tsv", sep='\t', index=False)

# Estimated expression (TPM) for each gene.
assert (TPMS[['Gene IDs', 'Gene names']].groupby('Gene IDs').nunique()['Gene names'] == 1).all()
TPMS = TPMS.drop(columns=['Group ID', 'Group name'])
TPMS = TPMS.groupby(['Gene IDs', 'Gene names'], as_index=False).sum()

TPMS.to_csv(ld.SUPPLEMENTARY_TABLES / "Supplementary Table 10: Gene TPMs.tsv", sep='\t', index=False)

# Summary of differentially expressed transcripts mapping to protein-coding / non-coding genes.
columns = ['ID', 'name', *[(ifn, 'mock', col) for ifn in ld.IFNS for col in ['log2FoldChange', 'pvalue', 'padj']]]
rename = {(ifn, 'mock', col): f"{ifn}: {col}" for ifn in ld.IFNS for col in ['log2FoldChange', 'pvalue', 'padj']}
for partition, table_id in [('protein_coding', '11'), ('non_coding', '12')]:
    summary = SUMMARY.loc[SUMMARY['partition'] == partition, columns].copy()
    summary = summary.sort_values(by='ID').rename(columns={'ID': 'Group ID', 'name': 'Group name'} | rename)
    summary.to_csv(
        ld.SUPPLEMENTARY_TABLES / f"Supplementary Table {table_id}: Differentially expressed {partition} groups.tsv",
        sep='\t', index=False
    )
