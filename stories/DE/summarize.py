from collections import defaultdict
from pathlib import Path

import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import false_discovery_control

import ld
from stories import terminus

GROUPS = terminus.GROUPS.load()


def postprocess(table: Path):
    df = pd.read_csv(table, index_col=0)

    # Drop untested records
    untested = df['log2FoldChange'].isna()
    df = df[~untested].copy()
    assert df.isna().sum().sum() == 0, f"Missing values in {df}"

    # Annotate group biotypes and names
    df['partition'] = df.index.map(lambda x: GROUPS[x].partition)

    # Group by partition and re-apply the FDR correction
    for partition, sub in df.groupby('partition'):
        sub = sub.drop(columns=['partition'])
        sub['padj'] = false_discovery_control(sub['pvalue'], method='bh')

        sub['category'] = 'Not significant'
        sub.loc[
            (sub['padj'] <= ld.thresholds.padj) & (sub['log2FoldChange'] >= ld.thresholds.log2fc), 'category'
        ] = 'Significant up'
        sub.loc[
            (sub['padj'] <= ld.thresholds.padj) & (sub['log2FoldChange'] <= -ld.thresholds.log2fc), 'category'
        ] = 'Significant down'

        saveto = ld.DESeq2.partitions / partition / table.name
        saveto.parent.mkdir(parents=True, exist_ok=True)
        sub.to_csv(saveto, index=True)


Parallel(n_jobs=-1, backend='threading')(
    delayed(postprocess)(table) for table in ld.DESeq2.tests.glob("*_vs_*.csv.gz")
)

# Postprocess the rlog/vsd tables
for table, saveto in (
        (ld.DESeq2.tests / "rld.csv.gz", ld.DESeq2.rld),
        (ld.DESeq2.tests / "vsd.csv.gz", ld.DESeq2.vsd),
):
    df = pd.read_csv(table, index_col=0)
    df = df.rename(columns={k: k.split('/')[-2] for k in df.columns})
    df['Name'] = df.index.map(lambda x: GROUPS[x].name)
    df.index.name = 'ID'
    df.to_csv(saveto, index=True)

# Make a summary of all statistical comparisons
summary = {"protein_coding": [], "non_coding": [], "artifacts": []}
for part in ld.DESeq2.partitions.iterdir():
    for path in part.glob("*.csv.gz"):
        df = pd.read_csv(path, index_col=0)

        target, control = path.name[:-7].split("_vs_")
        df = df.rename(columns={
            "category": (target, control, "category"),
            "log2FoldChange": (target, control, "log2FoldChange"),
            "padj": (target, control, "padj"), "baseMean": (target, control, "baseMean"),
            "lfcSE": (target, control, "lfcSE"), "pvalue": (target, control, "pvalue")
        }).reset_index(names='ID')
        df['partition'] = part.name
        df = df.set_index(['ID', 'partition'])
        summary[part.name].append(df)

# Merge all statistical comparisons into a single DataFrame
summary = {k: pd.concat(v, axis=1, join='outer') for k, v in summary.items()}
summary = pd.concat(summary.values(), axis=0, ignore_index=False)
assert summary.isna().sum().sum() == 0, "There are NaN values in the statistical results"
summary = summary.reset_index(level='partition', drop=False)

# Load & calculate median TPM per condition
TPMS = pd.read_csv(terminus.TPMS, sep='\t')

conditions = defaultdict(list)
for col in TPMS.select_dtypes(include=float).columns:
    conditions[col.split('-')[1]].append(col)

for condition, cols in conditions.items():
    TPMS[condition] = TPMS[cols].median(axis=1)
    TPMS = TPMS.drop(columns=cols)

TPMS = TPMS.set_index('ID')
TPMS = TPMS.rename(columns={k: (k, 'TPM') for k in TPMS.columns})

# Add median TPM to the summary
summary = summary.reset_index()
summary = pd.merge(summary, TPMS, left_on='ID', right_index=True, how='left')
assert summary.isna().sum().sum() == 0, "There are NaN values after merging TPM data"
summary = summary.set_index(['ID'])

# Assign biotypes/names to each transcript group
summary['name'] = summary.index.map(lambda x: GROUPS[x].name)
summary['type'] = summary.index.map(lambda x: GROUPS[x].type)

# Derive tags for each transcript group
alltags = []
for _, row in summary.iterrows():
    tags = set()

    if row["mock", "TPM"] >= ld.thresholds.min_background_tpm:
        tags.add(("mock", "Background"))

    for ifn in ld.IFNS:
        # Annotated IFN responses
        for control in ['mock']:
            if (ifn, control, 'category') in row and row[ifn, control, 'category'] != 'Not significant':
                tags.add((ifn, control, row[ifn, control, 'category']))
        # Annotated background RNAs
        if (
                (ifn, "mock", "Significant down") not in tags and
                (ifn, "mock", "Significant up") not in tags and
                (row[ifn, "TPM"] >= ld.thresholds.min_background_tpm)
        ):
            tags.add((ifn, "Background"))

    # Mark the total number of IFNs where the gene is up- or downregulated relative to the mock
    cnts = defaultdict(int)
    for tag in tags:
        if len(tag) == 3 and tag[1] == "mock":
            cnts[tag[2]] += 1
    if "Significant up" in cnts and "Significant down" in cnts:
        report = {x for x in tags if len(x) == 3 and x[1] == "mock"}
        print(f"Both up and down regulation detected for {row['name']}: {report}")

    for k, v in cnts.items():
        tags.add((f"IFN-{v}", k))

    # Mark background RNAs that are background in all conditions
    total = 0
    for condition in ["mock", *ld.IFNS]:
        if (condition, "Background") in tags:
            tags.remove((condition, "Background"))
            total += 1
    if total == len(ld.IFNS) + 1:
        tags.add("Background")

    alltags.append(tags)
summary['tags'] = alltags

summary.to_pickle(ld.DESeq2.summary, protocol=-1)
summary.to_csv(ld.DESeq2.summary.with_suffix('.csv.gz'), compression='gzip')
