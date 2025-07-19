from collections import defaultdict
from pathlib import Path

import pandas as pd
from joblib import Parallel, delayed

import ld
from stories import DE, cCRE, terminus

TAGS = pd.read_pickle(ld.TAGS)

# Sanity check - background RNAs must not be differentially expressed in any of the comparisons
for tags in TAGS['tags']:
    if 'Background' in tags:
        assert not any(x[2].startswith("Significant") for x in tags if isinstance(x, tuple) and x[1] == 'mock'), tags

# Transpose the tags and get a mapping from each tag to the genes
categories = defaultdict(set)
for ind, tags in TAGS[['ID', 'tags']].itertuples(index=False, name=None):
    for tag in tags:
        categories[tag].add(ind)

# Print stats
for k, v in categories.items():
    print(f"{k}: {len(v)}")

# Derive comparisons
# 1. To the background
# 2. Between ifns
# 3. Between ifn(K)s
# 4. scRNA-seq specific genes
comparisons = {}

for ifn in DE.IFNS:
    for direction in "up", "down":
        # Relative to the background
        comparisons[f"{ifn}_{direction}-vs-bckg"] = ((ifn, "mock", f"Significant {direction}"), "Background")
        # Relative to the other ifns
        for otherifn in DE.IFNS:
            if otherifn == ifn:
                continue
            comparisons[f"{ifn}_{direction}-vs-{otherifn}_{direction}"] = (
                (ifn, "mock", f"Significant {direction}"), (otherifn, "mock", f"Significant {direction}")
            )

for K in range(1, 6):
    for direction in "up", "down":
        # Relative to the background
        comparisons[f"IFN-{K}_{direction}-vs-bckg"] = ((f"IFN-{K}", f"Significant {direction}"), "Background")
        # Relative to the other ifns
        for otherK in range(1, K):
            assert otherK != K
            comparisons[f"IFN-{K}_{direction}-vs-IFN-{otherK}_{direction}"] = (
                (f"IFN-{K}", f"Significant {direction}"), (f"IFN-{otherK}", f"Significant {direction}")
            )

# scRNA-seq categories
comparisons["monocytes-and-lymphocytes-vs-monocytes"] = (("scRNA-seq", "Monocytes & Lymphocytes"),
                                                         ("scRNA-seq", "Monocyte-specific"))
comparisons["monocytes-and-lymphocytes-vs-lymphocytes"] = (("scRNA-seq", "Monocytes & Lymphocytes"),
                                                           ("scRNA-seq", "Lymphocyte-specific"))
comparisons["monocytes-vs-lymphocytes"] = (("scRNA-seq", "Monocyte-specific"), ("scRNA-seq", "Lymphocyte-specific"))
comparisons["lymphocytes-vs-monocytes"] = (("scRNA-seq", "Lymphocyte-specific"), ("scRNA-seq", "Monocyte-specific"))

for k, v in comparisons.items():
    print(f"{k} -> {v}")
    assert len(v) == 2, f"Invalid comparison {k}: {v}"
    for v in v:
        assert v in categories, f"Invalid category {v} in comparison {k}"

# Load the matching between cCREs and transcripts
sequences = cCRE.sequences()[['roi-type', 'seqid', 'roi-start', 'roi-end', 'sequence']]
overlaps = cCRE.overlaps()

sequences = sequences.merge(overlaps, on=['roi-type', 'seqid', 'roi-start', 'roi-end'], how='outer')
assert sequences.isna().sum().sum() == 0, "Some sequences are missing in the overlaps data"

# Match transcripts to groups
tx2group = pd.read_csv(terminus.TX2GROUP, sep="\t")[['transcript_id', 'group', 'name']]
sequences = sequences.merge(tx2group, left_on='Transcript ID', right_on='transcript_id', how='right')
sequences = sequences[['group', 'name', 'roi-type', 'sequence']].rename(columns={'group': 'ID'})
assert sequences.isna().sum().sum() == 0, "Some sequences are missing in the tx2group data"


def job(title: str, roi: str, target: set[str], background: set[str], sequences: pd.DataFrame, saveto: Path):
    # This is not necessary. E.g., upregulated RNAs could be upregulated by multiple IFNs
    # assert len(target & background) == 0, f"Target and background overlap in {name} ({roi})"

    # Subset the data
    mask = (sequences["ID"].isin(target | background)) & (sequences["roi-type"] == roi)
    df = sequences[mask].drop(columns='roi-type').copy()

    # Report missing sequences
    missing = (target | background) - set(df["ID"])
    if missing:
        print(f"{title} ({roi}): {len(missing)} missing sequences: {', '.join(sorted(missing)[:5])}")

    # Squish regions with the same sequence
    names, ids = defaultdict(set), defaultdict(set)
    for ind, name, sequence in df[['ID', 'name', 'sequence']].itertuples(index=False, name=None):
        names[sequence].add(name)
        ids[sequence].add(ind)
    names = {k: tuple(sorted(v)) for k, v in names.items()}
    ids = {k: tuple(sorted(v)) for k, v in ids.items()}

    df["name"] = df['sequence'].map(names)
    df["ID"] = df['sequence'].map(ids)
    df = df.drop_duplicates()

    before, unidirectional, contradictory = len(df), [], []
    for _, group in df.groupby("sequence", as_index=False):
        istarget = any(ind in target for x in group["ID"] for ind in x)
        isbackground = any(ind in background for x in group["ID"] for ind in x)
        assert istarget or isbackground, f"Neither target nor background in {title} ({roi}): {group['ID'].unique()}"
        if istarget and isbackground:
            contradictory.append(group)
        else:
            group["category"] = "target" if istarget else "background"
            unidirectional.append(group)

    if contradictory:
        contradictory = pd.concat(contradictory, ignore_index=True)
        names = {name for group in contradictory['name'] for name in group}
        print(
            f"{title} ({roi}): {len(contradictory)} contradictory sequences "
            f"({len(contradictory) / (len(contradictory) + len(unidirectional)):.2%}): {list(names)[:5]}..."
        )

    df = pd.concat(unidirectional, ignore_index=True)

    # Skip if there are < 25 sequences in each test category
    cnts = df["category"].value_counts()
    if any(cnts * ld.streme.hofract < 25) or len(cnts) != 2:
        return f"{title} ({roi}): skipped"

    # Shuffle the sequences
    df = df.sample(frac=1, random_state=39)

    # Save the result
    report = f"{title} ({roi})"
    saveto.mkdir(parents=True, exist_ok=True)
    for cat, fname in ("target", ld.streme.target), ("background", ld.streme.background):
        subset = df.loc[df["category"] == cat]
        report += f"\n\t{cat}: {len(subset)}"
        with open(saveto / fname, "w") as stream:
            for name, sequence in subset[["name", "sequence"]].itertuples(index=False, name=None):
                stream.write(f">{name}\n{sequence}\n")
    return report


# For whatever reason, python struggles to terminate the joblib threads.
reports = Parallel(n_jobs=-1, backend='sequential', pre_dispatch='all')(
    delayed(job)(name, roi, categories[target], categories[background], sequences, ld.streme.comparisons / roi / name)
    for name, (target, background) in comparisons.items() for roi in ["PLS"]  # , "pELS", "DNase-H3K4me3"]
)

for report in reports:
    print(report, "\n")
