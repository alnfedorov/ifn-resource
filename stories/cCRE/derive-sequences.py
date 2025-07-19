import pandas as pd
from biobit import io

import ld
import utils.fasta
from assemblies import GRCh38

# Load all overlaps of cCREs with transcripts
df = pd.read_pickle(ld.cCRE.overlaps.pkl)
cCREs = df[["roi-type", "seqid", "roi-start", "roi-end"]].drop_duplicates()

# Normalize the length of all elements
mask, nstart, nend = [], [], []
sizes = GRCh38.seqid.sizes()
for etype, seqid, start, end in cCREs[["roi-type", "seqid", "roi-start", "roi-end"]].itertuples(index=False):
    fitto = ld.sequences.lengths[etype]
    length = end - start

    if length > fitto:
        mask.append(False)
    elif length < fitto:
        padleft = (fitto - length) // 2
        padright = fitto - length - padleft

        start = max(0, start - padleft)
        end = min(sizes[seqid], end + padright)

        if end - start != fitto:
            mask.append(False)
        else:
            mask.append(True)
            nstart.append(start)
            nend.append(end)
    else:
        mask.append(True)
        nstart.append(start)
        nend.append(end)

print(f"Dropped {len(cCREs) - sum(mask)} cCREs that could not be normalized to the desired length")
cCREs.loc[mask, "roi-norm-start"] = nstart
cCREs.loc[mask, "roi-norm-end"] = nend
cCREs = cCREs[mask].astype({"roi-norm-start": int, "roi-norm-end": int}).copy()
assert all(cCREs["roi-norm-end"] - cCREs["roi-norm-start"] == cCREs["roi-type"].map(ld.sequences.lengths))

# Fetch all sequences
reader = io.fasta.IndexedReader(GRCh38.fasta)
cCREs["sequence"] = [
    utils.fasta.fetch(reader, seqid, (start, end))
    for seqid, start, end in cCREs[["seqid", "roi-norm-start", "roi-norm-end"]].itertuples(index=False)
]

ld.sequences.saveto.parent.mkdir(exist_ok=True, parents=True)
cCREs.to_pickle(ld.sequences.saveto)
