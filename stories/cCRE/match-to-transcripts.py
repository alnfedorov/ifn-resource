from collections import defaultdict

import pandas as pd
from biobit.collections import interval_tree
from biobit.core.loc import Orientation, Interval
from biobit.io.bed import Bed6
from joblib import Parallel, delayed

import ld
import utils
from assemblies import GRCh38

GENCODE = GRCh38.gencode.load()

tss_index = defaultdict(interval_tree.BitsBuilder)
transcripts_index = defaultdict(interval_tree.BitsBuilder)
all_transcripts = []

for rna in GENCODE.rnas.values():
    if not utils.rnas.is_within_universe(rna):
        continue

    record = {
        "Transcript ID": rna.ind, "seqid": rna.loc.seqid, "strand": rna.loc.strand,
        "start": rna.loc.start, "end": rna.loc.end,
    }
    tss = rna.loc.start if rna.loc.strand == "+" else rna.loc.end
    tss_index[rna.loc.seqid].add((tss - 1, tss + 1), record)
    transcripts_index[rna.loc.seqid].add((rna.loc.start, rna.loc.end), record)

    all_transcripts.append(rna)

tss_index = interval_tree.Forest(tss_index)
transcripts_index = interval_tree.Forest(transcripts_index)

# Load all cCREs from ENCODE and select only those that are relevant for the analysis
df = pd.read_csv(GRCh38.cCRE, sep='\t', names=["seqid", "start", "end", "ind", "name", "ccre"])
df['ccre'] = df['ccre'].apply(lambda x: set(x.split(',')))

ROIs = {}
for roi, saveto, color in [
    ("PLS", ld.cCRE.PLS, "255,0,0"),
    ("pELS", ld.cCRE.pELS, "0,255,0"),
    ("DNase-H3K4me3", ld.cCRE.DNase_H3K4me3, "0,0,255")
]:
    ccre = df[df['ccre'].apply(lambda x: roi in x)]
    ccre = ccre[['seqid', 'start', 'end', 'name']].copy()

    # Save as BED
    bed = [
        Bed6(seqid, (start, end), name, 0, Orientation.Dual)
        for seqid, start, end, name in ccre.itertuples(index=False)
    ]
    utils.bed.tbindex(bed, Bed6, saveto)

    # Save the ROIs partitioned by contig
    for seqid, partition in ccre.groupby('seqid'):
        ROIs[seqid, roi] = partition


# Map ROIs to transcripts
def job(roi: str, seqid: str, offset: int, coordinates: pd.DataFrame, index: interval_tree.Bits[dict] | None):
    if not index:
        print(f"No index for {roi} on {seqid}, skipping")
        return []
    assert (coordinates['seqid'] == seqid).all(), f"{coordinates['seqid'].unique()} != {seqid}"

    queries = [Interval(start - offset, end + offset) for start, end in zip(coordinates['start'], coordinates['end'])]
    overlaps = index.batch_intersect_intervals(queries)

    results = []
    for (_, start, end, name), (_, overlap) in zip(coordinates.itertuples(index=False), overlaps):
        for rna in overlap:
            assert rna["seqid"] == seqid, f"{rna['seqid']} != {seqid}"

            results.append({
                "Transcript ID": rna["Transcript ID"],
                "seqid": seqid,
                "rna-start": rna["start"], "rna-end": rna["end"], "rna-strand": rna["strand"],
                "roi-type": roi, "roi-start": start, "roi-end": end, "roi-name": name
            })
    return results


indices = {"PLS": tss_index, "pELS": tss_index, "DNase-H3K4me3": transcripts_index}
results = Parallel(n_jobs=-1, backend='threading', pre_dispatch='all')(
    delayed(job)(roi, seqid, ld.cCRE.overlaps.max_distances[roi], coordinates, indices[roi].get(seqid))
    for (seqid, roi), coordinates in ROIs.items()
)
results = pd.concat([pd.DataFrame(x) for x in results]).drop_duplicates()
results["imputed"] = False

# Impute promoters for transcripts without matched ENCODE PLS
records = []
has_pls = set(results.loc[results['roi-type'] == "PLS", "Transcript ID"])
for rna in all_transcripts:
    if rna.ind in has_pls:
        continue

    if rna.loc.strand == "+":
        start = rna.loc.start - ld.sequences.lengths["PLS"] // 2
    else:
        assert rna.loc.strand == '-'
        start = rna.loc.end - ld.sequences.lengths["PLS"] // 2
    end = start + ld.sequences.lengths["PLS"]
    assert end - start == ld.sequences.lengths["PLS"] and start < end

    records.append({
        "Transcript ID": rna.ind,
        "seqid": rna.loc.seqid, "rna-start": rna.loc.start, "rna-end": rna.loc.end, "rna-strand": rna.loc.strand,
        "roi-type": "PLS", "roi-start": start, "roi-end": end, "roi-name": f"Imputed[{rna.ind}]", "imputed": True
    })
print(f"Total parsed transcripts: {len(all_transcripts):,}")
print(f"\tImputed PLS: {len(records):,}")
print(f"\tENCODE PLS: {len(has_pls):,}")
imputed = pd.DataFrame(records)

# Combine the data
assert (results.columns == imputed.columns).all(), set(results.columns).symmetric_difference(imputed.columns)
results = pd.concat([results, imputed])

# Save as pkl
ld.cCRE.overlaps.saveto.mkdir(exist_ok=True, parents=True)
results.to_pickle(ld.cCRE.overlaps.pkl)

# Save as BED
bed = []
for tid, roi, seqid, start, end, imputed in results[
    ["Transcript ID", "roi-type", "seqid", "roi-start", "roi-end", "imputed"]
].itertuples(index=False, name=None):
    imputed = "Impute" if imputed else "cCRE"
    bed.append(Bed6(
        seqid, (start, end), f"{roi}-{tid}[{imputed}]", 0, Orientation.Dual,
    ))
utils.bed.tbindex(bed, Bed6, ld.cCRE.overlaps.bed)
