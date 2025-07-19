import pandas as pd
from biobit import io
from joblib import Parallel, delayed

import ld
from assemblies import GRCh38
from stories import cCRE
from utils import motifs, fasta


def screen(seqid: str, start: int, end: int, allmotifs: motifs.ZeroOrderMotifsCollection):
    reader = io.fasta.IndexedReader(GRCh38.fasta)
    forward = fasta.fetch(reader, seqid, (start, end), strand='+')
    revcomp = fasta.fetch(reader, seqid, (start, end), strand='-')

    # Calculate motif response score for each promoter as a maximum of its forward and reverse complement scores
    scores = motifs.score(forward.upper(), revcomp.upper(), allmotifs)

    record = {"seqid": seqid, "roi-norm-start": start, "roi-norm-end": end}
    for motif, score in zip(allmotifs.motifs, scores):
        record[motif.ind, motif.target] = score
    return record


# Load all motifs and calculate PWMs
with open(ld.jaspar.nonredundant) as f:
    jaspar = motifs.parse.jaspar(f)

pwms = tuple(pfm.to_pwm_hocomoco() for pfm in jaspar.motifs)
database = motifs.ZeroOrderMotifsCollection("ACGT", motifs=pwms)

# Load promoter sequences
sequences = cCRE.sequences()
sequences = sequences[sequences['roi-type'] == 'PLS']
regions = sequences[['seqid', 'roi-norm-start', 'roi-norm-end']].drop_duplicates()

# Calculate per-motif scores for each promoter
print(f"Calculating per-motif scores for {len(regions)} promoters...")
records = Parallel(n_jobs=-1, verbose=100, pre_dispatch='all', batch_size=1024)(
    delayed(screen)(seqid, start, end, database)
    for seqid, start, end in regions.itertuples(index=False, name=None)
)
df = pd.DataFrame(records)

# Save the results
ld.response.scores.parent.mkdir(parents=True, exist_ok=True)
df.to_pickle(ld.response.scores, protocol=-1)
