import gzip
from collections import defaultdict

import pandas as pd
from biobit.core.loc import Strand

import utils.rnas
from assemblies import GRCh38
from stories import nextflow

# Always keep high-quality transcripts with these tags
ALLOWED_TAGS = ('MANE Select', 'MANE Plus Clinical')
ALLOWED_TSL = (1,)

# Otherwise, keep transcripts where all splicing junctions are supported by at least X reads
dfs = []
for fname in (nextflow.FULL_RESULTS / "star_salmon" / "log").glob("*.SJ.out.tab"):
    df = pd.read_csv(fname, sep="\t", header=None, names=[
        "seqid", "start", "end", "strand", "intron_motif", "annotated",
        "unique_reads", "multimapping_reads", "max_spliced_alignment_overhang"
    ])
    df['start'] -= 1  # Convert to 0-based
    dfs.append(df)
df = pd.concat(dfs, ignore_index=True)

df = df.groupby(['seqid', 'start', 'end', 'strand', 'intron_motif', 'annotated'], as_index=False).agg({
    'unique_reads': 'sum',
    'multimapping_reads': 'sum',
    'max_spliced_alignment_overhang': 'max'
})

# Remove poorly supported splice junctions
mask = (
        (df['max_spliced_alignment_overhang'] >= 20) &
        ((df['unique_reads'] >= 5) | (df['multimapping_reads'] >= 20)) &
        (df['strand'] != 0)  # Exclude undefined strands
)
df = df[mask].copy()

# Make the index of 'good' splice junctions
df['strand'] = df['strand'].map({1: Strand.Forward, 2: Strand.Reverse})
allowed_junctions = defaultdict(set)
for seqid, strand, start, end in zip(df['seqid'], df['strand'], df['start'], df['end']):
    allowed_junctions[seqid, strand].add((start, end))

# Make a list of allowed RNA and Gene IDs
allowed_genes, allowed_rnas = set(), set()
gencode = GRCh38.gencode.load()

summary = defaultdict(int)
for rna in gencode.rnas.values():
    # Skip RNAs that are not within the overall transcriptome universe
    if not utils.rnas.is_within_universe(rna):
        continue

    if any(tag in rna.attrs.tags for tag in ALLOWED_TAGS) or rna.attrs.TSL in ALLOWED_TSL:
        allowed_rnas.add(rna.ind)
        allowed_genes.add(rna.gene)
        summary['Allowed by tags'] += 1
        continue

    junctions = allowed_junctions[rna.loc.seqid, rna.loc.strand]
    allowed = True
    for nxt, prv in zip(rna.exons[1:], rna.exons[:-1]):
        assert prv.end < nxt.start, f"Exons {prv} and {nxt} are not in order for RNA {rna.ind}"
        if (prv.end, nxt.start) not in junctions:
            allowed = False
            break
    if allowed:
        # If the RNA is not allowed by tags, but all junctions are supported, then we should keep it
        allowed_rnas.add(rna.ind)
        allowed_genes.add(rna.gene)
        summary['Allowed by junctions'] += 1
    else:
        summary['Rejected by junctions'] += 1

print("RNA counts:")
for key, value in summary.items():
    print(f"\t{key}: {value:,}")

with (
    open(nextflow.RESOURCES / "filtered-annotation.gtf", "w") as saveto,
    gzip.open(GRCh38.gencode.gtf, 'rt') as gtf
):
    for line in gtf:
        if line.startswith("#"):
            saveto.write(line)
            continue

        _, source, feature, _, *_, attributes = line.split("\t")
        attributes = dict(x.split(maxsplit=1) for x in attributes.strip().split(";") if x)

        match feature:
            case "gene":
                assert "gene_id" in attributes and \
                       attributes['gene_id'][0] == '"' and \
                       attributes['gene_id'][-1] == '"'

                gid = attributes['gene_id'][1:-1]
                if gid in allowed_genes:
                    saveto.write(line)
            case _:
                assert "gene_id" in attributes and \
                       attributes['gene_id'][0] == '"' and \
                       attributes['gene_id'][-1] == '"'
                assert "transcript_id" in attributes and \
                       attributes['transcript_id'][0] == '"' and \
                       attributes['transcript_id'][-1] == '"'

                gid = attributes['gene_id'][1:-1]
                tid = attributes['transcript_id'][1:-1]
                if tid in allowed_rnas:
                    assert gid in allowed_genes
                    saveto.write(line)
