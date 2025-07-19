import gzip

import utils
from assemblies import GRCh38
from stories import nextflow

allowed_genes, allowed_rnas = set(), set()
gencode = GRCh38.gencode.load()
for rna in gencode.rnas.values():
    if utils.rnas.is_within_universe(rna):
        allowed_rnas.add(rna.ind)
        allowed_genes.add(rna.gene)

with (
    open(nextflow.RESOURCES / "full-annotation.gtf", "w") as saveto,
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
