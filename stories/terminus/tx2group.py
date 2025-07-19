from collections import defaultdict

import pandas as pd
from biobit.toolkit.annotome import Annotome
from biobit.toolkit.annotome.transcriptome import RNA

import ld
import utils
from assemblies import GRCh38


def resolve_biotypes(rnas: list[RNA[GRCh38.gencode.AttrRNA]]) -> GRCh38.gencode.RNAType | None:
    if len(rnas) == 0:
        return None
    elif len(rnas) == 1:
        return rnas[0].attrs.type

    biotypes = {rna.attrs.type for rna in rnas}
    if "protein_coding" in biotypes:
        # Pretend that the main RNA biotype is protein_coding if it is present
        return "protein_coding"
    elif biotypes.issubset({
        "transcribed_unprocessed_pseudogene", "transcribed_processed_pseudogene", "translated_processed_pseudogene",
        "processed_pseudogene", "unprocessed_pseudogene", "transcribed_unitary_pseudogene", "lncRNA"
    }):
        return "pseudogene"
    elif biotypes.issubset({
        'retained_intron', 'nonsense_mediated_decay', 'non_stop_decay', 'protein_coding_CDS_not_defined'
    }):
        return "protein_coding_CDS_not_defined"
    elif biotypes == {'lncRNA', 'retained_intron'}:
        return 'retained_intron'
    elif biotypes == {'TEC', 'lncRNA'}:
        return 'lncRNA'
    elif biotypes == {'artifact', 'processed_pseudogene'}:
        return "pseudogene"
    else:
        return None


def compose_name(gencode: Annotome, rnas: list[RNA[GRCh38.gencode.AttrRNA]]) -> str:
    # Name should be short and take the form of Gene[tind-1/tind-2/...]|Gene[tind-1/tind-2/...]|etc
    per_gene = defaultdict(list)
    for rna in rnas:
        per_gene[rna.gene].append(rna.attrs.name)

    names = []
    for gene, rna_names in per_gene.items():
        gname = gencode.genes[gene].attrs.name
        if gname.startswith("ENSG"):
            # Gene name is an Ensembl ID - simply concatenate transcript names
            names.append("/".join(sorted(rna_names)))
            continue
        # Gene name is a proper name - use it as a prefix
        postfixes = []
        for rna_name in rna_names:
            if rna_name.startswith(gname):
                postfixes.append(rna_name[len(gname):].lstrip("-"))
            else:
                postfixes.append(rna_name)
        postfixes = sorted(set(postfixes))
        if len(postfixes) >= 5:
            names.append(f"{gname}[N={len(postfixes)}]")
        else:
            names.append(f"{gname}[" + "/".join(postfixes) + "]")
    return "|".join(sorted(names))


# Parse all predicted groups (must be identical across all outputs)
lines = set()
for file in ld.terminus.output.glob("*/clusters.txt"):
    with open(file, 'r') as stream:
        lines.add(stream.read())

assert len(lines) == 1

# Parse collapses identified by the Terminus
groups = {}
for line in lines.pop().split("\n"):
    if not line:
        continue
    name, *transcripts = line.split(",")
    assert name.startswith("NewTr")
    for tid in transcripts:
        assert tid not in groups
        groups[tid] = name

# Parse all RNAs that were used to create groups (e.g., belong to the filtered RNA universe)
RNAs = set()
for file in ld.SALMON.glob("*/quant.sf"):
    salmon = pd.read_csv(file, sep="\t", usecols=["Name"])
    RNAs.update(salmon["Name"])

# Fill-in singleton groups and make backward mapping to genes
gencode = GRCh38.gencode.load()
all_rnas_per_gene = defaultdict(set)
for rna in RNAs:
    rna = gencode.rnas[rna]
    all_rnas_per_gene[rna.gene].add(rna.ind)
    if rna.ind not in groups:  # If not in groups, it is a singleton = a group by itself
        groups[rna.ind] = rna.ind

# Revmap groups to transcript IDs
revmap = defaultdict(list)
for tid, name in groups.items():
    revmap[name].append(tid)

sizes = defaultdict(int)
result = {}
for tids in revmap.values():
    sizes[len(tids)] += 1
    rnas = [gencode.rnas[tid] for tid in tids]

    # Singletons
    if len(rnas) == 1:
        rna = rnas.pop()
        result[rna.ind] = utils.rnas.RNAGroup(rna.ind, rna.attrs.name, rna.attrs.type, (rna,))
        continue

    ind = "|".join(sorted(tids))
    biotypes = {x.attrs.type for x in rnas}
    name = compose_name(gencode, rnas)
    assert ind not in result

    # Are all these RNAs coming from the same gene?
    genes = {rna.gene for rna in rnas}
    if len(genes) == 1:
        gene = gencode.genes[genes.pop()]

        # Collapse the whole gene
        if len(tids) == len(all_rnas_per_gene[gene.ind]):
            assert gene.ind not in result
            result[gene.ind] = utils.rnas.RNAGroup(gene.ind, name, gene.attrs.type, tuple(rnas))
            continue

        # Collapse if the gene biotype is among the transcript biotypes
        if gene.attrs.type in biotypes:
            result[ind] = utils.rnas.RNAGroup(ind, name, gene.attrs.type, tuple(rnas))
            continue

    # Manually resolved cases
    match name:
        case "7SK[.1-201]|RN7SK[201]":
            assert 'snRNA' in biotypes
            result[ind] = utils.rnas.RNAGroup(ind, '7SK', 'snRNA', tuple(rnas))
            continue
        case "ENST00000797550|RMRP-201":
            assert 'ribozyme' in biotypes
            result[ind] = utils.rnas.RNAGroup(ind, name, 'ribozyme', tuple(rnas))
            continue

    if len(biotypes) == 1:
        result[ind] = utils.rnas.RNAGroup(ind, name, biotypes.pop(), tuple(rnas))
        continue

    # Fall-back to the resolution based on the most well-annotated transcript
    # MANE Select, MANE Plus Clinical
    welldefined = [rna for rna in rnas if rna.attrs.tags & {'MANE Select', 'MANE Plus Clinical'}]
    resolved = resolve_biotypes(welldefined)
    if resolved is not None:
        result[ind] = utils.rnas.RNAGroup(ind, name, resolved, tuple(rnas))
        continue

    # TSL - transcript support level
    welldefined = sorted([rna for rna in rnas if rna.attrs.TSL in {1, 2}], key=lambda x: x.attrs.TSL)
    welldefined = [rna for rna in welldefined if rna.attrs.TSL == welldefined[0].attrs.TSL]
    resolved = resolve_biotypes(welldefined)
    if resolved is not None:
        result[ind] = utils.rnas.RNAGroup(ind, name, resolved, tuple(rnas))
        continue

    # GENCODE primary
    welldefined = [rna for rna in rnas if rna.attrs.tags & {'GENCODE primary'}]
    resolved = resolve_biotypes(welldefined)
    if resolved is not None:
        result[ind] = utils.rnas.RNAGroup(ind, name, resolved, tuple(rnas))
        continue

    # All other cases
    resolved = resolve_biotypes(rnas)
    if resolved is not None:
        result[ind] = utils.rnas.RNAGroup(ind, name, resolved, tuple(rnas))
        continue

    raise ValueError(f"Unresolved biotypes: {name} {biotypes}")

print("Final size of transcript groups:")
total = sum(sizes.values())
for k, size in sorted(sizes.items()):
    print(f"\t{k} -> {size}({size / total:.2%})")

# Save groups
ld.terminus.groups.dump(groups)
ld.GROUPS.dump(result)

# Make a tx2group mapping
df = [
    {"transcript_id": rna.ind, "group": clind, "name": group.name, "type": group.type}
    for clind, group in result.items() for rna in group.members
]
tx2group = pd.DataFrame(df)
tx2group.to_csv(ld.TX2GROUP, sep="\t", index=False)
