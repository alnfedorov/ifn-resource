from dataclasses import dataclass

from biobit.toolkit.annotome.transcriptome import RNA

from assemblies.GRCh38.gencode import AttrRNA, RNAType


def is_within_universe(rna: RNA[AttrRNA]) -> bool:
    # https://www.gencodegenes.org/pages/faq.html [What does level 1, 2 or 3 mean in the GTF/GFF3?]
    # Level 1 - validated Pseudogene loci that were jointly predicted by the Yale Pseudopipe and UCSC Retrofinder
    #           pipelines as well as by Havana manual annotation; other transcripts that were verified experimentally
    #           by RT-PCR and sequencing through the GENCODE experimental pipeline.
    # Level 2 - manual annotation Havana manual annotation (and Ensembl annotation where it is identical to Havana).
    # Level 3 - automated annotation. Ensembl loci where they are different from the Havana annotation or where no
    #           Havana annotation can be found.
    return rna.attrs.level == 1 or rna.attrs.level == 2


def is_high_quality(rna: RNA[AttrRNA]) -> bool:
    return is_within_universe(rna) and (
            rna.attrs.TSL == 1 or rna.attrs.TSL == 2 or
            'MANE Select' in rna.attrs.tags or 'MANE Plus Clinical' in rna.attrs.tags
    )


@dataclass(frozen=True, slots=True)
class RNAGroup:
    ind: str
    name: str
    type: RNAType
    members: tuple[RNA[AttrRNA], ...]

    def __post_init__(self):
        assert all(is_within_universe(rna) for rna in self.members)
        assert self.type in RNAType.__args__, self.type

    @property
    def partition(self) -> str | None:
        if self.type in {
            'protein_coding', 'protein_coding_LoF', 'non_stop_decay', 'nonsense_mediated_decay',
            'protein_coding_CDS_not_defined', 'retained_intron',
            'IG_C_gene', 'TR_D_gene', 'TR_J_gene', 'IG_V_gene', 'IG_J_gene', 'TR_V_gene', 'IG_D_gene',
            'TR_C_gene',
        }:
            return "protein_coding"
        elif self.type in {
            # Pseudogenes that may:
            # * Contain functional polyadenylation sites
            # * Acquire poly-A tails during RNA processing, which are retained after reverse transcription
            'IG_pseudogene', 'IG_C_pseudogene', 'IG_V_pseudogene', 'IG_J_pseudogene',
            'TR_J_pseudogene', 'TR_V_pseudogene',
            'unitary_pseudogene', 'transcribed_unitary_pseudogene',
            'pseudogene', 'processed_pseudogene', 'unprocessed_pseudogene', 'transcribed_unprocessed_pseudogene',
            'translated_processed_pseudogene', 'transcribed_processed_pseudogene', 'processed_transcript',
            # Non-coding RNAs that might be polyadenylated
            'ribozyme', 'lncRNA',
        }:
            return "non_coding"
        else:
            # Non-coding RNAs that are not polyadenylated and likely to be artifacts of the purification process
            assert self.type in {
                'rRNA', 'misc_RNA', 'Mt_tRNA', 'vault_RNA', 'TEC', 'miRNA', 'snoRNA',
                'artifact', 'sRNA', 'rRNA_pseudogene', 'snRNA', 'scRNA', 'Mt_rRNA', 'scaRNA',
            }
            return "artifacts"
