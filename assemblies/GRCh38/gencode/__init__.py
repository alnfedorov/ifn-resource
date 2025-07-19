from pathlib import Path
from typing import Literal

from attr import define, field
from biobit.toolkit import annotome as at

ROOT = Path(__file__).parent

gff3 = ROOT / "gencode.v48.primary_assembly.annotation.gff3.gz"
gtf = ROOT / "gencode.v48.primary_assembly.annotation.gtf.gz"
index = ROOT / "GRCh38.gencode.v48.annotome.pkl"

GeneSource = Literal['ENSEMBL', 'HAVANA']

GeneType = Literal[
    'translated_processed_pseudogene', 'snoRNA', 'snRNA', 'TR_C_gene', 'pseudogene', 'IG_C_gene', 'IG_C_pseudogene',
    'Mt_tRNA', 'transcribed_processed_pseudogene', 'transcribed_unprocessed_pseudogene', 'IG_pseudogene', 'TR_D_gene',
    'IG_D_gene', 'rRNA_pseudogene', 'scaRNA', 'unprocessed_pseudogene', 'IG_J_pseudogene', 'TR_V_pseudogene', 'miRNA',
    'processed_pseudogene', 'TR_V_gene', 'artifact', 'IG_J_gene', 'Mt_rRNA', 'IG_V_pseudogene', 'rRNA', 'TR_J_gene',
    'misc_RNA', 'lncRNA', 'vault_RNA', 'unitary_pseudogene', 'IG_V_gene', 'transcribed_unitary_pseudogene', 'sRNA',
    'TEC', 'protein_coding', 'ribozyme', 'TR_J_pseudogene'
]
GeneLevel = Literal[1, 2, 3]


@define(hash=True, slots=True, frozen=True, eq=True, order=True, repr=True, str=True)
class AttrGene:
    source: GeneSource
    level: GeneLevel
    name: str
    type: GeneType

    def __attrs_post_init__(self):
        if self.type not in GeneType.__args__:
            raise ValueError(f"Invalid biotype: {self.biotype}")
        if self.source not in GeneSource.__args__:
            raise ValueError(f"Invalid source: {self.source}")


RNASource = Literal['ENSEMBL', 'HAVANA']

RNAType = Literal[
    'translated_processed_pseudogene', 'retained_intron', 'snoRNA', 'processed_transcript', 'snRNA', 'TR_C_gene',
    'pseudogene', 'IG_C_gene', 'IG_C_pseudogene', 'Mt_tRNA', 'transcribed_processed_pseudogene',
    'transcribed_unprocessed_pseudogene', 'IG_pseudogene', 'TR_D_gene', 'IG_D_gene', 'non_stop_decay',
    'rRNA_pseudogene', 'scaRNA', 'unprocessed_pseudogene', 'IG_J_pseudogene', 'TR_V_pseudogene', 'miRNA',
    'processed_pseudogene', 'TR_V_gene', 'protein_coding_LoF', 'nonsense_mediated_decay', 'artifact', 'IG_J_gene',
    'Mt_rRNA', 'IG_V_pseudogene', 'rRNA', 'protein_coding_CDS_not_defined', 'TR_J_gene', 'misc_RNA', 'lncRNA',
    'vault_RNA', 'unitary_pseudogene', 'IG_V_gene', 'transcribed_unitary_pseudogene', 'sRNA', 'TEC', 'protein_coding',
    'ribozyme', 'TR_J_pseudogene'
]
RNATag = Literal[
    'Ensembl canonical', 'GENCODE basic', 'GENCODE primary', 'MANE Select', 'MANE Plus Clinical',
    'appris_principal_1', 'retained_intron_CDS', '5_nested_supported_extension', '454_RNA_Seq_supported',
    'non_canonical_polymorphism', 'mRNA_start_NF', 'inferred_exon_combination', 'bicistronic',
    'NMD_likely_if_extended', 'alternative_3_UTR', 'non_canonical_TEC', 'non_canonical_genome_sequence_error',
    '5_standard_supported_extension', 'non_canonical_U12', 'non_canonical_other',
    'NAGNAG_splice_site', '3_standard_supported_extension', 'low_sequence_quality', 'not_best_in_genome_evidence',
    'upstream_ATG', 'dotter_confirmed', 'retained_intron_first', 'exp_conf', 'appris_principal_4',
    '3_nested_supported_extension', 'appris_alternative_2', 'mRNA_end_NF', 'appris_principal_5', 'pseudo_consens',
    'appris_principal_2', 'non_ATG_start', 'non_canonical_conserved', 'RP_supported_TIS', 'cds_start_NF',
    'upstream_uORF', 'non_submitted_evidence', 'CAGE_supported_TSS', 'appris_alternative_1', 'stop_codon_readthrough',
    'seleno', 'RNA_Seq_supported_only', 'appris_principal_3', 'inferred_transcript_model',
    'retained_intron_final', 'RNA_Seq_supported_partial', 'CCDS', 'NMD_exception', 'nested_454_RNA_Seq_supported',
    'sequence_error', 'alternative_5_UTR', 'downstream_ATG', 'cds_end_NF', 'readthrough_transcript',
    'not_organism_supported', 'overlapping_uORF', 'TAGENE'
]
RNATSL = Literal[1, 2, 3, 4, 5, None]
RNALevel = Literal[1, 2, 3]


@define(hash=True, slots=True, frozen=True, eq=True, order=True, repr=True, str=True)
class AttrRNA:
    source: RNASource
    level: RNALevel
    name: str
    type: RNAType
    tags: frozenset[RNATag]
    TSL: RNATSL
    CDS: frozenset[str]

    def __attrs_post_init__(self):
        if self.source not in RNASource.__args__:
            raise ValueError(f"Invalid source: {self.source}")
        if self.type not in RNAType.__args__:
            raise ValueError(f"Invalid biotype: {self.type}")
        if any(x not in RNATag.__args__ for x in self.tags):
            raise ValueError(f"Invalid tags: {self.tags}")
        if self.TSL not in RNATSL.__args__:
            raise ValueError(f"Invalid TSL: {self.TSL}")


CDSSource = Literal['ENSEMBL', 'HAVANA']


@define(hash=True, slots=True, frozen=True, eq=True, order=True, repr=True, str=True)
class AttrCDS:
    source: CDSSource
    transcripts: frozenset[str] = field(converter=lambda x: frozenset(x))

    def __attrs_post_init__(self):
        if self.source not in CDSSource.__args__:
            raise ValueError(f"Invalid source: {self.source}")


def load() -> at.Annotome[AttrGene, AttrRNA, AttrCDS]:
    return at.read_pkl(index.as_posix())
