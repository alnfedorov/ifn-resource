from pathlib import Path

from . import seqid, gencode

ROOT = Path(__file__).parent

name = "GRCh38"
organism = "Homo sapiens"

fasta = ROOT / "GRCh38.primary_assembly.genome.fa.bgz"
chromsizes = ROOT / "GRCh38.primary_assembly.genome.sizes"

cCRE = ROOT / "GRCh38-cCREs.bed.gz"  # ENCODE cCRE version 3
