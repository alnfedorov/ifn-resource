from pathlib import Path
from typing import assert_never

from biobit.core.loc import IntoStrand, Strand, IntoInterval
from biobit.io.fasta import IndexedReader, Reader

COMPLEMENTARITY = str.maketrans("ACGT", "TGCA")


def fetch(reader: IndexedReader, seqid: str, interval: IntoInterval, strand: IntoStrand = "+") -> str:
    seq = reader.fetch(seqid, interval).upper()
    match Strand(strand):
        case Strand.Forward:
            return seq
        case Strand.Reverse:
            return seq[::-1].translate(COMPLEMENTARITY)
        case _:
            assert_never(strand)


def read(path: Path | str) -> dict[str, str]:
    records = Reader(path).read_to_end()
    return {r.id: r.seq for r in records}
