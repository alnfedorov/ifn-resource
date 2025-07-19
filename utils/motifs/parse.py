import os
from io import TextIOBase

import numpy as np

from .motif import ZeroOrderMotifsCollection, PositionFrequencyMatrix


def jaspar(file: os.PathLike[str] | TextIOBase) -> ZeroOrderMotifsCollection[PositionFrequencyMatrix]:
    if not isinstance(file, TextIOBase):
        with open(file) as stream:
            return jaspar(stream)
    assert isinstance(file, TextIOBase)

    # Load JASPAR motifs
    lines = file.readlines()
    if len(lines) % 5 != 0:
        raise ValueError(f"Expected multiple of 5 lines for JASPAR motifs, but got: {len(lines)}")

    motifs = []
    for ind in range(0, len(lines), 5):
        title, A, C, G, T = lines[ind:ind + 5]
        assert title.startswith(">"), f"Unexpected title: {title}"

        pfm = []
        for nuc, ntype in zip([A, C, G, T], ["A", "C", "G", "T"]):
            nuc = nuc.strip().replace("[", "").replace("]", "").split()
            assert nuc[0] == ntype, f"Unexpected line: {nuc}"
            nuc = nuc[1:]
            nuc = [float(x) for x in nuc]
            pfm.append(nuc)

        # Meta information
        meta = title.strip()[1:].split("\t")
        if len(meta) == 1:
            ind, target = meta[0], ""
        elif len(meta) == 2:
            ind, target = meta
        else:
            raise ValueError(f"Unexpected meta information: {meta}")

        # Position frequency matrix
        pcm = PositionFrequencyMatrix(ind, target, np.array(pfm, dtype=np.float32))
        motifs.append(pcm)

    return ZeroOrderMotifsCollection("ACGT", motifs=motifs)
