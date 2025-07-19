import numpy as np
import numpy.typing as npt

from utils.motifs.motif import ZeroOrderMotifsCollection, ZeroOrderMotif


def _OHE_DNA() -> npt.NDArray[np.float32]:
    mapping = (
        ('A', [0]),
        ('C', [1]),
        ('G', [2]),
        ('T', [3]),
        ('U', [3]),
        ('M', [0, 1]),
        ('R', [0, 2]),
        ('W', [0, 3]),
        ('S', [1, 2]),
        ('Y', [1, 3]),
        ('K', [2, 3]),
        ('V', [0, 1, 2]),
        ('H', [0, 1, 3]),
        ('D', [0, 2, 3]),
        ('B', [1, 2, 3]),
        ('X', [0, 1, 2, 3]),
        ('N', [0, 1, 2, 3]),
    )
    ohe = np.zeros((4, 255), dtype=np.float32)
    for letter, indices in mapping:
        if not isinstance(letter, str) or len(letter) != 1:
            raise ValueError('Incorrect mapping format')
        upper, lower = ord(letter.upper()), ord(letter.lower())
        weight = 1 / len(indices)
        for ind in indices:
            ohe[ind, upper] = weight
            ohe[ind, lower] = weight
    return ohe


def score[T: ZeroOrderMotif](
        forward: str, revcomp: str, motifs: ZeroOrderMotifsCollection[T]
) -> list[float]:
    # Check that each motif has the same alphabet size
    sizes = set(motif.nletters() for motif in motifs.motifs)
    if len(sizes) != 1 or sizes != {len(motifs.alphabet)}:
        raise ValueError(
            f"All motifs must have the same alphabet size, but got: {sizes} for expected alphabet: {motifs.alphabet}"
        )

    if motifs.alphabet != "ACGT":
        raise ValueError("Only DNA alphabet (ACGT) is supported")

    # One-hot-encoding
    encoding = _OHE_DNA()
    ohe = []
    for seq in forward, revcomp:
        seqarray = np.fromiter(seq.encode("ASCII"), dtype=np.uint8, count=len(seq))
        ohe.append(encoding[:, seqarray])

    # Prepare a buffer for results
    fwdbuffer = np.empty(len(forward), dtype=np.float32)
    revbuffer = np.empty(len(revcomp), dtype=np.float32)

    # Score each motif
    results = []
    for motif in motifs.motifs:
        revpwm = motif.matrix[:, ::-1]
        size = len(forward) - len(motif) + 1

        # Score forward & reverse complement sequences
        for encoded, buffer in zip(ohe, (fwdbuffer, revbuffer)):
            response = buffer[:size]
            response[:] = 0

            # Convolution with the PWM nucleotide-wise
            for enc, rpwm in zip(encoded, revpwm):
                response += np.convolve(enc, rpwm, mode="valid")
        results.append(
            float(max(fwdbuffer[:size].max(), revbuffer[:size].max()))
        )
    return results
