from abc import ABC
from typing import Any

import numpy as np
import numpy.typing as npt
from attrs import define, field


@define(slots=True, frozen=True)
class ZeroOrderMotif(ABC):
    ind: str
    target: str
    matrix: npt.NDArray[np.float32]

    def limits(self) -> tuple[float, float]:
        return self.matrix.min(axis=0).sum(), self.matrix.max(axis=0).sum()

    def nletters(self) -> int:
        return self.matrix.shape[0]

    def __len__(self) -> int:
        return self.matrix.shape[1]


@define(slots=True, frozen=True)
class ZeroOrderMotifsCollection[T:ZeroOrderMotif]:
    alphabet: str
    attributes: dict[str, Any] = field(factory=dict)
    motifs: tuple[T, ...] = field(factory=tuple)

    def __len__(self) -> int:
        return len(self.motifs)


@define(slots=True, frozen=True)
class PositionWeightMatrix(ZeroOrderMotif):
    pass


@define(slots=True, frozen=True)
class PositionProbabilityMatrix(ZeroOrderMotif):
    def __attrs_post_init__(self):
        if not np.allclose(self.matrix.sum(axis=0), 1.0):
            raise ValueError(f"Columns must sum to 1.0, but got: {self.matrix.sum(axis=0)} for {self}")

    def to_pwm(self, bckfreq: float | npt.NDArray[np.float32] | None = None) -> PositionWeightMatrix:
        if bckfreq is None:
            bckfreq = 1 / self.matrix.shape[0]
        elif isinstance(bckfreq, np.ndarray) and bckfreq.size != self.matrix.shape[0]:
            raise ValueError(
                f"Background frequencies must have the same length as the matrix: "
                f"{self.matrix.shape[0]}, but got: {bckfreq.size} for {self}"
            )
        matrix = self.matrix.copy()
        matrix /= bckfreq
        np.log2(matrix, out=matrix)
        return PositionWeightMatrix(self.ind, self.target, matrix)


@define(slots=True, frozen=True)
class PositionFrequencyMatrix(ZeroOrderMotif):
    def to_ppm(self, pseudocnt: float = 1.0) -> PositionProbabilityMatrix:
        matrix = self.matrix.copy()
        matrix += pseudocnt
        return PositionProbabilityMatrix(self.ind, self.target, matrix / matrix.sum(axis=0))

    def to_pwm_hocomoco(self) -> PositionWeightMatrix:
        cnts = self.matrix.sum(axis=0)
        pseudocnt = np.log(cnts)
        matrix = np.log(
            (self.matrix + pseudocnt * 0.25) / ((cnts + pseudocnt) * 0.25)
        )
        return PositionWeightMatrix(self.ind, self.target, matrix)
