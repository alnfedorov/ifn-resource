from pathlib import Path

ROOT = Path(__file__).parent
RESULTS = ROOT / "results"
RESOURCES = ROOT / "resources"


class cCRE:
    saveto = RESULTS / "cCRE"

    PLS = saveto / "PLS.bed.gz"
    pELS = saveto / "pELS.bed.gz"
    DNase_H3K4me3 = saveto / "DNase-H3K4me3.bed.gz"

    class overlaps:
        max_distances = {"PLS": 500, "pELS": 2_000, "DNase-H3K4me3": 0}
        saveto = RESULTS / "cCRE"
        pkl = saveto / "overlaps.pkl"
        bed = saveto / "overlaps.bed.gz"


class sequences:
    lengths = {"PLS": 350, "pELS": 350, "DNase-H3K4me3": 350}
    saveto = RESULTS / "sequences.pkl"
