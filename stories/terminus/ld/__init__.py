from pathlib import Path

from utils import PklData
from utils.rnas import RNAGroup

ROOT = Path(__file__).parent

SALMON = ROOT / "salmon"
RESULTS = ROOT / "results"


class terminus:
    output = RESULTS / "terminus"
    groups: PklData[dict[str, str]] = PklData(RESULTS / "terminus.pkl")
    executable = ROOT / "bin" / "terminus"


TPMS = RESULTS / "tpms.merged.tsv.gz"
READS = RESULTS / "reads.merged.tsv.gz"

TX2GROUP = RESULTS / "tx2group.tsv"
GROUP2GENE = RESULTS / "group2gene.tsv"
GROUPS: PklData[dict[str, RNAGroup]] = PklData(RESULTS / "clusters.pkl")
