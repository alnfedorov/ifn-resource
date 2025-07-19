import subprocess
from pathlib import Path

from biobit.io.bed import AnyBed


def tbindex(bed: list[AnyBed], Bed: type[AnyBed], saveto: Path):
    assert saveto.suffixes[-1] in {".bgz", ".gz", ".bgzip"}
    saveto.parent.mkdir(parents=True, exist_ok=True)

    bed = sorted(bed, key=lambda x: (x.seqid, x.interval.start))
    with Bed.Writer(saveto, compression="bgz") as writer:
        writer.write_records(bed)

    # Tabix index the file
    subprocess.run(f"tabix -f -p bed {saveto}", shell=True, check=True)
