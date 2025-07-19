from pathlib import Path
from subprocess import check_call
from typing import Iterable

import pandas as pd
from joblib import Parallel, delayed

import ld
from stories import nextflow, terminus

DESEQ2 = Path(__file__).with_suffix(".R")
assert DESEQ2.exists()

SAVETO = ld.DESeq2.tests
SAVETO.mkdir(exist_ok=True, parents=True)


def run_test(samples: Path, baseline: str, comparisons: Iterable[str]):
    params = [
        samples, terminus.TX2GROUP, baseline,
        str(ld.thresholds.log2fc), str(ld.thresholds.padj), SAVETO, "$".join(comparisons)
    ]
    print(f"Running DESeq2 with parameters: {params}")
    check_call(["Rscript", DESEQ2, *params])


# Create the samples table
records = []
for exp in nextflow.seqproj.experiments:
    quant = nextflow.FILTERED_RESULTS / "salmon" / str(exp.ind + "+" + exp.sample.attributes["title"]) / "quant.sf"
    assert quant.exists(), f"Missing salmon results for {exp.sample.attributes['title']}"
    attributes = exp.sample.attributes
    records.append({
        "condition": attributes["condition"], "donor": attributes["donor"], "title": attributes["title"], "file": quant
    })
samples = SAVETO / "samples.csv"
pd.DataFrame(records).to_csv(samples, index=False)

comparisons = {
    "mock": ["IFNa2a", "IFNa10", "IFNb", "IFNa1", "IFNo"],
    # "IFNb": ["IFNa2a", "IFNa10", "IFNa1", "IFNo"],
    # "IFNa2a": ["IFNa10", "IFNa1", "IFNo"],
    # "IFNa10": ["IFNa1", "IFNo"],
    # "IFNa1": ["IFNo"],
}

Parallel(n_jobs=-1)(delayed(run_test)(samples, baseline, cmps) for baseline, cmps in comparisons.items())
