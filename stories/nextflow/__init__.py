from pathlib import Path

from biobit.toolkit import seqproj, nfcore

ROOT = Path(__file__).parent
RESOURCES = ROOT / "resources"
FULL_RESULTS = ROOT / "full-results"
FILTERED_RESULTS = ROOT / "filtered-results"

seqproj = nfcore.rnaseq.parse.into_seqproj(
    seqproj.adapter.yaml.load(ROOT / 'seq-project.yaml'),
    FULL_RESULTS,
    seqexp2descriptor=lambda exp: nfcore.rnaseq.descriptor.from_seqexp(
        exp, title_builder=lambda x: x.sample.attributes['title']
    )
)
