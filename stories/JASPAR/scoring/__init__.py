import pandas as pd

from . import ld


def motifs() -> pd.DataFrame:
    return pd.read_pickle(ld.response.per_motif)


def clusters() -> pd.DataFrame:
    return pd.read_pickle(ld.response.per_cluster)
