import pandas as pd

from . import ld


def overlaps() -> pd.DataFrame:
    return pd.read_pickle(ld.cCRE.overlaps.pkl)


def sequences():
    return pd.read_pickle(ld.sequences.saveto)
