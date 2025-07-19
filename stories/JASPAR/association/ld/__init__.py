from pathlib import Path

ROOT = Path(__file__).parent
RESULTS = ROOT / "results"

TXGROUP_SUMMARY = RESULTS / "txgroup_summary.pkl"
STAT_TESTS = RESULTS / "stat_tests.pkl"


class thresholds:
    min_tpm = 10.0  # Minimum TPM for a transcript to be considered expressed
    zscore = (-1, 2)  # (lower, upper) bounds for z-score
    qvalue = 0.01
