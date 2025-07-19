#!/bin/bash
set -euxo pipefail

# See: https://github.com/COMBINE-lab/terminus/issues/5

# Path to the terminus executable
TERMINUS_BIN="ld/bin/terminus"

# Clean up old results and re-create the target directory
rm -rf ld/results/terminus
mkdir -p ld/results/terminus

# Run terminus group on each salmon directory in parallel
for fld in ld/salmon/*; do
    "$TERMINUS_BIN" group --dir "$fld" --min-spread 0.05 -o ld/results/terminus &
done

# Wait for all background processes to finish
wait

# Run terminus collapse on the results to get the final output.
# Threshold is set to 16% (any 3 samples in 18 experiments).
"$TERMINUS_BIN" collapse -c 0.16 -t "$(nproc)" -d ld/salmon/* -o ld/results/terminus
