import pandas as pd

import ld

SUMMARY = pd.read_pickle(ld.DESeq2.summary)

########################################################################################################################
# Separate reporting
########################################################################################################################
columns = [(ifn, 'mock', 'category') for ifn in ld.IFNS]

print('Total number of DETs up/down regulated by IFNs:')
for direction in ['up', 'down']:
    for partition in ['protein_coding', 'non_coding']:
        # Total number DETs up/down regulated by all IFNs
        mask = (
                (SUMMARY['partition'] == partition) &
                (SUMMARY[columns] == f'Significant {direction}').all(axis=1)
        )
        print(f"\t{partition} DETs {direction}-regulated by ALL IFNs: {mask.sum()}")

        # Total number DETs up/down regulated by any IFNs
        mask = (
                (SUMMARY['partition'] == partition) &
                (SUMMARY[columns] == f'Significant {direction}').any(axis=1)
        )
        print(f"\t{partition} DETs {direction}-regulated by ANY IFN: {mask.sum()}")

# IFN-specific DETs
print('Total number of IFN-specific DETs:')
for ifn in ld.IFNS:
    for partition in ['protein_coding', 'non_coding']:
        exclude_tags = {
            (other_ifn, 'mock', f'Significant {direction}')
            for direction in ["up", "down"] for other_ifn in ld.IFNS
            if other_ifn != ifn
        }
        for direction in ['up', 'down']:
            mask = (
                    (SUMMARY['partition'] == partition) &
                    SUMMARY['tags'].apply(
                        lambda tags: (ifn, 'mock', f"Significant {direction}") in tags and len(tags & exclude_tags) == 0
                    )
            )
            print(f"\t{ifn} {direction}-regulated {partition} transcripts: {mask.sum()}")
