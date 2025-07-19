import pandas as pd

import ld
from stories import terminus

GROUPS = terminus.GROUPS.load()

rld = pd.read_csv(ld.DESeq2.rld, index_col=0)

# Derive row data
row_data = rld[['Name']]
row_data['Partition'] = row_data.index.map(lambda x: GROUPS[x].partition)
row_data['Ensembl IDs'] = row_data.index.map(lambda x: ";".join(x.ind for x in GROUPS[x].members))

row_data = row_data[row_data['Partition'] != 'artifacts'].copy()

# Keep only relevant rows in the RLD
rld = rld.loc[row_data.index].drop(columns=['Name'])

# Calculate Z-scores
zscore = rld.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

# Load TPMs and select relevant rows
tpms = pd.read_csv(terminus.TPMS, sep="\t", index_col=0)
tpms = tpms.loc[row_data.index]

# Derive column data
col_data = pd.DataFrame(index=rld.columns)
col_data['Donor'] = col_data.index.map(lambda x: x.split('+')[1].split('-')[0])
col_data['Treatment'] = col_data.index.map(lambda x: x.split('+')[1].split('-')[1])
col_data = col_data.sort_values(by=['Treatment', 'Donor'])

# Rename the columns to match the R standard
mapping = {x: x.split('+')[1].replace('-', '_') for x in col_data.index}
col_data.index = col_data.index.map(mapping)

# Reorder the columns to follow a desired order
col_data['Order'] = col_data['Treatment'].map({"mock": 0, "IFNa1": 1, "IFNa2a": 2, "IFNa10": 3, "IFNo": 4, "IFNb": 5})
col_data = col_data.sort_values(by=['Order', 'Donor']).drop(columns=['Order'])

# Reorder columns in all tables to match the column data
rld = rld.rename(columns=mapping)[col_data.index]
tpms = tpms.rename(columns=mapping)[col_data.index]
zscore = zscore.rename(columns=mapping)[col_data.index]

# Re-index the data to satisfy iSEE requirements
assert row_data['Name'].is_unique, "Row names must be unique for iSEE"
mapping = row_data['Name'].to_dict()

row_data = row_data.reset_index().set_index('Name')
rld.index = rld.index.map(mapping)
zscore.index = zscore.index.map(mapping)
tpms.index = tpms.index.map(mapping)

# Save the tables
ld.iSEE.root.mkdir(parents=True, exist_ok=True)

rld.to_csv(ld.iSEE.rld, index=True)
zscore.to_csv(ld.iSEE.zscore, index=True)
tpms.to_csv(ld.iSEE.tpm, index=True)

row_data.to_csv(ld.iSEE.row_data, index=True)
col_data.to_csv(ld.iSEE.col_data, index=True)
