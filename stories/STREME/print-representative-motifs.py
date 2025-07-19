import xml.etree.ElementTree as ET

import pandas as pd

import ld

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

# Parse all STREME results
records = []
for streme in ld.streme.comparisons.glob(f"*/*/streme/streme.xml"):
    comparison = streme.parent.parent.name
    region = streme.parent.parent.parent.name

    with open(streme, "r") as f:
        et = ET.parse(f)
        train_positives = int(et.find('model/train_positives').attrib['count'])
        train_negatives = int(et.find('model/train_negatives').attrib['count'])
        test_positives = int(et.find('model/test_positives').attrib['count'])
        test_negatives = int(et.find('model/test_negatives').attrib['count'])

        for motif in et.findall('motifs/motif'):
            records.append({
                'comparison': comparison,
                'region': region,
                'train_positives': train_positives,
                'train_negatives': train_negatives,
                'test_positives': test_positives,
                'test_negatives': test_negatives,
                **motif.attrib
            })

df = pd.DataFrame(records).drop(columns=[
    'alt', 'width', 'initial_width', 'score_threshold', 'is_palindromic', 'elapsed_time', 'site_distr', 'site_hist',
    'train_dtc', 'train_bernoulli', 'test_dtc', 'test_bernoulli',
    'test_log_evalue', 'test_log_pvalue', 'train_log_pvalue',
    'total_sites', 'max_sites', 'npassing'  # Not sure what these are
]).astype({
    'test_evalue': float, 'test_pvalue': float, 'train_pvalue': float,
    'train_pos_count': int, 'train_neg_count': int, 'test_pos_count': int, 'test_neg_count': int,
})

# Normalize and select relevant columns
df[['target', 'background']] = df['comparison'].str.split('-vs-', expand=True)
df['id'] = df['id'].str.split('-', expand=True)[1]

df['train_pos_fraction'] = df['train_pos_count'] / df['train_positives']
df['train_neg_fraction'] = df['train_neg_count'] / df['train_negatives']
df['train_enrichment'] = df['train_pos_fraction'] / df['train_neg_fraction']

df['test_pos_fraction'] = df['test_pos_count'] / df['test_positives']
df['test_neg_fraction'] = df['test_neg_count'] / df['test_negatives']
df['test_enrichment'] = df['test_pos_fraction'] / df['test_neg_fraction']

df = df[[
    'target', 'background', 'region', 'id', 'test_evalue', 'test_pvalue',
    'train_enrichment', 'test_enrichment',
    'train_pos_fraction', 'train_neg_fraction',
    'train_positives', 'train_negatives',
    'test_pos_fraction', 'test_neg_fraction',
    'test_positives', 'test_negatives'
]].sort_values(by=['test_evalue'])

# print(f'Total rows: {len(df)}')
# print(df.head(50))

# Remove all low quality motifs
high_quality = df[
    (df['test_evalue'] <= 0.01) &
    (df['test_enrichment'] >= 4) & (df['train_enrichment'] >= 4)
    ]
print(f'Total rows: {len(high_quality)}')
print(high_quality)
