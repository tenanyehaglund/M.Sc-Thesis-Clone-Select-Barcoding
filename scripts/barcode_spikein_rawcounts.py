import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import re

# --- User parameters ---
data_dir = '/content/drive/MyDrive/Barcode Analysis (hot Tenanye and Hye Jeong)'
target_r1 = 'AGAGTCCGTCTCTCAGA'
target_r2 = 'TCTGAGAGACGGACTCT'
file_pattern = os.path.join(data_dir, '*_barcodes_R1_R2_counts.csv')

# --- Load and parse files ---
records = []
for path in glob.glob(file_pattern):
    fname = os.path.basename(path)
    # Extract day number (e.g., Day0, Day2, Day4, etc.)
    m_day = re.search(r'Day(\d+)', fname)
    day = int(m_day.group(1)) if m_day else None
    # Extract replicate label (e.g., '1', '2', 'NR')
    m_rep = re.search(r'Day\d+_([^_]+)_barcodes', fname)
    rep = m_rep.group(1) if m_rep else ''
    # Read the counts CSV
    df_counts = pd.read_csv(path)
    # Find the spike-in count
    subset = df_counts[
        (df_counts['barcode_R1'] == target_r1) &
        (df_counts['barcode_R2'] == target_r2)
    ]
    count = int(subset['count'].iloc[0]) if not subset.empty else 0
    records.append({'day': day, 'replicate': rep, 'count': count})

# Build DataFrame and sort
df = pd.DataFrame(records)
df = df.sort_values(['day', 'replicate'])

# --- Plot per-replicate over time ---
plt.figure(figsize=(8, 5))
for rep, grp in df.groupby('replicate'):
    plt.plot(grp['day'], grp['count'], marker='o', label=f'Rep {rep}')
# Plot the mean across replicates
df_mean = df.groupby('day')['count'].mean().reset_index()
plt.plot(df_mean['day'], df_mean['count'],
         linestyle='--', marker='x', label='Mean')

plt.xlabel('Day')
plt.ylabel('Spike-in Count')
plt.title('Spike-in Barcode Abundance Over Time')
plt.legend()
plt.grid(True)
plt.show()
