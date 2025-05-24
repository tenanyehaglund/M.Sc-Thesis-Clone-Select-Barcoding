# 0) Install deps
!pip install --quiet pandas numpy scipy matplotlib seaborn scikit-bio

# 2) Imports
import os, glob, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal
from skbio.stats.distance import DistanceMatrix, permdisp
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import MDS

# 3) Paths
PROJECT_ROOT = '/content/drive/MyDrive/Barcode Analysis (hot Tenanye and Hye Jeong)'
RAW_PATTERN  = os.path.join(PROJECT_ROOT, '**', '*_barcodes_R1_R2_counts.csv')
OUT_DIR      = os.path.join(PROJECT_ROOT, 'Normalized_Paired_End_Processed')
SPIKE_R1     = 'AGAGTCCGTCTCTCAGA'
SPIKE_R2     = 'TCTGAGAGACGGACTCT'

os.makedirs(OUT_DIR, exist_ok=True)

# ──────────────────────────────────────────────────────────────────────────────
# A) Recursively find all raw count files
# ──────────────────────────────────────────────────────────────────────────────
files = glob.glob(RAW_PATTERN, recursive=True)
if not files:
    raise FileNotFoundError(f"No files found under {PROJECT_ROOT} matching '*_barcodes_R1_R2_counts.csv'")
print(f"✔ Found {len(files)} raw CSVs")

# Load and collect features
samples = []
feats = set()
for fn in files:
    name = os.path.basename(fn)
    day  = int(re.search(r'Day(\d+)', name).group(1))
    rep  = re.search(r'Day\d+_([^_]+)_barcodes', name).group(1)
    df   = pd.read_csv(fn)
    df['feat'] = df['barcode_R1'] + '_' + df['barcode_R2']
    feats |= set(df['feat'])
    samples.append({'day':day,'rep':rep,'df':df})

feats = sorted(feats)

# Build raw abundance matrix & metadata
raw_mat = []
meta    = []
for s in samples:
    cnt = dict(zip(s['df']['feat'], s['df']['count']))
    row = np.array([cnt.get(f,0) for f in feats], dtype=int)
    raw_mat.append(row)
    meta.append((s['day'], s['rep']))
raw_mat = np.vstack(raw_mat)
depths  = raw_mat.sum(axis=1)

# ──────────────────────────────────────────────────────────────────────────────
# B) Rarefy to the minimum depth
# ──────────────────────────────────────────────────────────────────────────────
min_depth = int(depths.min())
def rarefy(counts, depth):
    pool = np.repeat(np.arange(len(counts)), counts)
    if len(pool) < depth:
        return counts.copy()
    pick = np.random.choice(pool, depth, replace=False)
    return np.bincount(pick, minlength=len(counts))

np.random.seed(0)
norm_mat = np.vstack([rarefy(r, min_depth) for r in raw_mat])

# Save normalized matrix
pd.DataFrame(norm_mat, columns=feats).to_csv(
    os.path.join(OUT_DIR, 'normalized_counts_matrix.csv'), index=False
)

# ──────────────────────────────────────────────────────────────────────────────
# C) Diversity metrics on normalized data
# ──────────────────────────────────────────────────────────────────────────────
def compute_div(counts):
    nz = counts[counts>0]
    tot = nz.sum()
    rich = nz.size
    p    = nz/tot if tot>0 else []
    sh   = -np.sum(p*np.log(p)) if rich>0 else 0
    si   = np.sum(p*p) if rich>0 else 0
    ev   = sh/np.log(rich) if rich>1 else 0
    return rich, sh, si, ev

metrics = np.array([compute_div(r) for r in norm_mat])
df_div = pd.DataFrame(metrics, columns=['richness','shannon','simpson','evenness'])
df_div['day'] = [d for d,_ in meta]
df_div['rep'] = [r for _,r in meta]
df_div.to_csv(os.path.join(OUT_DIR,'diversity_normalized.csv'), index=False)

# Boxplots & KW
with open(os.path.join(OUT_DIR,'diversity_stats_normalized.txt'),'w') as fout:
    fig, axes = plt.subplots(2,2, figsize=(10,8))
    for ax, m in zip(axes.flatten(), df_div.columns[:4]):
        sns.boxplot(data=df_div, x='day', y=m, ax=ax)
        sns.swarmplot(data=df_div, x='day', y=m, color='k', size=3, ax=ax)
        H,p = kruskal(*[grp[m].values for _,grp in df_div.groupby('day')])
        fout.write(f"{m}: KW H={H:.3f}, p={p:.3e}\n")
        ax.set_title(f"{m} (p={p:.3e})")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR,'diversity_boxplots_normalized.png'), dpi=300)
    plt.close()
print("✔ Diversity (normalized) complete")

# ──────────────────────────────────────────────────────────────────────────────
# D) Gini index
# ──────────────────────────────────────────────────────────────────────────────
def gini(x):
    x = np.sort(x.astype(float))
    cum = np.cumsum(x)
    if cum[-1]==0: return 0
    n = len(x)
    return (n+1 - 2*(cum/cum[-1]).sum())/n

gini_vals = [gini(r) for r in norm_mat]
df_gini = pd.DataFrame({
    'day':[d for d,_ in meta],
    'rep':[r for _,r in meta],
    'gini':gini_vals})
df_gini.to_csv(os.path.join(OUT_DIR,'gini_normalized.csv'), index=False)

H,p = kruskal(*[grp['gini'].values for _,grp in df_gini.groupby('day')])
with open(os.path.join(OUT_DIR,'gini_stats_normalized.txt'),'w') as f:
    f.write(f"Gini KW H={H:.3f}, p={p:.3e}\n")

plt.figure(figsize=(4,3))
sns.boxplot(data=df_gini, x='day', y='gini')
sns.swarmplot(data=df_gini, x='day', y='gini', color='k', size=3)
plt.title(f"Gini (p={p:.3e})")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR,'gini_boxplot_normalized.png'), dpi=300)
plt.close()
print("✔ Gini (normalized) complete")

# ──────────────────────────────────────────────────────────────────────────────
# E) Spike-in relative abundance
# ──────────────────────────────────────────────────────────────────────────────
feat_key = f"{SPIKE_R1}_{SPIKE_R2}"
if feat_key not in feats:
    raise ValueError("Spike-in feature not found among raw data")
idx = feats.index(feat_key)

abs_ct = norm_mat[:, idx]
rel_ab = abs_ct / norm_mat.sum(axis=1)
df_s = pd.DataFrame({
    'day':     [d for d,_ in meta],
    'rep':     [r for _,r in meta],
    'absolute':abs_ct,
    'relative':rel_ab
})
df_s.to_csv(os.path.join(OUT_DIR,'spikein_abundance_normalized.csv'), index=False)

plt.figure(figsize=(6,4))
sns.lineplot(data=df_s, x='day', y='relative', hue='rep', marker='o')
plt.title("Spike-in Relative Abundance")
plt.ylabel("Relative abundance")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR,'spikein_rel_abundance.png'), dpi=300)
plt.close()
print("✔ Spike-in relative abundance complete")

print("All normalized analyses complete. Outputs at:", OUT_DIR)
