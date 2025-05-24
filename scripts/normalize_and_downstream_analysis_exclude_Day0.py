# 0) Install deps
!pip install --quiet pandas numpy scipy matplotlib seaborn scikit-bio scikit-posthocs

# 2) Imports
import os, glob, re, itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal
import scikit_posthocs as sp
from skbio.stats.distance import DistanceMatrix, permdisp
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import MDS

# 3) Paths
PROJECT_ROOT = '/content/drive/MyDrive/Barcode Analysis (hot Tenanye and Hye Jeong)'
RAW_PATTERN  = os.path.join(PROJECT_ROOT, '**', '*_barcodes_R1_R2_counts.csv')
OUT_DIR      = os.path.join(PROJECT_ROOT, 'Paired_End_Processed_Ex_Day0')
SPIKE_R1     = 'AGAGTCCGTCTCTCAGA'
SPIKE_R2     = 'TCTGAGAGACGGACTCT'

os.makedirs(OUT_DIR, exist_ok=True)

# ──────────────────────────────────────────────────────────────────────────────
# A) Find all raw count files, exclude Day 0
# ──────────────────────────────────────────────────────────────────────────────
files = glob.glob(RAW_PATTERN, recursive=True)
if not files:
    raise FileNotFoundError(f"No files found under {PROJECT_ROOT}")
print(f"✔ Found {len(files)} raw CSVs (including Day 0)")

samples = []
feats   = set()
for fn in files:
    name = os.path.basename(fn)
    m = re.search(r'Day(\d+)', name)
    if not m:
        continue
    day = int(m.group(1))
    if day == 0:
        continue
    rep = re.search(r'Day\d+_([^_]+)_barcodes', name).group(1)
    df  = pd.read_csv(fn)
    df['feat'] = df['barcode_R1'] + '_' + df['barcode_R2']
    feats |= set(df['feat'])
    samples.append({'day': day, 'rep': rep, 'df': df})

days_used = sorted({s['day'] for s in samples})
print(f"✔ Using {len(samples)} samples from Days {days_used}")

feats = sorted(feats)

# Build raw abundance matrix & metadata
raw_mat = []
meta    = []
for s in samples:
    cnt = dict(zip(s['df']['feat'], s['df']['count']))
    row = np.array([cnt.get(f, 0) for f in feats], dtype=int)
    raw_mat.append(row)
    meta.append((s['day'], s['rep']))
raw_mat = np.vstack(raw_mat)
depths  = raw_mat.sum(axis=1)

# ──────────────────────────────────────────────────────────────────────────────
# B) Rarefy all samples to the minimum depth
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

pd.DataFrame(norm_mat, columns=feats) \
  .to_csv(os.path.join(OUT_DIR, 'normalized_counts_matrix.csv'), index=False)

# ──────────────────────────────────────────────────────────────────────────────
# C) Diversity metrics & statistics
# ──────────────────────────────────────────────────────────────────────────────
def compute_div(counts):
    nz   = counts[counts > 0]
    tot  = nz.sum()
    rich = nz.size
    p    = nz/tot if tot else []
    sh   = -np.sum(p * np.log(p)) if rich else 0
    si   = np.sum(p * p) if rich else 0
    ev   = sh / np.log(rich) if rich > 1 else 0
    return rich, sh, si, ev

metrics = np.array([compute_div(r) for r in norm_mat])
df_div  = pd.DataFrame(metrics, columns=['richness','shannon','simpson','evenness'])
df_div['day'] = [d for d,_ in meta]
df_div['rep'] = [r for _,r in meta]
df_div.to_csv(os.path.join(OUT_DIR, 'diversity_normalized.csv'), index=False)

with open(os.path.join(OUT_DIR, 'diversity_stats_normalized.txt'), 'w') as fout:
    fig, axes = plt.subplots(2, 2, figsize=(12,10))
    for ax, metric in zip(axes.flatten(), df_div.columns[:4]):
        # Kruskal-Wallis
        groups = [g[metric].values for _,g in df_div.groupby('day')]
        H, p_kw = kruskal(*groups)
        fout.write(f"{metric}: KW H={H:.3f}, p={p_kw:.3e}\n")
        
        # Dunn post-hoc
        dunn = sp.posthoc_dunn(df_div, val_col=metric, group_col='day', p_adjust='bonferroni')
        dunn.to_csv(os.path.join(OUT_DIR, f'dunn_{metric}.csv'))
        fout.write(f"\nDunn {metric} (bonferroni):\n{dunn}\n\n")
        
        # Plot
        sns.boxplot(data=df_div, x='day', y=metric, order=days_used, ax=ax)
        sns.swarmplot(data=df_div, x='day', y=metric, order=days_used,
                      color='k', size=3, ax=ax)
        ax.set_title(f"{metric.capitalize()} (KW p={p_kw:.3e})")
        
        # Annotate pairwise
        y_max = df_div[metric].max()
        y_min = df_div[metric].min()
        yrange = y_max - y_min
        offset = yrange * 0.1
        line_h = yrange * 0.02
        i = 0
        for d1, d2 in itertools.combinations(days_used, 2):
            p_val = dunn.loc[d1, d2]
            if p_val < 0.05:
                sig = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*'
                x1, x2 = days_used.index(d1), days_used.index(d2)
                y = y_max + offset * i
                ax.plot([x1, x1, x2, x2],
                        [y, y+line_h, y+line_h, y],
                        color='k')
                ax.text((x1+x2)/2, y+line_h*1.1, sig,
                        ha='center', va='bottom')
                i += 1

    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, 'diversity_boxplots_normalized.png'), dpi=300)
    plt.close()

print("✔ Diversity (normalized, ex Day 0) complete")

# ──────────────────────────────────────────────────────────────────────────────
# D) Gini index
# ──────────────────────────────────────────────────────────────────────────────
def gini(x):
    x = np.sort(x.astype(float))
    cum = np.cumsum(x)
    return 0 if cum[-1]==0 else (len(x)+1 - 2*(cum/cum[-1]).sum())/len(x)

df_gini = pd.DataFrame({
    'day': [d for d,_ in meta],
    'rep': [r for _,r in meta],
    'gini': [gini(r) for r in norm_mat]
})
df_gini.to_csv(os.path.join(OUT_DIR, 'gini_normalized.csv'), index=False)

H_g, p_g = kruskal(*[g['gini'].values for _,g in df_gini.groupby('day')])
dunn_g = sp.posthoc_dunn(df_gini, val_col='gini', group_col='day', p_adjust='bonferroni')
dunn_g.to_csv(os.path.join(OUT_DIR, 'dunn_gini.csv'))

with open(os.path.join(OUT_DIR, 'gini_stats_normalized.txt'), 'w') as f:
    f.write(f"Gini: KW H={H_g:.3f}, p={p_g:.3e}\n\nDunn (bonferroni):\n{dunn_g}\n")

plt.figure(figsize=(6,4))
sns.boxplot(data=df_gini, x='day', y='gini', order=days_used)
sns.swarmplot(data=df_gini, x='day', y='gini', order=days_used,
              color='k', size=3)
plt.title(f"Gini (KW p={p_g:.3e})")

# annotate
ymax = df_gini['gini'].max(); ymin = df_gini['gini'].min()
yr   = ymax - ymin; off = yr*0.1; lh=yr*0.02; i=0
for d1,d2 in itertools.combinations(days_used,2):
    pval = dunn_g.loc[d1,d2]
    if pval < 0.05:
        sig = '***' if pval<0.001 else '**' if pval<0.01 else '*'
        x1,x2 = days_used.index(d1), days_used.index(d2)
        y = ymax + off*i
        plt.plot([x1,x1,x2,x2],[y,y+lh,y+lh,y], color='k')
        plt.text((x1+x2)/2, y+lh*1.1, sig,
                 ha='center',va='bottom')
        i += 1

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, 'gini_boxplot_normalized.png'), dpi=300)
plt.close()
print("✔ Gini (normalized, ex Day 0) complete")

# ──────────────────────────────────────────────────────────────────────────────
# E) Spike-in relative abundance
# ──────────────────────────────────────────────────────────────────────────────
feat_key = f"{SPIKE_R1}_{SPIKE_R2}"
if feat_key not in feats:
    raise ValueError("Spike-in feature not found")
idx = feats.index(feat_key)

abs_ct = norm_mat[:, idx]
rel_ab = abs_ct / norm_mat.sum(axis=1)
df_s   = pd.DataFrame({
    'day':     [d for d,_ in meta],
    'rep':     [r for _,r in meta],
    'absolute': abs_ct,
    'relative': rel_ab
})
df_s.to_csv(os.path.join(OUT_DIR, 'spikein_abundance_normalized.csv'), index=False)

plt.figure(figsize=(6,4))
sns.lineplot(data=df_s, x='day', y='relative',
             hue='rep', marker='o', sort=True)
plt.title("Spike-in Relative Abundance (ex Day 0)")
plt.ylabel("Relative abundance")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, 'spikein_rel_abundance.png'), dpi=300)
plt.close()
print("✔ Spike-in relative abundance (ex Day 0) complete")

print("All analyses excluding Day 0 completed—outputs in:", OUT_DIR)
