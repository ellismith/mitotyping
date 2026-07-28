"""
04_visualization.py
===================
Generate all figures for the U01 macaque brain aging mitotyping analysis.

Reads from:
    results/all_mitopps.csv
    results/mitotyping_age_results_clean.csv
    results/mitotyping_age_results_sensitivity.csv

Writes to:
    figures/   (300dpi PNGs, auto-named from plot titles)

Run:
    conda run -n mito_env_test python3 04_visualization.py

All figures are saved via save=True in mitotyping_plots. To regenerate a
subset, comment out the sections you don't need.
"""

import matplotlib
matplotlib.use('Agg')   # non-interactive backend for cluster use

import sys
import os
sys.path.insert(0, '/scratch/easmit31/mitotyping')

import importlib
import mitotyping_plots as mp
importlib.reload(mp)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ── paths ─────────────────────────────────────────────────────────────────────
RESULTS_DIR = '/scratch/easmit31/mitotyping/results/'
FIGURES_DIR = '/scratch/easmit31/mitotyping/figures/'
os.makedirs(FIGURES_DIR, exist_ok=True)

# ── load data ─────────────────────────────────────────────────────────────────
print("Loading data...")
combined_mitopps = pd.read_csv(f'{RESULTS_DIR}all_mitopps.csv')
results_df       = pd.read_csv(f'{RESULTS_DIR}mitotyping_age_results_clean.csv')
results_sens     = pd.read_csv(f'{RESULTS_DIR}mitotyping_age_results_sensitivity.csv')

cell_types = sorted(combined_mitopps['cell_type'].unique())
print(f"Cell types: {cell_types}")
print(f"Sig hits: {(results_df['fdr_age'] < 0.05).sum()}")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — Overview heatmaps (age effects, all cell types)
# ══════════════════════════════════════════════════════════════════════════════
print("\n[1/6] Overview age effect heatmaps...")

# main figure — all cell types, louvains averaged, vmax=0.03
mp.plot_mitotyping(results_df, mode='age', average_louvains=True,
                   only_sig=True, vmax=0.03, save=True)
plt.close('all')

# L1 level summary — compact 7-row overview
mp.plot_mitotyping(results_df, mode='age', level='L1',
                   average_louvains=True, only_sig=False, save=True)
plt.close('all')

# by region — cols = regions averaged across cell types
mp.plot_mitotyping(results_df, mode='age', by_region=True,
                   only_sig=False, vmax=0.03, save=True)
plt.close('all')

# clustered rows
mp.plot_mitotyping(results_df, mode='age', average_louvains=True,
                   only_sig=True, vmax=0.03, cluster_rows=True, save=True)
plt.close('all')


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — Per-region heatmaps (faceted, cols = cell types)
# ══════════════════════════════════════════════════════════════════════════════
print("\n[2/6] Per-region faceted heatmaps...")

# significant only
mp.plot_mitotyping(results_df, mode='age', facet_region=True,
                   only_sig=True, vmax=0.03, save=True)
plt.close('all')

# all pathways (no sig filter) — shows null regions too
mp.plot_mitotyping(results_df, mode='age', facet_region=True,
                   only_sig=False, vmax=0.03, save=True)
plt.close('all')


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — Per-cell-type heatmaps
# ══════════════════════════════════════════════════════════════════════════════
print("\n[3/6] Per-cell-type heatmaps...")

for ct in cell_types:
    print(f"  {ct}...")

    # mitoPPS — all louvains x regions
    mp.plot_mitotyping(combined_mitopps, mode='mitopps', cell_type=ct,
                       level='L3', save=True)
    plt.close('all')

    # mitoPPS — louvains averaged, one col per region
    mp.plot_mitotyping(combined_mitopps, mode='mitopps', cell_type=ct,
                       level='L3', average_louvains=True, save=True)
    plt.close('all')

    # mitoPPS — per animal ordered by age
    mp.plot_mitotyping(combined_mitopps, mode='mitopps', cell_type=ct,
                       level='L3', per_animal=True, save=True)
    plt.close('all')

    # age effects — louvains averaged
    sub = results_df[results_df['cell_type'] == ct]
    mp.plot_mitotyping(sub, mode='age', cell_type=ct,
                       level='L3', average_louvains=True,
                       only_sig=False, vmax=0.03, save=True)
    plt.close('all')

    # age effects — faceted by region (one panel per region)
    mp.plot_mitotyping(sub, mode='age', cell_type=ct,
                       level='L3', facet_region=True,
                       only_sig=False, vmax=0.03, save=True)
    plt.close('all')


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — Summary plots
# ══════════════════════════════════════════════════════════════════════════════
print("\n[4/6] Summary plots...")

direction, consistency_mat = mp.summarize_mitotyping(results_df)

# direction bar chart — all cell types including 0-hit ones
all_cts = results_df['cell_type'].dropna().unique()
direction_all = (results_df[results_df['fdr_age'] < 0.05]
                 .groupby('cell_type')
                 .apply(lambda x: pd.Series({
                     'n_sig':          len(x),
                     'n_increasing':   (x['beta_age'] > 0).sum(),
                     'n_decreasing':   (x['beta_age'] < 0).sum(),
                     'pct_increasing': round((x['beta_age'] > 0).mean() * 100, 1),
                 }))
                 .reindex(all_cts).fillna(0)
                 .reset_index().sort_values('n_sig', ascending=False))

mp.plot_direction(direction_all, save=True)
plt.close('all')

mp.plot_consistency(consistency_mat, save=True)
plt.close('all')

mp.plot_consistency(consistency_mat, min_frac=0.5, save=True)
plt.close('all')

# dotplot — significant hits by region
mp.plot_significant_dotplot(results_df, results_sens, save=True)
plt.close('all')


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — Aging pathways (top hits per cell type + scatter vs age)
# ══════════════════════════════════════════════════════════════════════════════
print("\n[5/6] Aging pathway plots...")

for ct in cell_types:
    n_sig = (results_df[results_df['cell_type'] == ct]['fdr_age'] < 0.05).sum()
    print(f"  {ct} ({n_sig} sig hits)...")
    mp.plot_aging_pathways(results_df, combined_mitopps,
                           cell_type=ct, save=True)
    plt.close('all')


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6 — Data availability summary
# ══════════════════════════════════════════════════════════════════════════════
print("\n[6/6] Data availability figure...")

# coverage table — requires h5ad access, run separately if needed
# from mitotyping_analysis import run_coverage_summary
# run_coverage_summary(save=True)

# distribution of mitoPPS scores
META    = {'animal_id','louvain','region','age','sex','n_cells','cell_type'}
pw_cols = [c for c in combined_mitopps.columns if c not in META]

fig, ax = plt.subplots(figsize=(7, 4))
vals = combined_mitopps[pw_cols].values.flatten()
vals = vals[~np.isnan(vals)]
ax.hist(vals, bins=100, color='steelblue', edgecolor='none', log=True)
ax.axvline(np.percentile(vals, 99), color='red', linewidth=1.5,
           label=f"99th pct = {np.percentile(vals, 99):.2f}")
ax.set_xlabel('mitoPPS value')
ax.set_ylabel('count (log scale)')
ax.set_title('mitoPPS score distribution')
ax.legend(fontsize=8)
plt.tight_layout()
fig.savefig(f'{FIGURES_DIR}mitoPPS_distribution.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"  Saved: mitoPPS_distribution.png")


print(f"\nDone. All figures saved to {FIGURES_DIR}")
