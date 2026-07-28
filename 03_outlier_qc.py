"""
03_outlier_qc.py
QC pass over the age-regression results: find combos with suspiciously
uniform betas (possible artifacts), flag outlier animal-combo pairs
(mean mitoPPS > 5x within-combo median), rerun regression with them
dropped (sensitivity analysis), and save all diagnostic plots to disk.

Usage:
    python 03_outlier_qc.py
"""

import os

import matplotlib
matplotlib.use("Agg")  # no display available outside a notebook
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mitotyping_analysis import MITOPPS_PATH, REGRESSION_PATH, OUTDIR
import mitotyping_outliers as mto

QC_DIR = f"{OUTDIR}/qc_plots"
os.makedirs(QC_DIR, exist_ok=True)


def savefig(name):
    path = f"{QC_DIR}/{name}.png"
    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  saved {path}")


def suspicious_combos(results_df):
    combo_stats = mto.compute_combo_stats(results_df)

    print("Most suspicious combos (lowest CV = most uniform across pathways):")
    print(combo_stats.nsmallest(20, 'cv')[
        ['cell_type', 'louvain', 'region', 'mean_beta', 'std_beta', 'cv', 'pct_positive', 'n_sig']
    ])

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(combo_stats['cv'].dropna(), bins=50, color='steelblue', edgecolor='none')
    ax.set_xlabel('coefficient of variation of beta_age across pathways')
    ax.set_ylabel('n combos')
    ax.set_title('How uniform are age effects across pathways within each combo?\n'
                 'Low CV = all pathways moving together (suspicious)')
    savefig('combo_cv_distribution')

    return combo_stats


def inspect_worst_combo(combo_stats, combined_mitopps):
    worst = combo_stats[combo_stats['pct_positive'] > 0.95].nsmallest(1, 'cv').iloc[0]
    print(f"Examining: {worst['cell_type']} | louvain {worst['louvain']} | {worst['region']}")
    print(f"Mean beta: {worst['mean_beta']:.4f}, CV: {worst['cv']:.3f}")

    sub = combined_mitopps[
        (combined_mitopps['cell_type'] == worst['cell_type']) &
        (combined_mitopps['louvain'].astype(str) == str(worst['louvain'])) &
        (combined_mitopps['region'] == worst['region'])
    ].copy()

    from mitotyping_analysis import META
    pw_cols = [c for c in combined_mitopps.columns if c not in META]

    print(f"\nN animals: {sub['animal_id'].nunique()}")
    print(f"N cells per sample: {sub['n_cells'].describe()}")
    print("\nmitoPPS summary across all pathways:")
    print(sub[pw_cols].stack().describe())

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for ax, pw in zip(axes, pw_cols[:3]):
        ax.scatter(sub['age'], sub[pw], s=30, alpha=0.7)
        z = np.polyfit(sub['age'].astype(float), sub[pw].astype(float), 1)
        x = np.linspace(sub['age'].min(), sub['age'].max(), 100)
        ax.plot(x, np.polyval(z, x), 'r-', linewidth=1.5)
        ax.set_xlabel('age')
        ax.set_ylabel('mitoPPS')
        ax.set_title(pw[:30])
    plt.suptitle(f"{worst['cell_type']} | louvain {worst['louvain']} | {worst['region']}", fontsize=10)
    savefig('worst_combo_scatter')

    mean_mitopps_per_animal = sub[pw_cols].mean(axis=1)
    print("\nMean mitoPPS per animal:")
    print(pd.DataFrame({
        'animal_id': sub['animal_id'].values,
        'age': sub['age'].values,
        'n_cells': sub['n_cells'].values,
        'mean_mitopps': mean_mitopps_per_animal.values,
    }).sort_values('mean_mitopps', ascending=False).head(10))

    return worst


def find_and_check_outliers(combined_mitopps, results_df, threshold=5):
    outlier_pairs = mto.find_outlier_pairs(combined_mitopps, threshold=threshold)
    print(f"Outlier animal-combo pairs found: {len(outlier_pairs)}")
    for pair in sorted(outlier_pairs):
        print(f"  {pair}")

    combined_mitopps_sens, results_sens = mto.run_sensitivity(combined_mitopps, outlier_pairs)

    print("Significant hits - original:", (results_df['fdr_age'] < 0.05).sum())
    print("Significant hits - sensitivity:", (results_sens['fdr_age'] < 0.05).sum())

    sig_orig = set(zip(results_df[results_df['fdr_age'] < 0.05]['cell_type'],
                        results_df[results_df['fdr_age'] < 0.05]['pathway']))
    sig_sens = set(zip(results_sens[results_sens['fdr_age'] < 0.05]['cell_type'],
                        results_sens[results_sens['fdr_age'] < 0.05]['pathway']))

    print(f"\nHits in original only: {len(sig_orig - sig_sens)}")
    print(f"Hits in sensitivity only: {len(sig_sens - sig_orig)}")
    print(f"Hits in both: {len(sig_orig & sig_sens)}")

    return outlier_pairs, combined_mitopps_sens, results_sens


def plot_outlier_diagnostics(combined_mitopps, outlier_pairs):
    plot_df = mto.build_outlier_table(combined_mitopps, outlier_pairs)

    # ── 3-panel overview: distribution, fold-vs-age, rate by cell type ──
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    ax = axes[0]
    ax.hist(plot_df[~plot_df['is_outlier']]['mean_mitopps'],
            bins=100, color='steelblue', alpha=0.7, label='normal', log=True)
    ax.hist(plot_df[plot_df['is_outlier']]['mean_mitopps'],
            bins=30, color='red', alpha=0.8, label='outlier', log=True)
    ax.axvline(5, color='k', linewidth=1, linestyle='--', label='5x median threshold approx')
    ax.set_xlabel('mean mitoPPS per animal-combo')
    ax.set_ylabel('count (log scale)')
    ax.set_title('Distribution of mean mitoPPS\nred = flagged outliers')
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.scatter(plot_df[~plot_df['is_outlier']]['fold_median'],
               plot_df[~plot_df['is_outlier']]['age'],
               color='steelblue', s=3, alpha=0.3, label='normal')
    ax.scatter(plot_df[plot_df['is_outlier']]['fold_median'],
               plot_df[plot_df['is_outlier']]['age'],
               color='red', s=30, alpha=0.9, label='outlier', zorder=5)
    ax.axvline(5, color='k', linewidth=1, linestyle='--', label='threshold')
    ax.set_xlabel('fold over combo median')
    ax.set_ylabel('animal age')
    ax.set_title('Outliers vs age\nare outliers concentrated at certain ages?')
    ax.legend(fontsize=8)

    ax = axes[2]
    outlier_counts = (plot_df[plot_df['is_outlier']]
                      .groupby('cell_type').size()
                      .sort_values(ascending=True))
    normal_counts = (plot_df[~plot_df['is_outlier']]
                     .groupby('cell_type').size()
                     .reindex(outlier_counts.index))
    pct_outlier = (outlier_counts / (outlier_counts + normal_counts) * 100).round(1)

    ax.barh(range(len(outlier_counts)), pct_outlier.values, color='red', alpha=0.7)
    ax.set_yticks(range(len(outlier_counts)))
    ax.set_yticklabels(outlier_counts.index, fontsize=8)
    ax.set_xlabel('% of animal-combo pairs flagged as outliers')
    ax.set_title('Outlier rate by cell type')

    plt.suptitle(f'mitoPPS outlier animal-combo pairs (n={len(outlier_pairs)})\n'
                 '5x within-combo median threshold', fontsize=11, y=1.02)
    savefig('outlier_overview')

    print("\nOutliers per cell type:")
    print(plot_df[plot_df['is_outlier']]['cell_type'].value_counts())
    print("\nOutliers per region:")
    print(plot_df[plot_df['is_outlier']]['region'].value_counts())
    print(f"\nOutlier animals: {sorted(plot_df[plot_df['is_outlier']]['animal_id'].unique())}")

    # ── fold-over-median distribution, global ──
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(plot_df['fold_median'].clip(0, 20), bins=100,
            color='steelblue', edgecolor='none', log=True)
    ax.axvline(5, color='red', linewidth=1.5, linestyle='--', label='current threshold (5x)')
    ax.set_xlabel('fold over within-combo median')
    ax.set_ylabel('count (log scale)')
    ax.set_title('Distribution of fold-over-median\nwhere is the natural break?')
    ax.legend()
    savefig('fold_median_distribution')

    print(pd.Series(plot_df['fold_median']).describe(percentiles=[.9, .95, .99, .999]))

    # ── per-cell-type grid ──
    cell_types = sorted(plot_df['cell_type'].unique())
    n_cts = len(cell_types)
    n_cols = 3
    n_rows = int(np.ceil(n_cts / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 5, n_rows * 3))
    axes = axes.flatten()

    for ax, ct in zip(axes, cell_types):
        sub = plot_df[plot_df['cell_type'] == ct]['fold_median'].clip(0, 20)
        ax.hist(sub, bins=50, color='steelblue', edgecolor='none', log=True)

        p99 = sub.quantile(0.99)
        p999 = sub.quantile(0.999)
        ax.axvline(p99, color='orange', linewidth=1.5, linestyle='--', label=f'99th={p99:.1f}x')
        ax.axvline(p999, color='red', linewidth=1.5, linestyle='--', label=f'99.9th={p999:.1f}x')
        ax.axvline(5, color='gray', linewidth=1, linestyle=':', label='5x (global)')

        n_out = (sub > 5).sum()
        ax.set_title(f"{ct}\n(n_outliers>5x: {n_out})", fontsize=9)
        ax.set_xlabel('fold over median', fontsize=7)
        ax.set_ylabel('count (log)', fontsize=7)
        ax.legend(fontsize=6, frameon=False)

    for ax in axes[n_cts:]:
        ax.set_visible(False)

    plt.suptitle('Fold-over-median distribution per cell type\n'
                 'orange=99th pct, red=99.9th pct, gray=5x global threshold',
                 fontsize=11, y=1.01)
    savefig('fold_median_by_cell_type')

    print("\n99th and 99.9th percentiles of fold-over-median per cell type:")
    print(plot_df.groupby('cell_type')['fold_median'].quantile([0.99, 0.999]).unstack().round(2))

    return plot_df


if __name__ == "__main__":
    combined_mitopps = pd.read_csv(MITOPPS_PATH)
    results_df = pd.read_csv(REGRESSION_PATH)

    combo_stats = suspicious_combos(results_df)
    inspect_worst_combo(combo_stats, combined_mitopps)
    outlier_pairs, combined_mitopps_sens, results_sens = find_and_check_outliers(
        combined_mitopps, results_df)
    plot_outlier_diagnostics(combined_mitopps, outlier_pairs)

    # persist sensitivity results for downstream use
    results_sens.to_csv(f"{OUTDIR}/mitotyping_age_results_sensitivity.csv", index=False)
    print(f"\nSaved sensitivity regression results to "
          f"{OUTDIR}/mitotyping_age_results_sensitivity.csv")
