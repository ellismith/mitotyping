"""
mitotyping_outliers.py
Outlier / QC identification for the mitotyping pipeline.

Factored out of mitotyping_u01_2.ipynb. Pure analysis — no plotting.
The notebook keeps all figures; it just calls these to get the numbers.

Functions
---------
compute_combo_stats  : per-combo uniformity of age betas across pathways
                       (low CV / pct_positive near 0 or 1 = suspicious artifact)
find_outlier_pairs   : animal-combo pairs with mean mitoPPS > threshold x combo median
build_outlier_table  : long table of every animal-combo with outlier flag + fold-over-median
run_sensitivity      : drop outlier pairs, rerun age regression, return sens data + results
"""

import numpy as np
import pandas as pd

from mitotyping_analysis import META, MIN_ANIMALS, run_age_regression


def compute_combo_stats(results_df):
    """
    Summarize how uniform the age effect is across pathways within each
    (cell_type x louvain x region) combo.

    Real biology = pathway-specific effects (high CV, mixed directions).
    Artifact      = all pathways move together (low CV, pct_positive ~ 0 or 1).

    Returns a dataframe with mean_beta, std_beta, cv, pct_positive, n_sig, n_pathways.
    """
    combo_stats = (
        results_df
        .groupby(['cell_type', 'louvain', 'region'])
        .apply(lambda x: pd.Series({
            'mean_beta':    x['beta_age'].mean(),
            'std_beta':     x['beta_age'].std(),
            'cv':           x['beta_age'].std() / abs(x['beta_age'].mean())
                            if x['beta_age'].mean() != 0 else np.nan,
            'pct_positive': (x['beta_age'] > 0).mean(),
            'n_sig':        (x['fdr_age'] < 0.05).sum(),
            'n_pathways':   len(x),
        }))
        .reset_index()
    )
    return combo_stats


def find_outlier_pairs(combined_mitopps, threshold=5):
    """
    Systematically flag outlier (animal, cell_type, louvain, region) pairs where
    the animal's mean mitoPPS exceeds `threshold` x the within-combo median.

    Returns a set of (animal_id, cell_type, str(louvain), region) tuples.
    """
    pw_cols = [c for c in combined_mitopps.columns if c not in META]
    outlier_pairs = set()

    for keys, d in combined_mitopps.groupby(['cell_type', 'louvain', 'region']):
        d = d.copy()
        d['mean_mitopps'] = d[pw_cols].mean(axis=1)
        median = d['mean_mitopps'].median()
        outliers = d[d['mean_mitopps'] > threshold * median]
        for _, row in outliers.iterrows():
            outlier_pairs.add((row['animal_id'], keys[0], str(keys[1]), keys[2]))

    return outlier_pairs


def build_outlier_table(combined_mitopps, outlier_pairs):
    """
    Build a long per-(animal x combo) table with mean mitoPPS, combo median,
    fold-over-median, and an is_outlier flag. Used to drive the outlier figures.
    """
    pw_cols = [c for c in combined_mitopps.columns if c not in META]
    records = []

    for keys, d in combined_mitopps.groupby(['cell_type', 'louvain', 'region']):
        d = d.copy()
        d['mean_mitopps'] = d[pw_cols].mean(axis=1)
        median = d['mean_mitopps'].median()
        for _, row in d.iterrows():
            is_outlier = (row['animal_id'], keys[0], str(keys[1]), keys[2]) in outlier_pairs
            records.append({
                'cell_type':    keys[0],
                'louvain':      keys[1],
                'region':       keys[2],
                'animal_id':    row['animal_id'],
                'age':          row['age'],
                'mean_mitopps': row['mean_mitopps'],
                'combo_median': median,
                'fold_median':  row['mean_mitopps'] / (median + 1e-10),
                'is_outlier':   is_outlier,
            })

    return pd.DataFrame(records)


def run_sensitivity(combined_mitopps, outlier_pairs, min_animals=MIN_ANIMALS):
    """
    Sensitivity analysis: drop the flagged outlier pairs and rerun the age
    regression. Returns (combined_mitopps_sens, results_sens).
    """
    combined_mitopps_sens = combined_mitopps[
        ~combined_mitopps.apply(
            lambda r: (r['animal_id'], r['cell_type'], str(r['louvain']), r['region'])
                      in outlier_pairs, axis=1)
    ].copy()

    results_sens = run_age_regression(combined_mitopps_sens, min_animals=min_animals)
    return combined_mitopps_sens, results_sens
