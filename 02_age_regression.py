"""
02_age_regression.py
Load combined_mitopps.csv (from 01_run_pipeline.py) and run the
mitoPPS ~ age + sex OLS regression per (cell_type x louvain x region) combo.
Writes results to REGRESSION_PATH and prints summary stats.

Usage:
    python 02_age_regression.py
"""

import os

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from mitotyping_analysis import (
    MITOPPS_PATH, REGRESSION_PATH,
    run_age_regression,
)


def load_or_run_regression(combined_mitopps):
    if not os.path.exists(REGRESSION_PATH):
        results_df = run_age_regression(combined_mitopps)
        results_df.to_csv(REGRESSION_PATH, index=False)

    results_df = pd.read_csv(REGRESSION_PATH)
    print(f"Loaded: {len(results_df)} tests")
    print(f"Significant (FDR<0.05): {(results_df['fdr_age'] < 0.05).sum()}")
    print(f"Cell types in results: {sorted(results_df['cell_type'].dropna().unique())}")
    return results_df


def global_fdr_correction(results_df):
    """Optional: global BH-FDR across all tests, for comparison against the
    within-combo FDR already computed in run_age_regression."""
    mask = ~results_df['pval_age'].isna()
    results_df['fdr_age_global'] = np.nan
    results_df.loc[mask, 'fdr_age_global'] = multipletests(
        results_df.loc[mask, 'pval_age'], method='fdr_bh')[1]

    print(f"Significant - within-combo FDR<0.05: {(results_df['fdr_age'] < 0.05).sum()}")
    print(f"Significant - global FDR<0.05:       {(results_df['fdr_age_global'] < 0.05).sum()}")

    sig_global = results_df[results_df['fdr_age_global'] < 0.05]
    print("\nSurviving hits by cell type:")
    print(sig_global['cell_type'].value_counts())
    print("\nSurviving hits by pathway:")
    print(sig_global['pathway'].value_counts().head(20))
    return results_df


def print_significant_summary(results_df):
    sig = results_df[results_df['fdr_age'] < 0.05][
        ['cell_type', 'louvain', 'region', 'pathway', 'beta_age', 'fdr_age', 'n_animals']
    ].sort_values('fdr_age')
    print(f"\nTotal significant (within-combo FDR<0.05): {len(sig)}")
    print(sig.to_string())

    print("\nRegions with significant hits:")
    print(sig['region'].value_counts())
    print("\nRegions with NO significant hits:")
    all_regions = set(results_df['region'].unique())
    sig_regions = set(sig['region'].unique())
    print(all_regions - sig_regions)
    return sig


if __name__ == "__main__":
    combined_mitopps = pd.read_csv(MITOPPS_PATH)
    results_df = load_or_run_regression(combined_mitopps)
    results_df = global_fdr_correction(results_df)
    results_df.to_csv(REGRESSION_PATH, index=False)  # persist fdr_age_global too
    print_significant_summary(results_df)
