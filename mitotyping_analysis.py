"""
mitotyping_analysis.py
Core mitotyping analysis for the U01 macaque brain-aging snRNA-seq dataset.

Config + pipeline functions factored out of mitotyping_u01_2.ipynb.
Edit the config block below to change paths/filters — the notebook imports
these names, so this file is the single source of truth for settings.

Functions
---------
run_mitotyping      : pseudobulk -> MitoPathway scores -> mitoPPS, per cell type
run_age_regression  : mitoPPS ~ age + sex OLS per (cell_type x louvain x region)
flag_bad_combos     : QC for combos with over-correlated pathway mitoPPS
"""

import os
import gc
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import issparse
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import multipletests

# shared metadata column set (non-pathway columns)
META = {'animal_id', 'louvain', 'region', 'age', 'sex', 'n_cells', 'cell_type'}

# ── standard cell-type abbreviations (matches other pipelines) ────────────
# Keys = abbreviation (already alphabetical), values = full name used in
# raw data / results dataframes.
CELL_TYPE_ABBR = {
    'Astrocyte':        'AST',
    'Basket':           'BC',
    'Cerebellar':       'CER',
    'Ependymal':        'EPEN',
    'Glutamatergic':    'EXC',
    'GABA':             'INH',
    'Midbrain':         'MBN',
    'Microglia':        'MGL',
    'MSN':              'MSN',
    'Oligodendrocyte':  'OLIG',
    'OPC':              'OPC',
    'Vascular':         'VASC',
}
CELL_TYPE_ORDER = sorted(CELL_TYPE_ABBR.values())


def apply_cell_type_abbr(df):
    """
    Relabel a dataframe's 'cell_type' column to standard abbreviations.
    Note: mitotyping_plots.cell_type_colors now uses abbreviated keys too,
    so call this before passing data into any mp.plot_* function, not after.
    """
    df = df.copy()
    df['cell_type'] = df['cell_type'].map(CELL_TYPE_ABBR).fillna(df['cell_type'])
    return df

# ==========================================================================
# config
# ==========================================================================
# ── data paths ────────────────────────────────────────────────────────────
DATA_DIR       = "/data/CEM/smacklab/data/bican/u01_cell_class/intermediate_cell_class"
MITOCARTA_PATH = "/scratch/easmit31/mitotyping/Human.MitoCarta3.0.xls"
OUTDIR         = "/scratch/easmit31/mitotyping/results/"
ORTHOLOG_CSV   = "/scratch/easmit31/data/human-macaque-orthologs/ensembl113_mmul10_macaque_human.csv"

# ── filters ───────────────────────────────────────────────────────────────
MIN_CELLS         = 100  # min cells per (animal × louvain × region) pseudobulk group
MIN_ANIMALS       = 10   # min unique animals per (louvain × region) combo for regression
MIN_AGE           = 1    # exclude animals younger than 1 year
MIN_PATHWAY_GENES = 5    # min detected genes to include a pathway in scoring

# ── OPC/Oligo split ───────────────────────────────────────────────────────
# within the opc-olig h5ad, louvains 12 and 13 are OPCs, all others are oligodendrocytes
OPC_LOUVAINS = ['12', '13']

# ── output paths ──────────────────────────────────────────────────────────
MITOPPS_PATH    = f"{OUTDIR}/all_mitopps.csv"
REGRESSION_PATH = f"{OUTDIR}/mitotyping_age_results_clean.csv"

# ── cell type → h5ad filename mapping ────────────────────────────────────
CELL_TYPES = {
    'Astrocyte':     'Res1_astrocytes_update.h5ad',
    'Basket':        'Res1_basket-cells_update.h5ad',
    'Cerebellar':    'Res1_cerebellar-neurons_subset.h5ad',
    'Ependymal':     'Res1_ependymal-cells_new.h5ad',
    'GABA':          'Res1_GABAergic-neurons_subset.h5ad',
    'Glutamatergic': 'Res1_glutamatergic-neurons_update.h5ad',
    'MSN':           'Res1_medium-spiny-neurons_subset.h5ad',
    'Microglia':     'Res1_microglia_new.h5ad',
    'Midbrain':      'Res1_midbrain-neurons_update.h5ad',
    'OPC_Oligo':     'Res1_opc-olig_subset.h5ad',
    'Vascular':      'Res1_vascular-cells_subset.h5ad',
}

os.makedirs(OUTDIR, exist_ok=True)


# ==========================================================================
# pipeline
# ==========================================================================
def run_mitotyping(cell_type, fname, ens2human, mito_genes, pathways):
    """
    Run full mitotyping pipeline for one cell type.

    Steps:
      1. Load h5ad
      2. Map Ensembl IDs → human gene symbols
      3. Subset to ~964 mitochondrial genes (drops full 29k matrix to save memory)
      4. Pseudobulk: mean expression per (animal × louvain × region)
      5. Filter: ≥ MIN_CELLS per group, ≥ MIN_ANIMALS per louvain-region combo
      6. MitoPathway scores: mean expression across pathway member genes
      7. mitoPPS: normalize each pathway against all others to remove
         the effect of overall mitochondrial abundance

    OPC_Oligo is automatically split:
      louvains 12, 13 → 'OPC'
      all others      → 'Oligodendrocyte'

    Note: .X is already log-normalized in these h5ads — no additional
    normalization is applied.

    Returns mitopps_df or None if skipped.
    Skips processing if output CSV already exists.
    """
    ct_safe  = cell_type.replace('/', '_')
    out_path = f"{OUTDIR}/{ct_safe}_mitopps.csv"

    if os.path.exists(out_path):
        print(f"  already done — loading {out_path}")
        return pd.read_csv(out_path)

    # ── load ──────────────────────────────────────────────────────────
    adata = ad.read_h5ad(f"{DATA_DIR}/{fname}")
    adata = adata[adata.obs['age'].astype(float) >= MIN_AGE].copy()
    print(f"  cells after age filter: {adata.n_obs}")

    # ── ortholog mapping ──────────────────────────────────────────────
    # rename var_names from Ensembl IDs to human gene symbols
    adata.var_names = pd.Series(adata.var_names).map(ens2human).fillna(
        pd.Series(adata.var_names)).values

    # ── subset to mito genes and free full matrix ─────────────────────
    # key memory optimization: extract only ~964 mito gene columns
    # then delete the full 29k × n_cells matrix immediately
    # NOTE: .X is already log-normalized — no normalize_total/log1p needed
    mito_genes_in_data = sorted(mito_genes & set(adata.var_names))
    assert len(mito_genes_in_data) > 500, \
        f"only {len(mito_genes_in_data)} mito genes found — check ortholog mapping"

    gene_idx = np.where(pd.Index(adata.var_names).isin(mito_genes_in_data))[0]
    X   = adata.X[:, gene_idx]
    if issparse(X): X = X.toarray()
    obs = adata.obs[['animal_id', 'louvain', 'region', 'age', 'sex']].copy()
    del adata; gc.collect()  # free full matrix — saves several GB
    print(f"  mito genes: {len(mito_genes_in_data)}, X: {X.shape}")

    # ── pseudobulk ────────────────────────────────────────────────────
    # mean expression per (animal × louvain × region) group
    obs['group'] = (obs['animal_id'].astype(str) + '__' +
                    obs['louvain'].astype(str)    + '__' +
                    obs['region'].astype(str))

    rows, meta_rows = [], []
    for grp, idx in obs.groupby('group').groups.items():
        if len(idx) < MIN_CELLS: continue
        animal, louvain, region = grp.split('__', 2)
        rows.append(X[obs.index.get_indexer(idx)].mean(axis=0))

        # split OPC_Oligo by louvain cluster label
        ct_label = ('OPC'              if cell_type == 'OPC_Oligo' and louvain in OPC_LOUVAINS
                    else 'Oligodendrocyte' if cell_type == 'OPC_Oligo'
                    else cell_type)

        meta_rows.append({'animal_id': animal, 'louvain': louvain,
                          'region': region, 'age': obs.loc[idx[0], 'age'],
                          'sex': obs.loc[idx[0], 'sex'], 'n_cells': len(idx),
                          'cell_type': ct_label})

    del X, obs; gc.collect()

    if not meta_rows:
        print("  no groups passed MIN_CELLS — skipping"); return None

    meta       = pd.DataFrame(meta_rows)
    pseudobulk = pd.concat(
        [meta, pd.DataFrame(np.vstack(rows), columns=mito_genes_in_data)], axis=1)
    del rows, meta_rows; gc.collect()

    # drop louvain-region combos with too few animals for reliable regression
    n_animals  = pseudobulk.groupby(['louvain','region'])['animal_id'].nunique()
    valid      = n_animals[n_animals >= MIN_ANIMALS].index
    keep       = pseudobulk.apply(lambda r: (r['louvain'], r['region']) in valid, axis=1)
    pseudobulk = pseudobulk[keep].reset_index(drop=True)
    print(f"  pseudobulk: {len(pseudobulk)} samples, {len(valid)} combos")
    if len(pseudobulk) == 0:
        print("  nothing passed filters — skipping"); return None

    # ── MitoPathway scores ────────────────────────────────────────────
    # for each Level-3 pathway: mean expression across detected member genes
    meta_cols    = ['animal_id','louvain','region','age','sex','n_cells','cell_type']
    pathways_use = pathways[pathways['level'] == 'Level3'].copy()
    pathways_use['matched_genes'] = pathways_use['gene_list'].apply(
        lambda gl: sorted(set(gl) & set(mito_genes_in_data)))
    pathways_use = pathways_use[
        pathways_use['matched_genes'].apply(len) >= MIN_PATHWAY_GENES]

    scores_dict = {}
    for _, row in pathways_use.iterrows():
        genes = [g for g in row['matched_genes'] if g in pseudobulk.columns]
        if genes:
            scores_dict[row['MitoPathway']] = pseudobulk[genes].mean(axis=1).values

    pathway_scores = pd.concat([
        pseudobulk[meta_cols].reset_index(drop=True),
        pd.DataFrame(scores_dict)
    ], axis=1)
    del pseudobulk; gc.collect()

    # ── mitoPPS ───────────────────────────────────────────────────────
    # normalize each pathway against all others to remove overall abundance effect
    # for sample Si and pathway Pi:
    #   mitoPPS(Si, Pi) = mean over Pn≠Pi of [score(Si,Pi)/score(Si,Pn)]
    #                                        / [mean over Sx of score(Sx,Pi)/score(Sx,Pn)]
    # value=1: average priority | >1: prioritized | <1: deprioritized
    pw_cols       = [c for c in pathway_scores.columns if c not in meta_cols]
    scores_matrix = pathway_scores[pw_cols].values.astype(float)
    n_s, n_p      = scores_matrix.shape
    print(f"  mitoPPS: {n_s} samples × {n_p} pathways...")

    mitopps_mat = np.zeros((n_s, n_p))
    for i in range(n_p):
        ratios = []
        for j in range(n_p):
            if i == j: continue
            r_ij = scores_matrix[:, i] / (scores_matrix[:, j] + 1e-10)
            ratios.append(r_ij / (r_ij.mean() + 1e-10))
        mitopps_mat[:, i] = np.mean(ratios, axis=0)

    mitopps_df = pd.concat([
        pathway_scores[meta_cols].reset_index(drop=True),
        pd.DataFrame(mitopps_mat, columns=pw_cols)
    ], axis=1)
    del scores_matrix, mitopps_mat; gc.collect()
    print(f"  mitoPPS range: {mitopps_df[pw_cols].values.min():.3f} – "
          f"{mitopps_df[pw_cols].values.max():.3f}")

    mitopps_df.to_csv(out_path, index=False)
    print(f"  saved: {ct_safe}_mitopps.csv")
    return mitopps_df


def run_age_regression(combined_mitopps, min_animals=MIN_ANIMALS):
    """
    Linear regression of mitoPPS ~ age + sex.

    Fits one model per pathway per (cell_type × louvain × region) combination.
    BH-FDR correction is applied across pathways within each combination.

    Parameters
    ----------
    combined_mitopps : concatenated output from run_mitotyping() for all cell types
    min_animals      : minimum unique animals required to run regression

    Returns
    -------
    results_df : long-format table with beta_age, se_age, pval_age, fdr_age
    """
    META    = {'animal_id','louvain','region','age','sex','n_cells','cell_type'}
    pw_cols = [c for c in combined_mitopps.columns if c not in META]
    groups  = ['cell_type','louvain','region'] if 'cell_type' in combined_mitopps.columns \
              else ['louvain','region']
    results = []

    for keys, d in combined_mitopps.groupby(groups):
        if d['animal_id'].nunique() < min_animals: continue
        ct, lou, reg = (keys if len(groups) == 3 else (None, keys[0], keys[1]))

        pvals, betas, ses = [], [], []
        for pw in pw_cols:
            try:
                fit = ols(f"Q('{pw}') ~ age + C(sex)", data=d).fit()
                betas.append(fit.params['age'])
                ses.append(fit.bse['age'])
                pvals.append(fit.pvalues['age'])
            except:
                betas.append(np.nan); ses.append(np.nan); pvals.append(np.nan)

        # BH-FDR across pathways within this combo
        mask = ~np.isnan(pvals)
        fdr  = np.full(len(pvals), np.nan)
        if mask.sum() > 0:
            fdr[mask] = multipletests(np.array(pvals)[mask], method='fdr_bh')[1]

        for k, pw in enumerate(pw_cols):
            results.append({'cell_type': ct, 'louvain': lou, 'region': reg,
                            'pathway': pw, 'beta_age': betas[k], 'se_age': ses[k],
                            'pval_age': pvals[k], 'fdr_age': fdr[k],
                            'n_animals': d['animal_id'].nunique()})

    results_df = pd.DataFrame(results)
    print(f"Total tests: {len(results_df)}")
    print(f"Significant (FDR<0.05): {(results_df['fdr_age'] < 0.05).sum()}")
    return results_df


def flag_bad_combos(combined_mitopps, corr_threshold=0.5):
    """
    QC: identify louvain × region × cell_type combos where mitoPPS values have high correlation 
    across all pathways for a given sample 

    Parameters
    ----------
    combined_mitopps : output of run pipeline loop
    corr_threshold   : mean pairwise correlation threshold (default 0.5)

    Returns
    -------
    bad_df : dataframe of flagged combos sorted by correlation
    """
    META    = {'animal_id','louvain','region','age','sex','n_cells','cell_type'}
    pw_cols = [c for c in combined_mitopps.columns if c not in META]
    bad     = []

    for keys, d in combined_mitopps.groupby(['cell_type','louvain','region']):
        corr     = d[pw_cols].corr().values
        mean_cor = corr[corr < 1].mean()
        if mean_cor > corr_threshold:
            bad.append({'cell_type': keys[0], 'louvain': keys[1], 'region': keys[2],
                        'mean_pathway_corr': round(mean_cor, 3),
                        'n_animals': d['animal_id'].nunique()})

    bad_df = pd.DataFrame(bad).sort_values('mean_pathway_corr', ascending=False) \
             if bad else pd.DataFrame()
    print(f"Bad combos (mean pathway corr > {corr_threshold}): {len(bad_df)}")
    if len(bad_df) > 0:
        print(bad_df.head(20))
    return bad_df
