"""
01_run_pipeline.py
Load MitoCarta + ortholog map + pathway hierarchy, then run the mitotyping
pipeline (pseudobulk -> MitoPathway scores -> mitoPPS) for every cell type
in CELL_TYPES. Writes combined_mitopps.csv to OUTDIR.

Usage:
    python 01_run_pipeline.py
"""

import gc
import os

import anndata as ad
import pandas as pd

from mitotyping_analysis import (
    DATA_DIR, MITOCARTA_PATH, OUTDIR, ORTHOLOG_CSV,
    MITOPPS_PATH, CELL_TYPES,
    run_mitotyping,
)


def preview_one_h5ad(fname="Res1_astrocytes_update.h5ad"):
    """Quick sanity check on one h5ad before running the full loop."""
    adata = ad.read_h5ad(f"{DATA_DIR}/{fname}", backed='r')
    print("Shape (cells x genes):", adata.shape)
    print("Obs columns:", adata.obs.columns.tolist())
    print("Layers:", list(adata.layers.keys()))
    print("N animals:", adata.obs['animal_id'].nunique())
    print("N regions:", adata.obs['region'].nunique())
    print("N louvains:", adata.obs['louvain'].nunique())
    print(f"Age range: {adata.obs['age'].min():.1f} - {adata.obs['age'].max():.1f}")
    print("First 5 var_names (should be Ensembl IDs):", list(adata.var_names[:5]))
    adata.file.close()
    del adata


def load_ortholog_map():
    mart = pd.read_csv(ORTHOLOG_CSV)
    mart = mart[mart['Human homology type'] == 'ortholog_one2one']
    ens2human = dict(zip(mart['Gene stable ID'], mart['Human gene name']))
    print(f"Orthologs loaded: {len(ens2human)} one-to-one mappings")
    return ens2human


def load_mito_genes():
    mito = pd.read_excel(MITOCARTA_PATH, sheet_name="A Human MitoCarta3.0")
    mito_genes = set(mito['Symbol'].dropna())
    print(f"MitoCarta genes: {len(mito_genes)}")
    return mito_genes


def load_pathways():
    pathways = pd.read_excel(MITOCARTA_PATH, sheet_name="C MitoPathways")
    pathways = pathways.dropna(subset=['Genes']).reset_index(drop=True)
    pathways['gene_list'] = pathways['Genes'].str.split(', ')

    hierarchy = pathways['MitoPathways Hierarchy'].str.split(' > ', expand=True)
    hierarchy.columns = ['L1', 'L2', 'L3'][:hierarchy.shape[1]]
    pathways = pd.concat([pathways, hierarchy], axis=1)

    def get_level(row):
        if row['MitoPathway'] == row.get('L1'):
            return 'Level1'
        if row['MitoPathway'] == row.get('L2'):
            return 'Level2'
        if row['MitoPathway'] == row.get('L3'):
            return 'Level3'
        return 'Unknown'

    pathways['level'] = pathways.apply(get_level, axis=1)
    print(f"Pathways loaded: {len(pathways)}")
    print("Pathway counts per level:")
    print(pathways['level'].value_counts())
    return pathways


def run_all_cell_types(ens2human, mito_genes, pathways):
    all_mitopps = []
    for cell_type, fname in CELL_TYPES.items():
        ct_name = cell_type.replace('/', '_')
        out_path = f"{OUTDIR}/{ct_name}_mitopps.csv"

        if os.path.exists(out_path):
            print(f"{cell_type}: loading existing CSV")
            all_mitopps.append(pd.read_csv(out_path))
        else:
            print(f"\n{'='*60}\n{cell_type}\n{'='*60}")
            result = run_mitotyping(cell_type, fname, ens2human, mito_genes, pathways)
            if result is not None:
                all_mitopps.append(result)
            gc.collect()

    combined_mitopps = pd.concat(all_mitopps, ignore_index=True)
    combined_mitopps.to_csv(MITOPPS_PATH, index=False)
    print(f"\nDone. {len(combined_mitopps)} pseudobulk samples")
    print(f"Cell types: {sorted(combined_mitopps['cell_type'].unique())}")
    return combined_mitopps


if __name__ == "__main__":
    preview_one_h5ad()
    ens2human = load_ortholog_map()
    mito_genes = load_mito_genes()
    pathways = load_pathways()
    combined_mitopps = run_all_cell_types(ens2human, mito_genes, pathways)
    print(f"\nSaved combined mitoPPS table to {MITOPPS_PATH}")
