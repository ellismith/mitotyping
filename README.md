# Mitotyping Analysis — U01 Macaque Brain Aging

Quantifying mitochondrial pathway prioritization (mitoPPS) and its relationship with aging across cell types and brain regions in the U01 rhesus macaque snRNA-seq dataset.

Based on [Monzel & Picard (2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11844627/).

---

## Methods summary

- **Data**: ~53 rhesus macaques, 11 brain regions, 11 cell types; snRNA-seq h5ads at `/data/CEM/smacklab/data/bican/u01_cell_class/intermediate_cell_class/`
- **Gene mapping**: macaque Ensembl IDs → human gene symbols via one-to-one orthologs; enables matching to MitoCarta 3.0
- **Pseudobulk**: mean log-normalized expression per animal × louvain × region; filters: ≥100 cells/group, ≥10 animals/combo; yields **9,024 pseudobulk samples**
- **MitoPathway scores**: mean expression across member genes; 69 Level-3 pathways retained (≥5 detected genes; 34 L3 dropped, 46 L1/L2 excluded by design)
- **mitoPPS**: pairwise ratio normalization removes overall mitochondrial abundance effect; value=1 is average, >1 prioritized, <1 deprioritized
- **QC**: 41 outlier animal × subtype × region pairs identified (mean mitoPPS >5× within-combo median); sensitivity analysis confirmed 0 significant hits lost on removal
- **Age regression**: OLS `mitoPPS ~ age + sex` per pathway per combination; **21,459 tests** total
- **Multiple testing**: BH-FDR within each cell type × louvain × region combination; **78 hits at FDR<0.05**; 2 survive global FDR correction

---

## Key results

- **78 significant pathway-age associations** (within-combo FDR<0.05)
- **2 survive global FDR correction**: CI subunits and Fe-S-containing proteins, both in Cerebellar lCb
- **Most hits by cell type**: Cerebellar (30), Glutamatergic (20), Basket (11), Astrocyte (7)
- **OPC**: 0 significant hits — biological finding, not artifact
- **Consistent decreases with age**: CI/CIV subunits, mtDNA repair, mitochondrial ribosome
- **Consistent increases**: Selenoproteins, NAD biosynthesis, ROS/glutathione metabolism

---

## Cell type abbreviations

| Abbreviation | Full name |
|---|---|
| AST | Astrocyte |
| BC | Basket cells |
| CER | Cerebellar neurons |
| EXC | Glutamatergic neurons |
| INH | GABAergic neurons |
| MBN | Midbrain neurons |
| MGL | Microglia |
| MSN | Medium spiny neurons |
| OPC | OPC |
| OLIG | Oligodendrocyte |
| VASC | Vascular cells |

---

## Repository structure

```
mitotyping/
├── mitotyping_analysis.py     # core pipeline: run_mitotyping(), run_age_regression()
├── mitotyping_plots.py        # all visualization functions
├── mitotyping_u01.py          # main analysis notebook (jupytext --to notebook)
├── mitotyping_outliers.py     # outlier QC functions
├── generate_explorer.py       # generates interactive HTML explorer from results
├── explorer_template.html     # HTML/JS template for explorer (edit this, not the output)
├── run_all.py                 # orchestrates full pipeline
├── 01_run_pipeline.py         # step 1: pseudobulk + mitoPPS
├── 02_age_regression.py       # step 2: OLS regression
├── 03_outlier_qc.py           # step 3: QC and sensitivity analysis
├── 04_visualization.py        # step 4: generate all figures
├── replot_aging.py            # regenerate aging pathway figures
├── plot_per_celltype.py       # per cell type heatmaps
├── plot_by_region_abbr.py     # per region heatmaps
├── plot_data_summary.py       # data availability figure
├── pathway_gene_table.py      # MitoCarta gene detection summary
└── README.md
```

---

## Usage

### Environment setup

```bash
module load mamba/latest
mamba create -n mito_env_mamba python=3.10 -y
source activate mito_env_mamba
pip install scanpy anndata pandas numpy scipy openpyxl xlrd statsmodels seaborn matplotlib jupytext scikit-learn
```

### Run full pipeline

```bash
conda run -n mito_env_test python3 run_all.py
```

### Convert notebook

```bash
jupytext --to notebook mitotyping_u01.py
```

### Generate interactive explorer

```bash
conda run -n mito_env_test python3 generate_explorer.py
# output: /scratch/easmit31/mitotyping/figures/mitotyping_explorer_MMDD_HHMM.html
# open in any browser — no internet required
```

---

## Data paths (ASU Sol cluster)

| Item | Path |
|---|---|
| h5ad files | `/data/CEM/smacklab/data/bican/u01_cell_class/intermediate_cell_class/` |
| Ortholog mapping | `/scratch/easmit31/data/human-macaque-orthologs/ensembl113_mmul10_macaque_human.csv` |
| MitoCarta 3.0 | `/scratch/easmit31/mitotyping/Human.MitoCarta3.0.xls` |
| Results | `/scratch/easmit31/mitotyping/results/` |
| Figures | `/scratch/easmit31/mitotyping/figures/` |

---

## Results files

| File | Description |
|---|---|
| `results/all_mitopps.csv` | 9,024 pseudobulk samples × 69 mitoPPS values + metadata |
| `results/mitotyping_age_results_clean.csv` | 21,459 regression tests (beta, SE, p-value, FDR) |
| `results/mitotyping_age_results_sensitivity.csv` | Same after removing 41 outlier pairs (79 sig hits) |
| `results/pathway_gene_counts.csv` | MitoCarta pathway gene detection summary |

---

## Reference

Monzel AS, Picard M. (2025). Mitochondrial pathway prioritization scores reveal cell-type-specific mitotypes. *PNAS*. https://pmc.ncbi.nlm.nih.gov/articles/PMC11844627/
