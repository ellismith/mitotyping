# Mitotyping Analysis — U01 Macaque Brain Aging

Quantifying mitochondrial pathway prioritization (mitoPPS) and its relationship with aging across cell types and brain regions in the U01 rhesus macaque snRNA-seq dataset.

Based on the mitotyping framework from [Monzel & Picard (2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11844627/).

---

## Overview

Mitochondria are not generic energy factories — each cell type maintains a distinct mitochondrial program (a **mitotype**) reflecting its metabolic priorities. This pipeline quantifies mitotypes from single-nucleus RNA-seq data and asks which mitochondrial pathways change with age, in which cell types, and in which brain regions.

**Dataset:** U01 macaque brain aging cohort
- ~53 rhesus macaques, ages 2–21 years
- 11 brain regions (ACC, CN, dlPFC, EC, HIP, IPP, lCb, M1, MB, mdTN, NAc)
- 11 cell types (astrocytes, basket cells, cerebellar neurons, GABA neurons, glutamatergic neurons, medium spiny neurons, microglia, midbrain neurons, OPC, oligodendrocytes, vascular cells)

---

## Pipeline

### 1. Load and map genes
Single-nucleus h5ad files contain macaque Ensembl IDs. These are mapped to human gene symbols using one-to-one orthologs, enabling matching to MitoCarta 3.0.

### 2. Pseudobulk
Cells are collapsed to mean expression per **(animal × louvain subtype × brain region)** group.

Filters:
- ≥ 100 cells per group
- ≥ 10 unique animals per louvain × region combination
- Animals ≥ 1 year old

### 3. MitoPathway scores
For each of 149 Level-3 MitoCarta pathways, pathway score = mean log-normalized expression of detected member genes (≥ 5 genes required).

### 4. mitoPPS
Normalizes each pathway score against all others to remove the effect of overall mitochondrial abundance.

\`\`\`
mitoPPS(Si, Pi) = mean over Pn≠Pi of:
    [score(Si,Pi) / score(Si,Pn)] / [mean over Sx of score(Sx,Pi)/score(Sx,Pn)]
\`\`\`

- **mitoPPS = 1** → average priority
- **mitoPPS > 1** → pathway prioritized relative to others
- **mitoPPS < 1** → pathway deprioritized

### 5. Age regression
\`\`\`
mitoPPS ~ age + sex
\`\`\`
OLS regression per pathway per combination. BH-FDR correction within each louvain × region × cell type combination across 149 pathways.

### 6. QC
Outlier animal × subtype × region combinations identified where mitoPPS values exceed 5× the within-combo median. Sensitivity analysis confirms findings are robust to their exclusion.

---

## Key findings

- **78 significant pathway-age associations** (within-combo FDR < 0.05)
- **2 survive global FDR correction**: CI subunits and Fe-S-containing proteins in cerebellar neurons
- **Consistent decreases with age**: CI subunits, CIV subunits, mtDNA repair, mitochondrial ribosome
- **Consistent increases with age**: Selenoproteins, NAD biosynthesis, ROS/glutathione metabolism
- **Most hits**: cerebellar neurons (30), basket cells (11), glutamatergic neurons (20)
- **OPCs**: no significant effects at FDR < 0.05

---

## Repository structure

\`\`\`
mitotyping/
├── mitotyping_u01.py          # main pipeline notebook (convert with jupytext)
├── mitotyping_plots.py        # all visualization functions
├── mitotyping_regression.py   # standalone age regression script
└── README.md
\`\`\`

---

## Requirements

\`\`\`bash
module load mamba/latest
mamba create -n mito_env_mamba python=3.10 -y
source activate mito_env_mamba
pip install scanpy anndata pandas numpy scipy openpyxl xlrd statsmodels seaborn matplotlib jupytext scikit-learn
\`\`\`

Convert notebook:
\`\`\`bash
jupytext --to notebook mitotyping_u01.py
\`\`\`

---

## Data

- **h5ads**: \`/data/CEM/smacklab/data/bican/u01_cell_class/intermediate_cell_class/\`
- **Orthologs**: \`/scratch/easmit31/data/human-macaque-orthologs/ensembl113_mmul10_macaque_human.csv\`
- **MitoCarta 3.0**: [Download here](https://www.broadinstitute.org/mitocarta/mitocarta30-inventory-mammalian-mitochondrial-proteins-and-pathways)

---

## References

Monzel AS, Picard M. (2025). Mitochondrial pathway prioritization scores reveal cell-type-specific mitotypes. *PNAS*. https://pmc.ncbi.nlm.nih.gov/articles/PMC11844627/
