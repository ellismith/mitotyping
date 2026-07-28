# Mitotyping Analysis Guide
### U01 Macaque Brain Aging — mitochondrial pathway prioritization

This guide walks through the full analysis pipeline, explains what each script does,
where results are saved, and how to use each visualization function.
It's written for someone who may not be deeply familiar with the code,
but wants to understand the logic and know where to look when things need changing.

---

## What is mitotyping?

Mitochondria are not generic — each cell type has a distinct **mitotype**: a characteristic
pattern of which mitochondrial pathways it invests in. An astrocyte invests heavily in fatty
acid oxidation; a neuron prioritizes OXPHOS. These patterns shift with aging.

The **mitoPPS** (mitochondrial pathway prioritization score) captures this specialization.
Instead of asking "how much is this pathway expressed?", it asks "is this pathway
*more expressed than expected* relative to all other mitochondrial pathways in this cell?"

A mitoPPS of **1.0** = average priority.
A mitoPPS **> 1** = this pathway is prioritized (invested in more than others).
A mitoPPS **< 1** = this pathway is deprioritized.

The method comes from [Monzel & Picard (2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11844627/).

---

## Dataset at a glance

| Feature | Detail |
|---|---|
| Species | Rhesus macaque (*Macaca mulatta*) |
| N animals | ~53 (ages 2–21 years) |
| Brain regions | 11 (ACC, CN, dlPFC, EC, HIP, IPP, lCb, M1, MB, mdTN, NAc) |
| Cell types | 11 (see abbreviation table below) |
| Pathways analyzed | 69 Level-3 MitoCarta pathways |
| Total regression tests | 21,459 |
| Significant hits (FDR<0.05) | 78 |

### Cell type abbreviations

| Abbreviation | Full name | Notes |
|---|---|---|
| AST | Astrocyte | All 11 regions |
| BC | Basket cells | lCb only |
| CER | Cerebellar neurons | lCb only |
| EXC | Glutamatergic neurons | Most regions except lCb, MB |
| INH | GABAergic neurons | Most regions except lCb, MB |
| MBN | Midbrain neurons | MB only |
| MGL | Microglia | All 11 regions |
| MSN | Medium spiny neurons | CN and NAc only |
| OPC | OPC | All 11 regions |
| OLIG | Oligodendrocyte | All 11 regions |
| VASC | Vascular cells | Most regions |

Cell type × region combinations with fewer than 100 cells per animal or fewer
than 10 animals are excluded. **The only non-biological exclusion** is EXC in CN
(only 4/53 animals have ≥100 glutamatergic neurons there).

---

## Repository structure

```
mitotyping/
├── mitotyping_analysis.py     # pipeline: pseudobulk, pathway scoring, mitoPPS, regression
├── mitotyping_plots.py        # all visualization functions
├── mitotyping_u01.py          # main analysis notebook (run as .py or convert with jupytext)
├── mitotyping_outliers.py     # outlier detection and QC
├── generate_explorer.py       # builds the interactive HTML explorer
├── explorer_template.html     # HTML/CSS/JS template (edit this, not the generated output)
├── run_all.py                 # runs the full pipeline end to end
├── 01_run_pipeline.py         # step 1: pseudobulk + mitoPPS
├── 02_age_regression.py       # step 2: OLS regression
├── 03_outlier_qc.py           # step 3: QC and sensitivity analysis
├── 04_visualization.py        # step 4: all figure generation
└── README.md
```

---

## Environment setup

Run once in terminal before starting:

```bash
module load mamba/latest
mamba create -n mito_env_mamba python=3.10 -y
source activate mito_env_mamba
pip install scanpy anndata pandas numpy scipy openpyxl xlrd statsmodels \
            seaborn matplotlib jupytext scikit-learn
```

To convert the `.py` notebook to Jupyter format:

```bash
jupytext --to notebook mitotyping_u01.py
```

---

## Data paths

All paths live in `mitotyping_analysis.py` — **change them there**, not scattered across scripts.

```python
DATA_DIR      = "/data/CEM/smacklab/data/bican/u01_cell_class/intermediate_cell_class/"
MITOCARTA_PATH = "/scratch/easmit31/mitotyping/Human.MitoCarta3.0.xls"
ORTHOLOG_CSV  = "/scratch/easmit31/data/human-macaque-orthologs/ensembl113_mmul10_macaque_human.csv"
OUTDIR        = "/scratch/easmit31/mitotyping/results/"
```

---

## Pipeline walkthrough

The pipeline runs in four steps. Each step saves a CSV to `results/`.

### Step 1 — Pseudobulk

**Script:** `01_run_pipeline.py` → calls `mta.run_mitotyping()`

**What it does:**

For each cell type h5ad file:
1. Reads only metadata in `backed='r'` mode (no RAM cost)
2. Maps macaque Ensembl IDs → human gene symbols using one-to-one orthologs
3. Subsets to ~964 MitoCarta genes, discards the rest immediately
4. **Pseudobulks**: collapses all cells → mean log-normalized expression
   per `(animal × louvain × region)` group
5. Applies filters:

```python
MIN_CELLS   = 100    # minimum cells per pseudobulk group
MIN_ANIMALS = 10     # minimum animals per louvain×region combo
MIN_AGE     = 1      # exclude animals younger than 1 year
```

6. Scores each of the 69 qualifying MitoCarta Level-3 pathways as
   the mean expression of its member genes
7. Computes **mitoPPS** via pairwise ratio normalization

**The mitoPPS formula** (simplified):

```
mitoPPS(sample S, pathway P) =
    mean over all other pathways Q:
        [score(S, P) / score(S, Q)]     ← this sample's ratio
        ─────────────────────────────────
        [mean over all samples of that ratio]  ← expected ratio
```

This removes the effect of overall mitochondrial abundance — a sample
with globally high mitochondrial expression doesn't get artificially high
scores for all pathways.

**Where cell types are defined:**

In `mitotyping_analysis.py`, the `CELL_TYPES` dict maps each cell type name
to its h5ad filename. **This is where you add or remove cell types.**

```python
CELL_TYPES = {
    'Astrocyte':      'Res1_astrocytes_update.h5ad',
    'Basket':         'Res1_basket-cells_update.h5ad',
    'Cerebellar':     'Res1_cerebellar-neurons_subset.h5ad',
    # ... etc
}
```

The OPC/Oligodendrocyte split is also defined here:

```python
OPC_LOUVAINS = {'12', '13'}   # louvains 12+13 → OPC, all others → Oligodendrocyte
```

**Output:** `results/all_mitopps.csv`

```
Shape: 9,024 rows × 76 columns
Rows:  one per (animal × louvain × region × cell type) pseudobulk sample
Cols:  animal_id, louvain, region, age, sex, n_cells, cell_type, [69 pathway scores]
```

---

### Step 2 — Age regression

**Script:** `02_age_regression.py` → calls `mta.run_age_regression()`

**What it does:**

For each (cell type × louvain × region) combination, fits:

```
mitoPPS ~ age + C(sex)
```

using OLS (ordinary least squares). This gives a `beta_age` coefficient — how
much the mitoPPS of that pathway changes per year of age — and a p-value.

**Multiple testing correction:** BH-FDR is applied *within* each
(cell type × louvain × region) combination, across the 69 pathways.
This means each combination has its own corrected p-values.

**Output:** `results/mitotyping_age_results_clean.csv`

```
Shape: 21,459 rows × 10 columns
Cols:  cell_type, louvain, region, pathway, beta_age, se_age, pval_age, fdr_age, n_animals
```

---

### Step 3 — Outlier QC

**Script:** `03_outlier_qc.py` → calls `mto.flag_bad_combos()`

**What it does:**

Identifies animal × combo pairs where the mean mitoPPS across all pathways
exceeds 5× the within-combo median. These are cases where the mitoPPS
normalization likely failed (e.g. one pathway near zero makes all ratios
blow up).

**Result:** 41 outlier pairs found. Sensitivity analysis (removing all 41)
showed **0 significant hits lost** and **1 gained**. We keep the raw data
for the primary analysis and report the sensitivity results separately.

**Output:**
- `results/mitotyping_age_results_sensitivity.csv` — regression without outliers
- Outlier report printed to stdout

---

### Step 4 — Visualization

**Script:** `04_visualization.py`

Calls all plotting functions in `mitotyping_plots.py` and saves figures to
`figures/`. See the plotting section below for details on each function.

---

## Plotting functions

All functions live in `mitotyping_plots.py`. Import once:

```python
import importlib
import mitotyping_plots as mp
importlib.reload(mp)   # pick up any edits without restarting
```

All functions accept `save=True` to write a 300dpi PNG to `figures/`.
Filenames are auto-generated from the plot title.

---

### `mp.plot_mitotyping()` — the main heatmap

This is the workhorse. It handles both mitoPPS values and age effects,
at any pathway level, with many layout options.

```python
mp.plot_mitotyping(
    df,                        # combined_mitopps (mode=mitopps) or results_df (mode=age)
    mode        = "mitopps",   # "mitopps" | "age"
    level       = "L3",        # "L1" | "L2" | "L3"  — pathway aggregation level
    cell_type   = None,        # filter to one cell type, e.g. "Astrocyte"
    category    = None,        # filter to one L1 category, e.g. "OXPHOS"
    average_louvains = False,  # one column per cell_type×region (averages across louvains)
    per_animal  = False,       # one column per animal, ordered youngest→oldest
    by_region   = False,       # one column per region (age mode)
    facet_region = False,      # one heatmap per region (age mode)
    only_sig    = True,        # show only pathways with ≥1 significant hit (age mode)
    fdr_threshold = 0.05,      # significance threshold for dots
    cluster_rows = False,      # hierarchically cluster pathway rows
    cluster_cols = False,      # hierarchically cluster columns
    vmax        = None,        # color scale max (e.g. 0.03); None = auto from data
    save        = False,       # save PNG to figures/
)
```

**How the heatmap is constructed** (simplified):

```python
# 1. filter the dataframe to the requested cell type / category
df = df[df['cell_type'] == cell_type]

# 2. build a label for each column — e.g. "Astrocyte, ACC"
df['combo'] = df['cell_type'] + ', ' + df['region']

# 3. pivot: rows = pathways, columns = combos, values = mitoPPS or beta_age
mat = df.pivot_table(index='pathway', columns='combo',
                     values='beta_age', aggfunc='mean')

# 4. optionally aggregate to L1/L2
mat.index = [L1_map.get(p, p) for p in mat.index]
mat = mat.groupby(level=0).mean()

# 5. optionally cluster rows
row_order = leaves_list(linkage(pdist(mat.values), method='ward'))
mat = mat.iloc[row_order]

# 6. draw with matplotlib — color scale centered at 0 for age effects
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
ax.imshow(mat.values, cmap='RdBu_r', norm=norm)

# 7. add significance dots
for each cell where fdr < threshold:
    ax.plot(col, row, 'k.', markersize=3)
```

**Common usage patterns:**

```python
# overview: all cell types, age effects, louvains averaged, vmax=0.03
mp.plot_mitotyping(results_df, mode='age', average_louvains=True,
                   only_sig=True, vmax=0.03, save=True)

# one heatmap per region, cols = cell types
mp.plot_mitotyping(results_df, mode='age', facet_region=True,
                   only_sig=False, vmax=0.03, save=True)

# mitoPPS for one cell type, per animal ordered by age
mp.plot_mitotyping(combined_mitopps, cell_type='Astrocyte',
                   per_animal=True, save=True)

# OXPHOS pathways only
mp.plot_mitotyping(results_df, mode='age', category='OXPHOS',
                   average_louvains=True, save=True)
```

**Where regions and cell types are looped over:**

In `04_visualization.py`, the per-cell-type plots are generated by:

```python
cell_types = sorted(combined_mitopps['cell_type'].unique())
for ct in cell_types:
    mp.plot_mitotyping(combined_mitopps, cell_type=ct, per_animal=True, save=True)
```

To add or exclude a cell type, edit the `CELL_TYPES` dict in `mitotyping_analysis.py`.
The loop automatically picks up whatever is in the data.

---

### `mp.plot_aging_pathways()` — top aging pathways + scatter

Two plots for one cell type at a time:
1. Heatmap of the top aging pathways with animals ordered youngest→oldest
2. Scatter plots of individual pathway scores vs age with regression lines

```python
mp.plot_aging_pathways(
    results_df,          # regression results
    combined_mitopps,    # raw mitoPPS values (for scatter dots)
    cell_type = "Astrocyte",
    top_n     = 10,      # n most negative + n most positive betas
    fdr_threshold = 0.05,
    save      = False,
)
```

**How pathway selection works:**
- If any pathways survive FDR<0.05 for this cell type, those are shown
- If none survive (e.g. Microglia), falls back to the `top_n` pathways
  by most extreme effect size (positive and negative)
- Scatter plots show the first 6 pathways selected

**Colorscale:** centered at mitoPPS=1 (biological neutral). White = average
priority, red = prioritized, blue = deprioritized.

---

### `mp.summarize_mitotyping()` — direction + consistency

```python
direction, consistency_mat = mp.summarize_mitotyping(results_df, fdr_threshold=0.05)
```

Returns:
- **direction**: DataFrame with n_sig, n_increasing, n_decreasing per cell type
- **consistency_mat**: pathways × cell types matrix of fraction of regions significant

---

### `mp.plot_direction()` — direction bar chart

```python
mp.plot_direction(direction, save=True)
```

Horizontal bar chart: each cell type has a red bar (pathways increasing with age)
and a blue bar (pathways decreasing). Shows both magnitude and direction of aging effects.

---

### `mp.plot_consistency()` — consistency heatmap

```python
mp.plot_consistency(consistency_mat, save=True)

# show only pathways significant in ≥50% of regions for at least one cell type
mp.plot_consistency(consistency_mat, min_frac=0.5, save=True)
```

Heatmap where each value = fraction of regions where that pathway is significant
for that cell type. High values = robust finding across regions.

---

### `mp.plot_significant_dotplot()` — significant hits overview

```python
mp.plot_significant_dotplot(results_df, results_sens, fdr_threshold=0.05, save=True)
```

Dot plot grouped by region → pathway. Each dot = one significant
(cell type × louvain × region × pathway) hit. Dot size = -log10(FDR),
color = cell type. Stars = hits that only appear in the sensitivity analysis
(without outliers).

---

### `mp.plot_mitopps_pca()` — PCA of mitoPPS

```python
coords, pca_obj, pca_df = mp.plot_mitopps_pca(
    combined_mitopps,
    color_by = 'cell_type',   # 'cell_type' | 'region' | 'age'
    pc_x = 1, pc_y = 2,       # which PCs to plot
    save = False,
)
```

Each dot is one pseudobulk sample. Tight clusters by cell type confirm that
cell types have distinct mitotypes. Color by age to see if there's a
directional shift.

---

### `mp.plot_pathway_rank()` — pathway ranking

```python
mp.plot_pathway_rank(
    combined_mitopps,
    cell_type_a = "Astrocyte",
    cell_type_b = None,        # None = compare against all other cell types averaged
    top_n       = 20,
    save        = False,
)
```

Ranked bar chart: which pathways does cell type A invest in more/less than
cell type B? Bars colored by L1 category.

---

### `mp.plot_pathway_bivariate()` — two-pathway scatter

```python
mp.plot_pathway_bivariate(
    combined_mitopps,
    pathway_x = "CI subunits",
    pathway_y = "Fatty acid oxidation",
    color_by  = 'cell_type',
    save      = False,
)
```

Scatter of two pathway scores across all pseudobulk samples.

---

### `mp.plot_pathway_bubbles()` — bubble plot

```python
mp.plot_pathway_bubbles(
    combined_mitopps,
    pathways_to_show = ["CI subunits", "TCA cycle", "mtDNA repair", ...],
    n_top = 10,   # if pathways_to_show=None, show top n by variance
    save  = False,
)
```

Compact matrix view: pathways × cell types, bubble size = z-score magnitude,
color = direction.

---

## Results files

| File | Description |
|---|---|
| `results/all_mitopps.csv` | 9,024 pseudobulk samples × 69 mitoPPS values + metadata |
| `results/mitotyping_age_results_clean.csv` | 21,459 regression tests — the main results |
| `results/mitotyping_age_results_sensitivity.csv` | Same after removing 41 outlier pairs |
| `results/pathway_gene_counts.csv` | Which pathways were included/excluded and why |

---

## Interactive explorer

The HTML explorer embeds all data and runs in any browser — no internet required.

**To regenerate after updating results:**

```bash
conda run -n mito_env_test python3 generate_explorer.py
# saves: figures/mitotyping_explorer_MMDD_HHMM.html
```

**To change layout or controls:**
Edit `explorer_template.html` — this is the source. **Do not edit** the generated
HTML output directly; it will be overwritten next time you regenerate.

**To download:** use the OOD file browser at `ood04.sol.rc.asu.edu`,
navigate to `figures/`, and click Download.

---

## Key decisions and their rationale

| Decision | Rationale |
|---|---|
| Use Level-3 pathways only | L1/L2 are groupings of L3 — scoring them would double-count |
| MIN_PATHWAY_GENES=5 | Fewer than 5 genes → mean too noisy to be reliable |
| MIN_CELLS=100 | Ensures pseudobulk mean is stable |
| MIN_ANIMALS=10 | Need enough samples to run a regression |
| BH-FDR within combo, not global | Global correction across 21k tests is too conservative |
| Keep outliers in primary analysis | Sensitivity analysis confirmed 0 hits lost on removal |
| `.X` not re-normalized | Data is already log-normalized in the h5ads |

---

## Common questions

**Why are some cell types missing from certain regions?**
Biologically, some cell types are region-specific (basket cells only in cerebellum,
midbrain neurons only in MB). Others are excluded because not enough animals
have ≥100 cells in that region. See `results/pathway_gene_counts.csv` and the
data availability figure.

**Why does OPC have 0 significant hits?**
OPC is present in all 11 regions with good sample sizes. The null result is real —
OPC mitotypes appear stable with age in this dataset.

**What does a negative beta_age mean?**
The pathway becomes less prioritized with age — the cell type invests less in
it relative to other pathways as animals get older.

**How do I filter to look at one region?**
```python
# filter before plotting
sub = results_df[results_df['region'] == 'lCb']
mp.plot_mitotyping(sub, mode='age', average_louvains=True, save=True)
```

**How do I search for a pathway by keyword?**
```python
pw_cols = [c for c in combined_mitopps.columns
           if c not in {'animal_id','louvain','region','age','sex','n_cells','cell_type'}]
matches = [p for p in pw_cols if 'tca' in p.lower()]
print(matches)
```

**How do I add a new cell type?**
1. Add its h5ad filename to `CELL_TYPES` in `mitotyping_analysis.py`
2. Re-run `01_run_pipeline.py` to generate new mitoPPS data
3. Re-run `02_age_regression.py` for regression results
4. Re-run `04_visualization.py` for figures
5. Re-run `generate_explorer.py` for the interactive explorer
