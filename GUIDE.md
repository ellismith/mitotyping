# Mitotyping Analysis Guide
### Macaque Brain Aging — mitochondrial pathway prioritization

## What is mitotyping?

Mitochondria aren't just energy factories — they're specialized. Each cell type has
a distinct **mitotype**: a characteristic pattern of which mitochondrial pathways it
invests in most heavily, that reflects its metabolic priorities. These differences are encoded in 
gene expression — specifically in which mitochondrial pathway genes a cell transcribes more or less of.

The **mitoPPS** (mitochondrial pathway prioritization score) captures this specialization.
Rather than asking "how much is this pathway expressed?", it asks: *"is this pathway
more expressed than expected, relative to all other mitochondrial pathways in this
cell at this moment?"*

A mitoPPS of **1.0** = average priority (neither up nor down).
A mitoPPS **> 1** = this pathway is prioritized — the cell invests in it more than
in other pathways on average.
A mitoPPS **< 1** = this pathway is deprioritized relative to others.

This normalization is important because it removes the effect of overall mitochondrial
mass or activity — a cell that just happens to have more mitochondria overall doesn't
get artificially high scores for every pathway.

The method comes from [Monzel & Picard (2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11844627/).

MitoCarta 3.0 is a curated database of 1,136 human mitochondrial genes organized into 149 functional 
pathways across 7 broad categories. This gives us a vocabulary for describing mitotypes quantitatively.

---

## Dataset at a glance

| Feature | Detail |
|---|---|
| Species | Rhesus macaque (*Macaca mulatta*) |
| N animals | 53 (ages 3.2–20.8 years) |
| Brain regions | 11 (ACC, CN, dlPFC, EC, HIP, IPP, lCb, M1, MB, mdTN, NAc) |
| Cell types analyzed | 11 |
| Pseudobulk samples | 9,024 (animal × louvain × region groups passing QC) |
| MitoCarta pathways analyzed | 69 Level-3 pathways |
| Total regression tests | 21,459 |
| Significant hits (FDR < 0.05) | 78 |

### Cell type abbreviations and region coverage

The pipeline uses the full cell type names internally. Abbreviations appear in
plots and the interactive explorer.

| Full name (in data) | Abbreviation | Regions with data | N animals (range) |
|---|---|---|---|
| Astrocyte | AST | All 11 regions | 34–50 |
| Basket | BC | lCb only | 43 |
| Cerebellar | CER | lCb only | 46 |
| Glutamatergic | EXC | ACC, dlPFC, EC, HIP, IPP, M1, mdTN | 47–52 |
| GABA | INH | ACC, CN, IPP, M1, NAc, dlPFC, mdTN | 22–52 |
| Midbrain | MBN | MB only | 44 |
| Microglia | MGL | ACC, HIP, IPP, MB, mdTN, NAc | 25–34 |
| MSN | MSN | CN, NAc only | 41–50 |
| OPC | OPC | All 11 regions | 31–51 |
| Oligodendrocyte | OLIG | All 11 regions | 35–52 |
| Vascular | VASC | ACC, IPP, M1, MB, dlPFC, lCb, mdTN | 11–32 |

---

## Repository structure

```
mitotyping/
├── mitotyping_analysis.py     # pipeline functions: pseudobulk, scoring, mitoPPS, regression
├── mitotyping_plots.py        # all visualization functions
├── mitotyping_u01.py          # main analysis notebook (convert with jupytext)
├── mitotyping_outliers.py     # outlier detection and QC
├── generate_explorer.py       # builds the self-contained interactive HTML explorer
├── explorer_template.html     # HTML/CSS/JS template — edit this, not the output file
├── run_all.py                 # runs the full pipeline end to end
├── 01_run_pipeline.py         # step 1: pseudobulk + mitoPPS
├── 02_age_regression.py       # step 2: OLS age regression
├── 03_outlier_qc.py           # step 3: outlier detection + sensitivity analysis
├── 04_visualization.py        # step 4: generate all figures
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

---

## Data paths

All paths are defined in one place — `mitotyping_analysis.py`. Change them there,
not scattered across scripts.

```python
DATA_DIR       = "/data/CEM/smacklab/data/bican/u01_cell_class/intermediate_cell_class/"
MITOCARTA_PATH = "/scratch/easmit31/mitotyping/Human.MitoCarta3.0.xls"
ORTHOLOG_CSV   = "/scratch/easmit31/data/human-macaque-orthologs/ensembl113_mmul10_macaque_human.csv"
OUTDIR         = "/scratch/easmit31/mitotyping/results/"
```

---

## Pipeline walkthrough

### Step 1 — Pseudobulk and mitoPPS scoring

**Script:** `01_run_pipeline.py` → calls `mta.run_mitotyping()`

**What is pseudobulking?**
Single-nucleus RNA-seq gives you one expression profile per *cell*. But for
regression analysis, you need one profile per *sample* — a stable average that
represents a cell population rather than the noise of a single cell. Pseudobulking
collapses all cells of the same type in the same animal and region into a single
mean expression vector. This is the standard approach for bulk-like analysis of
single-cell data.

**What the function does, step by step:**

1. Reads each h5ad file using `backed='r'` (reads only metadata into memory —
   the actual count matrix stays on disk until needed, saving RAM)
2. Maps macaque Ensembl gene IDs → human gene symbols using one-to-one orthologs.
   This is necessary because MitoCarta uses human gene symbols.
3. Subsets to the ~964 MitoCarta genes detected in the data; discards all other genes
4. Collapses cells to mean log-normalized expression per **(animal × louvain × region)**
5. Applies three filters:

```python
MIN_CELLS   = 100    # an animal must have ≥100 cells in this louvain×region group
                     # to contribute a pseudobulk sample
MIN_ANIMALS = 10     # after applying MIN_CELLS, at least 10 animals must remain
                     # for this louvain×region combo to be included at all
MIN_AGE     = 1      # exclude animals younger than 1 year (juvenile animals
                     # may have very different developmental profiles)
```

6. Scores each MitoCarta Level-3 pathway as the **mean expression** of its detected
   member genes in that pseudobulk sample
7. Computes **mitoPPS** via pairwise ratio normalization (see formula below)

**Why louvain clusters?**
Louvain clusters are subclusters within each cell type — e.g. within astrocytes,
louvain 0 might represent a subtype concentrated in cortex while louvain 5 is more
prevalent in subcortical regions. By pseudobulking at the louvain level rather than
the cell type level, we preserve subtype-specific variation. The tradeoff is that
per-louvain sample sizes are smaller, which is why some combos fail MIN_ANIMALS=10.

**The mitoPPS formula:**

```
mitoPPS(sample S, pathway P) =
    mean over all other pathways Q ≠ P:
        [score(S, P) / score(S, Q)]           ← ratio in this sample
        ──────────────────────────────────────
        [mean over all samples of score(S,P)/score(S,Q)]  ← expected ratio
```

Intuitively: mitoPPS asks whether pathway P is expressed *more relative to all
other pathways* in this sample than it typically is across all samples. A ratio
above 1 means yes — this pathway is being prioritized. Think of it like a fold-change:
1 = same as the reference average.
2 = twice as prioritized as average.
0.5 = half as prioritized as average
This cancels out any differences in overall mitochondrial activity between samples.

**Why 69 pathways and not all 149?**

MitoCarta organizes pathways into a three-level hierarchy:
- **L1** (7 pathways): very broad categories like "OXPHOS" and "Metabolism"
- **L2** (39 pathways): intermediate groupings like "Complex I" and "Fatty acid metabolism"
- **L3** (103 pathways): the most specific level, e.g. "CI subunits", "mtDNA repair"

We score **L3 only**. L1 and L2 are just groupings of L3 pathways — scoring them
would mean double-counting the same genes. Of the 103 L3 pathways, 34 are dropped
because fewer than 5 of their member genes are detected in the macaque data
(`MIN_PATHWAY_GENES = 5`). The most common reasons: the pathway is tiny (CII
subunits has only 4 genes total), macaque doesn't have one-to-one orthologs for
those genes, or the genes are too lowly expressed in nuclei to show up in snRNA-seq.

**Where cell types are defined** — in `mitotyping_analysis.py`:

```python
CELL_TYPES = {
    'Astrocyte':      'Res1_astrocytes_update.h5ad',
    'Basket':         'Res1_basket-cells_update.h5ad',
    'Cerebellar':     'Res1_cerebellar-neurons_subset.h5ad',
    'GABA':           'Res1_GABAergic-neurons_subset.h5ad',
    'Glutamatergic':  'Res1_glutamatergic-neurons_update.h5ad',
    'MSN':            'Res1_medium-spiny-neurons_subset.h5ad',
    'Microglia':      'Res1_microglia_new.h5ad',
    'Midbrain':       'Res1_midbrain-neurons_update.h5ad',
    'OPC':            'Res1_opc-olig_subset.h5ad',
    'Oligodendrocyte':'Res1_opc-olig_subset.h5ad',  # same file, split by louvain
    'Vascular':       'Res1_vascular-cells_subset.h5ad',
}
```

OPC and Oligodendrocyte share one h5ad file. The split is by louvain cluster:

```python
OPC_LOUVAINS = {'12', '13'}   # louvains 12 and 13 → OPC
                               # louvains 0,2,3,4,5,6,7,8,9,10,11 → Oligodendrocyte
```

**Output:** `results/all_mitopps.csv`

```
Shape:   9,024 rows × 76 columns
One row per (animal × louvain × region × cell type) pseudobulk sample

Key columns:
  animal_id   unique animal identifier (e.g. "56K", "2F2")
  cell_type   cell type name as in CELL_TYPES dict (e.g. "Astrocyte", "GABA")
  louvain     louvain cluster number (integer)
  region      brain region abbreviation (e.g. "ACC", "lCb")
  age         animal age in years (float)
  sex         animal sex ("M" or "F")
  n_cells     number of cells that contributed to this pseudobulk sample

  [69 pathway columns]
              one column per MitoCarta L3 pathway
              values are mitoPPS scores (float, centered around 1.0)
              e.g. "CI subunits", "TCA cycle", "mtDNA repair"
```

---

### Step 2 — Age regression

**Script:** `02_age_regression.py` → calls `mta.run_age_regression()`

**What it does:**

For each (cell type × louvain × region) combination, fits an OLS (ordinary least
squares) linear regression for each of the 69 pathways:

```
mitoPPS ~ age + C(sex)
```

`C(sex)` treats sex as a categorical covariate — it controls for sex differences
in mitoPPS without modeling them explicitly. The key output is `beta_age`: how much
the mitoPPS of that pathway changes per year of age, after accounting for sex.

**Multiple testing correction:**
BH-FDR (Benjamini-Hochberg false discovery rate) is applied *within* each
(cell type × louvain × region) combination, across the 69 pathways. This means
each combination gets its own set of corrected p-values. We don't apply one
global correction across all 21,459 tests, because that would be overly
conservative given the strong biological structure (most pathways won't change
in most cell types).

**Output:** `results/mitotyping_age_results_clean.csv`

```
Shape:   21,459 rows × 10 columns
One row per (pathway × louvain × region × cell type) regression

Key columns:
  cell_type    cell type (e.g. "Cerebellar", "Astrocyte")
  louvain      louvain cluster number
  region       brain region (e.g. "lCb", "ACC")
  pathway      MitoCarta L3 pathway name (e.g. "CI subunits")
  beta_age     regression coefficient for age — how much mitoPPS changes per year
               negative = pathway deprioritized with age
               positive = pathway prioritized with age
               typical range: -0.02 to +0.02
  se_age       standard error of beta_age
  pval_age     raw p-value from OLS
  fdr_age      BH-FDR corrected p-value (corrected within each combo)
               values < 0.05 are considered significant
  n_animals    number of animals in this louvain×region combo
```

---

### Visualization

**Script:** `mitotyping_plots.py`
Saves figures to `figures/`. This is the single place where all plot calls live — if you want to regenerate
specific figures or change parameters, edit this file.

---

## Plotting functions

All functions are in `mitotyping_plots.py`. To use interactively:

```python
import importlib
import mitotyping_plots as mp
importlib.reload(mp)   # picks up any edits without restarting Python
```

All functions accept `save=True` to write a 300dpi PNG to `figures/`.
Filenames are auto-generated from the plot title.

---

### `mp.plot_mitotyping()` — the main heatmap

The workhorse visualization. Handles both raw mitoPPS values and age effect betas,
at any pathway hierarchy level, with many layout options.

**Data needed:** `combined_mitopps` for `mode='mitopps'`, or `results_df` for `mode='age'`

```python
mp.plot_mitotyping(
    df,
    mode             = "mitopps",  # "mitopps" = show raw scores
                                   # "age"     = show age effect (beta_age)
    level            = "L3",       # "L1" | "L2" | "L3"
                                   # L1/L2 average the L3 scores within each group
    cell_type        = None,       # e.g. "Astrocyte" — filter to one cell type
    category         = None,       # e.g. "OXPHOS" — filter to one L1 category
    average_louvains = False,      # True = one column per cell_type×region
                                   #        (averages across all louvain clusters)
                                   # False = one column per louvain×region combo
    per_animal       = False,      # True = one column per animal, youngest→oldest
    by_region        = False,      # True = one column per region (age mode)
    facet_region     = False,      # True = one separate heatmap per region
    only_sig         = True,       # True = only show pathways with ≥1 significant hit
                                   # (only applies in age mode)
    fdr_threshold    = 0.05,       # significance threshold for the dots overlay
    cluster_rows     = False,      # hierarchically cluster pathway rows
    cluster_cols     = False,      # hierarchically cluster columns
    vmax             = None,       # clips color scale to ±vmax (e.g. 0.03)
                                   # None = auto-scale from data
    save             = False,
)
```

**How the heatmap is constructed** (simplified):

```python
# 1. filter to requested cell type / L1 category
df = df[df['cell_type'] == cell_type]

# 2. build a label for each column
df['combo'] = df['cell_type'] + ', ' + df['region']   # e.g. "Astrocyte, ACC"

# 3. pivot table: rows = pathways, columns = combos
#    each cell = mean mitoPPS (or mean beta_age) for that pathway × combo
mat = df.pivot_table(index='pathway', columns='combo',
                     values='beta_age', aggfunc='mean')

# 4. optionally aggregate from L3 → L1 or L2
mat.index = [L1_map.get(p, p) for p in mat.index]
mat = mat.groupby(level=0).mean()

# 5. optionally cluster rows using Ward linkage
#    (groups pathways with similar age-effect profiles together)
row_order = leaves_list(linkage(pdist(mat.values), method='ward'))
mat = mat.iloc[row_order]

# 6. draw — for age mode, color scale centered at 0
#    (red = increases with age, blue = decreases, white = no change)
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
ax.imshow(mat.values, cmap='RdBu_r', norm=norm)

# 7. overlay dots for significant cells (fdr < threshold)
for each cell where fdr < threshold:
    ax.plot(col, row, 'k.', markersize=3)
```

**Common usage:**

```python
# all cell types, age effects, averaged across louvains
mp.plot_mitotyping(results_df, mode='age', average_louvains=True,
                   only_sig=True, vmax=0.03, save=True)

# one heatmap per region — good for seeing which regions have the most hits
mp.plot_mitotyping(results_df, mode='age', facet_region=True,
                   only_sig=False, vmax=0.03, save=True)

# raw mitoPPS for one cell type, one column per animal ordered by age
mp.plot_mitotyping(combined_mitopps, cell_type='Astrocyte',
                   per_animal=True, save=True)

# OXPHOS pathways only
mp.plot_mitotyping(results_df, mode='age', category='OXPHOS',
                   average_louvains=True, save=True)
```

**The cell type loop is** in `04_visualization.py`:

```python
cell_types = sorted(combined_mitopps['cell_type'].unique())
for ct in cell_types:
    mp.plot_mitotyping(combined_mitopps, cell_type=ct, per_animal=True, save=True)
```

---

### `mp.plot_aging_pathways()` — top aging pathways + scatter vs age

Two plots for one cell type: (1) a heatmap of the top aging pathways with animals
ordered youngest → oldest (so a color gradient left to right indicates an age effect),
and (2) scatter plots of individual pathway scores vs age with regression lines.

**Data needed:** both `results_df` and `combined_mitopps`

```python
mp.plot_aging_pathways(
    results_df,
    combined_mitopps,
    cell_type     = "Astrocyte",
    top_n         = 10,       # show top N most negative AND top N most positive
                              # betas — so up to 2×top_n pathways total
    fdr_threshold = 0.05,
    save          = False,
)
```

**Pathway selection:** if any pathways survive FDR<0.05 for this cell type, those
are shown. If none survive (Microglia, OPC, Vascular have 0 hits), falls back to
the `top_n` most extreme pathways by raw effect size.

**Colorscale:** centered at mitoPPS=1 (the biological neutral point). White = average
priority, red = prioritized, blue = deprioritized. This is different from the age
effects heatmap which is centered at 0 (no change).

---

### `mp.summarize_mitotyping()` — direction + consistency summary

```python
direction, consistency_mat = mp.summarize_mitotyping(results_df, fdr_threshold=0.05)
```

Returns two objects:
- **`direction`**: a DataFrame with one row per cell type showing how many pathways
  are significantly increasing vs decreasing with age
- **`consistency_mat`**: a pathways × cell types matrix where each value is the
  fraction of regions where that pathway is significant for that cell type.
  A value of 1.0 means the pathway is significant in every region; 0.5 means
  half the regions.

---

### `mp.plot_direction()` — direction bar chart

```python
mp.plot_direction(direction, save=True)
```

Horizontal bar chart: red bars = pathways increasing with age, blue = decreasing.
Shows both the number and direction of significant aging effects per cell type.

---

### `mp.plot_consistency()` — consistency heatmap

```python
mp.plot_consistency(consistency_mat, save=True)

# filter to pathways significant in ≥50% of regions for at least one cell type
# — these are the most robust, region-consistent findings
mp.plot_consistency(consistency_mat, min_frac=0.5, save=True)
```

---

### `mp.plot_significant_dotplot()` — significant hits overview

```python
mp.plot_significant_dotplot(results_df, results_sens, fdr_threshold=0.05, save=True)
```

One dot per significant (cell type × louvain × region × pathway) hit.
Grouped by region → pathway. Dot size = −log10(FDR) (bigger = more significant),
color = cell type. Stars mark hits that only appear in the sensitivity analysis
(i.e. only significant after removing outliers).

---

### `mp.plot_mitopps_pca()` — PCA of mitoPPS space

```python
coords, pca_obj, pca_df = mp.plot_mitopps_pca(
    combined_mitopps,
    color_by = 'cell_type',   # 'cell_type' | 'region' | 'age'
    pc_x = 1, pc_y = 2,
    save = False,
)
```

PCA (principal component analysis) reduces the 69-dimensional mitoPPS space to
2 dimensions for visualization. Each dot = one pseudobulk sample. If dots cluster
tightly by cell type, it confirms that cell types have distinct mitotypes.

---

### `mp.plot_pathway_rank()` — pathway ranking between cell types

```python
mp.plot_pathway_rank(
    combined_mitopps,
    cell_type_a = "Astrocyte",
    cell_type_b = None,    # None = compare against all other cell types averaged
    top_n       = 20,
    save        = False,
)
```

Ranked bar chart: which pathways does cell type A invest in more or less than
cell type B? Bars are colored by L1 category.

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

Scatter of two pathway scores against each other across all pseudobulk samples.
Useful for seeing whether two pathways are co-prioritized within or across cell types.

---

### `mp.plot_pathway_bubbles()` — bubble plot

```python
mp.plot_pathway_bubbles(
    combined_mitopps,
    pathways_to_show = ["CI subunits", "TCA cycle", "mtDNA repair", ...],
    n_top = 10,    # if pathways_to_show=None, uses top n by variance
    save  = False,
)
```

Compact matrix: pathways × cell types. Bubble size = z-score magnitude (how far
from the mean), color = direction (red = prioritized, blue = deprioritized).

---

## Results files

| File | Description |
|---|---|
| `results/all_mitopps.csv` | 9,024 pseudobulk samples × 69 mitoPPS scores + metadata |
| `results/mitotyping_age_results_clean.csv` | 21,459 regression tests — primary results |
| `results/mitotyping_age_results_sensitivity.csv` | Same, outliers removed (21,321 tests, 79 sig) |
| `results/pathway_gene_counts.csv` | Which pathways included/excluded and why |

**Note on outlier QC:** the outlier detection identified 41 suspicious animal×combo
pairs, but we did **not** remove them from the primary analysis. The sensitivity
file exists to show the results are robust to their inclusion.

---

## Interactive explorer

The HTML explorer embeds all data and runs entirely in the browser — no server,
no internet connection needed. It's a single self-contained file you can email or
share via Slack.

**To regenerate after updating results:**

```bash
conda run -n mito_env_test python3 generate_explorer.py
# output: figures/mitotyping_explorer_MMDD_HHMM.html
```

**To change the layout or controls:** edit `explorer_template.html`. This is the
source file with all the HTML, CSS, and JavaScript. The generator script reads this
template and injects the real data into it. **Never edit the generated HTML output
directly** — it gets overwritten every time you regenerate.

**To download:** use the OOD file browser at `ood04.sol.rc.asu.edu`, navigate to
`figures/`, and click Download on the latest timestamped file.

---

## Key analysis decisions

| Decision | Rationale |
|---|---|
| Level-3 pathways only | L1/L2 are groupings of L3 — scoring them would double-count the same genes |
| MIN_PATHWAY_GENES=5 | Fewer than 5 genes → the pathway mean is too noisy to be reliable |
| MIN_CELLS=100 | Pseudobulk means are only stable when averaging enough cells |
| MIN_ANIMALS=10 | Need enough biological replicates (animals) to run a regression |
| Per-louvain pseudobulk | Preserves subtype-level variation; tradeoff is smaller sample sizes per combo |
| BH-FDR within combo | Global correction across 21k tests would be far too conservative |
| Outliers kept in primary | Sensitivity analysis confirmed 0 hits lost on removal |
| `.X` not re-normalized | The h5ads already contain log-normalized counts |

---

## Common questions

**Why do some cell types appear in fewer regions?**
Most absences are biological — basket cells are a cerebellar interneuron type,
midbrain neurons only exist in the midbrain, MSN are striatal. EXC is additionally
absent from CN and NAc because those brain areas have very few glutamatergic neurons
(median ~10 cells per animal in CN), so the per-louvain filters exclude them.

**Why does OPC have 0 significant hits?**
OPC is present in all 11 regions with reasonable animal counts (31–51 per region).
The null result is real — OPC mitotypes appear stable with age in this dataset.

**What does a negative beta_age mean?**
The pathway becomes less prioritized relative to other mitochondrial pathways as
animals age. For example, CI subunits in Cerebellar neurons has beta_age = -0.007,
meaning mitoPPS declines by ~0.007 per year — about a 14% decline over 20 years.

**How do I filter to one region?**
```python
sub = results_df[results_df['region'] == 'lCb']
mp.plot_mitotyping(sub, mode='age', average_louvains=True, save=True)
```

**How do I search for a pathway by keyword?**
```python
pw_cols = [c for c in combined_mitopps.columns
           if c not in {'animal_id','louvain','region','age','sex','n_cells','cell_type'}]
matches = [p for p in pw_cols if 'tca' in p.lower()]
print(matches)
# ['TCA cycle', 'TCA-associated']
```

**How do I add a new cell type?**
1. Add its h5ad filename to `CELL_TYPES` in `mitotyping_analysis.py`
2. Re-run `01_run_pipeline.py` to generate mitoPPS data
3. Re-run `02_age_regression.py` for regression results
4. Re-run `04_visualization.py` for figures
5. Re-run `generate_explorer.py` to update the interactive explorer
