"""
generate_explorer.py
Loads real mitotyping data and injects it into explorer_template.html.
Edit explorer_template.html for layout/style/JS changes.
Edit this file for data loading/aggregation changes.

Usage:
    conda run -n mito_env_test python3 generate_explorer.py
Output:
    /scratch/easmit31/mitotyping/figures/mitotyping_explorer.html
"""
import json
import pandas as pd
import numpy as np
import os

RESULTS_DIR  = "/scratch/easmit31/mitotyping/results/"
FIGURES_DIR  = "/scratch/easmit31/mitotyping/figures/"
MITOCARTA    = "/scratch/easmit31/mitotyping/Human.MitoCarta3.0.xls"
SCRIPT_PATH  = "/scratch/easmit31/mitotyping/mitotyping_plots.py"
TEMPLATE     = "/scratch/easmit31/mitotyping/explorer_template.html"
from datetime import datetime
OUTPUT       = f"{FIGURES_DIR}mitotyping_explorer_{datetime.now().strftime('%m%d_%H%M')}.html"
os.makedirs(FIGURES_DIR, exist_ok=True)

def safe_round(x, nd=4):
    """Round a value, but convert NaN/Inf to None (-> JSON null) instead of
    letting a bare NaN token leak into the exported JS. A raw NaN there breaks
    the heatmap's color function (cfn) and silently aborts the whole render
    partway through — this is what was producing the blank/half-drawn
    heatmaps, not a data-density issue."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    if np.isnan(v) or np.isinf(v):
        return None
    return round(v, nd)

# ── load mitoPPS ──────────────────────────────────────────────────────────────
print("Loading mitoPPS...")
mito = pd.read_csv(f"{RESULTS_DIR}all_mitopps.csv")
mito['louvain'] = mito['louvain'].astype(str)
META = {'animal_id','louvain','region','age','sex','n_cells','cell_type'}
pw_cols = [c for c in mito.columns if c not in META]

# ── load age results ──────────────────────────────────────────────────────────
print("Loading age results...")
results = pd.read_csv(f"{RESULTS_DIR}mitotyping_age_results_clean.csv")
results['louvain'] = results['louvain'].astype(str)

# ── aggregate 1: mean per cell_type x region (averaged across louvains/animals)
print("Aggregating: cell type x region...")
mitopps_avg = {}
for (ct, reg), grp in mito.groupby(['cell_type','region']):
    if ct not in mitopps_avg: mitopps_avg[ct] = {}
    mitopps_avg[ct][reg] = {pw: safe_round(grp[pw].mean(), 4) for pw in pw_cols}

# ── aggregate 2: mean per cell_type x animal (ordered by age) ────────────────
# NOTE: this pools each animal's cells across ALL regions for that cell type.
# It's used by the "Per animal" heatmap view (which is itself region-agnostic).
print("Aggregating: per animal...")
mitopps_animal = {}
animal_ages = {}
for (ct, animal_id), grp in mito.groupby(['cell_type','animal_id']):
    if ct not in mitopps_animal: mitopps_animal[ct] = {}
    if ct not in animal_ages: animal_ages[ct] = {}
    mitopps_animal[ct][animal_id] = {pw: safe_round(grp[pw].mean(), 4) for pw in pw_cols}
    animal_ages[ct][animal_id] = safe_round(grp['age'].iloc[0], 2)

# sort animals by age per cell type
animal_order = {}
for ct in animal_ages:
    animal_order[ct] = sorted(animal_ages[ct].keys(), key=lambda a: animal_ages[ct][a])

# ── aggregate 2b: mean per cell_type x region x animal (ordered by age) ──────
# Region-resolved counterpart of aggregate 2, keyed "cell_type|region" to match
# MITOPPS_AVG / AGE. This is what the heatmap hover-scatter should use whenever
# the hovered cell corresponds to a specific cell type AND a specific region
# (e.g. the "per louvain" view, or the avg view with a region filter / the
# cell-type-x-region "separate" column mode) — otherwise the scatter pools
# animals across every region for that cell type, which inflates n far beyond
# the n shown for that specific louvain x region combo.
print("Aggregating: per animal x region...")
mitopps_animal_reg = {}
animal_ages_reg = {}
for (ct, reg, animal_id), grp in mito.groupby(['cell_type','region','animal_id']):
    key = f"{ct}|{reg}"
    if key not in mitopps_animal_reg: mitopps_animal_reg[key] = {}
    if key not in animal_ages_reg: animal_ages_reg[key] = {}
    mitopps_animal_reg[key][animal_id] = {pw: safe_round(grp[pw].mean(), 4) for pw in pw_cols}
    animal_ages_reg[key][animal_id] = safe_round(grp['age'].iloc[0], 2)

animal_order_reg = {}
for key in animal_ages_reg:
    animal_order_reg[key] = sorted(animal_ages_reg[key].keys(), key=lambda a: animal_ages_reg[key][a])

# ── aggregate 3: mean per cell_type x louvain x region ───────────────────────
print("Aggregating: per louvain...")
mitopps_louvain = {}
for (ct, lou, reg), grp in mito.groupby(['cell_type','louvain','region']):
    if ct not in mitopps_louvain: mitopps_louvain[ct] = {}
    key = f"{lou},{reg}"
    mitopps_louvain[ct][key] = {pw: safe_round(grp[pw].mean(), 4) for pw in pw_cols}

# louvain x region column order per cell type (numeric sort on louvain)
louvain_order = {}
for ct in mitopps_louvain:
    keys = list(mitopps_louvain[ct].keys())
    louvain_order[ct] = sorted(keys, key=lambda k: (k.split(',')[1], int(k.split(',')[0]) if k.split(',')[0].isdigit() else 0))

# ── aggregate 4: age effects ──────────────────────────────────────────────────
print("Aggregating age effects...")
age_data = {}
for (ct, reg, pw), grp in results.groupby(['cell_type','region','pathway']):
    if ct not in age_data: age_data[ct] = {}
    if reg not in age_data[ct]: age_data[ct][reg] = {}
    age_data[ct][reg][pw] = {
        'b': safe_round(grp['beta_age'].mean(), 5),
        'f': safe_round(grp['fdr_age'].min(),   5),
    }

# ── n animals per cell_type x region ─────────────────────────────────────────
n_animals = {}
for (ct, reg), grp in mito.groupby(['cell_type','region']):
    if ct not in n_animals: n_animals[ct] = {}
    n_animals[ct][reg] = int(grp['animal_id'].nunique())

# ── MitoCarta hierarchy ───────────────────────────────────────────────────────
print("Loading MitoCarta...")
mc = pd.read_excel(MITOCARTA, sheet_name="C MitoPathways")
mc = mc.dropna(subset=['MitoPathways Hierarchy','Genes'])
mc['L1'] = mc['MitoPathways Hierarchy'].str.split(' > ').str[0]
mc['L2'] = mc['MitoPathways Hierarchy'].str.split(' > ').str[1]

def get_level(row):
    if row['MitoPathway'] == row['L1']: return 'L1'
    if row['MitoPathway'] == row['L2']: return 'L2'
    return 'L3'
mc['level'] = mc.apply(get_level, axis=1)

pw_to_l1 = dict(zip(mc['MitoPathway'], mc['L1']))
pw_to_l2 = dict(zip(mc['MitoPathway'], mc['L2']))
l1_cats  = sorted(mc['L1'].dropna().unique().tolist())
l2_by_l1 = (mc.groupby('L1')['L2']
               .unique()
               .apply(lambda x: sorted(set(v for v in x.tolist() if isinstance(v, str))))
               .to_dict())

# ── metadata ──────────────────────────────────────────────────────────────────
cell_types = sorted(mito['cell_type'].unique().tolist())

FULL_NAMES = {
    'Astrocyte':      'Astrocyte (AST)',
    'Basket':         'Basket cells (BC)',
    'Cerebellar':     'Cerebellar neurons (CER)',
    'Glutamatergic':  'Glutamatergic neurons (EXC)',
    'GABA':           'GABAergic neurons (INH)',
    'Midbrain':       'Midbrain neurons (MBN)',
    'Microglia':      'Microglia (MGL)',
    'MSN':            'Medium spiny neurons (MSN)',
    'OPC':            'OPC',
    'Oligodendrocyte':'Oligodendrocyte (OLIG)',
    'Vascular':       'Vascular cells (VASC)',
}
regions    = sorted(mito['region'].unique().tolist())
pathways   = pw_cols
n_sig      = int((results['fdr_age'] < 0.05).sum())
n_animals_total = int(mito['animal_id'].nunique())

print(f"Cell types: {cell_types}")
print(f"Regions: {regions}")
print(f"Pathways: {len(pathways)}")
print(f"Sig hits: {n_sig}")
print(f"Animals: {n_animals_total}")

# ── substitution dict ─────────────────────────────────────────────────────────
subs = {
    '__MITOPPS_AVG_JSON__':        json.dumps(mitopps_avg),
    '__MITOPPS_ANIMAL_JSON__':     json.dumps(mitopps_animal),
    '__MITOPPS_ANIMAL_REG_JSON__': json.dumps(mitopps_animal_reg),
    '__MITOPPS_LOUVAIN_JSON__':    json.dumps(mitopps_louvain),
    '__ANIMAL_ORDER_JSON__':       json.dumps(animal_order),
    '__ANIMAL_ORDER_REG_JSON__':   json.dumps(animal_order_reg),
    '__ANIMAL_AGES_JSON__':        json.dumps(animal_ages),
    '__ANIMAL_AGES_REG_JSON__':    json.dumps(animal_ages_reg),
    '__LOUVAIN_ORDER_JSON__':      json.dumps(louvain_order),
    '__AGE_JSON__':                json.dumps(age_data),
    '__N_ANIM_JSON__':             json.dumps(n_animals),
    '__PW_TO_L1_JSON__':           json.dumps(pw_to_l1),
    '__PW_TO_L2_JSON__':           json.dumps(pw_to_l2),
    '__L2_BY_L1_JSON__':           json.dumps({k: list(v) for k,v in l2_by_l1.items()}),
    '__CT_JSON__':                 json.dumps(cell_types),
    '__REG_JSON__':                json.dumps(regions),
    '__PW_JSON__':                 json.dumps(pathways),
    '__L1_CATS_JSON__':            json.dumps(l1_cats),
    '__PW_L2_JSON__':              json.dumps(sorted({l2 for ls in l2_by_l1.values() for l2 in ls})),
    '__SCRIPT_PATH__':             SCRIPT_PATH,
    '__RESULTS_PATH__':            RESULTS_DIR,
    '__FIGURES_PATH__':            FIGURES_DIR,
    '__N_CELL_TYPES__':            str(len(cell_types)),
    '__N_REGIONS__':               str(len(regions)),
    '__N_PATHWAYS__':              str(len(pathways)),
    '__N_SIG__':                   str(n_sig),
    '__N_ANIMALS__':               str(n_animals_total),
    '__CT_OPTIONS__':              '\n'.join(f'<option value="{ct}">{FULL_NAMES.get(ct,ct)}</option>' for ct in cell_types),
    '__REG_OPTIONS__':             '\n'.join(f'<option value="{r}">{r}</option>' for r in regions),
    '__L1_OPTIONS__':              '\n'.join(f'<option value="{c}">{c}</option>' for c in l1_cats),
}

print(f"Loading template: {TEMPLATE}")
with open(TEMPLATE) as f:
    html = f.read()

for placeholder, value in subs.items():
    html = html.replace(placeholder, value)

with open(OUTPUT, 'w') as f:
    f.write(html)

print(f"\nSaved: {OUTPUT}")
print("Open in any browser — no internet required")
