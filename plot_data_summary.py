import matplotlib
matplotlib.use('Agg')
import pandas as pd
import anndata as ad
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import sys
sys.path.insert(0, '/scratch/easmit31/mitotyping')
import mitotyping_analysis as ma

DATA_DIR  = '/data/CEM/smacklab/data/bican/u01_cell_class/intermediate_cell_class/'
OUTDIR    = '/scratch/easmit31/mitotyping/figures/'
REGIONS   = ['ACC','CN','dlPFC','EC','HIP','IPP','lCb','M1','MB','mdTN','NAc']
MIN_CELLS   = 100
MIN_ANIMALS = 10
LOW_N       = 20  # threshold below which to flag as low n

# OPC louvains within opc-olig h5ad
OPC_LOUVAINS = {'12', '13'}

CELL_TYPES = {
    'AST':  ('Res1_astrocytes_update.h5ad',          None),
    'BC':   ('Res1_basket-cells_update.h5ad',         None),
    'CER':  ('Res1_cerebellar-neurons_subset.h5ad',   None),
    'INH':  ('Res1_GABAergic-neurons_subset.h5ad',    None),
    'EXC':  ('Res1_glutamatergic-neurons_update.h5ad',None),
    'MSN':  ('Res1_medium-spiny-neurons_subset.h5ad', None),
    'MGL':  ('Res1_microglia_new.h5ad',               None),
    'MBN':  ('Res1_midbrain-neurons_update.h5ad',     None),
    'OPC':  ('Res1_opc-olig_subset.h5ad',             'opc'),
    'OLIG': ('Res1_opc-olig_subset.h5ad',             'olig'),
    'VASC': ('Res1_vascular-cells_subset.h5ad',       None),
}

CT_ORDER = ['AST','BC','CER','EXC','INH','MBN','MGL','MSN','OPC','OLIG','VASC']

# ── build status matrix ───────────────────────────────────────────────────────
status_mat  = pd.DataFrame(0, index=CT_ORDER, columns=REGIONS)
n_pass_mat  = pd.DataFrame(0, index=CT_ORDER, columns=REGIONS)
n_total_mat = pd.DataFrame(0, index=CT_ORDER, columns=REGIONS)
reason_mat  = pd.DataFrame('', index=CT_ORDER, columns=REGIONS)

loaded = {}  # cache h5ads so opc-olig only loaded once
for ct, (fname, subset) in CELL_TYPES.items():
    print(f"Processing {ct}...")
    if fname not in loaded:
        adata = ad.read_h5ad(f'{DATA_DIR}/{fname}', backed='r')
        loaded[fname] = adata.obs[['animal_id','region','louvain']].copy()
        adata.file.close()

    obs = loaded[fname].copy()

    # filter to OPC or Oligo louvains if needed
    if subset == 'opc':
        obs = obs[obs['louvain'].astype(str).isin(OPC_LOUVAINS)]
    elif subset == 'olig':
        obs = obs[~obs['louvain'].astype(str).isin(OPC_LOUVAINS)]

    for reg in REGIONS:
        sub        = obs[obs['region'] == reg]
        n_total    = sub['animal_id'].nunique()
        per_animal = sub.groupby('animal_id').size()
        n_pass     = int((per_animal >= MIN_CELLS).sum())

        n_total_mat.loc[ct, reg] = n_total
        n_pass_mat.loc[ct, reg]  = n_pass

        if n_total == 0:
            status_mat.loc[ct, reg] = 0
            reason_mat.loc[ct, reg] = 'absent'
        elif n_pass < MIN_ANIMALS:
            status_mat.loc[ct, reg] = 1
            if n_pass == 0:
                reason_mat.loc[ct, reg] = f'sparse\n({n_total} animals\n<100 cells)'
            else:
                reason_mat.loc[ct, reg] = f'too few\n({n_pass}/{n_total}\npass ≥100)'
        else:
            status_mat.loc[ct, reg] = 2
            reason_mat.loc[ct, reg] = f'{n_pass}'

# ── plot ──────────────────────────────────────────────────────────────────────
colors = {0: '#f0ede8', 1: '#f4a261', 2: '#2196a6'}
low_color = '#a8d5db'  # lighter teal for low-n included combos
labels = {0: 'Biologically absent', 1: 'Excluded by filter', 2: 'Included in analysis'}

fig, ax = plt.subplots(figsize=(14, 8))

for i, ct in enumerate(CT_ORDER):
    for j, reg in enumerate(REGIONS):
        s      = status_mat.loc[ct, reg]
        n_pass = n_pass_mat.loc[ct, reg]
        reason = reason_mat.loc[ct, reg]

        # use lighter color for low-n included combos
        if s == 2 and n_pass < LOW_N:
            color = low_color
        else:
            color = colors[s]

        rect = plt.Rectangle([j-0.48, i-0.48], 0.96, 0.96,
                              facecolor=color, edgecolor='white', linewidth=1.5)
        ax.add_patch(rect)

        if s == 2:
            txt_color = 'white' if n_pass > 30 else '#1a1a1a'
            label = str(int(n_pass))
            if n_pass < LOW_N:
                label += '*'
                txt_color = '#1a4a50'
            ax.text(j, i, label, ha='center', va='center',
                    fontsize=10, fontweight='500', color=txt_color)
        elif s == 1:
            lines = reason.split('\n')
            for li, line in enumerate(lines):
                ax.text(j, i + (li - len(lines)/2 + 0.5)*0.28,
                        line, ha='center', va='center',
                        fontsize=6.5, color='#7a3000')
        else:
            ax.text(j, i, '—', ha='center', va='center',
                    fontsize=12, color='#c0b8b0')

ax.set_xlim(-0.5, len(REGIONS)-0.5)
ax.set_ylim(-0.5, len(CT_ORDER)-0.5)
ax.set_xticks(range(len(REGIONS)))
ax.set_xticklabels(REGIONS, fontsize=12)
ax.set_yticks(range(len(CT_ORDER)))
ax.set_yticklabels(CT_ORDER, fontsize=12)
ax.set_title('Data availability per cell type × region\n'
             'Numbers = animals included | * = low n (<20) | filters: ≥100 cells/animal, ≥10 animals/combo',
             fontsize=12, pad=12)

patches = [
    mpatches.Patch(color=colors[2],  label='Included (≥20 animals)'),
    mpatches.Patch(color=low_color,  label=f'Included but low n (<{LOW_N} animals)*'),
    mpatches.Patch(color=colors[1],  label='Excluded by filter'),
    mpatches.Patch(color=colors[0],  label='Biologically absent'),
]
ax.legend(handles=patches, loc='upper center', bbox_to_anchor=(0.5, -0.06),
          fontsize=10, frameon=False, ncol=4)

plt.tight_layout()
fpath = f'{OUTDIR}/data_availability_summary.png'
fig.savefig(fpath, dpi=300, bbox_inches='tight')
print(f'Saved: {fpath}')
