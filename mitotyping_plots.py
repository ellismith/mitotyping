import os
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import matplotlib.patches as mpatches


# ── color palettes matching your R colors ────────────────────────────────────
# NOTE: keys updated to standard abbreviations (AST, BC, CER, EXC, INH, MBN,
# MGL, MSN, OLIG, OPC, VASC) to match the relabeling done in
# mitotyping_analysis.apply_cell_type_abbr(). If you still have data using
# full cell-type names (Astrocyte, Glutamatergic, etc.), color lookups here
# will fall back to gray until that data is relabeled too.
cell_type_colors = {
    'MSN':   '#5e3c99',
    'CER':   '#44839f',
    'BC':    '#4e379a',
    'AST':   '#d53a50',
    'OLIG':  '#fb8d64',
    'OPC':   '#fbae64',
    'VASC':  '#d1538f',
    'MGL':   '#bb471c',
    'EPEN':  '#fc6238',
    'EXC':   '#3d9995',
    'INH':   '#4d6aa8',
    'MBN':   '#65bfb0',
}

region_colors = {
    'NAc':   '#148e34',
    'ACC':   '#5555b0',
    'IPP':   '#6c9cc0',
    'CN':    '#2caa4d',
    'HIP':   '#84db55',
    'lCb':   '#e6cd69',
    'mdTN':  '#117867',
    'dlPFC': '#345296',
    'M1':    '#467aa1',
    'EC':    '#532b70',
    'MB':    '#fc7462',
}

# ── MitoCarta — loaded once on import ────────────────────────────────────────
MITOCARTA_PATH = "/scratch/easmit31/mitotyping/Human.MitoCarta3.0.xls"
FIGURE_DIR     = "/scratch/easmit31/mitotyping/figures/"
os.makedirs(FIGURE_DIR, exist_ok=True)

pathways = pd.read_excel(MITOCARTA_PATH, sheet_name="C MitoPathways")
pathways = pathways.dropna(subset=['Genes']).reset_index(drop=True)
pathways['gene_list'] = pathways['Genes'].str.split(', ')
hierarchy = pathways['MitoPathways Hierarchy'].str.split(' > ', expand=True)
hierarchy.columns = ['L1','L2','L3'][:hierarchy.shape[1]]
pathways  = pd.concat([pathways, hierarchy], axis=1)

def _get_level(row):
    if row['MitoPathway'] == row.get('L1'): return 'Level1'
    if row['MitoPathway'] == row.get('L2'): return 'Level2'
    if row['MitoPathway'] == row.get('L3'): return 'Level3'
    return 'Unknown'

pathways['level'] = pathways.apply(_get_level, axis=1)
pathway_to_L1     = dict(zip(pathways['MitoPathway'], pathways['L1']))

META = {'animal_id','louvain','region','age','sex','n_cells','cell_type',
        'subtype','pathway','beta_age','se_age','pval_age','fdr_age','n_animals'}


# ── filename helper ───────────────────────────────────────────────────────────
def _title_to_filename(title):
    """Convert a plot title to a clean filename."""
    fname = title.replace('\n','_').replace(' | ','_').replace(' ','_')
    fname = ''.join(c for c in fname if c.isalnum() or c in '_-')
    return fname[:120] + '.png'


# ── internal helpers ──────────────────────────────────────────────────────────
def _draw_heatmap(mat, fdr_mat=None, title="", fdr_threshold=0.05,
                  cmap='viridis', cbar_label='', save=False, vmax=None):
    """Draw a heatmap with optional FDR dot overlay and optional save."""
    n_rows, n_cols = mat.shape
    fig, ax = plt.subplots(figsize=(max(14, n_cols*0.5), max(4, n_rows*0.35)))
    vals = mat.values.astype(float)
    if cmap == 'RdBu_r':
        vmax = np.nanpercentile(np.abs(vals[~np.isnan(vals)]), 95)
        norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    else:
        norm = None
    im = ax.imshow(vals, aspect='auto', cmap=cmap, norm=norm)
    if fdr_mat is not None:
        for r in range(n_rows):
            for c in range(n_cols):
                v = fdr_mat.values[r, c]
                if not np.isnan(v) and v < fdr_threshold:
                    ax.plot(c, r, 'k.', markersize=3, alpha=0.8)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(mat.columns, fontsize=12, rotation=45, ha='right')
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(mat.index, fontsize=12)
    ax.set_title(title, fontsize=14, pad=10)
    cbar = plt.colorbar(im, ax=ax, label=cbar_label, shrink=0.4) 
    cbar.set_label(cbar_label, fontsize=16)    # colorbar title/label — bumped
    cbar.ax.tick_params(labelsize=13)          # colorbar tick numbers — bumped
    
    plt.tight_layout()
    if save:
        fpath = os.path.join(FIGURE_DIR, _title_to_filename(title))
        fig.savefig(fpath, dpi=300, bbox_inches='tight')
        print(f"Saved: {fpath}")
    plt.show()


def _make_combo_label(df, groupby, average_louvains, cell_type):
    """Build x-axis grouping label."""
    if average_louvains:
        if 'cell_type' in df.columns and cell_type is None:
            return df['cell_type'] + '\n' + df['region']
        return df['region']
    if groupby == "louvain_region": return df['louvain'].astype(str) + '\n' + df['region']
    if groupby == "region":         return df['region']
    if groupby == "louvain":        return df['louvain'].astype(str)
    return df['louvain'].astype(str) + '\n' + df['region']


def _agg_level(mat, level):
    """Aggregate matrix rows to L1 or L2 by averaging across member pathways."""
    if level == 'L3': return mat
    lmap      = pathways.set_index('MitoPathway')[level].to_dict()
    mat       = mat.copy()
    mat.index = [lmap.get(p, p) for p in mat.index]
    return mat.groupby(level=0).mean()


def _sort_cols(mat, per_animal):
    """Sort columns: by first part alphabetically, then second part numerically."""
    if per_animal: return mat
    def sort_key(x):
        parts = x.split(', ')
        if len(parts) == 1:
            # single value e.g. region name or cell type
            p0 = parts[0]
            return (p0.lower(), 0, p0)
        p0, p1 = parts[0], parts[1]
        # p0 could be louvain number or cell type depending on groupby
        # p1 could be region or louvain
        p0_num = int(p0) if p0.isdigit() else 0
        p1_num = int(p1) if p1.isdigit() else 0
        p0_str = p0 if not p0.isdigit() else ''
        p1_str = p1 if not p1.isdigit() else ''
        # if p0 is a number (louvain) sort by region first then louvain
        if p0.isdigit():
            return (p1_str.lower(), p0_num, p0_str)
        # if p0 is a string (cell type) sort by cell type then region
        return (p0_str.lower(), p1_num, p1_str)
    return mat[sorted(mat.columns, key=sort_key)]


def _cluster(mat, fdr_mat=None, cluster_rows=False, cluster_cols=False):
    """Hierarchically cluster rows and/or columns using ward linkage."""
    mat_filled = mat.fillna(mat.mean(axis=1).mean())
    if cluster_rows and len(mat) > 2:
        row_order = leaves_list(linkage(pdist(mat_filled.values), method='ward'))
        mat = mat.iloc[row_order]
        if fdr_mat is not None: fdr_mat = fdr_mat.iloc[row_order]
    if cluster_cols and mat.shape[1] > 2:
        col_order = leaves_list(linkage(pdist(mat_filled.values.T), method='ward'))
        mat = mat.iloc[:, col_order]
        if fdr_mat is not None: fdr_mat = fdr_mat.iloc[:, col_order]
    return mat, fdr_mat


# ── main plot function ────────────────────────────────────────────────────────
def plot_mitotyping(df, mode="mitopps", level="L3", category=None,
                    cell_type=None, groupby="louvain_region",
                    average_louvains=False, per_animal=False, louvain=None,
                    fdr_threshold=0.05, only_sig=True,
                    by_region=False, facet_region=False,
                    cluster_rows=False, cluster_cols=False,
                    save=False, vmax=None):
    """
    Unified heatmap for all mitotyping visualizations.

    Parameters
    ----------
    df               : combined_mitopps (mode=mitopps/scores) or results_df (mode=age)
    mode             : "mitopps" | "scores" | "age"
    level            : "L1" | "L2" | "L3"
    category         : filter to one L1 category e.g. "OXPHOS"
    cell_type        : filter to one cell type, or None for all
    groupby          : "louvain_region" | "louvain" | "region"
    average_louvains : one col per cell_type x region
    per_animal       : one col per animal ordered youngest to oldest
    louvain          : filter to one louvain cluster
    fdr_threshold    : FDR threshold for significance dots (mode=age)
    only_sig         : only show pathways with >= 1 hit (mode=age)
    by_region        : cols = regions (mode=age)
    facet_region     : one heatmap per region (mode=age)
    cluster_rows     : hierarchically cluster pathway rows
    cluster_cols     : hierarchically cluster columns
    save             : save PNG to FIGURE_DIR with auto-generated filename

    Examples
    --------
    plot_mitotyping(combined_mitopps)
    plot_mitotyping(combined_mitopps, level="L1", average_louvains=True)
    plot_mitotyping(combined_mitopps, cell_type="Astrocyte", per_animal=True)
    plot_mitotyping(results_df, mode="age", average_louvains=True)
    plot_mitotyping(results_df, mode="age", category="OXPHOS", save=True)
    plot_mitotyping(results_df, mode="age", facet_region=True)
    plot_mitotyping(results_df, mode="age", cluster_rows=True, cluster_cols=True)
    """
    df = df.copy()
    if 'louvain' in df.columns: df['louvain'] = df['louvain'].astype(str)

    # filters
    if cell_type is not None and 'cell_type' in df.columns:
        df = df[df['cell_type'] == cell_type]
        print(f"Cell type: {cell_type} ({len(df)} rows)")
    if louvain is not None and 'louvain' in df.columns:
        df = df[df['louvain'] == str(louvain)]
        print(f"Louvain: {louvain} ({len(df)} rows)")
    if len(df) == 0: return print("No data after filtering.")

    fdr_mat = None

    # ── age mode ──────────────────────────────────────────────────────
    if mode == "age":
        if category:
            df = df[df['pathway'].map(pathway_to_L1) == category]
            print(f"Category '{category}': {df['pathway'].nunique()} pathways")

        # one heatmap per region
        if facet_region:
            for region in sorted(df['region'].unique()):
                sub = df[df['region'] == region].copy()
                sub['combo'] = sub['cell_type'] if 'cell_type' in sub.columns else sub['louvain']
                mat     = sub.pivot_table(index='pathway', columns='combo',
                                          values='beta_age', aggfunc='mean')
                fdr_mat = sub.pivot_table(index='pathway', columns='combo',
                                          values='fdr_age',  aggfunc='min')
                if only_sig:
                    sig = (fdr_mat < fdr_threshold).any(axis=1)
                    mat, fdr_mat = mat[sig], fdr_mat[sig]
                if len(mat) == 0: continue
                mat, fdr_mat    = _agg_level(mat, level), _agg_level(fdr_mat, level)
                mat, fdr_mat    = _cluster(mat, fdr_mat, cluster_rows, cluster_cols)
                mat             = mat[sorted(mat.columns)]
                fdr_mat         = fdr_mat[sorted(fdr_mat.columns)]
                mat.columns     = [c.replace('\n',', ') for c in mat.columns]
                fdr_mat.columns = mat.columns
                parts = [f"age | {region}", f"level={level}"]
                if category: parts.append(category)
                parts.append(f"FDR<{fdr_threshold}")
                _draw_heatmap(mat, fdr_mat=fdr_mat, title=" | ".join(parts),
                              fdr_threshold=fdr_threshold,
                              cmap='RdBu_r', cbar_label='beta (age effect)', save=save)
            return

        # build combo label
        if by_region:
            df['combo'] = df['region']
        elif average_louvains:
            df['combo'] = (df['cell_type'] + '\n' + df['region']
                           if 'cell_type' in df.columns and cell_type is None
                           else df['region'])
        else:
            df['combo'] = _make_combo_label(df, groupby, average_louvains, cell_type)

        mat     = df.pivot_table(index='pathway', columns='combo',
                                  values='beta_age', aggfunc='mean')
        fdr_mat = df.pivot_table(index='pathway', columns='combo',
                                  values='fdr_age',  aggfunc='min')

        if only_sig:
            sig = (fdr_mat < fdr_threshold).any(axis=1)
            mat, fdr_mat = mat[sig], fdr_mat[sig]
            print(f"Significant pathways (FDR<{fdr_threshold}): {sig.sum()}")
        if len(mat) == 0: return print("No significant pathways.")

        mat, fdr_mat = _agg_level(mat, level), _agg_level(fdr_mat, level)
        mat, fdr_mat = _cluster(mat, fdr_mat, cluster_rows, cluster_cols)

        if not cluster_cols:
            mat     = _sort_cols(mat, per_animal=False)
            fdr_mat = fdr_mat.reindex(columns=mat.columns)

        mat.columns     = [c.replace('\n',', ') for c in mat.columns]
        fdr_mat.columns = mat.columns

        parts = ["age", f"level={level}"]
        if cell_type:                   parts.append(cell_type)
        elif 'cell_type' in df.columns: parts.append("all cell types")
        if category:         parts.append(category)
        if by_region:        parts.append("by region")
        if average_louvains: parts.append("louvains averaged")
        if cluster_rows:     parts.append("rows clustered")
        if cluster_cols:     parts.append("cols clustered")
        parts.append(f"FDR<{fdr_threshold}")
        if vmax is not None: parts.append(f"vmax={vmax}")
        _draw_heatmap(mat, fdr_mat=fdr_mat, title=" | ".join(parts),
                      fdr_threshold=fdr_threshold,
                      cmap='RdBu_r', cbar_label='beta (age effect)', save=save, vmax=vmax)
        return

    # ── mitopps / scores mode ─────────────────────────────────────────
    pw_cols = [c for c in df.columns if c not in META]
    if category: pw_cols = [p for p in pw_cols if pathway_to_L1.get(p) == category]
    if not pw_cols: return print("No pathways to plot.")

    if per_animal:
        animal_meta = df.groupby('animal_id')['age'].first().reset_index().sort_values('age')
        mat         = df.groupby('animal_id')[pw_cols].mean().T
        mat         = mat[animal_meta['animal_id'].values]
        age_map     = dict(zip(animal_meta['animal_id'], animal_meta['age'].round(1)))
        mat.columns = [f"{a}\n({age_map[a]}y)" for a in mat.columns]
    else:
        df['combo'] = _make_combo_label(df, groupby, average_louvains, cell_type)
        mat = df[pw_cols + ['combo']].groupby('combo')[pw_cols].mean().T

    mat = _agg_level(mat, level)
    mat, _ = _cluster(mat, None, cluster_rows, cluster_cols and not per_animal)

    if not cluster_cols and not per_animal:
        mat = _sort_cols(mat, per_animal)

    mat.columns = [c.replace('\n',', ') for c in mat.columns]

    parts = [mode, f"level={level}"]
    if cell_type:                   parts.append(cell_type)
    elif 'cell_type' in df.columns: parts.append("all cell types")
    if category:         parts.append(category)
    if per_animal:       parts.append("per animal, age ordered")
    if average_louvains: parts.append("louvains averaged")
    if cluster_rows:     parts.append("rows clustered")
    if cluster_cols:     parts.append("cols clustered")
    _draw_heatmap(mat, title=" | ".join(parts), cmap='viridis', cbar_label=mode, save=save)


# ── summary functions ─────────────────────────────────────────────────────────
def summarize_mitotyping(results_df, fdr_threshold=0.05):
    """
    Summarize age regression results.

    Returns
    -------
    direction       : per cell type — n_sig, n_increasing, n_decreasing
    consistency_mat : pathways x cell types — fraction of regions significant
    """
    if 'cell_type' not in results_df.columns:
        results_df = results_df.copy()
        results_df['cell_type'] = 'unknown'

    sig = results_df[results_df['fdr_age'] < fdr_threshold].copy()
    if len(sig) == 0:
        return print(f"No significant results at FDR<{fdr_threshold}")

    direction = (sig.groupby('cell_type')
                 .apply(lambda x: pd.Series({
                     'n_sig':          len(x),
                     'n_increasing':   (x['beta_age'] > 0).sum(),
                     'n_decreasing':   (x['beta_age'] < 0).sum(),
                     'pct_increasing': round((x['beta_age'] > 0).mean() * 100, 1),
                 }))
                 .reset_index()
                 .sort_values('n_sig', ascending=False))

    print("=== Direction summary ===")
    print(direction.to_string(index=False))

    n_regions     = results_df.groupby(['cell_type','pathway'])['region'].nunique()
    n_sig_regions = sig.groupby(['cell_type','pathway'])['region'].nunique()
    consistency   = (n_sig_regions / n_regions).rename('frac_sig').reset_index()

    consistency_mat = (consistency
                       .pivot_table(index='pathway', columns='cell_type', values='frac_sig')
                       .fillna(0))
    consistency_mat = consistency_mat.loc[
        consistency_mat.mean(axis=1).sort_values(ascending=False).index]

    print(f"\n=== Consistency matrix: {consistency_mat.shape} ===")
    print(f"Pathways significant in ≥1 region for ≥1 cell type: "
          f"{(consistency_mat > 0).any(axis=1).sum()}")

    return direction, consistency_mat


def plot_direction(direction_summary, save=False):
    """Bar chart of significant hits per cell type, split by direction."""
    df  = direction_summary.sort_values('n_sig', ascending=True)
    fig, ax = plt.subplots(figsize=(6, max(4, len(df) * 0.4)))
    ax.barh(df['cell_type'],  df['n_increasing'], color='#d73027', label='increasing')
    ax.barh(df['cell_type'], -df['n_decreasing'], color='#4575b4', label='decreasing')
    ax.axvline(0, color='k', linewidth=0.8)
    ax.set_xlabel('n significant pathways')
    ax.set_title('Age-associated mitoPPS changes\nby cell type and direction', fontsize=10)
    ax.legend(frameon=False)
    plt.tight_layout()
    if save:
        fpath = os.path.join(FIGURE_DIR, 'direction_summary.png')
        fig.savefig(fpath, dpi=300, bbox_inches='tight')
        print(f"Saved: {fpath}")
    plt.show()


def plot_consistency(consistency_mat, min_frac=0.0, save=False):
    """
    Heatmap of pathway consistency across cell types.
    Value = fraction of regions where pathway is significant.
    """
    mat = consistency_mat.copy()
    if min_frac > 0:
        keep = (mat >= min_frac).any(axis=1)
        mat  = mat[keep]
        print(f"Pathways with frac_sig >= {min_frac}: {keep.sum()}")
    if len(mat) == 0: return print("Nothing to plot.")
    mat = mat[sorted(mat.columns)]
    n_rows, n_cols = mat.shape
    fig, ax = plt.subplots(figsize=(max(8, n_cols*0.6), max(4, n_rows*0.35)))
    im = ax.imshow(mat.values, aspect='auto', cmap='Reds', vmin=0, vmax=1)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(mat.columns, fontsize=8, rotation=45, ha='right')
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(mat.index, fontsize=7)
    title = f"Pathway consistency (min_frac={min_frac})" if min_frac > 0 else "Pathway consistency"
    ax.set_title(f"{title}\nvalue = fraction of regions significant", fontsize=10, pad=10)
    plt.colorbar(im, ax=ax, label='fraction of regions (FDR<0.05)', shrink=0.4)
    plt.tight_layout()
    if save:
        fpath = os.path.join(FIGURE_DIR, f'consistency_minfrac{min_frac}.png')
        fig.savefig(fpath, dpi=300, bbox_inches='tight')
        print(f"Saved: {fpath}")
    plt.show()


# ── paper-equivalent plots ────────────────────────────────────────────────────
def plot_pathway_rank(combined_mitopps, cell_type_a, cell_type_b=None, top_n=20, save=False):
    """
    Rank pathways by mean mitoPPS difference between two cell types.
    Equivalent to Fig 4A.
    """
    pw_cols = [c for c in combined_mitopps.columns if c not in META]
    means   = combined_mitopps.groupby('cell_type')[pw_cols].mean()

    if cell_type_a not in means.index:
        return print(f"'{cell_type_a}' not found. Available: {list(means.index)}")

    a    = means.loc[cell_type_a]
    b    = means.loc[cell_type_b] if cell_type_b and cell_type_b in means.index \
           else means.drop(cell_type_a).mean()
    diff = (a - b).sort_values()

    unique_cats = sorted(set(pathway_to_L1.get(p,'Other') for p in diff.index))
    palette     = plt.cm.Set2(np.linspace(0, 1, len(unique_cats)))
    cat_colors  = dict(zip(unique_cats, palette))
    colors      = [cat_colors.get(pathway_to_L1.get(p,'Other'),'gray') for p in diff.index]

    label_b = cell_type_b if cell_type_b else "all others"
    title   = f"Pathway ranking | {cell_type_a} vs {label_b}"

    fig, ax = plt.subplots(figsize=(7, max(6, len(diff)*0.18)))
    ax.barh(range(len(diff)), diff.values, color=colors, height=0.8)
    ax.set_yticks(range(len(diff)))
    ax.set_yticklabels(diff.index, fontsize=6)
    ax.axvline(0, color='k', linewidth=0.8)
    ax.ticklabel_format(axis='x', useOffset=False, style='plain')
    ax.set_xlabel(f'mean mitoPPS: {cell_type_a} minus {label_b}', fontsize=8, labelpad=12)
    ax.set_title(title, fontsize=10)
    patches = [mpatches.Patch(color=c, label=cat) for cat, c in cat_colors.items()]
    ax.legend(handles=patches, fontsize=6, frameon=False,
              bbox_to_anchor=(1.01,1), loc='upper left')
    plt.tight_layout()
    if save:
        fig.savefig(os.path.join(FIGURE_DIR, _title_to_filename(title)),
                    dpi=300, bbox_inches='tight')
    plt.show()


def plot_pathway_bivariate(combined_mitopps, pathway_x, pathway_y,
                            color_by='cell_type', save=False):
    """Scatter of two pathway scores. Equivalent to Fig 4B/C."""
    if pathway_x not in combined_mitopps.columns or pathway_y not in combined_mitopps.columns:
        return print("Pathway not found. Check spelling.")

    df      = combined_mitopps[[pathway_x, pathway_y, color_by]].dropna()
    groups  = sorted(df[color_by].unique())
    palette = plt.cm.tab20(np.linspace(0, 1, len(groups)))
    cmap    = dict(zip(groups, palette))
    title   = f"Bivariate | {pathway_x} vs {pathway_y}"

    fig, ax = plt.subplots(figsize=(6, 5))
    for grp, sub in df.groupby(color_by):
        ax.scatter(sub[pathway_x], sub[pathway_y],
                   color=cmap[grp], label=grp, alpha=0.5, s=15, linewidths=0)
    lim = [min(df[pathway_x].min(), df[pathway_y].min()),
           max(df[pathway_x].max(), df[pathway_y].max())]
    ax.plot(lim, lim, 'k--', linewidth=0.8, alpha=0.4)
    ax.set_xlabel(pathway_x, fontsize=9)
    ax.set_ylabel(pathway_y, fontsize=9)
    ax.set_title(title, fontsize=9)
    ax.legend(fontsize=6, frameon=False, bbox_to_anchor=(1.01,1), loc='upper left')
    plt.tight_layout()
    if save:
        fig.savefig(os.path.join(FIGURE_DIR, _title_to_filename(title)),
                    dpi=300, bbox_inches='tight')
    plt.show()


def plot_radar(combined_mitopps, cell_types=None, top_n=10, save=False):
    """Radar chart of top pathways per cell type. Equivalent to Fig 4D."""
    pw_cols = [c for c in combined_mitopps.columns if c not in META]
    means   = combined_mitopps.groupby('cell_type')[pw_cols].mean()
    if cell_types: means = means.loc[cell_types]

    zscores   = (means - means.mean()) / means.std()
    top_paths = sorted(set(p for ct in zscores.index
                           for p in zscores.loc[ct].nlargest(top_n).index))
    data      = zscores[top_paths]
    n_vars    = len(top_paths)
    angles    = np.linspace(0, 2*np.pi, n_vars, endpoint=False).tolist() + [0]
    palette   = plt.cm.Set2(np.linspace(0, 1, len(means.index)))

    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
    for i, ct in enumerate(data.index):
        vals = data.loc[ct].tolist() + [data.loc[ct].tolist()[0]]
        ax.plot(angles, vals, color=palette[i], linewidth=1.5, label=ct)
        ax.fill(angles, vals, color=palette[i], alpha=0.1)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(top_paths, fontsize=6)
    ax.set_title("Mitotype radar chart", fontsize=10, pad=20)
    ax.legend(fontsize=7, frameon=False, bbox_to_anchor=(1.3, 1.1))
    plt.tight_layout()
    if save:
        fig.savefig(os.path.join(FIGURE_DIR, 'radar_chart.png'),
                    dpi=300, bbox_inches='tight')
    plt.show()


def plot_mitopps_pca(combined_mitopps, color_by='cell_type', n_components=3,
                     pc_x=1, pc_y=2, age_trajectory=False, save=False):
    """PCA of mitoPPS scores. Equivalent to Fig 5B-D."""
    pw_cols = [c for c in combined_mitopps.columns if c not in META]
    df      = combined_mitopps.dropna(subset=pw_cols)
    X_std   = StandardScaler().fit_transform(df[pw_cols].values)
    pca     = PCA(n_components=n_components)
    coords  = pca.fit_transform(X_std)
    var_ex  = pca.explained_variance_ratio_ * 100
    pc_xi, pc_yi = pc_x-1, pc_y-1
    title   = f"mitoPPS PCA | colored by {color_by}"

    fig, ax = plt.subplots(figsize=(7, 6))
    if color_by == 'age':
        sc = ax.scatter(coords[:,pc_xi], coords[:,pc_yi],
                        c=df['age'].values, cmap='plasma', s=15, alpha=0.7, linewidths=0)
        plt.colorbar(sc, ax=ax, label='age (years)', shrink=0.6)
    else:
        groups  = sorted(df[color_by].unique())
        palette = plt.cm.tab20(np.linspace(0, 1, len(groups)))
        cmap    = dict(zip(groups, palette))
        for grp in groups:
            mask = df[color_by] == grp
            ax.scatter(coords[mask,pc_xi], coords[mask,pc_yi],
                       color=cmap[grp], label=grp, s=15, alpha=0.6, linewidths=0)
        if age_trajectory:
            for grp in groups:
                mask  = (df[color_by] == grp).values
                order = np.argsort(df['age'].values[mask])
                ax.plot(coords[mask,pc_xi][order], coords[mask,pc_yi][order],
                        color=cmap[grp], linewidth=0.8, alpha=0.4)
        ax.legend(fontsize=7, frameon=False, bbox_to_anchor=(1.01,1), loc='upper left')

    ax.set_xlabel(f'PC{pc_x} ({var_ex[pc_xi]:.1f}%)', fontsize=9)
    ax.set_ylabel(f'PC{pc_y} ({var_ex[pc_yi]:.1f}%)', fontsize=9)
    ax.set_title(title, fontsize=10)
    ax.axhline(0, color='k', linewidth=0.5, alpha=0.3)
    ax.axvline(0, color='k', linewidth=0.5, alpha=0.3)
    plt.tight_layout()
    if save:
        fig.savefig(os.path.join(FIGURE_DIR, _title_to_filename(title)),
                    dpi=300, bbox_inches='tight')
    plt.show()
    return coords, pca, df


def plot_mitopps_pca_trajectory(combined_mitopps, cell_type=None, save=False):
    """PCA averaged to animal level with age trajectory arrow. Equivalent to Fig 6C."""
    pw_cols = [c for c in combined_mitopps.columns if c not in META]
    df      = combined_mitopps.copy()
    if cell_type: df = df[df['cell_type'] == cell_type]

    animal_df = df.groupby('animal_id').agg(
        age=('age','first'), **{pw: (pw,'mean') for pw in pw_cols}
    ).reset_index().sort_values('age')

    X      = StandardScaler().fit_transform(animal_df[pw_cols].values)
    pca    = PCA(n_components=2)
    coords = pca.fit_transform(X)
    var_ex = pca.explained_variance_ratio_ * 100
    ages   = animal_df['age'].values
    title  = f"mitoPPS aging trajectory{' | ' + cell_type if cell_type else ''}"

    fig, ax = plt.subplots(figsize=(7, 6))
    sc = ax.scatter(coords[:,0], coords[:,1], c=ages, cmap='plasma',
                    s=80, zorder=3, linewidths=0.5, edgecolors='white')

    young = coords[ages <= np.percentile(ages, 33)]
    old   = coords[ages >= np.percentile(ages, 67)]
    y_cen, o_cen = young.mean(axis=0), old.mean(axis=0)
    ax.annotate('', xy=o_cen, xytext=y_cen,
                arrowprops=dict(arrowstyle='->', color='black', lw=2, mutation_scale=20))
    ax.scatter(*y_cen, color='blue',   s=100, zorder=5, label='young centroid')
    ax.scatter(*o_cen, color='yellow', s=100, zorder=5, label='old centroid',
               edgecolors='k', linewidths=0.5)

    plt.colorbar(sc, ax=ax, label='age (years)', shrink=0.6)
    ax.set_xlabel(f'PC1 ({var_ex[0]:.1f}%)', fontsize=9)
    ax.set_ylabel(f'PC2 ({var_ex[1]:.1f}%)', fontsize=9)
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=7, frameon=False)
    ax.axhline(0, color='k', linewidth=0.5, alpha=0.3)
    ax.axvline(0, color='k', linewidth=0.5, alpha=0.3)
    plt.tight_layout()
    if save:
        fig.savefig(os.path.join(FIGURE_DIR, _title_to_filename(title)),
                    dpi=300, bbox_inches='tight')
    plt.show()
    print(f"Young→old centroid distance: {np.linalg.norm(o_cen - y_cen):.3f}")


def plot_pca_loadings(pca, pw_cols, top_n=10, n_pcs=3, save=False):
    """PCA pathway loadings bar chart. Equivalent to Fig 5E."""
    fig, axes = plt.subplots(1, n_pcs, figsize=(n_pcs*5, 6))
    for i, ax in enumerate(axes):
        loadings = pd.Series(pca.components_[i], index=pw_cols)
        top      = pd.concat([loadings.nlargest(top_n), loadings.nsmallest(top_n)])
        colors   = ['#d73027' if v > 0 else '#4575b4' for v in top.values]
        ax.barh(range(len(top)), top.values, color=colors)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top.index, fontsize=7)
        ax.axvline(0, color='k', linewidth=0.8)
        ax.set_title(f'PC{i+1} ({pca.explained_variance_ratio_[i]*100:.1f}% var)', fontsize=9)
        ax.set_xlabel('loading', fontsize=8)
    plt.suptitle('PCA pathway loadings', fontsize=11, y=1.02)
    plt.tight_layout()
    if save:
        fig.savefig(os.path.join(FIGURE_DIR, 'pca_loadings.png'),
                    dpi=300, bbox_inches='tight')
    plt.show()


def plot_pathway_bubbles(combined_mitopps, pathways_to_show=None, cell_types=None,
                          n_top=10, save=False):
    """Bubble plot of pathways x cell types. Equivalent to Fig 5F."""
    pw_cols = [c for c in combined_mitopps.columns if c not in META]
    means   = combined_mitopps.groupby('cell_type')[pw_cols].mean()
    if cell_types: means = means.loc[[c for c in cell_types if c in means.index]]
    if pathways_to_show is None:
        pathways_to_show = means.std().nlargest(n_top).index.tolist()

    data  = means[pathways_to_show].T
    zdata = data.subtract(data.mean(axis=1), axis=0).divide(data.std(axis=1), axis=0).clip(-3, 3)
    title = "mitoPPS bubble plot"

    n_rows, n_cols = len(pathways_to_show), len(means.index)
    fig, ax = plt.subplots(figsize=(max(6, n_cols*0.7), max(4, n_rows*0.4)))
    for i, pw in enumerate(pathways_to_show):
        for j, ct in enumerate(means.index):
            val   = zdata.loc[pw, ct] if ct in zdata.columns else 0
            size  = min(abs(val)*200, 500)
            color = '#d73027' if val > 0 else '#4575b4'
            ax.scatter(j, i, s=size, color=color, alpha=0.7, linewidths=0)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(means.index, fontsize=8, rotation=45, ha='right')
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(pathways_to_show, fontsize=7)
    ax.set_title(f"{title}\nsize=|z-score|, red=prioritized, blue=deprioritized", fontsize=9)
    ax.grid(True, alpha=0.2)
    plt.tight_layout()
    if save:
        fig.savefig(os.path.join(FIGURE_DIR, _title_to_filename(title)),
                    dpi=300, bbox_inches='tight')
    plt.show()


def plot_aging_pathways(results_df, combined_mitopps, cell_type,
                         top_n=10, fdr_threshold=0.05, save=False):
    """
    Top aging pathways heatmap + scatter vs age. Equivalent to Fig 6D/E.
    """
    res = results_df[results_df['cell_type'] == cell_type].copy()
    res = res.groupby('pathway').agg(
        beta_age=('beta_age','mean'), fdr_age=('fdr_age','min')
    ).reset_index()

    sig = res[res['fdr_age'] < fdr_threshold].sort_values('beta_age')
    top = pd.concat([sig.head(top_n), sig.tail(top_n)]).drop_duplicates()
    if len(top) == 0:
        top = res.sort_values('beta_age').pipe(
            lambda x: pd.concat([x.head(top_n), x.tail(top_n)])).drop_duplicates()
        print(f"No significant pathways — showing top {top_n} by effect size")

    pw_cols     = top['pathway'].tolist()
    ct_data     = combined_mitopps[combined_mitopps['cell_type'] == cell_type].copy()
    animal_means = ct_data.groupby('animal_id').agg(
        age=('age','first'), **{pw: (pw,'mean') for pw in pw_cols}
    ).reset_index().sort_values('age')

    # ── heatmap ───────────────────────────────────────────────────────
    mat         = animal_means[pw_cols].T
    mat.columns = [f"{str(a)} ({float(age):.0f}y)"
                   for a, age in zip(animal_means['animal_id'], animal_means['age'])]
    n_rows, n_cols = mat.shape
    fig, ax = plt.subplots(figsize=(max(10, n_cols*0.3), max(4, n_rows*0.4)))
    vals = mat.values.astype(float)
    im   = ax.imshow(vals, aspect='auto', cmap='viridis')
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(mat.columns, fontsize=8, rotation=45, ha='right')
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(mat.index, fontsize=10)
    title_hm = f"{cell_type} — top aging pathways ordered by age"
    ax.set_title(title_hm, fontsize=13, pad=10)
    plt.colorbar(im, ax=ax, label='mitoPPS', shrink=0.4)
    plt.tight_layout()
    if save:
        fig.savefig(os.path.join(FIGURE_DIR, _title_to_filename(title_hm)),
                    dpi=300, bbox_inches='tight')
    plt.show()

    # ── scatter plots ─────────────────────────────────────────────────
    plot_pws = pw_cols[:6]
    n_c, n_r = 3, int(np.ceil(len(plot_pws)/3))
    fig, axes = plt.subplots(n_r, n_c, figsize=(n_c*4, n_r*3))
    axes = axes.flatten() if n_r*n_c > 1 else [axes]
    for ax, pw in zip(axes, plot_pws):
        sub = animal_means[['age', pw]].dropna()
        x, y = sub['age'].astype(float).values, sub[pw].astype(float).values
        ax.scatter(x, y, s=20, alpha=0.7, color='steelblue')
        xfit = np.linspace(x.min(), x.max(), 100)
        ax.plot(xfit, np.polyval(np.polyfit(x, y, 1), xfit), 'r-', linewidth=1.5)
        beta = top[top['pathway']==pw]['beta_age'].values[0]
        fdr  = top[top['pathway']==pw]['fdr_age'].values[0]
        ax.set_title(f"{pw[:35]}\nbeta={beta:.4f}, FDR={fdr:.3f}", fontsize=10)
        ax.set_xlabel('age (years)', fontsize=10)
        ax.set_ylabel('mitoPPS', fontsize=10)
    for ax in axes[len(plot_pws):]: ax.set_visible(False)
    title_sc = f"{cell_type} — mitoPPS vs age"
    plt.suptitle(title_sc, fontsize=13, y=1.01)
    plt.tight_layout()
    if save:
        fig.savefig(os.path.join(FIGURE_DIR, _title_to_filename(title_sc)),
                    dpi=300, bbox_inches='tight')
    plt.show()


def plot_mitopps_pca_trajectory(combined_mitopps, cell_type=None):
    """
    PCA of mitoPPS averaged to animal level, with age trajectory.
    Each dot = one animal. Color = age.
    """
    META    = {'animal_id','louvain','region','age','sex','n_cells','cell_type'}
    pw_cols = [c for c in combined_mitopps.columns if c not in META]

    df = combined_mitopps.copy()
    if cell_type:
        df = df[df['cell_type'] == cell_type]

    # average to one row per animal
    animal_df = df.groupby('animal_id').agg(
        age=('age','first'),
        **{pw: (pw,'mean') for pw in pw_cols}
    ).reset_index().sort_values('age')

    X      = StandardScaler().fit_transform(animal_df[pw_cols].values)
    pca    = PCA(n_components=2)
    coords = pca.fit_transform(X)
    var_ex = pca.explained_variance_ratio_ * 100

    fig, ax = plt.subplots(figsize=(7, 6))

    # color by age
    sc = ax.scatter(coords[:,0], coords[:,1],
                    c=animal_df['age'].values, cmap='plasma',
                    s=60, zorder=3, linewidths=0)

    # connect by age order
    ax.plot(coords[:,0], coords[:,1], '-', color='gray',
            linewidth=0.8, alpha=0.4, zorder=2)

    plt.colorbar(sc, ax=ax, label='age (years)', shrink=0.6)
    ax.set_xlabel(f'PC1 ({var_ex[0]:.1f}%)', fontsize=9)
    ax.set_ylabel(f'PC2 ({var_ex[1]:.1f}%)', fontsize=9)
    title = f"mitoPPS aging trajectory"
    if cell_type: title += f" | {cell_type}"
    ax.set_title(title, fontsize=10)
    ax.axhline(0, color='k', linewidth=0.5, alpha=0.3)
    ax.axvline(0, color='k', linewidth=0.5, alpha=0.3)
    plt.tight_layout()
    plt.show()


def plot_significant_dotplot(results_df, results_sens=None, fdr_threshold=0.05, save=False):
    """
    Dot plot of significant age effects.
    Circles = main analysis
    Stars = only significant without outliers
    """
    sig = results_df[results_df["fdr_age"] < fdr_threshold].copy()

    agg = sig.groupby(["cell_type", "region", "pathway"]).agg(
        beta_age=("beta_age", "mean"),
        fdr_age=("fdr_age", "min"),
        n_louvains=("louvain", "nunique") if "louvain" in sig.columns else ("beta_age", "size")
    ).reset_index()

    agg["log_fdr"] = -np.log10(agg["fdr_age"])

    agg_sens = None
    new_hits = set()

    if results_sens is not None:
        sig_sens = results_sens[results_sens["fdr_age"] < fdr_threshold].copy()

        agg_sens = sig_sens.groupby(["cell_type", "region", "pathway"]).agg(
            beta_age=("beta_age", "mean"),
            fdr_age=("fdr_age", "min"),
            n_louvains=("louvain", "nunique") if "louvain" in sig_sens.columns else ("beta_age", "size")
        ).reset_index()

        agg_sens["log_fdr"] = -np.log10(agg_sens["fdr_age"])

        orig_keys = set(zip(agg.cell_type, agg.region, agg.pathway))
        sens_keys = set(zip(agg_sens.cell_type, agg_sens.region, agg_sens.pathway))
        new_hits = sens_keys - orig_keys

        print(f"Main hits: {len(orig_keys)}")
        print(f"Sensitivity hits: {len(sens_keys)}")
        print(f"Sensitivity-only: {len(new_hits)}")

    region_order = ["ACC", "CN", "dlPFC", "EC", "HIP", "IPP", "lCb", "M1", "MB", "mdTN", "NAc"]
    all_regions = set(agg.region)
    if agg_sens is not None:
        all_regions |= set(agg_sens.region)
    region_order = [r for r in region_order if r in all_regions]

    y_pos = {}
    y_ticks = []
    y_dividers = []
    gap = 1
    current_y = 0

    for reg in region_order:
        reg_data = agg[agg.region == reg]
        if agg_sens is not None:
            reg_data = pd.concat([reg_data, agg_sens[agg_sens.region == reg]])

        reg_data = reg_data.drop_duplicates(["pathway"])
        if len(reg_data) == 0:
            continue

        pw_order = reg_data.groupby("pathway")["beta_age"].mean().sort_values().index

        for pw in pw_order:
            y_pos[(reg, pw)] = current_y
            y_ticks.append((current_y, pw, 6.5, "#333333"))
            current_y += 1

        y_ticks.append((current_y - len(pw_order) / 2, reg, 8, region_colors.get(reg, "black")))
        y_dividers.append(current_y - .5 + gap / 2)
        current_y += gap

    total_y = current_y
    fig, ax = plt.subplots(figsize=(10, max(8, total_y * .28)))

    ct_list = sorted(set(agg.cell_type) | (set() if agg_sens is None else set(agg_sens.cell_type)))

    for _, row in agg.iterrows():
        if (row.region, row.pathway) not in y_pos:
            continue

        y = y_pos[(row.region, row.pathway)]
        color = cell_type_colors.get(row.cell_type, "gray")
        size = row.log_fdr * 40
        jitter = (ct_list.index(row.cell_type) - len(ct_list) / 2) * 0.0005
        x = row.beta_age + jitter

        ax.scatter(x, y, s=size, color=color, alpha=.85, linewidths=0, zorder=3)

        if row.n_louvains > 1:
            ax.annotate(
                str(int(row.n_louvains)),
                (x, y),
                textcoords="offset points",
                xytext=(0, 4),
                ha="center",
                va="bottom",
                fontsize=6,
                color="black",
                zorder=4
            )

    if agg_sens is not None:
        for _, row in agg_sens.iterrows():
            key = (row.cell_type, row.region, row.pathway)
            if key not in new_hits:
                continue

            y = y_pos[(row.region, row.pathway)]
            color = cell_type_colors.get(row.cell_type, "gray")
            size = row.log_fdr * 55
            jitter = (ct_list.index(row.cell_type) - len(ct_list) / 2) * 0.0005
            x = row.beta_age + jitter

            ax.scatter(
                x, y,
                marker="*",
                s=size,
                color=color,
                edgecolor="black",
                linewidth=.6,
                zorder=5
            )

            if row.n_louvains > 1:
                ax.annotate(
                    str(int(row.n_louvains)),
                    (x, y),
                    textcoords="offset points",
                    xytext=(0, 4),
                    ha="center",
                    va="bottom",
                    fontsize=6,
                    color="black",
                    zorder=6
                )

    ax.axvline(0, color="k", lw=.8)

    pw_ticks = [t for t in y_ticks if t[1] not in region_order]
    reg_ticks = [t for t in y_ticks if t[1] in region_order]

    ax.set_yticks([y for y, _, _, _ in pw_ticks])
    ax.set_yticklabels([lab for _, lab, _, _ in pw_ticks], fontsize=6.5)
    ax.set_ylim(-.5, total_y)

    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks([y for y, _, _, _ in reg_ticks])
    ax2.set_yticklabels([lab for _, lab, _, _ in reg_ticks], fontsize=9)

    for tick, (_, _, _, col) in zip(ax2.get_yticklabels(), reg_ticks):
        tick.set_color(col)
        tick.set_fontweight("bold")

    for yd in y_dividers[:-1]:
        ax.axhline(yd, color="gray", lw=1, alpha=.3)

    ax.set_xlabel("mean beta_age (averaged across louvains)")
    ax.set_title(
        f"Significant age effects (FDR < {fdr_threshold})\n"
        f"grouped by region → pathway | size = -log10(FDR)",
        fontsize=10
    )

    handles = [
        plt.scatter([], [], color=cell_type_colors.get(ct, "gray"), s=50, label=ct)
        for ct in ct_list
    ]

    if len(new_hits):
        handles.append(
            plt.scatter([], [], marker="*", color="gray", edgecolor="black", s=90, label="only significant without outliers")
        )

    ax.legend(handles=handles, fontsize=7, frameon=False, bbox_to_anchor=(1.15, 1), loc="upper left")

    plt.tight_layout()

    if save:
        plt.savefig(os.path.join(FIGURE_DIR, "significant_dotplot.png"), dpi=300, bbox_inches="tight")

    plt.show()
