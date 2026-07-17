import matplotlib
matplotlib.use('Agg')
import sys
sys.path.insert(0, '/scratch/easmit31/mitotyping')
import importlib
import mitotyping_plots as mp
import mitotyping_analysis as ma
importlib.reload(mp)
importlib.reload(ma)
import pandas as pd

results_df = pd.read_csv('/scratch/easmit31/mitotyping/results/mitotyping_age_results_clean.csv')

# apply cell type abbreviations
results_abbr = ma.apply_cell_type_abbr(results_df)

print("Cell types after abbreviation:", sorted(results_abbr['cell_type'].unique()))
print("Significant hits:", (results_abbr['fdr_age'] < 0.05).sum())

# also update cell_type_colors in mp to use abbreviated keys
mp.cell_type_colors = {
    ma.CELL_TYPE_ABBR.get(k, k): v
    for k, v in mp.cell_type_colors.items()
}

print("\nPlotting all regions with abbreviated cell type labels...")
mp.plot_mitotyping(
    results_abbr,
    mode='age',
    facet_region=True,
    only_sig=False,
    fdr_threshold=0.05,
    vmax=0.03,
    save=True
)

print("\nPlotting significant only...")
mp.plot_mitotyping(
    results_abbr,
    mode='age',
    facet_region=True,
    only_sig=True,
    fdr_threshold=0.05,
    vmax=0.03,
    save=True
)

print(f"\nDone. Figures in {mp.FIGURE_DIR}")
