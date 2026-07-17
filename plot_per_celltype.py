import matplotlib
matplotlib.use('Agg')
import sys
sys.path.insert(0, '/scratch/easmit31/mitotyping')
import importlib
import mitotyping_plots as mp
importlib.reload(mp)
import pandas as pd
import matplotlib.pyplot as plt

results_df       = pd.read_csv('/scratch/easmit31/mitotyping/results/mitotyping_age_results_clean.csv')
combined_mitopps = pd.read_csv('/scratch/easmit31/mitotyping/results/all_mitopps.csv')

cell_types = sorted(combined_mitopps['cell_type'].unique())

for ct in cell_types:
    print(f"\n{'='*40}\n{ct}\n{'='*40}")

    # mitoPPS heatmap — all louvains x regions
    print(f"  mitoPPS (all louvains)...")
    mp.plot_mitotyping(combined_mitopps, mode='mitopps', cell_type=ct,
                       level='L3', save=True)
    plt.close('all')

    # mitoPPS heatmap — averaged across louvains
    print(f"  mitoPPS (louvains averaged)...")
    mp.plot_mitotyping(combined_mitopps, mode='mitopps', cell_type=ct,
                       level='L3', average_louvains=True, save=True)
    plt.close('all')

    # age effects heatmap
    sub = results_df[results_df['cell_type'] == ct]
    n_sig = (sub['fdr_age'] < 0.05).sum()
    print(f"  age effects ({n_sig} sig hits)...")
    mp.plot_mitotyping(sub, mode='age', cell_type=ct,
                       level='L3', average_louvains=True,
                       only_sig=False, vmax=0.03, save=True)
    plt.close('all')

    # per animal
    print(f"  per animal...")
    mp.plot_mitotyping(combined_mitopps, mode='mitopps', cell_type=ct,
                       level='L3', per_animal=True, save=True)
    plt.close('all')

print(f"\nDone. Figures in {mp.FIGURE_DIR}")

# rerun just the facet by region plots
print("\nPlotting age effects faceted by region per cell type...")
for ct in cell_types:
    sub = results_df[results_df['cell_type'] == ct]
    n_sig = (sub['fdr_age'] < 0.05).sum()
    print(f"  {ct} ({n_sig} sig hits)...")
    mp.plot_mitotyping(sub, mode='age', cell_type=ct,
                       level='L3', facet_region=True,
                       only_sig=False, vmax=0.03, save=True)
    plt.close('all')

print(f"Done. Figures in {mp.FIGURE_DIR}")
