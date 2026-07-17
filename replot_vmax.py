import matplotlib
matplotlib.use('Agg')
import sys
sys.path.insert(0, '/scratch/easmit31/mitotyping')
import importlib
import mitotyping_plots as mp
importlib.reload(mp)
import pandas as pd
results_df = pd.read_csv('/scratch/easmit31/mitotyping/results/mitotyping_age_results_clean.csv')
for vmax in [None, 0.03, 0.01]:
    print(f"\nPlotting vmax={vmax}...")
    mp.plot_mitotyping(results_df, mode='age', average_louvains=True,
                       save=True, vmax=vmax)
    mp.plot_mitotyping(results_df, mode='age', level='L1',
                       average_louvains=True, save=True, vmax=vmax)
    mp.plot_mitotyping(results_df, mode='age', by_region=True,
                       save=True, vmax=vmax)
print(f"\nDone. Figures in {mp.FIGURE_DIR}")
