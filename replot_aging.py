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

for ct in sorted(results_df['cell_type'].unique()):
    print(f'Plotting {ct}...')
    mp.plot_aging_pathways(results_df, combined_mitopps, cell_type=ct, save=True)
    plt.close('all')

print('Done.')
