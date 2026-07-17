path = '/scratch/easmit31/mitotyping/mitotyping_plots.py'
txt  = open(path).read()

# fix heatmap tick label sizes in plot_aging_pathways
old1 = '''    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(mat.columns, fontsize=5, rotation=45, ha='right')
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(mat.index, fontsize=7)
    title_hm = f"{cell_type} — top aging pathways ordered by age"
    ax.set_title(title_hm, fontsize=10, pad=10)'''

new1 = '''    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(mat.columns, fontsize=8, rotation=45, ha='right')
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(mat.index, fontsize=10)
    title_hm = f"{cell_type} — top aging pathways ordered by age"
    ax.set_title(title_hm, fontsize=13, pad=10)'''

# fix scatter plot font sizes
old2 = '''        ax.set_title(f"{pw[:35]}\\nbeta={beta:.4f}, FDR={fdr:.3f}", fontsize=7)
        ax.set_xlabel('age (years)', fontsize=7)
        ax.set_ylabel('mitoPPS', fontsize=7)'''

new2 = '''        ax.set_title(f"{pw[:35]}\\nbeta={beta:.4f}, FDR={fdr:.3f}", fontsize=10)
        ax.set_xlabel('age (years)', fontsize=10)
        ax.set_ylabel('mitoPPS', fontsize=10)'''

old3 = '''    title_sc = f"{cell_type} — mitoPPS vs age"
    plt.suptitle(title_sc, fontsize=10, y=1.01)'''

new3 = '''    title_sc = f"{cell_type} — mitoPPS vs age"
    plt.suptitle(title_sc, fontsize=13, y=1.01)'''

fixes = [(old1,new1,'heatmap tick labels'),(old2,new2,'scatter labels'),(old3,new3,'scatter title')]
for old, new, label in fixes:
    if old in txt:
        txt = txt.replace(old, new)
        print(f"Fixed: {label}")
    else:
        print(f"Pattern not found: {label}")

open(path, 'w').write(txt)
print("Done.")
