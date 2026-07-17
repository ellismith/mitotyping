path = '/scratch/easmit31/mitotyping/mitotyping_plots.py'
txt  = open(path).read()

old = '''    vals = mat.values.astype(float)
    vmax = np.nanpercentile(np.abs(vals[~np.isnan(vals)]), 95)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    im   = ax.imshow(vals, aspect='auto', cmap='RdBu_r', norm=norm)'''

new = '''    vals = mat.values.astype(float)
    vmax_val = np.nanpercentile(np.abs(vals[~np.isnan(vals)] - 1), 95)
    norm = TwoSlopeNorm(vmin=1-vmax_val, vcenter=1, vmax=1+vmax_val)
    im   = ax.imshow(vals, aspect='auto', cmap='RdBu_r', norm=norm)'''

if old in txt:
    txt = txt.replace(old, new)
    open(path, 'w').write(txt)
    print("Fixed: colorscale now centered at mitoPPS=1")
else:
    print("Pattern not found — print the relevant lines:")
    for i, line in enumerate(txt.split('\n')):
        if 'vmax' in line and 'aging' not in line.lower():
            print(f"  {i}: {line}")
