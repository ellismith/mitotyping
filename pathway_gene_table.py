import pandas as pd
import numpy as np
import anndata as ad

MITOCARTA = "/scratch/easmit31/mitotyping/Human.MitoCarta3.0.xls"
ORTHOLOG  = "/scratch/easmit31/data/human-macaque-orthologs/ensembl113_mmul10_macaque_human.csv"
RESULTS   = "/scratch/easmit31/mitotyping/results/"
DATA_DIR  = "/data/CEM/smacklab/data/bican/u01_cell_class/intermediate_cell_class/"

print("Loading orthologs...")
mart = pd.read_csv(ORTHOLOG)
mart = mart[mart['Human homology type'] == 'ortholog_one2one']
ens2human = dict(zip(mart['Gene stable ID'], mart['Human gene name']))
print(f"Orthologs: {len(ens2human)}")

print("Loading gene names from h5ad...")
adata = ad.read_h5ad(f'{DATA_DIR}/Res1_astrocytes_update.h5ad', backed='r')
var_names = pd.Series(adata.var_names).map(ens2human).fillna(pd.Series(adata.var_names)).values
detected_genes = set(var_names)
adata.file.close()
print(f"Detected genes: {len(detected_genes)}")

print("Loading MitoCarta...")
mc = pd.read_excel(MITOCARTA, sheet_name="C MitoPathways")
mc = mc.dropna(subset=['Genes','MitoPathways Hierarchy'])
mc['gene_list'] = mc['Genes'].str.split(', ')
mc['L1'] = mc['MitoPathways Hierarchy'].str.split(' > ').str[0]
mc['L2'] = mc['MitoPathways Hierarchy'].str.split(' > ').str[1]

def get_level(row):
    if row['MitoPathway'] == row['L1']: return 'L1'
    if row['MitoPathway'] == row['L2']: return 'L2'
    return 'L3'
mc['level'] = mc.apply(get_level, axis=1)

print("Computing gene counts per pathway...")
rows = []
for _, row in mc.iterrows():
    genes          = row['gene_list']
    n_total        = len(genes)
    n_detected     = len([g for g in genes if g in detected_genes])
    in_analysis    = row['level'] == 'L3' and n_detected >= 5
    if row['level'] in ['L1','L2']:
        reason = 'L1/L2 — excluded by design'
    elif n_detected < 5:
        reason = f'< 5 detected ({n_detected})'
    else:
        reason = 'included'
    rows.append({
        'pathway':        row['MitoPathway'],
        'level':          row['level'],
        'L1':             row['L1'],
        'n_total_genes':  n_total,
        'n_detected':     n_detected,
        'in_analysis':    in_analysis,
        'reason':         reason,
    })

df = pd.DataFrame(rows).sort_values(['level','n_detected'], ascending=[True,False])
out = f"{RESULTS}pathway_gene_counts.csv"
df.to_csv(out, index=False)
print(f"\nSaved: {out}")

print(f"\n=== Summary ===")
print(df.groupby('reason').size().sort_values(ascending=False))

print(f"\n=== L3 excluded (< 5 detected genes) ===")
excl = df[(df['level']=='L3')&(df['n_detected']<5)][['pathway','n_total_genes','n_detected']]
print(excl.to_string(index=False))

print(f"\n=== Included L3 pathways ({df[df['reason']=='included'].shape[0]}) ===")
incl = df[df['reason']=='included'][['pathway','L1','n_total_genes','n_detected']]
print(incl.to_string(index=False))
