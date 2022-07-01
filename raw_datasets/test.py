"""Test File"""

import data_cleaners as dc

# devo
data = dc.pd.read_csv("~/rosetta-antibody-ddgs/raw_datasets/interface_data_use.csv")
# "1YY9", "3GBN", removed for now since they aren't done
pdbs = ["1DQJ", "1MHP", "1MLC", "1N8Z", "1VFB", "4FQY"]

per_res_score = dc.pd.DataFrame({
    "#PDB": pdbs,
    "sidechain_constr": dc.np.zeros(6),
    "constr": dc.np.zeros(6),
    "harm": dc.np.zeros(6),
    "unconstr": dc.np.zeros(6),
    "unrelaxed": dc.np.zeros(6)
})

rmsds = dc.pd.DataFrame({
    "#PDB": pdbs,
    "sidechain_constr": dc.np.ones(6),
    "constr": dc.np.ones(6),
    "harm": dc.np.ones(6),
    "unconstr": dc.np.ones(6),
    "unrelaxed": dc.np.ones(6)
})

df1 = per_res_score.melt(id_vars=["#PDB"], var_name='method', value_name='score')
df2 = rmsds.melt(id_vars=['#PDB'], var_name='method', value_name='RMSD')

df1 = df1.set_index(
    ["#PDB", df1.groupby(["#PDB"]).cumcount()])
df2 = df2.set_index(
    ["#PDB", df1.groupby(["#PDB"]).cumcount()])


df3 = (dc.pd.concat([df1, df2], axis=1)
         .sort_index()
         .reset_index())
df3 = df3.T.drop_duplicates().T
print(df3)