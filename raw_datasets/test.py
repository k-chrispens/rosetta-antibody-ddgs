"""Test File"""

import data_cleaners as dc

sipmab = dc.pd.read_excel(
    "~/rosetta-antibody-ddgs/raw_datasets/SiPMAB_subsetted.xlsx", sheet_name="Sheet1")
ab_bind = dc.pd.read_excel("~/rosetta-antibody-ddgs/raw_datasets/AB_bind_raw.xlsx",
                           sheet_name="AB-Bind_experimental_data.csv", usecols=[0, 1, 2, 3, 4, 5, 11, 12, 13])

# devo
ab_unique = ab_bind["#PDB"].unique()
sipmab_unique = sipmab["#PDB"].unique()

new_df = ab_bind

for i in sipmab_unique:
    if i in ab_unique:
        continue
    else:
        new_df = dc.pd.concat([new_df, sipmab.loc[sipmab["#PDB"] == i]])