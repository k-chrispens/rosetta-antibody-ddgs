"""Test File"""

import data_cleaners as dc

skempi = dc.pd.read_csv(
    "~/rosetta-antibody-ddgs/raw_datasets/skempi_v2_cleaned.csv")

# devo
print(skempi.head())

skempi["#Pdb"] = skempi["#Pdb"].apply(
    lambda x: dc.re.sub(r"(\w{4}).*", r"\1", x))
skempi["Mutation(s)_PDB"] = skempi["Mutation(s)_PDB"].apply(
    lambda x: dc.re.sub(r"(\w)(\w)(\w+)", r"\2:\1\3", dc.re.sub(r",", r";", x)))
prot1 = skempi["Protein 1"].str.contains(
    "fab|mab", case=False) | skempi["Protein 1"].str.contains("antibody|Fv", case=False)
prot2 = skempi["Protein 2"].str.contains(
    "fab|mab", case=False) | skempi["Protein 2"].str.contains("antibody|Fv", case=False)

abs = prot1 | prot2
print(sum(abs))
skempi = skempi[abs]

skempi = skempi.assign(ddG = lambda x: dc.ddg_from_kd(
    x["Affinity_mut_parsed"], x["Temperature"], x["Affinity_wt_parsed"]))

print(skempi["ddG"])
