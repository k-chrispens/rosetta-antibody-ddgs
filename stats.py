"""Running statistics on full dataset to get some of the characteristics.
Author: Karson Chrispens
Requires ab_proj env."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re

sns.set_context("notebook", font_scale = 0.6)

data = pd.read_csv("~/rosetta-antibody-ddgs/raw_datasets/interface_data_use.csv")

### Plotting average LD per PDB

averages = pd.DataFrame()
averages["#PDB"] = data["#PDB"].unique()
averages_list = []

for pdb in averages["#PDB"]:
    lds = data[data["#PDB"] == pdb]["LD"]
    averages_list.append(sum(lds) / len(lds))

averages["Average LD"] = averages_list
barplt = sns.barplot(x="Average LD", y="#PDB", data = averages)
plt.savefig("./images/average_LD_per_PDB.png")
plt.clf()

### Average LD per source

averages = pd.DataFrame()
averages["Source"] = data["Source"].unique()
averages_list = []

for source in averages["Source"]:
    lds = data[data["Source"] == source]["LD"]
    averages_list.append(sum(lds) / len(lds))

averages["Average LD"] = averages_list
barplt = sns.barplot(x="Source", y="Average LD", data=averages)
plt.savefig("./images/average_LD_per_source.png")
plt.clf()


### Plotting number of mutations to each AA in the whole dataset

aas = ['C', 'D', 'S', 'Q', 'K',
     'I', 'P', 'T', 'F', 'N',
     'G', 'H', 'L', 'R', 'W',
     'A', 'V', 'E', 'Y', 'M']
muts_to_aa = []

for aa in aas:
    # count regex matches in mutations
    all_mut_aa = data["Mutations"].apply(lambda x: re.findall(
        fr"\w:\w\d+{aa}", x))
    num_mut_per_data = all_mut_aa.apply(len)
    total_mut = sum(num_mut_per_data)
    muts_to_aa.append(total_mut)

data_muts_to_aa = pd.DataFrame(data = {"Amino Acids": aas, "Mutations to Amino Acid": muts_to_aa})
sns.barplot(x = "Amino Acids", y = "Mutations to Amino Acid", data = data_muts_to_aa)
plt.savefig("./images/mutations_per_aa.png")
plt.clf()

### Plotting mutations per PDB

mut_pdbs = pd.DataFrame()
mut_pdbs["#PDB"] = data["#PDB"].unique()
mut_pdbs_list = []

for pdb in mut_pdbs["#PDB"]:
    lds = data[data["#PDB"] == pdb]["LD"]
    mut_pdbs_list.append(sum(lds))

mut_pdbs["Number of Mutations"] = mut_pdbs_list
barplt = sns.barplot(x="Number of Mutations", y="#PDB", data = mut_pdbs)
plt.savefig("./images/mutations_per_PDB.png")
plt.clf()

### Mutations per PDB without Phillips

mut_pdbs = mut_pdbs[~(mut_pdbs["#PDB"].str.contains("3GBN") | mut_pdbs["#PDB"].str.contains("4FQY"))]
barplt = sns.barplot(x="Number of Mutations", y="#PDB", data=mut_pdbs)
plt.savefig("./images/mutations_per_PDB_no_phillips.png")
plt.clf()

### Average ddG for interface vs non-interface mutations

data = pd.read_csv("./raw_datasets/interface_data_use.csv")
non_int = data[data["Interface?"] == False]
interface = data[data["Interface?"] == True]

non_int_avg = sum(non_int["ddG(kcal/mol)"]) / len(non_int)
interface_avg = sum(interface["ddG(kcal/mol)"]) / len(interface)
avgs = pd.DataFrame({
    "Mutation Position": ["Interface", "Non-Interface"],
    "Average ddG": [interface_avg, non_int_avg]
})
barplt = sns.barplot(x = "Mutation Position", y = "Average ddG", data=avgs)
plt.savefig("./images/avg_ddG_interfaces.png")
plt.clf()

### ddG histogram for interface vs non-interface

non_int_ddgs = non_int["ddG(kcal/mol)"]
interface_ddgs = interface["ddG(kcal/mol)"]

df = pd.DataFrame({
    "Non-Interface": non_int_ddgs,
    "Interface": interface_ddgs
})

histplt = sns.histplot(data = df, kde = True, element="step")
plt.xlabel("ΔΔG (kcal/mol)")
plt.savefig("./images/ddG_hist_interfaces.png")
plt.clf()