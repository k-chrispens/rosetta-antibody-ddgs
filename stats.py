""""""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import *
from patsy import dmatrices
import statsmodels.api as sm
import seaborn as sns
import re

sns.set_context("poster", font_scale = 1.2)

data = pd.read_csv("~/rosetta-antibody-ddgs/raw_datasets/use_this_data.csv")

### Plotting average LD per PDB NOTE: for interface mutations


averages = pd.DataFrame()
averages["#PDB"] = data["#PDB"].unique()
averages_list = []
interface_data = data[data["Interface?"] == True]

for pdb in averages["#PDB"]:
    lds = interface_data[interface_data["#PDB"] == pdb]["LD"]
    averages_list.append(sum(lds) / len(lds))

averages["Average LD"] = averages_list
fig, ax = plt.subplots(figsize=(9, 13))
barplt = sns.barplot(x="Average LD", y="#PDB", data = averages)
plt.savefig("./images/average_LD_per_PDB.png")
plt.clf()

### Average LD per source

averages = pd.DataFrame()
averages["Source"] = data["Source"].unique()
averages_list = []

for source in averages["Source"]:
    lds = interface_data[interface_data["Source"] == source]["LD"]
    averages_list.append(sum(lds) / len(lds))

averages["Average LD"] = averages_list
fig, ax = plt.subplots(figsize=(13, 9))
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
fig, ax = plt.subplots(figsize=(11, 9))
sns.barplot(x = "Amino Acids", y = "Mutations to Amino Acid", data = data_muts_to_aa, palette="crest")
plt.title("Mutations per Amino Acid")
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

# data = pd.read_csv("./raw_datasets/interface_data_use.csv")
non_int = data[data["Interface?"] == False]
interface = data[data["Interface?"] == True]

non_int_avg = sum(non_int["ddG(kcal/mol)"]) / len(non_int)
interface_avg = sum(interface["ddG(kcal/mol)"]) / len(interface)
avgs = pd.DataFrame({
    "Mutation Position": ["Interface", "Non-Interface"],
    "Average ddG": [interface_avg, non_int_avg]
})
barplt = sns.barplot(x = "Mutation Position", y = "Average ddG", data=avgs)
plt.title("Average ΔΔG for Interface vs Non-Interface Mutations")
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
plt.title("ΔΔG of Interface and Non-Interface Mutations")
plt.xlabel("ΔΔG (kcal/mol)")
plt.savefig("./images/ddG_hist_interfaces.png")
plt.clf()

### ddG histogram for all mutations

fig, ax = plt.subplots(figsize=(11, 9))
histplt = sns.histplot(data = data["ddG(kcal/mol)"], kde = True, element="step", palette=colors.Colormap("crest"))
plt.title("ΔΔG of Mutations in Dataset")
plt.xlabel("ΔΔG (kcal/mol)")
plt.savefig("./images/ddG_hist.png")
plt.clf()

### Epistatic 1st order coeffs

coeff_dict = {
    "S29F": 0.3474748118778095,
    "N30S": 0.7850383407363816,
    "N31S": 0.6615118758001257,
    "S52I": 5.751941190773612,
    "S56T": 0.9459365035490297,
    "T57A": 0.6456928920891415,
    "A58N": -0.5386655036571883,
    "S70T": 0.049680007275027735,
    "I73K": 2.1620856562256816,
    "F74S": 2.366870250624797,
    "S75T": 0.1521897924114612,
    "N76S": 0.20765234633346258,
    "N82AS": 0.0618895688006494,
    "T83R": 0.14655872622630178,
    "F91Y": -0.0647305093343003,
    "S100BY": 0.28464140588828324,
}

coeffs = pd.DataFrame({
    "Coefficients": coeff_dict.values(),
    "Mutations": coeff_dict.keys()
})

barplt = sns.barplot(x="Mutations", y="Coefficients", data=coeffs, palette="rocket")
sns.despine(bottom=True)
barplt.axhline(y=0, color="black")
barplt.set_xlabel("Mutations (CR9114)")
barplt.set_ylabel("Model Coefficients (kcal/mol)")
plt.savefig("./images/coeffs_1_9114.png")
plt.clf()

coeff_dict = {
    "P28T": 1.3125487578097275,
    "R30S": 1.1271353818461578,
    "T57A": 0.05041079565338836,
    "K58N": 0.25192585367828907,
    "P61Q": 0.030367383609466614,
    "D73E": 0.08223430957471058,
    "F74S": 1.1464864109761814,
    "A75T": -0.10433763817940633,
    "G76S": 0.4610491870945383,
    "V78A": 0.29463292741194436,
    "V100L": 0.1561891795098677
}

coeffs = pd.DataFrame({
    "Coefficients": coeff_dict.values(),
    "Mutations": coeff_dict.keys()
})

barplt = sns.barplot(x="Mutations", y="Coefficients",
                     data=coeffs, palette="vlag")
sns.despine(bottom=True)
barplt.axhline(y=0, color="black")
barplt.set_xlabel("Mutations (CR6261)")
barplt.set_ylabel("Model Coefficients (kcal/mol)")
plt.savefig("./images/coeffs_1_6261.png")
plt.clf()

fig, ax = plt.subplots(figsize=(1, 5))
coeffs = pd.read_csv("./raw_datasets/9114_2order_1.csv", index_col="Mutation")
heatmap = sns.heatmap(data=coeffs, cmap="vlag", annot=True)
plt.title("Non-interacting\nCoefficients")
plt.savefig("./images/coeffs_2_1_9114_heatmap.png", bbox_inches="tight")
plt.clf()

coeffs = pd.read_csv("./raw_datasets/6261_2order_1.csv", index_col="Mutation")
heatmap = sns.heatmap(data=coeffs, cmap="vlag", annot=True)
plt.title("Non-interacting\nCoefficients")
plt.savefig("./images/coeffs_2_1_6261_heatmap.png", bbox_inches="tight")
plt.clf()


### Epistatic 2nd order coeffs

data = pd.read_csv("./raw_datasets/6261_2order.csv")
idx = np.unique(data[['Mut1','Mut2']].values.ravel())
data1 = data.pivot_table('Coefficient','Mut1','Mut2').reindex(index=idx, columns=idx)
data2 = data.pivot_table('Coefficient','Mut2','Mut1').reindex(index=idx, columns=idx)
data = data1.combine_first(data2)
data = data.fillna(1)
fig, ax = plt.subplots(figsize=(10, 5))
plot = sns.heatmap(data=data, cmap='vlag', center=0, annot=True)
plot.set_title("CR6261 2nd Order Coefficients (kcal/mol)")
plt.ylabel("")
plt.savefig("./images/coeffs_2_6261.png")
plt.clf()

data = pd.read_csv("./raw_datasets/9114_2order.csv")
idx = np.unique(data[['Mut1','Mut2']].values.ravel())
data1 = data.pivot_table('Coefficient','Mut1','Mut2').reindex(index=idx, columns=idx)
data2 = data.pivot_table('Coefficient','Mut2','Mut1').reindex(index=idx, columns=idx)
data = data1.combine_first(data2)
data = data.fillna(1)
fig, ax = plt.subplots(figsize=(10, 5))
plot = sns.heatmap(data=data, cmap='vlag', center=0, annot=True)
plot.set_title("CR9114 2nd Order Coefficients (kcal/mol)")
plt.ylabel("")
plt.savefig("./images/coeffs_2_9114.png")
plt.clf()

## Correlation Plots

ddgs = pd.read_excel("./FLEX_RUNS.xlsx", "Sheet1")
ddgs = ddgs.dropna(subset="8 50 r s")
print(ddgs.head())
y, X = dmatrices("Q('8 50 r s') ~ Q('ddG(kcal/mol)')",
                 data=ddgs, return_type='dataframe')
model = sm.OLS(y, X)
results = model.fit()
print(results.summary())
fig, ax = plt.subplots(figsize=(13, 10))
sns.regplot(y="8 50 r s", x="ddG(kcal/mol)", data=ddgs, color="darkcyan", truncate=False, scatter_kws={"s": 20})
plt.title("50 Structures, 5000 Backrub Steps")
plt.xlabel("Experimental ΔΔG (kcal/mol)")
plt.ylabel("Predicted ΔΔG (kcal/mol)")
plt.savefig("./images/current_best_corr.png")
plt.clf()

ddgs = pd.read_csv("./analysis_output/aligned_no_gam.csv")
y, X = dmatrices("Q('total_score') ~ Q('ddG(kcal/mol)')",
                 data=ddgs, return_type='dataframe')
model = sm.OLS(y, X)
results = model.fit()
print(results.summary())
fig, ax = plt.subplots(figsize=(13, 10))
sns.regplot(y="total_score", x="ddG(kcal/mol)", data=ddgs, color="darkcyan", truncate=False, scatter_kws={"s": 20})
plt.title("10 Structures, 10 Backrub Steps")
plt.xlabel("Experimental ΔΔG (kcal/mol)")
plt.ylabel("Predicted ΔΔG (kcal/mol)")
plt.savefig("./images/10_10_no_gam_corr.png")
plt.clf()

ddgs = pd.read_csv("./analysis_output/aligned_gam.csv")
y, X = dmatrices("Q('total_score') ~ Q('ddG(kcal/mol)')",
                 data=ddgs, return_type='dataframe')
model = sm.OLS(y, X)
results = model.fit()
print(results.summary())
fig, ax = plt.subplots(figsize=(13, 10))
sns.regplot(y="total_score", x="ddG(kcal/mol)", data=ddgs, color="darkcyan", truncate=False, scatter_kws={"s": 20})
plt.title("10 Structures, 10 Backrub Steps, GAM Reweight")
plt.xlabel("Experimental ΔΔG (kcal/mol)")
plt.ylabel("Predicted ΔΔG (kcal/mol)")
plt.savefig("./images/10_10_gam_corr.png")
plt.clf()

# Numbers of mutations in certain bins

data = pd.read_csv("./raw_datasets/use_this_data.csv")

df = pd.DataFrame({
    "Types": ["Single Mutation to Alanine", "Multiple Mutations", "Multiple Mutations (None Alanine)", "Small to Large Mutations"],
    "Counts": [data["Single_mut_A"].sum(), data["Multiple_mut"].sum(), data["Multiple_mut_no_A"].sum(), data["Small_to_large"].sum()]
})

fig, ax = plt.subplots(figsize=(13, 13))
mut_types = sns.barplot(x="Types", y="Counts", data=df, palette="crest")
plt.xticks(rotation=45, ha='right', rotation_mode="anchor")
plt.title("Types of Mutations in Dataset")
mut_types.set(xlabel=None)
plt.tight_layout()
plt.savefig("./images/mut_types.png")
plt.clf()