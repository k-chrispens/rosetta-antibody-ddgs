"""RMSD Graph Calculations using PyMOL RMSDs
Author: Karson Chrispens"""

from pyrosetta import *
from pyrosetta.rosetta.core.scoring import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
init()

pdbs = ["1DQJ", "1MHP", "1MLC", "1N8Z", "1VFB", "1YY9", "3GBN",
        "4FQY"]
crtsc_scores = []
crtsc_nobb_scores = []
harm_scores = []
unrelaxed_scores = []
no_constraints_scores = []
harm025_scores = []
sfxn = get_fa_scorefxn()

data = pd.read_csv("./PDBs/rmsds.csv")

rmsds = pd.DataFrame(
    columns=["#PDB", "harm", "harm_025", "unconst", "crtsc", "clean"])

for pdb in data["#PDB"].unique():
    data_pdb = data[data["#PDB"] == pdb]
    for type in data["type"].unique():
        rmsd_data = pd.DataFrame({
            f"{type}": data_pdb[data_pdb["type"] == type]["RMSD"],
            "#PDB": pdb
        })
        rmsds = pd.concat([rmsds, rmsd_data])

rmsds = rmsds.set_index("#PDB").groupby(level=0).transform(
    lambda x: sorted(x, key=lambda k: pd.isna(k))).dropna(axis=0, how="all").reset_index()


for pdb in pdbs:
    pose = pose_from_pdb(f"./PDBs/{pdb}.pdb")
    crtsc = pose_from_pdb(f"./PDBs/{pdb}_crtsc.pdb")
    crtsc_nobb = pose_from_pdb(f"./PDBs/{pdb}_clean.pdb")
    harm = pose_from_pdb(f"./PDBs/{pdb}_harm.pdb")
    unconst = pose_from_pdb(f"./PDBs/{pdb}_unconst.pdb")
    harm025 = pose_from_pdb(f"./PDBs/{pdb}_harm_025.pdb")
    crtsc_scores.append(sfxn.score(crtsc) / crtsc.total_residue())
    crtsc_nobb_scores.append(sfxn.score(crtsc_nobb) /
                             crtsc_nobb.total_residue())
    harm_scores.append(sfxn.score(harm) / harm.total_residue())
    unrelaxed_scores.append(sfxn.score(pose) / pose.total_residue())
    no_constraints_scores.append(sfxn.score(unconst) / unconst.total_residue())
    harm025_scores.append(sfxn.score(harm025) / unconst.total_residue())

per_res_score = pd.DataFrame({
    "#PDB": pdbs,
    "clean": crtsc_nobb_scores,
    "crtsc": crtsc_scores,
    "harm": harm_scores,
    "unconstr": no_constraints_scores,
    "unrelaxed": unrelaxed_scores,
    "harm_025": harm025_scores
})

df1 = per_res_score.melt(
    id_vars=["#PDB"], var_name='Method', value_name='Per Residue Score (REU)')
df2 = rmsds.melt(id_vars=['#PDB'], var_name='Method', value_name='RMSD (Å)')

df1 = df1.set_index(
    ["#PDB", df1.groupby(["#PDB"]).cumcount()])
df2 = df2.set_index(
    ["#PDB", df1.groupby(["#PDB"]).cumcount()])


df3 = (pd.concat([df1, df2], axis=1)
         .sort_index()
         .reset_index())
df3 = df3.T.drop_duplicates().T
print(df3)

plot = sns.scatterplot(x="RMSD (Å)", y="Per Residue Score (REU)",
                       hue="Method", style="#PDB", data=df3)
plt.legend(bbox_to_anchor=(1, 1), loc='upper left', borderaxespad=0)
plt.savefig("./images/rmsd_pymol_plot.png", bbox_inches='tight')
