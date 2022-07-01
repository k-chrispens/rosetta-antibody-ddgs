"""RMSD Graph Calculations
Author: Karson Chrispens"""

from pyrosetta import *
from pyrosetta.rosetta.core.scoring import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

pdbs = ["1DQJ", "1MHP", "1MLC", "1N8Z", "1VFB", "4FQY"] # "1YY9", "3GBN", removed for now since they aren't done
crtsc_scores = []
crtsc_nobb_scores = []
harm_scores = []
unrelaxed_scores = []
no_constraints_scores = []
crtsc_rmsds = []
crtsc_nobb_rmsds = []
harm_rmsds = []
unrelaxed_rmsds = []
no_constraints_rmsds = []
sfxn = get_fa_scorefxn()

for pdb in pdbs:
    pose = pose_from_pdb(f"./PDBs/{pdb}.pdb")
    crtsc = pose_from_pdb(f"./PDBs/{pdb}_crtsc.pdb")
    crtsc_nobb = pose_from_pdb(f"./PDBs/{pdb}_clean.pdb")
    harm = pose_from_pdb(f"./PDBs/{pdb}_harm.pdb")
    unconst = pose_from_pdb(f"./PDBs/{pdb}_unconst.pdb")
    crtsc_scores.append(sfxn.score(crtsc) / crtsc.total_residue())
    crtsc_nobb_scores.append(sfxn.score(crtsc_nobb) / crtsc_nobb.total_residue())
    harm_scores.append(sfxn.score(harm) / harm.total_residue())
    unrelaxed_scores.append(sfxn.score(pose) / pose.total_residue())
    no_constraints_scores.append(sfxn.score(unconst) / unconst.total_residue())
    crtsc_rmsds.append(all_atom_rmsd(pose, crtsc))
    crtsc_nobb_rmsds.append(all_atom_rmsd(pose, crtsc_nobb))
    harm_rmsds.append(all_atom_rmsd(pose, harm))
    unrelaxed_rmsds.append(0)
    no_constraints_rmsds.append(all_atom_rmsd(pose, unconst))

per_res_score = pd.DataFrame({
    "#PDB": pdbs,
    "sidechain_constr_scores": crtsc_nobb_scores,
    "constr_scores": crtsc_scores,
    "harm_scores": harm_scores,
    "unconstr_scores": no_constraints_scores,
    "unrelaxed_scores": unrelaxed_scores
})

rmsds = pd.DataFrame({
    "#PDB": pdbs,
    "sidechain_constr": crtsc_nobb_rmsds,
    "constr": crtsc_rmsds,
    "harm": harm_rmsds,
    "unconstr": no_constraints_rmsds,
    "unrelaxed": unrelaxed_rmsds
})


df = pd.concat([rmsds, per_res_score], axis=1)

df1 = df.melt(id_vars=["#PDB"],
              value_vars=['sidechain_constr_scores', "constr_scores",
                          "harm_scores", "unconstr_scores", "unrelaxed_scores"],
              var_name='#PDB', value_name='Per_Res_Score')
df2 = df.melt(id_vars=['#PDB'],
              value_vars=['sidechain_constr', "constr",
                          "harm", "unconstr", "unrelaxed"],
              var_name='#PDB', value_name='RMSD')

df1 = df1.set_index(
    ["#PDB", df1.groupby(["#PDB"]).cumcount()])
df2 = df2.set_index(
    ["#PDB", df1.groupby(["#PDB"]).cumcount()])


df3 = (pd.concat([df1, df2], axis=1)
         .sort_index(level=2)
         .reset_index(level=2, drop=True)
         .reset_index())

plot = sns.scatterplot(x="RMSD", y="Per_Res_Score", data=df3)
plt.savefig("./images/rmsd_plot.png")
