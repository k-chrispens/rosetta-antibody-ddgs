"""RMSD Graph Calculations
Author: Karson Chrispens"""

from pyrosetta import *
from rosetta.core.scoring import *
import numpy as np
import pandas as pd
import seaborn as sns

pdbs = ["1DQJ", "1MHP", "1MLC", "1N8Z", "1VFB", "4FQY"] # "1YY9", "3GBN", removed for now since they aren't done
crtsc_scores = np.array()
crtsc_nobb_scores = np.array()
harm_scores = np.array()
unrelaxed_scores = np.array()
no_constraints_scores = np.array()
crtsc_rmsds = np.array()
crtsc_nobb_rmsds = np.array()
harm_rmsds = np.array()
unrelaxed_rmsds = np.array()
no_constraints_rmsds = np.array()
sfxn = get_fa_scorefxn()

for pdb in pdbs:
    pose = pose_from_pdb(f"./PDBs/{pdb}.pdb")
    crtsc = pose_from_pdb(f"./PDBs/{pdb}_crtsc.pdb")
    crtsc_nobb = pose_from_pdb(f"./PDBs/{pdb}_clean.pdb")
    harm = pose_from_pdb(f"./PDBs/{pdb}_harm.pdb")
    unconst = pose_from_pdb(f"./PDBs/{pdb}_unconst.pdb")
    crtsc_scores = np.append(
        crtsc_scores, sfxn.score(crtsc) / crtsc.total_residue())
    crtsc_nobb_scores = np.append(crtsc_nobb_scores, sfxn.score(
        crtsc_nobb) / crtsc_nobb.total_residue())
    harm_scores = np.append(
        harm_scores, sfxn.score(harm) / harm.total_residue())
    unrelaxed_scores = np.append(
        unrelaxed_scores, sfxn.score(pose) / pose.total_residue())
    no_constraints_scores = np.append(
        no_constraints_scores, sfxn.score(unconst) / unconst.total_residue())
    crtsc_rmsds = np.append(
        crtsc_rmsds, all_atom_rmsd(pose, crtsc))
    crtsc_nobb_rmsds = np.append(
        crtsc_nobb_rmsds, all_atom_rmsd(pose, crtsc_nobb))
    harm_rmsds = np.append(
        harm_rmsds, all_atom_rmsd(pose, harm))
    unrelaxed_rmsds = np.append(
        unrelaxed_rmsds, 0)
    no_constraints_rmsds = np.append(
        no_constraints_rmsds, all_atom_rmsd(pose, unconst))

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
