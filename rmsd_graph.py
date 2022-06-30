"""RMSD Graph Calculations
Author: Karson Chrispens"""

from pyrosetta import *
from rosetta.core.scoring import *
import numpy as np
import pandas as pd
import seaborn as sns

pdbs = ["1DQJ", "1MHP", "1MLC", "1N8Z", "1VFB", "1YY9", "3GBN", "4FQY"]
crtsc_scores = np.array()
crtsc_nobb_scores = np.array()
harm_scores = np.array()
unrelaxed_scores = np.array()
no_constraints_scores = np.array()
sfxn = get_fa_scorefxn()

for pdb in pdbs:
    pose = pose_from_pdb(f"./PDBs/{pdb}.pdb")
    crtsc = pose_from_pdb(f"./PDBs/{pdb}_crtsc.pdb")
    crtsc_nobb = pose_from_pdb(f"./PDBs/{pdb}_clean.pdb")
    harm = pose_from_pdb(f"./PDBs/{pdb}_harm.pdb")
    unconst = pose_from_pdb(f"./PDBs/{pdb}_unconst.pdb")
    crtsc_scores = np.append(crtsc_scores, sfxn.score(crtsc) / crtsc.total_residue())
    crtsc_nobb_scores = np.append(crtsc_nobb_scores, sfxn.score(crtsc_nobb) / crtsc_nobb.total_residue())
    harm_scores = np.append(harm_scores, sfxn.score(harm) / harm.total_residue())
    unrelaxed_scores = np.append(unrelaxed_scores, sfxn.score(pose) / pose.total_residue())
    no_constraints_scores = np.append(no_constraints_scores, sfxn.score(unconst) / unconst.total_residue())

per_res_score = pd.DataFrame({
    "#PDB": pdbs,
    "sidechain_constr_scores": crtsc_nobb_scores,
    "constr_scores": crtsc_scores,
    "harm_scores": harm_scores,
    "unconstr_scores": no_constraints_scores,
    "unrelaxed_scores": unrelaxed_scores
})

rmsds = pd.read_csv("./raw_datasets/RMSD_graph_data.csv")
rmsds["unrelaxed"] = 0

df = pd.concat([rmsds, per_res_score], axis=1)

plot = sns.scatterplot()


