"""Harmonic Prerelax
Author: Karson Chrispens"""

from pyrosetta.rosetta.protocols.constraint_generator import *
import re
from pyrosetta.rosetta.core.select.movemap import *
from rosetta.core.pack.task import TaskFactory
from rosetta.protocols.relax import FastRelax
from rosetta.core.scoring import *
from pyrosetta import *
# import pandas as pd
init('-ex1 -ex2 -linmem_ig 10')  # add -ex1 -ex2

# So the process doesn't take as long
# data = pd.read_csv("./raw_datasets/interface_data_use.csv")
# pdbs = data["#PDB"].unique()

pdbs = ["1DQJ", "1MHP", "1MLC", "1N8Z", "1VFB", "1YY9", "3GBN", "4FQY"]
poses = []
for pdb in pdbs:
    pose = pose_from_pdb(f"./PDBs/{pdb}.pdb")
    poses.append(pose)

scorefxn = get_fa_scorefxn()
fr = FastRelax()
tf = TaskFactory()
tf.push_back(
    pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepacking())
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())

mmf = MoveMapFactory()
mmf.all_bb(True)
mmf.all_chi(True)

# Harmonic
apcg = AtomPairConstraintGenerator()
apcg.set_max_distance(9.0)
apcg.set_sd(0.25)  # changing the sd this time
apcg.set_ca_only(True)
apcg.set_use_harmonic_function(True)

add_csts = AddConstraints()
add_csts.add_generator(apcg)

harm_poses = []
for pose in poses:
    harm_pose = pose.clone()
    add_csts.apply(harm_pose)
    harm_poses.append(harm_pose)


scorefxn.set_weight(atom_pair_constraint, 1.0)
fr.set_scorefxn(scorefxn)
fr.set_task_factory(tf)
fr.set_movemap_factory(mmf)

for harm_pose in harm_poses:
    print("Before:", scorefxn(harm_pose))
    fr.apply(harm_pose)
    print("After:", scorefxn(harm_pose))
    # changed name for changed sd
    name = re.sub(r"(./PDBs/\w{4}).pdb",
                  r"\1_harm_025.pdb", harm_pose.pdb_info().name())
    harm_pose.dump_pdb(name)
