"""Harmonic Prerelax
Author: Karson Chrispens"""

from pyrosetta.rosetta.protocols.constraint_generator import *
import re
from pyrosetta.rosetta.core.select.movemap import *
from rosetta.core.pack.task import TaskFactory
from rosetta.protocols.relax import FastRelax
from rosetta.core.scoring import *
from pyrosetta import *
import pandas as pd
init('-ex1 -ex2 -linmem_ig 10')  # add -ex1 -ex2

data = pd.read_csv("./raw_datasets/interface_data_use.csv")
pdbs = data["#PDB"].unique()
# So the process doesn't take as long
pdbs = pdbs.sample(n=10, random_state=42)
poses = []
for pdb in pdbs:
    pose = pose_from_pdb(f"./PDBs/{pdb}.pdb")
    poses.append(pose)

scorefxn = get_fa_scorefxn()
unconst_sfxn = get_fa_scorefxn()
fr = FastRelax()
unconstrained_fr = FastRelax()
tf = TaskFactory()
tf.push_back(
    pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepacking())
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())

mmf = MoveMapFactory()
mmf.all_bb(False)
mmf.all_chi(True)

# Harmonic
apcg = AtomPairConstraintGenerator()
apcg.set_max_distance(9.0)
apcg.set_sd(0.5)
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
unconstrained_fr.set_scorefxn(unconst_sfxn) # Not sure if this makes a difference
unconstrained_fr.set_task_factory(tf)
unconstrained_fr.set_movemap_factory(mmf)

for pose in poses:
    print("Before:", scorefxn(pose))
    fr.apply(pose)
    print("After:", scorefxn(pose))
    name = re.sub(r"(./PDBs/\w{4}).pdb", r"\1_harm.pdb", pose.pdb_info().name())
    pose.dump_pdb(name)

