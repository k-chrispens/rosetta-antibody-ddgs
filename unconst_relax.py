"""Unconstrained Prerelax script
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

unconst_sfxn = get_fa_scorefxn()
unconstrained_fr = FastRelax()
tf = TaskFactory()
tf.push_back(
    pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepacking())
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())

mmf = MoveMapFactory()
mmf.all_bb(True)
mmf.all_chi(True)

unconstrained_fr.set_scorefxn(unconst_sfxn)
unconstrained_fr.set_task_factory(tf)
unconstrained_fr.set_movemap_factory(mmf)

for pose in poses:
    print("Before:", unconst_sfxn(pose))
    unconstrained_fr.apply(pose)
    print("After:", unconst_sfxn(pose))
    name = re.sub(r"(./PDBs/\w{4}).pdb",
                  r"\1_unconst.pdb", pose.pdb_info().name())
    pose.dump_pdb(name)
