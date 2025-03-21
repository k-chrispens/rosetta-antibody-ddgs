import re
from pyrosetta.rosetta.core.select.movemap import *
from rosetta.core.pack.task import TaskFactory
from rosetta.protocols.relax import FastRelax
from rosetta.core.scoring import *
from pyrosetta import *
import pandas as pd
init('-ex1 -ex2 -linmem_ig 10')  # add -ex1 -ex2

"""NOTE ON CONSTRAINTS: going to try overnight testing the relax.constrain_relax_to_start_coords(True),
but may want to actually loop through and do all the constraints for all the carbon atoms within 9, etc.
using a residue selector (to generate a .cst file?) then apply the cst file to the pose and score fxn.
NOTE: That didn't finish, it didn't even start fast relax. Almost got done with all the rotamers.
NOTE: Got done with 8 on Blanca over 24 hours."""


# data = pd.read_csv("./raw_datasets/interface_data_use.csv")
# pdbs = data["#PDB"].unique()
# poses = [] # NOTE COMMENTED OUT FOR NOW AS WE MAKE THE GRAPH
# for pdb in pdbs:
#     pose = pose_from_pdb(f"./PDBs/{pdb}.pdb")
#     poses.append(pose)

pdbs = ["1MLC", "1N8Z", "1VFB", "1YY9", "3GBN"] # "1DQJ", "1MHP", "4FQY" removed because I accidentally ran without the for loop at the end and then because blanca restarted or something
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

scorefxn.set_weight(atom_pair_constraint, 1.0)
fr.set_scorefxn(scorefxn)
fr.constrain_relax_to_start_coords(True)
fr.set_task_factory(tf)
fr.set_movemap_factory(mmf)

for pose in poses:
    print("Before:", scorefxn(pose))
    fr.apply(pose)
    print("After:", scorefxn(pose))
    name = re.sub(r"(./PDBs/\w{4}).pdb", r"\1_crtsc.pdb", pose.pdb_info().name())
    pose.dump_pdb(name)

# for pose in poses:
#     print("Before:", scorefxn(pose))
#     fr.apply(pose)
#     print("After:", scorefxn(pose))
#     pose.dump_pdb(f"{pose.pdb_info().name()}.clean")
