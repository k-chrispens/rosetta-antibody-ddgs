from pyrosetta.rosetta.core.select.movemap import *
from rosetta.core.pack.task import TaskFactory
from rosetta.protocols.relax import FastRelax
from pyrosetta import *
init('-ex1 -ex2 -linmem_ig 10')  # add -ex1 -ex2

scorefxn = get_fa_scorefxn()
pose = pose_from_pdb("./PDBs/3GBN.pdb")
fr = FastRelax()
tf = TaskFactory()
tf.push_back(
    pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepacking())
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())

mmf = MoveMapFactory()
mmf.all_bb(False)
mmf.all_chi(True)

fr.set_scorefxn(scorefxn)
fr.set_task_factory(tf)
fr.set_movemap_factory(mmf)
print("Before:", scorefxn(pose))
fr.apply(pose)
print("After:", scorefxn(pose))
pose.dump_pdb("3GBN_relax_test.pdb")
