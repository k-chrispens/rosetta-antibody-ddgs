"""Calculating ddGs. Adapted from Brian Petersen's new_ddg_unbind.py
Author: Karson Chrispens"""

import re
import time
import pandas as pd
import math
import matplotlib.pyplot as plt
import random
import numpy as np
from pyrosetta.rosetta.core.select.movemap import *
from pyrosetta.rosetta.core.select import *
from pyrosetta.rosetta.protocols import *
from pyrosetta.rosetta.core.pack.task import *
import sys
from pyrosetta import *
# NOTE: Why use soft_rep_design?
init("-ex1 -ex2 -linmem_ig 10 -use_input_sc -soft_rep_design -mute all")

data = pd.read_csv("./raw_datasets/interface_data_use.csv")

def pack_and_relax(pose, posi, amino, repack_range, scorefxn):

    mut_posi = []
    mut_posi.append(
        pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector())
    mut_posi[0].set_index(posi[0])
    for i in range(1, len(posi)):
      mut_posi.append(
          pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector())
      mut_posi[i].set_index(posi[i])
    # print(pyrosetta.rosetta.core.select.get_residues_from_subset(mut_posi.apply(pose)))
    if len(posi) == 1:
      comb_select = mut_posi[0]
    if len(posi) >= 2:
      comb_select = pyrosetta.rosetta.core.select.residue_selector.OrResidueSelector(
          mut_posi[0], mut_posi[1])
    if len(posi) >= 3:
      for i in range(2, len(posi)):
        comb_select = pyrosetta.rosetta.core.select.residue_selector.OrResidueSelector(
            comb_select, mut_posi[i])

    # Select Neighbor Position
    nbr_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
    nbr_selector.set_distance(repack_range)
    nbr_selector.set_focus_selector(comb_select)
    nbr_selector.set_include_focus_in_subset(True)
    # print(pyrosetta.rosetta.core.select.get_residues_from_subset(nbr_selector.apply(pose)))

    # Select No Design Area
    not_design = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(
        comb_select)
    # print(pyrosetta.rosetta.core.select.get_residues_from_subset(not_design.apply(pose)))

    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
    tf.push_back(
        pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
    tf.push_back(
        pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())

    # Disable Packing
    prevent_repacking_rlt = pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT()
    prevent_subset_repacking = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
        prevent_repacking_rlt, nbr_selector, True)
    tf.push_back(prevent_subset_repacking)

    # Disable design
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
        pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), not_design))

    # Enable design (change the residues to the mutated residues)
    for i in range(len(posi)):
      aa_to_design = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
      aa_to_design.aas_to_keep(amino[i])
      tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
          aa_to_design,  mut_posi[i]))

    mmf = MoveMapFactory()
    mmf.add_bb_action(mm_enable, nbr_selector)
    mmf.add_chi_action(mm_enable, nbr_selector)
    mm = mmf.create_movemap_from_pose(pose)

    packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(
        scorefxn)
    packer.task_factory(tf)
    minmover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
    minmover.score_function(scorefxn)
    minmover.movemap_factory(mmf)
    for _ in range(2):
        packer.apply(pose)
        minmover.apply(pose)

    # Fast Relax
    fr = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn_in = scorefxn, standard_repeats = 2)
    fr.constrain_relax_to_start_coords(True)
    fr.set_task_factory(tf)
    fr.set_movemap(mm)
    fr.apply(pose)


def unbind(pose):
    STEP_SIZE = 100
    JUMP = 2
    trans_mover = rigid.RigidBodyTransMover(pose, JUMP)
    trans_mover.step_size(STEP_SIZE)
    trans_mover.apply(pose)


def calc_ddg(pose, pos, wt, mut, repack_range, output_pdb=False):

    scorefxn = get_score_function()
    mutPose = pose.clone()
    original = pose.clone()
    unbound_mutPose = pose.clone()
    unbound_original = pose.clone()

    # Bound unmutated
    pack_and_relax(original, pos, wt, repack_range, scorefxn)
    if output_pdb:
        original.dump_pdb("1_bound_unmutated.pdb")
    bound_unmutated = scorefxn(original)

    # Bound mutated
    pack_and_relax(mutPose, pos, mut, repack_range, scorefxn)
    if output_pdb:
        mutPose.dump_pdb("2_bound_mutated.pdb")
    bound_mutated = scorefxn(mutPose)

    # Unbound unmutated
    unbind(unbound_original)
    pack_and_relax(unbound_original, pos, wt, repack_range, scorefxn)
    if output_pdb:
        unbound_original.dump_pdb("3_unbound_unmutated.pdb")
    unbound_unmutated = scorefxn(unbound_original)

    # Unbound mutated
    unbind(unbound_mutPose)
    pack_and_relax(unbound_mutPose, pos, mut, repack_range, scorefxn)
    if output_pdb:
        unbound_mutPose.dump_pdb("4_unbound_mutated.pdb")
    unbound_mutated = scorefxn(unbound_mutPose)

    print("unbound_unmutated", unbound_unmutated)
    print("bound_unmutated", bound_unmutated)
    print("unbound_mutated", unbound_mutated)
    print("bound_mutated", bound_mutated)
    ddG = (bound_mutated - unbound_mutated) - \
        (bound_unmutated - unbound_unmutated)
    return ddG


df = pd.DataFrame(columns=["Position", "WT_AA", "Mut_AA", "DDG"])
# FIXME REMAKE FOR ALL PDBS
muts = re.split(";", rando["Mutations"].values[0]) # FIXME LATER to take out rando
print(muts)
scorefxn = get_fa_scorefxn()

# for i in muts: # FIXME LATER for all mutations
pos = []
wt = []
mut = []
# pdb = rando["#PDB"].values[0] # FIXME LATER to take out rando
pose = pose_from_pdb(f"./PDBs/1YY9_harm.pdb") #FIXME LATER to do all pdbs, hardcoding for now. Delete samples/rando later too
repack_range=12
all = list(map(lambda x: re.sub(
    r"(\w):(\w)(\d+)(\w)", r"\1:\2:\3:\4", x), muts))
print(all)
for i in all:
    chain, start, posi, muta = re.split(":", i)
    pos.append(pose.pdb_info().pdb2pose(chain, int(posi)))
    wt.append(start)
    mut.append(muta)

start=time.time()
print("Mutations:", rando["Mutations"].values[0])
total=calc_ddg(pose, pos, wt, mut, repack_range, False)
print("DDG: ", total)
df=df.append({"Position": pos, "WT_AA": wt,
                "Mut_AA": mut, "DDG": total}, ignore_index=True)
end=time.time()
print("Total time:", end-start, "seconds")

print(df)
