"""Calculating ddGs.
REQUIRED - cmd line args: int repack_range, int packmin_ensemble_size, 
flag beta (scorefunction), str path (to output), flag cartesian, flag soft_rep, flag all_repack
Author: Karson Chrispens"""

import re
import time
import pandas as pd
import numpy as np
from pyrosetta.rosetta.core.select.movemap import *
from pyrosetta.rosetta.core.select import *
from pyrosetta.rosetta.protocols import *
from pyrosetta.rosetta.core.pack.task import *
from pyrosetta.rosetta.core.import_pose import *
import getopt, sys
from pyrosetta import *

args = sys.argv[1:]
options = "rpbo:csa"
long_options = ["repack_range=", "rounds_packmin=",
                "beta", "output_path=", "cartesian", "soft_rep", "all_repack"]
values_dict = {"r": 12, "p": 2, "b": False, "o": "./UNNAMED.csv", "c": False, "s": False, "a": False}

try:
    # Parsing argument
    arguments, values = getopt.getopt(args, options, long_options)
    
    print("args: ", args) # FIXME
    # checking each argument
    for currentArgument, currentValue in arguments:
        print("CHECK") # FIXME
        if currentArgument in ("-r", "--repack_range"):
            values_dict["r"] = currentValue
            print(f"Repack Range = {currentValue}")

        elif currentArgument in ("-p", "--rounds_packmin"):
            values_dict["p"] = currentValue
            print(f"Rounds of Pack and Minimization = {currentValue}")

        elif currentArgument in ("-b", "--beta"):
            values_dict["b"] = True
            print(f"Using beta score function")

        elif currentArgument in ("-o", "--output_path"):
            values_dict["o"] = currentValue
            print(f"Output path = {currentValue}")

        elif currentArgument in ("-c", "--cartesian"):
            values_dict["c"] = True
            print("Enabling cartesian minimization flag.")

        elif currentArgument in ("s", "--soft_rep"):
            values_dict["s"] = True
            print("Using soft_rep_design flag")

        elif currentArgument in ("a", "--all_repack"):
            values_dict["a"] = True
            print("Repacking all residues")

except getopt.error as err:
    # output error, and return with an error code
    print(str(err))

if values_dict["o"] == "./UNNAMED.csv":
    print("SOMETHING BREAKING")
    quit()

data = pd.read_csv("./raw_datasets/use_this_data.csv")

# INIT
# NOTE: Why use soft_rep_design?
if values_dict["b"] and values_dict["s"]:
    init("-beta -ex1 -ex2 -linmem_ig 10 -use_input_sc -soft_rep_design -mute all")
elif values_dict["b"]:
    init("-beta -ex1 -ex2 -linmem_ig 10 -use_input_sc -mute all")
elif values_dict["s"]:
    init("-ex1 -ex2 -linmem_ig 10 -use_input_sc -soft_rep_design -mute all")
else:
    init("-ex1 -ex2 -linmem_ig 10 -use_input_sc -mute all")


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

    # LOCAL OR GLOBAL REPACKING
    if not values_dict["a"]:
        # Disable Packing
        prevent_repacking_rlt = pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT()
        prevent_subset_repacking = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
        prevent_repacking_rlt, nbr_selector, True)
        tf.push_back(prevent_subset_repacking)
        # Prevent Design
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
            pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), not_design))
    else:
        # Prevent Design
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
            pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), not_design))

    # Enable design (change the residues to the mutated residues)
    for i in range(len(posi)):
        aa_to_design = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
        aa_to_design.aas_to_keep(amino[i])
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
            aa_to_design,  mut_posi[i]))

    # Minimizes residues within repack range of mutations
    mmf = MoveMapFactory()
    mmf.add_bb_action(mm_enable, nbr_selector)
    mmf.add_chi_action(mm_enable, nbr_selector)
    # mm = mmf.create_movemap_from_pose(pose) # ONLY NEEDED IF FAST RELAX

    packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(
        scorefxn)
    packer.task_factory(tf)
    minmover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
    minmover.score_function(scorefxn)
    minmover.movemap_factory(mmf)
    minmover.max_iter(2000)
    minmover.tolerance(0.00001)

    if values_dict["c"]:
        minmover.cartesian(True)
    
    packer.apply(pose)
    minmover.apply(pose)

    # Fast Relax FIXME
    # fr = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn_in = scorefxn, standard_repeats = 2)
    # fr.constrain_relax_to_start_coords(True)
    # fr.set_task_factory(tf)
    # fr.set_movemap(mm)
    # fr.apply(pose)


def unbind(pose, jump):
    STEP_SIZE = 100
    # JUMP WILL NEED TO ADJUST FOR EACH POSE.
    trans_mover = rigid.RigidBodyTransMover(pose, jump)
    trans_mover.step_size(STEP_SIZE)
    trans_mover.apply(pose)


def calc_ddg(pose, pos, wt, mut, repack_range, jump, output_pdb=False):

    # TESTING COPY VS CLONE
    mutPose = Pose()
    original = Pose()
    unbound_mutPose = Pose()
    unbound_original = Pose()
    mutPose.detached_copy(pose)
    original.detached_copy(pose)
    unbound_mutPose.detached_copy(pose)
    unbound_original.detached_copy(pose)

    # Bound unmutated
    pack_and_relax(original, pos, wt, repack_range, scorefxn)
    if output_pdb:
        original.dump_pdb("1_bound_unmutated.pdb")
    bound_unmutated = ddg_scorefxn(original)

    # Bound mutated
    pack_and_relax(mutPose, pos, mut, repack_range, scorefxn)
    if output_pdb:
        mutPose.dump_pdb("2_bound_mutated.pdb")
    bound_mutated = ddg_scorefxn(mutPose)

    # Unbound unmutated
    unbind(unbound_original, jump)
    pack_and_relax(unbound_original, pos, wt, repack_range, scorefxn)
    if output_pdb:
        unbound_original.dump_pdb("3_unbound_unmutated.pdb")
    unbound_unmutated = ddg_scorefxn(unbound_original)

    # Unbound mutated
    unbind(unbound_mutPose, jump)
    pack_and_relax(unbound_mutPose, pos, mut, repack_range, scorefxn)
    if output_pdb:
        unbound_mutPose.dump_pdb("4_unbound_mutated.pdb")
    unbound_mutated = ddg_scorefxn(unbound_mutPose)

    print("unbound_unmutated", unbound_unmutated)
    print("bound_unmutated", bound_unmutated)
    print("unbound_mutated", unbound_mutated)
    print("bound_mutated", bound_mutated)
    ddG = (bound_mutated - unbound_mutated) - \
        (bound_unmutated - unbound_unmutated)
    return ddG


pdbs = data["#PDB"].unique()
df = pd.DataFrame(columns=["#PDB", "Position", "WT_AA", "Mut_AA", "DDG"])

# need cartesian score function for minimization (if cart is chosen.)
if values_dict["b"] and values_dict["c"]:
    scorefxn = create_score_function("beta_nov16_cart.wts")
elif values_dict["c"]:
    scorefxn = create_score_function(
        "ref2015_cart.wts")
else:
    scorefxn = get_score_function()

# ddG score function should be regular, as the cart term can vary largely between structures.
ddg_scorefxn = get_score_function()

repack_range = int(values_dict["r"])  # Try 12, where did 8 come from?
# TO ALLOW PARALLEL RUNS AND TESTS: initial run was pdbs[:8], next run is pdbs[8:20], next after is [20:30],
# then [30:38]. These were generated based on approx times I wanted to let them run.
count = 0

for pdb in pdbs:
    points = data.loc[data["#PDB"] == pdb]
    points = points.loc[points["Interface?"] == True]
    pose = get_pdb_and_cleanup(f"./PDBs/{pdb}_all.pdb")
    for index, point in points.iterrows():
        muts = re.split(";", point["Mutations"])
        jump = point["Jump"]
        pos = []
        wt = []
        mut = []
        all = list(map(lambda x: re.sub(
            r"(\w):(\w)(\d+)(\w*)(\w)", r"\1:\2:\3:\5:\4", x), muts))
        print(all)
        for i in all:
            chain, start, posi, muta, ic = re.split(":", i)
            if ic:
                pos.append(pose.pdb_info().pdb2pose(chain, int(posi), ic))
            else:
                pos.append(pose.pdb_info().pdb2pose(chain, int(posi)))
            wt.append(start)
            mut.append(muta)

        start = time.time()
        print("Mutations:", point["Mutations"])
        total = 0
        for p in range(int(values_dict["p"])):
            total += calc_ddg(pose, pos, wt, mut, repack_range, jump, False)
        total = total / int(values_dict["p"])
        print("DDG: ", total)
        df = df.append({"#PDB": pdb, "Position": pos, "WT_AA": wt,
                        "Mut_AA": mut, "DDG": total}, ignore_index=True)
        end = time.time()
        print("Total time:", end-start, "seconds")
        print("Avg time per ensemble member:", (end-start)/int(values_dict["p"]), "seconds")
        count += 1
        if count % 20 == 0:
            df.to_csv(values_dict["o"], index=False)

df.to_csv(values_dict["o"], index=False)
