"""Calculating ddGs.
REQUIRED - cmd line args: int repack_range, int backrub_ensemble_size, 
flag beta (scorefunction), str path (to output), flag cartesian, 
flag soft_rep, flag all_repack, str pdbs (comma-separated)
Author: Karson Chrispens"""

from pyrosetta import *
from pyrosetta.rosetta.core.import_pose import *
from pyrosetta.rosetta.protocols import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.select.movemap import *
from pyrosetta.rosetta.core.scoring import *
import re
import time
import numpy as np
import pandas as pd
import getopt
import sys

args = sys.argv[1:]
options = "r:p:bo:csan:t:"
long_options = ["repack_range=", "backrub=",
                "beta", "output_path=", "cartesian", "soft_rep", "all_repack", "pdbs=", "steps="]
values_dict = {"r": 8, "p": 1, "b": False,
               "o": "./UNNAMED.csv", "c": False, "s": False, "a": False, "n": "all", "t": 5000}

try:
    # Parsing argument
    arguments, values = getopt.getopt(args, options, long_options)

    # checking each argument
    for currentArgument, currentValue in arguments:
        if currentArgument in ("-r", "--repack_range"):
            values_dict["r"] = currentValue
            print(f"Repack Range = {currentValue}")

        elif currentArgument in ("-p", "--backrub"):
            values_dict["p"] = currentValue
            print(f"Backrub Ensemble = {currentValue}")

        elif currentArgument in ("-b", "--beta"):
            values_dict["b"] = True
            print(f"Using beta score function")

        elif currentArgument in ("-o", "--output_path"):
            values_dict["o"] = currentValue
            print(f"Output path = {currentValue}")

        elif currentArgument in ("-c", "--cartesian"):
            values_dict["c"] = True
            print("Enabling cartesian minimization flag.")

        elif currentArgument in ("-s", "--soft_rep"):
            values_dict["s"] = True
            print("Using soft_rep_design flag")

        elif currentArgument in ("-a", "--all_repack"):
            values_dict["a"] = True
            print("Repacking all residues")
        
        elif currentArgument in ("-n", "--pdbs"):
            values_dict["n"] = re.split(",", currentValue)
            print("PDBs:", currentValue)

        elif currentArgument in ("-t", "--steps"):
            values_dict["t"] = currentValue
            print(f"Backrub Steps = {currentValue}")

except getopt.error as err:
    # output error, and return with an error code
    print(str(err))

data = pd.read_csv("./raw_datasets/use_this_data.csv")

# INIT
# NOTE: Why use soft_rep_design? Kellogg et al. 2011
if values_dict["b"] and values_dict["s"]:
    pyrosetta.init(
        "-beta -ex1 -ex2 -linmem_ig 10 -use_input_sc -soft_rep_design -mute all -backrub:mc_kt 1.2 -backrub:ntrials {} -nstruct 1".format(values_dict["t"]))
elif values_dict["b"]:
    pyrosetta.init("-beta -ex1 -ex2 -linmem_ig 10 -use_input_sc -mute all -backrub:mc_kt 1.2 -backrub:ntrials {} -nstruct 1".format(values_dict["t"]))
elif values_dict["s"]:
    pyrosetta.init("-ex1 -ex2 -linmem_ig 10 -use_input_sc -soft_rep_design -mute all -backrub:mc_kt 1.2 -backrub:ntrials {} -nstruct 1".format(values_dict["t"]))
else:
    # FIXME put backrub flags here  -mc_kt 1.2 -nstruct 50 (only for use with job distributor) -backrub:ntrials 35000
    pyrosetta.init(
        "-ex1 -ex2 -linmem_ig 10 -use_input_sc -mute all -backrub:mc_kt 1.2 -backrub:ntrials {} -nstruct 1".format(values_dict["t"]))
print(values_dict)
data = pd.read_csv("./raw_datasets/use_this_data.csv")


# TESTING (Seems like nbr_selector and scorefxn params are not required, but need to be required if this is a unit function)
def backrub_ensemble_gen(pose, nbr_selector, mmf, tf, scorefxn):

    start = time.time()
    backrubber = backrub.BackrubMover()
    backrubber.init_with_options()
    backrubber.set_movemap_factory(mmf)
    backrubber.set_min_atoms(3)
    backrubber.set_max_atoms(12)

    backrub_protocol = backrub.BackrubProtocol()
    backrub_protocol.set_backrub_mover(backrubber)
    backrub_protocol.set_taskfactory(tf)
    
    backrub_protocol.apply(pose)
    end = time.time()
    print("Backrub time: ", end-start, "seconds")

    # GMC = pyrosetta.rosetta.protocols.monte_carlo.GenericMonteCarloMover(mover=backrubber, scorefxn=scorefxn, maxtrials=50000, temperature=1.2, max_accepted_trials=)
    # GMC.set_mover(backrub)
    # GMC.set_scorefxn(scorefxn)
    # GMC.set_maxtrials(500)
    # GMC.set_temperature(1.0)
    # GMC.set_preapply(False)
    # GMC.set_recover_low(True)
    # GMC.apply(pose)

    # ensemble = []
    # for _ in range(values_dict["p"]):
    #     clone = Pose()
    #     clone.detached_copy(pose)
    #     backrubber.apply(clone)
    #     ensemble.append(clone)

    # return ensemble


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

    # Selecting only residues on current chain (restricts repacking to antibody only, assumes antigen is not as important)
    chain_selector = pyrosetta.rosetta.core.select.residue_selector.ChainSelector(
        pose.pdb_info().chain(posi[0]))
    for pos in posi:
        chain_selector = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector(
            chain_selector, pyrosetta.rosetta.core.select.residue_selector.ChainSelector(pose.pdb_info().chain(pos)))

    nbr_selector = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector(
        chain_selector, nbr_selector)

    # Select No Design Area
    not_design = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(
        comb_select)

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

    tf_mut = tf.clone()
    # Enable design (change the residues to the mutated residues)
    for i in range(len(posi)):
        aa_to_design = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
        aa_to_design.aas_to_keep(amino[i])
        tf_mut.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
            aa_to_design,  mut_posi[i]))

    # Minimizes residues within repack range of mutations
    mmf = MoveMapFactory()
    mmf.add_bb_action(mm_enable, nbr_selector)
    mmf.add_chi_action(mm_enable, nbr_selector)
    
    # mm = mmf.create_movemap_from_pose(pose) # ONLY NEEDED IF FAST RELAX

    # BACKRUBBING
    backrub_ensemble_gen(pose, nbr_selector, mmf, tf, scorefxn)
    mutPose = Pose()
    mutPose.detached_copy(pose)

    packer_mut = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(
        scorefxn)
    packer_wt = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(
        scorefxn)
    packer_mut.task_factory(tf_mut)
    packer_wt.task_factory(tf)

    # Want global minimization after backrub, this was turned on for minimization in flex ddg I think FIXME
    if values_dict["a"]:
        mmf.all_bb(True)
        mmf.all_chi(True)
    
    
    minmover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
    minmover.score_function(scorefxn)
    minmover.movemap_factory(mmf)
    minmover.max_iter(2000) # apparently 5000 was used in flex ddg FIXME
    minmover.tolerance(0.00001) # apparently 0.000001 was used in flex ddg FIXME
    minmover.abs_score_convergence_threshold(1.0)

    if values_dict["c"]:
        minmover.cartesian(True)

    packer_wt.apply(pose)
    packer_mut.apply(mutPose)
    minmover.apply(pose)
    minmover.apply(mutPose)

    return pose, mutPose


def unbind(pose, jump):
    STEP_SIZE = 100
    # JUMP NOTED FOR EACH POSE MANUALLY
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

    # Bound unmutated
    original, mutPose = pack_and_relax(mutPose, pos, mut, repack_range, scorefxn)    
    if output_pdb:
        original.dump_pdb("1_bound_unmutated.pdb")
    bound_unmutated = ddg_scorefxn(original)

    # Bound mutated
    if output_pdb:
        mutPose.dump_pdb("2_bound_mutated.pdb")
    bound_mutated = ddg_scorefxn(mutPose)
    rmsd_mutated = all_atom_rmsd(original, mutPose)
    
    unbound_original.detached_copy(original)
    unbound_mutPose.detached_copy(mutPose)

    # Unbound unmutated
    unbind(unbound_original, jump)
    # pack_and_relax(unbound_original, pos, wt, repack_range, scorefxn) # flex-ddg does NOT repack unbound poses
    if output_pdb:
        unbound_original.dump_pdb("3_unbound_unmutated.pdb")
    unbound_unmutated = ddg_scorefxn(unbound_original)

    # Unbound mutated
    unbind(unbound_mutPose, jump)
    # pack_and_relax(unbound_mutPose, pos, mut, repack_range, scorefxn) # flex-ddg does NOT repack unbound poses
    if output_pdb:
        unbound_mutPose.dump_pdb("4_unbound_mutated.pdb")
    unbound_mutated = ddg_scorefxn(unbound_mutPose)

    # print("unbound_unmutated", unbound_unmutated)
    # print("bound_unmutated", bound_unmutated)
    # print("unbound_mutated", unbound_mutated)
    # print("bound_mutated", bound_mutated)
    ddG = (bound_mutated - unbound_mutated) - \
        (bound_unmutated - unbound_unmutated)
    return ddG, rmsd_mutated


if values_dict["n"] == "all":
    pdbs = data["#PDB"].unique() # TESTING FIXME
else:
    pdbs = values_dict["n"] # TESTING FIXME
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

repack_range = int(values_dict["r"]) # 8 default is from flex-ddg

count = 0

for pdb in pdbs:
    points = data.loc[data["#PDB"] == pdb]
    points = points.loc[points["Interface?"] == True]
    points = points.loc[points["LD"] == 1]  # TESTING FIXME
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
        rmsd_total = 0
        for _ in range(int(values_dict["p"])):
            ddg, rmsd = calc_ddg(pose, pos, wt, mut, repack_range, jump, False)
            total += ddg
            rmsd_total += rmsd
        total = total / int(values_dict["p"])
        rmsd_total = rmsd_total / int(values_dict["p"])
        print("DDG: ", total)
        print("RMSD: ", rmsd_total)
        df = pd.concat([df, pd.DataFrame({"#PDB": pdb, "Position": pos, "WT_AA": wt,
                                          "Mut_AA": mut, "DDG": total, "RMSD": rmsd_total})], ignore_index=True, sort=True)
        end = time.time()
        print("Total time:", end-start, "seconds")
        print("Avg time per ensemble member:",
              (end-start)/int(values_dict["p"]), "seconds")
        count += 1
        if count % 2 == 0:
            df.to_csv(values_dict["o"], index=False)
            print("Wrote to csv.", flush=True)

df.to_csv(values_dict["o"], index=False)
