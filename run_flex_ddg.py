"""Taken from the flex-ddG tutorial
https://github.com/Kortemme-Lab/flex_ddG_tutorial/blob/master/run_example_2_saturation.py
and modified to fit the dataset I have"""

#!/usr/bin/python

import socket
import sys
import os
import subprocess
import numpy as np
import pandas as pd
import getopt
import re

use_multiprocessing = False # Maybe turn this on later if there's some pdbs that take especially long — also could delete some of the far residues.
if use_multiprocessing:
    import multiprocessing
    max_cpus = 4  # We might want to not run on the full number of cores, as Rosetta take about 2 Gb of memory per instance

args = sys.argv[1:]
options = "e:o:n:t:"
long_options = ["ensemble_size=",
                "output_path=", "pdbs=", "trials="]
values_dict = {"e": 1, "o": "UNNAMED", "n": "all", "t": 5000}

try:
    # Parsing argument
    arguments, values = getopt.getopt(args, options, long_options)

    # checking each argument
    for currentArgument, currentValue in arguments:
        
        if currentArgument in ("-e", "--ensemble_size"):
            values_dict["e"] = currentValue
            print(f"Backrub Ensemble = {currentValue}")

        elif currentArgument in ("-o", "--output_path"):
            values_dict["o"] = currentValue
            print(f"Output path = {currentValue}")

        elif currentArgument in ("-n", "--pdbs"):
            values_dict["n"] = re.split(",", currentValue)
            print("PDBs:", currentValue)

        elif currentArgument in ("-t", "--trials"):
            values_dict["t"] = currentValue
            print(f"Backrub Trials = {currentValue}")

except getopt.error as err:
    # output error, and return with an error code
    print(str(err))

###################################################################################################################################################################
# Important: The variables below are set to values that will make the run complete faster (as a tutorial example), but will not give scientifically valid results.
#            Please change them to the "normal" default values before a real run.
###################################################################################################################################################################

rosetta_scripts_path = os.path.abspath(
    "/projects/kach6913/rosetta.source.release-314/main/source/bin/rosetta_scripts.linuxgccrelease")
nstruct = int(values_dict['e'])  # Normally 35
max_minimization_iter = 5000  # Normally 5000
abs_score_convergence_thresh = 1.0  # Normally 1.0
number_backrub_trials = int(values_dict['t'])  # Normally 35000
# Can be whatever you want, if you would like to see results from earlier time points in the backrub trajectory. 7000 is a reasonable number, to give you three checkpoints for a 35000 step run, but you could also set it to 35000 for quickest run time (as the final minimization and packing steps will only need to be run one time).
backrub_trajectory_stride = int(values_dict['t'])
path_to_script = 'ddG_backrub_og.xml' # Use ddG_backrub.xml for REF2015 energy function.
# Getting data from the dataset
data = pd.read_csv("./raw_datasets/use_this_data.csv")
points = data.loc[data["Interface?"] == True]
points = points.loc[points["LD"] == 1] # Comment out to do all mutants

if not os.path.isfile(rosetta_scripts_path):
    print('ERROR: "rosetta_scripts_path" variable must be set to the location of the "rosetta_scripts" binary executable')
    print('This file might look something like: "rosetta_scripts.linuxgccrelease"')
    raise Exception('Rosetta scripts missing')


def run_flex_ddg_saturation(name, input_pdb_path, jump, mut_info, nstruct_i):
    output_directory = os.path.join(values_dict['o'], os.path.join(
        '%s' % (name), '%02d' % nstruct_i))
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)

    resfile_path = os.path.join(output_directory, 'mutate_%s.resfile' % (name))
    with open(resfile_path, 'w') as f:
        f.write('NATRO\nstart\n')
        for mut in mut_info:
            # chains, pos, mut
            mutation_chain, pos, mut_aa = mut
            f.write('%s %s PIKAA %s\n' %
                (pos[0], mutation_chain[0], mut_aa[0]))

    flex_ddg_args = [
        os.path.abspath(rosetta_scripts_path),
        "-s %s" % os.path.abspath(input_pdb_path),
        '-parser:protocol', os.path.abspath(path_to_script),
        '-parser:script_vars',
        'jump=' + str(jump),
        'mutate_resfile_relpath=' + os.path.abspath(resfile_path),
        'number_backrub_trials=%d' % number_backrub_trials,
        'max_minimization_iter=%d' % max_minimization_iter,
        'abs_score_convergence_thresh=%.1f' % abs_score_convergence_thresh,
        'backrub_trajectory_stride=%d' % backrub_trajectory_stride,
        '-restore_talaris_behavior',
        '-in:file:fullatom',
        '-ignore_unrecognized_res',
        '-ignore_zero_occupancy false',
        '-ex1',
        '-ex2',
        '-soft_rep_design'
    ]

    log_path = os.path.join(output_directory, 'rosetta.out')

    print('Running Rosetta with args:')
    print(' '.join(flex_ddg_args))
    print('Output logged to:', os.path.abspath(log_path))
    print()

    outfile = open(log_path, 'w')
    process = subprocess.Popen(flex_ddg_args, stdout=outfile,
                               stderr=subprocess.STDOUT, close_fds=True, cwd=output_directory)
    returncode = process.wait()
    outfile.close()


if __name__ == '__main__':

    for pdb in values_dict["n"]:
        path = f"./inputs/{pdb}_all.clean.pdb"
        points_pdb = points.loc[points["#PDB"] == pdb]
        for index, point in points_pdb.iterrows():
            for nstruct_i in range(1, nstruct + 1):
                residues_to_mutate = []
                name_muts = re.sub(";", "_", point["Mutations"])
                muts = re.split(";", point["Mutations"])
                jump = point["Jump"]
                pos = []
                chains = []
                mut = []
                all = list(map(lambda x: re.sub(
                    r"(\w):(\w)(\d+)(\w*)(\w)", r"\1:\2:\3:\5:\4", x), muts))

                for i in all:
                    chain, start, posi, muta, ic = re.split(":", i)
                    chains.append(chain)
                    pos.append(str(posi) + ic)
                    mut.append(muta)
                residues_to_mutate.append((chains, pos, mut))

                
                if use_multiprocessing:
                    pool = multiprocessing.Pool(processes=min(
                        max_cpus, multiprocessing.cpu_count()))

                if use_multiprocessing:
                    pool.apply_async(run_flex_ddg_saturation, args=args)
                else:
                    print(residues_to_mutate)
                    run_flex_ddg_saturation('{}_{}'.format(pdb, name_muts), path, jump, residues_to_mutate, nstruct_i)

                if use_multiprocessing:
                    pool.close()
                    pool.join()
