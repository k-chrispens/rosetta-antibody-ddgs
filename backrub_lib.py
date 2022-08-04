"""Taken from the flex-ddG tutorial
https://github.com/Kortemme-Lab/flex_ddG_tutorial/blob/master/run_example_2_saturation.py
and modified for only producing the backrub ensemble."""

#!/usr/bin/python

import socket
import sys
import os
import subprocess
import numpy as np
import pandas as pd
import getopt
import re

# Maybe turn this on later if there's some pdbs that take especially long — also could delete some of the far residues.
use_multiprocessing = False
if use_multiprocessing:
    import multiprocessing
    max_cpus = 2  # We might want to not run on the full number of cores, as Rosetta take about 2 Gb of memory per instance

args = sys.argv[1:]
options = "e:o:n:t:"
long_options = ["ensemble_size=",
                "output_path=", "pdbs=", "trials="]
values_dict = {"e": 1,
               "o": "UNNAMED", "n": "all", "t": 5000}

try:
    # Parsing argument
    arguments, values = getopt.getopt(args, options, long_options)

    # checking each argument
    for currentArgument, currentValue in arguments:
        
        if currentArgument in ("-e", "--ensemble_size"):
            values_dict["p"] = currentValue
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

nstructs = 50
rosetta_scripts_path = os.path.abspath(
    "/projects/kach6913/rosetta.source.release-314/main/source/bin/rosetta_scripts.linuxgccrelease")
max_minimization_iter = 5000
abs_score_convergence_thresh = 1.0
path_to_script = 'produce_backrub_lib.xml'
# Getting data from the dataset
data = pd.read_csv("./raw_datasets/use_this_data.csv")
points = data.loc[data["Interface?"] == True]  # TESTING

if not os.path.isfile(rosetta_scripts_path):
    print('ERROR: "rosetta_scripts_path" variable must be set to the location of the "rosetta_scripts" binary executable')
    print('This file might look something like: "rosetta_scripts.linuxgccrelease"')
    raise Exception('Rosetta scripts missing')


def make_ensemble_member(name, input_pdb_path, jump, mut_info):
    output_directory = os.path.join(values_dict['o'], os.path.join(
        '%s' % (name)))
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
        'max_minimization_iter=%d' % max_minimization_iter,
        'abs_score_convergence_thresh=%.1f' % abs_score_convergence_thresh,
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
        path = f"./inputs/{pdb}_all.pdb"
        points_pdb = points.loc[points["#PDB"] == pdb]  # TESTING

        for index, point in points_pdb.iterrows():
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
                pool.apply_async(make_ensemble_member, args=args)
            else:
                print(residues_to_mutate)
                make_ensemble_member('{}_{}'.format(
                    pdb, name_muts), path, jump, residues_to_mutate)

            if use_multiprocessing:
                pool.close()
                pool.join()
