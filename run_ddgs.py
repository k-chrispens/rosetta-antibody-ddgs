"""Taken from the flex-ddG tutorial
https://github.com/Kortemme-Lab/flex_ddG_tutorial/blob/master/run_example_2_saturation.py
and modified for ease of use"""

import socket
import sys
import os
import subprocess
import numpy as np
import pandas as pd
import getopt
import re

args = sys.argv[1:]
options = "e:o:t:i:j:p:"
long_options = ["ensemble_size=",
                "output_path=", "trials=", "input=", "jump=", "positions="]
values_dict = {"e": 10, "o": "UNNAMED", "t": 10000, "i": None, "j": 0, "p": []}

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

        elif currentArgument in ("-t", "--trials"):
            values_dict["t"] = currentValue
            print(f"Backrub Trials = {currentValue}")

        elif currentArgument in ("-i", "--input"):
            values_dict["i"] = currentValue
            print(f"Input PDB = {currentValue}")

        elif currentArgument in ("-j", "--jump"):
            values_dict["j"] = currentValue
            print(f"Jump = {currentValue}")

        elif currentArgument in ("-p", "--positions"):
            values_dict["p"] = re.split(",", currentValue)
            print(f"Positions = {currentValue}")
            
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
path_to_script = 'ddG_backrub_og.xml' # Use ddG_backrub.xml for REF2015 energy function. MUST ALSO DELETE "-restore_talaris_behavior"!

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
        # chains, pos, mut
        mutation_chain, pos, mut_aa = mut_info
        f.write('%s %s PIKAA %s\n' %
            (pos, mutation_chain, mut_aa))

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
        '-restore_talaris_behavior', # MUST DELETE FOR REF2015 FUNCTION
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

    path = values_dict["i"]
    jump = values_dict["j"]
    positions = values_dict["p"]
    for position in positions:
        position = re.sub(r"(\w):(\w)(\d+)(\w*)", r"\1:\2:\3:\4", position)

        for nstruct_i in range(1, nstruct + 1):
            residues_to_mutate = []
            for aa in "ACDEFGHIKLMNPQRSTVWY":
                chain, start, pos, ic = re.split(":", position)
                mut = aa
                if start == mut:
                    continue
                else:
                    pos = str(pos) + ic
                    residue_to_mutate = (chain, pos, mut)
                    pdb = re.sub(r"[.\w\/_]*\/(\w{4})[.\w\/_]*.pdb", r"\1", path)
                    print(residues_to_mutate)
                    run_flex_ddg_saturation('{}_{}{}{}'.format(pdb, start, pos, mut), path, jump, residue_to_mutate, nstruct_i)
