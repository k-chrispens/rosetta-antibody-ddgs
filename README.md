# Rosetta Antibody ΔΔGs

Code and datasets for an antibody binding ΔΔG protocol designed in Rosetta to improve simulation of antibody evolution. The protocol is that of [Barlow et al. 2018](https://pubs.acs.org/doi/10.1021/acs.jpcb.7b11367) ([GitHub](https://github.com/Kortemme-Lab/flex_ddG_tutorial)), with only a few parameter changes. 

## How to use

To obtain ddGs for a saturation mutagenesis type procedure, there are three main files to pay attention to: `predict_ddg.sh`, `ddgs_final.sh`, and `get_csv.sh`. 

For CU Boulder users, these are all set up for auto-submission of a job to Alpine on the amilan-ucb partition. A job will be submitted for each mutation position you enter. Please also remember to enter your email in the `--mail-user=` tag.

### Setup

Make sure you have installed the Anaconda environment used for this project contained in env.yaml.

```bash
conda env create --file=env.yaml
```

Then ensure that you have either cloned the whole repository or have downloaded the python and/or XML files necessary for the shell scripts to run — `ddG_backrub_og.xml`, `ddG_backrub.xml`, `calc_ddgs.py`, `analyze_ddGs.py`, `run_ddgs.py`, and the shell scripts themselves.

`harmonic_prerelax.py` is also helpful if your structure files have not been relaxed into Rosetta yet. To do this, run:

```bash
python harmonic_prerelax.py PATH_TO_PDB
```

Now, to run the PyRosetta version of the script, use `predict_ddg_py.sh`. To run the RosettaScripts version, use `predict_ddg.sh`. I have listed the process for `predict_ddg_py.sh` here. Both take the same command line arguments, but the RosettaScripts version requires that you set the path to your RosettaScripts executable in `run_ddgs.py`. The PyRosetta version has the REF2015 score function as default, while the RosettaScripts version defaults to the Talaris score function. The PyRosetta version does not yet output PDBs like the RosettaScripts version.

Running the script:

```bash
./predict_ddg_py.sh -p POSITION -o OUTPUT_PATH -i INPUT_PDB_PATH -j JUMP_NUMBER -n NUM_TASKS
```

Position arguments must be listed in the form (CHAIN):(WT_AA)(PDB_INDEX), e.g. H:D32 or L:S100B with an insertion code. Multiple positions can be entered at once, and should be entered as a bash array with spaces between positions, e.g. (H:D32 L:S100B ...)

The jump number will be used to properly unbind the antibody from the antigen, and can be found using the `find_jumps.ipynb` notebook.

The NUM_TASKS argument is the number of parallel jobs to run. However, the job set up by the `predict_ddg.sh` script will only run one job at a time, so this argument must be 1. This may be changed later.

An example script that can be run using files in this repository is:

```bash
./predict_ddg_py.sh -p H:D32 -o ./example -i ./inputs/1DQJ_all.clean.pdb -j 1 -n 1
```

## Questions and Issues
Please open an issue on GitHub if you have any questions or issues.
