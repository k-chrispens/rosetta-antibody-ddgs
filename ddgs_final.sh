#!/bin/bash

#SBATCH --job-name=rscripts-ddg
#SBATCH --output=rscripts_ddgcalc.%j.out
#SBATCH --time=24:00:00
#SBATCH --partition=amilan-ucb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=END
#SBATCH --mail-user=karson.chrispens@colorado.edu

source ~/.bashrc
which python
# RESET DEFAULT -t and -e to 10000 and 10
python run_ddgs.py -p $1 -o $2 -j $3 -i $4 -t 1000 -e 1
