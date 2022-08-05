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
python run_ddgs.py -n $1 -o $2 -t 10000 -e 10
