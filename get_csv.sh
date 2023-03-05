#!/bin/sh

#SBATCH --job-name=get_csv
#SBATCH --output=get_csv.%j.out
#SBATCH --time=24:00:00
#SBATCH --partition=amilan
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=END
#SBATCH --mail-user=karson.chrispens@colorado.edu

source ~/.bashrc
which python
python analyze_ddGs.py -o $1
