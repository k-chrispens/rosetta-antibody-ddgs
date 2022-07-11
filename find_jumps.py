"""Find jumps in fold trees for accurate modelling of unbinding
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
init()

data = pd.read_csv("./raw_datasets/use_this_data.csv")
pdbs = data["#PDB"].unique()

for pdb in pdbs:
    pose = pose_from_pdb("./PDBs/{pdb}_all.pdb")