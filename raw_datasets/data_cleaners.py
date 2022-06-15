"""
Functions for formatting and analyzing raw_datasets using pandas.

Used in the ab_proj virtual env.
"""

import numpy as np
import pandas as pd

def ab_bind_clean(file):
    """Cleaning and filtering for the AB-Bind database. 
    Filter out homology modeled PDBs, as well as unneeded data 
    (e.g. Rfree, pH, temperature)."""

def skempi_clean(file):
    """Cleaning and filtering for SKEMPI 2.0 database.
    Filter to include only relevant information: e.g. FIXME."""

def sipdab_clean(file):
    """Cleaning and filtering for SiPDAB database.
    Filter to include only relevant information: e.g. FIXME."""

def desai_clean(file):
    """Cleaning and filtering for Phillips et al. database.
    Filter to include only relevant information: e.g. FIXME"""