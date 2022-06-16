"""
Script for formatting and analyzing raw_datasets using pandas.

Used in the ab_proj virtual env. Requires functions from data_cleaners.py
"""

import data_cleaners as dc

ab_bind = dc.pd.read_excel("~/rosetta-antibody-ddgs/raw_datasets/AB_bind_raw.xlsx", sheet_name="AB-Bind_experimental_data.csv", usecols=[0, 1, 2, 3, 4, 5, 11, 12, 13])
# Using subsetted SiPDAB data that can be read into pandas more easily.
sipdab = dc.pd.read_excel("~/rosetta-antibody-ddgs/raw_datasets/SiPDAB_subsetted.xlsx", sheet_name="Sheet1")
# Using cleaned SKEMPI 2.0 data that has values for all affinities (though > and < values were also included) and contains antibodies only.
skempi = dc.pd.read_csv("~/rosetta_antibody_ddgs/raw_datasets/skempi_v2_cleaned.csv")
mason_etal = dc.pd.read_csv("~/rosetta_antibody_ddgs/raw_datasets/mason_etal_data.csv")
kiyoshi_etal = dc.pd.read_csv("~/rosetta_antibody_ddgs/raw_datasets/Kiyoshi_etal_data.csv")

# 