"""
Script for formatting and analyzing raw_datasets using pandas.

Used in the ab_proj virtual env. Requires functions from data_cleaners.py
"""

import data_cleaners as dc

# INITIAL IMPORT

ab_bind = dc.pd.read_excel("~/rosetta-antibody-ddgs/raw_datasets/AB_bind_minor.xlsx", sheet_name="AB-Bind_experimental_data.csv", usecols=[0, 1, 2, 3, 4, 5, 11, 12, 13])
# Using subsetted SiPDAB data that can be read into pandas more easily.
sipmab = dc.pd.read_excel("~/rosetta-antibody-ddgs/raw_datasets/SiPMAB_subsetted.xlsx", sheet_name="Sheet1")
# Using cleaned SKEMPI 2.0 data that has values for all affinities (though > and < values were also included) and contains antibodies only.
skempi = dc.pd.read_csv("~/rosetta-antibody-ddgs/raw_datasets/skempi_v2_cleaned.csv")
mason_etal = dc.pd.read_csv("~/rosetta-antibody-ddgs/raw_datasets/mason_etal_data.csv")
kiyoshi_etal = dc.pd.read_csv("~/rosetta-antibody-ddgs/raw_datasets/Kiyoshi_etal_data.csv")
# For missing logKd values, I averaged the logKd_x and logKd_y values.
phillips_etal_cr6261 = dc.pd.read_csv("~/rosetta-antibody-ddgs/raw_datasets/phillips_etal_CR6261_H1.csv")
phillips_etal_cr9114 = dc.pd.read_csv("~/rosetta-antibody-ddgs/raw_datasets/phillips_etal_CR9114_H3.csv")

# CLEANING

ab_bind = dc.ab_bind_clean(ab_bind)
sipmab = dc.sipmab_clean(sipmab)
skempi = dc.skempi_clean(skempi)
mason_etal = dc.mason_etal_clean(mason_etal)
kiyoshi_etal = dc.kiyoshi_clean(kiyoshi_etal)
phillips_etal_cr6261 = dc.phillips_clean(dc.return_mut_df(phillips_etal_cr6261, False), False)
phillips_etal_cr9114 = dc.phillips_clean(dc.return_mut_df(phillips_etal_cr9114, True), True)