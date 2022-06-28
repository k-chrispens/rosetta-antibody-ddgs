"""
Script for formatting and analyzing raw_datasets using pandas.
Author: Karson Chrispens
Used in the ab_proj virtual env. Requires functions from data_cleaners.py
"""

import data_cleaners as dc

# INITIAL IMPORT

ab_bind = dc.pd.read_excel("~/rosetta-antibody-ddgs/raw_datasets/AB_bind_minor.xlsx",
                           sheet_name="AB-Bind_experimental_data.csv",
                           usecols=[0, 1, 2, 3, 4, 5, 11, 12, 13])
# Using subsetted SiPDAB data that can be read into pandas more easily.
sipmab = dc.pd.read_excel(
    "~/rosetta-antibody-ddgs/raw_datasets/SiPMAB_subsetted.xlsx", sheet_name="Sheet1")
# Pulling cleaned SKEMPI 2.0 data that has values for all affinities
# (no > and < values were included) and contains antibodies only,
# as well as mutations in only the antibodies. This was obtained by
# scraping the SKEMPI 2.0 Website (skempiscrape.py), as the downloaded
# csv did not have all the data for unknown reasons.
skempi = dc.pd.read_csv(
    "~/rosetta-antibody-ddgs/raw_datasets/skempi_scraped.csv")
mason_etal = dc.pd.read_csv(
    "~/rosetta-antibody-ddgs/raw_datasets/mason_etal_data.csv")
kiyoshi_etal = dc.pd.read_csv(
    "~/rosetta-antibody-ddgs/raw_datasets/Kiyoshi_etal_data.csv")
# For missing logKd values, I averaged the logKd_x and logKd_y values.
phillips_etal_cr6261 = dc.pd.read_csv(
    "~/rosetta-antibody-ddgs/raw_datasets/phillips_etal_CR6261_H1.csv")
phillips_etal_cr9114 = dc.pd.read_csv(
    "~/rosetta-antibody-ddgs/raw_datasets/phillips_etal_CR9114_H3_filt.csv")

# CLEANING

ab_bind = dc.ab_bind_clean(ab_bind)
sipmab = dc.sipmab_clean(sipmab)
skempi = dc.skempi_rename(skempi)
mason_etal = dc.mason_etal_clean(mason_etal) # This overlaps, will need to do simple concat.
kiyoshi_etal = dc.kiyoshi_clean(kiyoshi_etal) # Already in SKEMPI 2.0, unfortunately.
phillips_etal_cr6261 = dc.return_mut_df(phillips_etal_cr6261, False)
phillips_etal_cr6261 = dc.phillips_clean(phillips_etal_cr6261, False)
phillips_etal_cr9114 = dc.return_mut_df(phillips_etal_cr9114, True)
phillips_etal_cr9114 = dc.phillips_clean(phillips_etal_cr9114, True)

ab_bind_skempi = dc.filter_overlap_and_combine(ab_bind, skempi)
ab_bind_skempi_sipmab = dc.filter_overlap_and_combine(ab_bind_skempi, sipmab)
ab_bind_skempi_sipmab_mason = dc.pd.concat([ab_bind_skempi_sipmab, mason_etal], join='inner')
ab_bind_skempi_sipmab_mason_phillips = dc.filter_overlap_and_combine(
    ab_bind_skempi_sipmab_mason, phillips_etal_cr6261)
all = dc.filter_overlap_and_combine(
    ab_bind_skempi_sipmab_mason_phillips, phillips_etal_cr9114)

try:
    all.to_csv('./raw_datasets/full_data.csv', index=False)
    print("Wrote file.")
except:
    print("Did not output.")