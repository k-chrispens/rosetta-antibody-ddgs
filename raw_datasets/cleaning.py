"""
Script for formatting and analyzing raw_datasets using pandas.

Used in the ab_proj virtual env. Requires functions from data_cleaners.py
"""

import data_cleaners as dc

ab_bind = dc.pd.read_excel("AB_bind_raw.xlsx", sheet_name="AB-Bind_experimental_data.csv", usecols=[0, 1, 2, 3, 4, 5, 11, 12])

print(ab_bind)