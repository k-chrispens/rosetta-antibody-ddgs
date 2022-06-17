"""Test File"""

import data_cleaners as dc

ab_bind = dc.pd.read_excel("~/rosetta-antibody-ddgs/raw_datasets/AB_bind_raw.xlsx",
                           sheet_name="AB-Bind_experimental_data.csv", usecols=[0, 1, 2, 3, 4, 5, 11, 12, 13])

print(ab_bind.head())
