"""Test File"""

import data_cleaners as dc

# devo
phillips_etal_cr6261 = dc.pd.read_csv(
    "~/rosetta-antibody-ddgs/raw_datasets/phillips_etal_CR6261_H1.csv")
print(phillips_etal_cr6261["variant"].is_unique)
phillips_etal_cr6261 = dc.return_mut_df(phillips_etal_cr6261, False)
print(phillips_etal_cr6261["Mutations"].is_unique)
# phillips_etal_cr6261 = dc.phillips_clean(phillips_etal_cr6261, False)
