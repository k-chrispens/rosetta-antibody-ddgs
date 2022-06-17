"""Test File"""

import data_cleaners as dc

phillips_etal = dc.pd.read_csv("~/rosetta-antibody-ddgs/raw_datasets/phillips_etal_CR6261_H1.csv")

testdf = dc.return_mut_df(phillips_etal, False)
# print(testdf)

# print(testdf.index[testdf["Mutations"].apply(len) == 11])

# print(dc.phillips_clean(testdf, False))
cleaneddf = dc.phillips_clean(testdf, False)
print(min(cleaneddf["ddG"]))
cleaneddf = cleaneddf.sort_values(by="ddG")
print(cleaneddf)

print(cleaneddf.describe())

