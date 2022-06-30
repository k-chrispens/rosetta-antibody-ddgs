"""Test File"""

import data_cleaners as dc

# devo
data = dc.pd.read_csv("~/rosetta-antibody-ddgs/raw_datasets/interface_data_use.csv")
rando = data.sample()
muts = dc.re.split(";", rando["Mutations"].values[0])
print(muts)
all = list(map(lambda x: dc.re.sub(r"(\w):(\w)(\d+)(\w)", r"\1:\2:\3:\4", x), muts))
print(all)
chain, wt, pos, mut = dc.re.split(":", all[0])
print(chain, wt, pos, mut)