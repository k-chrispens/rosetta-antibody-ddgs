import pandas as pd
import re
import data_cleaners as dc


def skempi_clean_test(file_df):
    """Cleaning and filtering for SKEMPI 2.0 database.
    Filter to include only relevant information. NOTE: Using 310K as temperature for all ddG calculations."""

    skempi = file_df.copy(True)
    # prot1 = skempi["Protein 1"].str.contains(
    #     "fab|mab", case=False) | skempi["Protein 1"].str.contains("antibody|Fv", case=False)
    # prot2 = skempi["Protein 2"].str.contains(
    #     "fab|mab", case=False) | skempi["Protein 2"].str.contains("antibody|Fv", case=False)
    undef = r"[><~nb]"
    undef_mut_affinity = skempi["Affinity Mut"].str.contains(undef)
    undef_wt_affinity = skempi["Affinity Wt"].str.contains(undef)
    undef_affinity = undef_mut_affinity | undef_wt_affinity

    skempi = skempi[~undef_affinity]

    skempi["PDB Code"] = skempi["PDB Code"].apply(
        lambda x: re.sub(r"(\w{4}).*", r"\1", x))
    skempi["Mutation Pdb"] = skempi["Mutation Pdb"].apply(
        lambda x: re.sub(r"(\w)(\w)(\w+)", r"\2:\1\3", re.sub(r",", r";", x)))

    skempi["LD"] = skempi["Mutation(s)_PDB"].apply(
        lambda x: x.count(";") + 1)
    # abs = prot1 | prot2
    # skempi = skempi[abs]

    skempi = skempi.assign(ddG=lambda x: dc.ddg_from_kd(
        x["Affinity Mut"], 310, x["Affinity Wt"]))
    skempi.rename(columns={"ddG": "ddG(kcal/mol)",
                  "PDB Code": "#PDB"}, inplace=True)

    return skempi


dfs = []
for page in range(1, 22):
    url = fr"https://life.bsc.es/pid/skempi2/database/browse/mutations?keywords=mutations.protein_1+contains+%22fab%22+or+mutations.protein_1+contains+%22fv%22+or+mutations.protein_2+contains+%22fv%22+or+mutations.protein_1+contains+%22mab%22+or+mutations.protein_2+contains+%22mab%22+or+mutations.protein_1+contains+%22antibody%22+or+mutations.protein_2+contains+%22antibody%22+or+mutations.protein_2+contains+%22fab%22&page={page}&pagination=50&pdbdb=rcsb"
    dfs.append(pd.read_html(url)[0])

concatdfs = [ skempi_clean_test(df) for df in dfs ]
results = pd.concat(concatdfs)

results.to_csv('skempi_scraped.csv')