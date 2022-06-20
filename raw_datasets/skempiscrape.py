import pandas as pd
import re
import data_cleaners as dc


def skempi_clean_test(file_df):
    """Cleaning and filtering for SKEMPI 2.0 database.
    Filter to include only relevant information. NOTE: Using 310K as temperature for all ddG calculations."""

    skempi = file_df.copy(True)

    undef = r"[><~nb]"  # filter out undefined entries
    undef_mut_affinity = skempi["Affinity Mut"].astype(str).str.contains(undef)
    undef_wt_affinity = skempi["Affinity Wt"].astype(str).str.contains(undef)
    undef_affinity = undef_mut_affinity | undef_wt_affinity

    skempi = skempi[~undef_affinity]

    skempi["PDB Code"] = skempi["PDB Code"].apply(
        lambda x: re.sub(r"(\w{4}).*", r"\1", x))
    skempi["Mutation Pdb"] = skempi["Mutation Pdb"].apply(
        lambda x: re.sub(r"(\w)(\w)(\w+)", r"\2:\1\3", re.sub(r",", r";", x)))

    skempi["LD"] = skempi["Mutation Pdb"].apply(
        lambda x: x.count(";") + 1)

    skempi["Affinity Mut"] = skempi["Affinity Mut"].astype(float)
    skempi["Affinity Wt"] = skempi["Affinity Wt"].astype(float)

    skempi = skempi.assign(ddG=lambda x: dc.ddg_from_kd(
        x["Affinity Mut"], 310, x["Affinity Wt"]))
    skempi.rename(columns={"ddG": "ddG(kcal/mol)",
                  "PDB Code": "#PDB"}, inplace=True)

    return skempi


dfs = []
for page in range(1, 22):
    url = fr"https://life.bsc.es/pid/skempi2/database/browse/mutations?keywords=mutations.protein_1+contains+%22fab%22+or+mutations.protein_1+contains+%22fv%22+or+mutations.protein_2+contains+%22fv%22+or+mutations.protein_1+contains+%22mab%22+or+mutations.protein_2+contains+%22mab%22+or+mutations.protein_1+contains+%22antibody%22+or+mutations.protein_2+contains+%22antibody%22+or+mutations.protein_2+contains+%22fab%22&page={page}&pagination=50&pdbdb=rcsb"
    dfs.append(pd.read_html(url)[0])

print("Done appending.")

concatdfs = []

for df in dfs:
    concatdfs.append(skempi_clean_test(df))

results = pd.concat(concatdfs)

print("Done concatenating.")

try:
    results.to_csv('skempi_scraped.csv')
    print("Wrote file.")
except:
    print("Did not output.")
