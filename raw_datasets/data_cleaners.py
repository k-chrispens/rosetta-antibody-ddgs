"""
Functions for formatting and analyzing raw_datasets using pandas.

Used in the ab_proj virtual env.
"""

import math
import numpy as np
import pandas as pd
import re


def ab_bind_clean(file_df):
    """Cleaning and filtering for the AB-Bind database. 
    Filter out homology modeled PDBs, as well as unneeded data 
    (e.g. Rfree, pH, temperature)."""
    ab_bind = file_df.copy(True)
    ab_bind = ab_bind.rename(columns={"PDB DOI": "Source"})

    # Remove homology modeled structures and structures with no heavy and light chain label
    homology = r"HM.+"
    notHL = r"[^HL]:"
    abcd = r"AB_CD" # to keep the one entry of two antibodies binding together
    to_keep = ab_bind["Partners(A_B)"].str.contains(abcd)
    keep = ab_bind[to_keep]
    to_remove = ab_bind['#PDB'].str.contains(homology)
    ab_bind = ab_bind[~to_remove]
    to_remove = ab_bind['Mutation'].str.contains(notHL)
    ab_bind = ab_bind[~to_remove]
    ab_bind = pd.concat([ab_bind, keep])

    ab_bind["Mutation"] = ab_bind["Mutation"].apply(
        lambda x: re.sub(r",", ";", x))

    ab_bind["LD"] = ab_bind["Mutation"].apply(lambda x: x.count(";") + 1)

    return ab_bind


def skempi_clean(file_df):
    """Cleaning and filtering for SKEMPI 2.0 database.
    Filter to include only relevant information."""

    skempi = file_df.copy(True)

    skempi["#Pdb"] = skempi["#Pdb"].apply(
        lambda x: re.sub(r"(\w{4}).*", r"\1", x))
    skempi["Mutation(s)_PDB"] = skempi["Mutation(s)_PDB"].apply(
        lambda x: re.sub(r"(\w)(\w)(\w+)", r"\2:\1\3", re.sub(r",", r";", x)))
    prot1 = skempi["Protein 1"].str.contains(
        "fab|mab", case=False) | skempi["Protein 1"].str.contains("antibody|Fv", case=False)
    prot2 = skempi["Protein 2"].str.contains(
        "fab|mab", case=False) | skempi["Protein 2"].str.contains("antibody|Fv", case=False)

    abs = prot1 | prot2
    skempi = skempi[abs]

    skempi = skempi.assign(ddG=lambda x: ddg_from_kd(
        x["Affinity_mut_parsed"], x["Temperature"], x["Affinity_wt_parsed"]))
    skempi.rename(columns={"ddG": "ddG(kcal/mol)"}, inplace=True)

    return skempi


def sipdab_clean(file_df):
    """Cleaning and filtering for SiPDAB database.
    Filter to include only relevant information: e.g. LD, PDB, onehot encoding, mutations, ddG."""


def ddg_from_kd(kd, temp, reference_kd):
    """Returns affinity ddG value from disassociation constant 
    using dG = RTlnK_d and a reference K_d."""
    return 8.31446261815324 / (4184) * temp * (np.log(kd) - np.log(reference_kd))


def phillips_clean(df: pd.DataFrame, cr9114: bool):
    """Filtering and calculations to fit the parsed Phillips et al. data 
    to the remainder of the data from other sources. 
    From PDB, T = 295K for CR6261, 293K for CR9114, but using body temp 310K 
    FIXME: Currently not sure what units are going in for ddG"""
    mut_df = df.copy(True)
    index = mut_df.index[mut_df['Mutations'].apply(
        len) == 16] if cr9114 else mut_df.index[mut_df['Mutations'].apply(len) == 11]

    reference = mut_df.iloc[index]["-logKD"]
    mut_df["-logKD"] = mut_df["-logKD"].apply(
        lambda x: ddg_from_kd(math.pow(10, -x), 310, math.pow(10, -reference)))

    mut_df.rename(columns={"-logKD": "ddG"}, inplace=True)
    return mut_df


def return_mut_df(file_df, cr9114: bool):
    """Function adapted from Brian Petersen to parse the 1hot data 
    from Phillips et al. to a list of mutations. Set cr9114 bool if working with CR9114 instead of CR6261."""
    if cr9114:
        mutations = ['F29S', 'S30N', 'S31N', 'I52S', 'T57S', 'A58T', 'N59A',
                     'T71S', 'K74I', 'S75F', 'T76S', 'S77N', 'S84N', 'R87T', 'Y95F', 'Y106S']
        pdb_id = "4FQY"
    else:
        mutations = ['T28P', 'S30R', 'A58T', 'N59K', 'Q62P',
                     'E74D', 'S75F', 'T76A', 'S77G', 'A79V', 'L104V']
        pdb_id = "3BGN"
    var_id = []
    var_mutations = []
    num_mutations = []

    for i in range(len(file_df['variant'])):
        # To deal with cases where 1hot length is not the same length as the total amount of mutations due to dropping off 0s
        one_hot = list(str(file_df['variant'].at[i]))
        if len(one_hot) != len(mutations):
            to_add = len(mutations) - len(one_hot)
            one_hot = ['0']*to_add + one_hot
        num_mutations.append(one_hot.count('1'))
        var_id.append(one_hot)
        var_mutations.append([mutations[j] for j in range(
            len(mutations)) if var_id[i][j] == '1'])

    log_kd = file_df['logKd']

    new_df = pd.DataFrame({
        "#PDB": pdb_id,
        "1hot": file_df['variant'],
        "Mutations": var_mutations,
        "-logKD": log_kd,
        "LD": num_mutations
    })
    # df of mutations and logKDs
    return new_df


def filter_overlap_and_combine(df1, df2):
    """Removes overlapping entries from two dataframes and 
    returns a combined dataset given that the two are combineable FIXME"""


def mason_etal_clean(filedf: pd.DataFrame):
    """Cleaning and filtering for Mason et al. data.
    Filter to include only relevant information: e.g. LD, PDB, onehot encoding, mutations, ddG.
    NOTE: The one hot encoding here does NOT account for mutations to
    different aas, only the position that was mutated in the original sequence."""

    newdf = filedf.copy(True)

    compare_to = "CSRWGGDGFYAMDYW"
    mut_start_pos = 96
    mutations = []
    one_hots = []
    cases = []

    for seq in newdf["Sequence"]:
        cases.append((compare_to, seq))

    for a, b in cases:
        muts = []
        one_hot = ""
        for i in range(len(a)):
            if a[i] == b[i]:
                one_hot += "0"
                continue
            else:
                one_hot += "1"
                mut = f"{a[i]}{i + mut_start_pos}{b[i]}"
                muts.append(mut)
        one_hots.append(one_hot)
        mutations.append(muts)

    newdf["1hot"] = one_hots
    newdf["Mutations"] = mutations
    newdf = newdf.drop(columns=["Sequence", "KD"])

    return newdf

def kiyoshi_clean(filedf: pd.DataFrame):
    """Cleaning and filtering for Kiyoshi et al. data.
    Filter to include only relevant information: e.g. LD, PDB, onehot encoding, mutations, ddG."""
    kiyoshi_etal = filedf.copy(True)
    kiyoshi_etal.drop(axis=1, index=0, inplace=True)

    kiyoshi_etal["Mutation"] = kiyoshi_etal["Mutation"].apply(
        lambda x: dc.re.sub(r"(\w)-(\w+)", r"\1:\2", x))

    return kiyoshi_etal