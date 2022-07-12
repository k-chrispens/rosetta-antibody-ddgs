"""
Functions for formatting and analyzing raw_datasets using pandas.
Author: Karson Chrispens
Used in the ab_proj virtual env.
"""

import math
from tokenize import String
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
    abcd = r"AB_CD"  # to keep the one entry of two antibodies binding together
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
    ab_bind.rename(columns={"Mutation": "Mutations"}, inplace=True)
    ab_bind["Source"] = "AB-Bind"  # Naming source for final dataset

    return ab_bind


def skempi_rename(file_df: pd.DataFrame):
    """Renaming a few things to be consistent with other dataframes"""

    skempi = file_df.copy(True)

    skempi.rename(columns={"Mutation Pdb": "Mutations"}, inplace=True)
    skempi["LD"] = skempi["Mutations"].apply(
        lambda x: x.count(";") + 1)
    skempi["Source"] = "SKEMPI 2.0"  # Naming source for final dataset

    return skempi


def sipmab_clean(file_df):
    """Cleaning and filtering for SiPMAB database.
    Filter to include only relevant information: e.g. LD, PDB, onehot encoding, mutations, ddG."""
    sipmab = file_df.copy(True)
    sipmab["Mutation"] = sipmab["Mutation"].apply(
        lambda x: re.sub(r"(\w)(\w+)", r"\1:\2", convert_3to1(x)))
    sipmab["LD"] = 1
    sipmab.rename(columns={"Mutation": "Mutations"}, inplace=True)
    sipmab["Source"] = "SiPMAB"  # Naming source for final dataset

    return sipmab


def ddg_from_kd(kd, temp, reference_kd):
    """Returns affinity ddG value from disassociation constant 
    using dG = RTlnK_d and a reference K_d."""
    return 8.31446261815324 / (4184) * temp * (np.log(kd) - np.log(reference_kd))


def phillips_clean(df: pd.DataFrame, cr9114: bool):
    """Filtering and calculations to fit the parsed Phillips et al. data 
    to the remainder of the data from other sources. 
    From PDB, T = 295K for CR6261, 293K for CR9114, but using body temp 310K."""
    mut_df = df.copy(True)
    mut_df = mut_df[mut_df['-logKD'] != 5]
    mut_df = mut_df[~mut_df['-logKD'].isnull()]
    index = mut_df.index[mut_df["LD"] == 0] # Switched mutations from germline to mutant to mutant to germline due to pdb
    reference = mut_df["-logKD"].iloc[index]
    mut_df["-logKD"] = mut_df["-logKD"].apply(
        lambda x: ddg_from_kd(math.pow(10, -x), 310, math.pow(10, -reference)))

    mut_df.rename(columns={"-logKD": "ddG(kcal/mol)"}, inplace=True)
    mut_df["Source"] = "Phillips et al. 2021"
    return mut_df


def return_mut_df(file_df, cr9114: bool):
    """Function adapted from Brian Petersen to parse the 1hot data 
    from Phillips et al. to a list of mutations. Set cr9114 bool if working with CR9114 instead of CR6261."""
    # Mutations are from PDB to germline
    if cr9114:
        mutations = ['S29F', 'N30S', 'N31S', 'S52I', 'S56T', 'T57A', 'A58N',
                     'S70T', 'I73K', 'F74S', 'S75T', 'N76S', 'N82AS', 'T83R', 'F91Y', 'S100BY']
        pdb_id = "4FQY"
    else:
        mutations = ['P28T', 'R30S', 'T57A', 'K58N', 'P61Q',
                     'D73E', 'F74S', 'A75T', 'G76S', 'V78A', 'V100L']
        pdb_id = "3GBN"
    var_mutations = []
    num_mutations = []

    for i in range(len(file_df['variant'])):
        muts = ""
        # To deal with cases where 1hot length is not the same length as the total amount of mutations due to dropping off 0s
        one_hot = list(str(file_df['variant'].at[i]))
        if len(one_hot) != len(mutations):
            to_add = len(mutations) - len(one_hot)
            one_hot = ['0']*to_add + one_hot
        for j in range(len(one_hot)):
            # Germline in dataset is all 0, we want that to be the full mutant
            if one_hot[j] == '0':
                muts += f"H:{mutations[j]};"
        muts = muts[:-1]
        if len(muts) == 0:
            num_mutations.append(0)
        else:
            num_mutations.append(muts.count(';') + 1)
        var_mutations.append(muts)

    log_kd = file_df['logKd']

    new_df = pd.DataFrame({
        "#PDB": pdb_id,
        "Mutations": var_mutations,
        "-logKD": log_kd,
        "LD": num_mutations
    })
    # df of mutations and logKDs
    return new_df


def filter_overlap_and_combine(df1: pd.DataFrame, df2: pd.DataFrame):
    """Removes overlapping entries from two dataframes and 
    returns a combined dataframe given that the two are combineable.
    NOTE: Prefers entries in df1 over df2 when overlap occurs. 
    The column with PDB IDs must be called \"#PDB\"."""
    df1_unique = df1["#PDB"].unique()
    df2_unique = df2["#PDB"].unique()

    new_df = df1

    for i in df2_unique:
        if i in df1_unique:
            continue
        else:
            new_df = pd.concat([new_df, df2.loc[df2["#PDB"] == i]], join="inner")
    return new_df


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
        muts = ""
        one_hot = ""
        lds = []
        ld = 0
        for i in range(len(a)):
            if a[i] == b[i]:
                one_hot += "0"
                continue
            else:
                ld += 1
                one_hot += "1"
                mut = f"H:{a[i]}{i + mut_start_pos}{b[i]};"
                muts = muts + mut
        muts = muts[:-1]
        lds.append(ld)
        one_hots.append(one_hot)
        mutations.append(muts)

    newdf["1hot"] = one_hots
    newdf["Mutations"] = mutations
    newdf = newdf.drop(columns=["Sequence", "KD"])
    newdf["Source"] = "Mason et al. 2021"

    newdf = newdf[newdf['1hot'] != "000000000000000"]
    return newdf


def kiyoshi_clean(filedf: pd.DataFrame):
    """Cleaning and filtering for Kiyoshi et al. data.
    Filter to include only relevant information: e.g. LD, PDB, onehot encoding, mutations, ddG."""
    kiyoshi_etal = filedf.copy(True)
    kiyoshi_etal.drop(axis=1, index=0, inplace=True)

    kiyoshi_etal["Mutation"] = kiyoshi_etal["Mutation"].apply(
        lambda x: re.sub(r"(\w)-(\w+)", r"\1:\2", x))
    kiyoshi_etal.rename(columns={"Mutation": "Mutations"}, inplace = True)
    kiyoshi_etal["Source"] = "Kiyoshi et al. 2014"

    return kiyoshi_etal


def convert_3to1(mutation: str):
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
         'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    for three in d.keys():
        if three in mutation:
            mutation = mutation.replace(three, d[three])

    return mutation
