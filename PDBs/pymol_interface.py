"""Script to create pymol function to quickly get interface between supplied chains
Author: Karson Chrispens"""

from pymol import cmd

one_letter = {'VAL': 'V', 'ILE': 'I', 'LEU': 'L', 'GLU': 'E', 'GLN': 'Q',
              'ASP': 'D', 'ASN': 'N', 'HIS': 'H', 'TRP': 'W', 'PHE': 'F', 'TYR': 'Y',
              'ARG': 'R', 'LYS': 'K', 'SER': 'S', 'THR': 'T', 'MET': 'M', 'ALA': 'A',
              'GLY': 'G', 'PRO': 'P', 'CYS': 'C'}


def list_res_at_interface(selection1, selection2):
    cmd.do("ls = []")
    # iterate both and name CA, ls.append(f"{one_letter[resn]}{resi}")
    cmd.select(f"({selection1} within 12.0 of {selection2}) and name CA")
    cmd.iterate("sele", r"ls.append(f'{one_letter[resn]}{resi}')")

    cmd.do("print ls")

cmd.extend("interface_res", list_res_at_interface)
