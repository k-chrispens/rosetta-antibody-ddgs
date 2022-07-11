"""
Will highlight mutations in the stem region of hemagglutinin for the given sequences. NOTE: MUST RENUMBER USING renumber.py FIRST
Author: Karson Chrispens
"""

from pymol import cmd

# Sequences from Phillips, slightly edited to align properly in PyMOL
new_caledonia_A = "ATYADTICIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDSHNGKLCLLKGIAPLQLGNCSVAGWILGNPECELLISKESWSYIVETPNPENGTCYPGYFADYEELREQLSSVSSFERFEIFPKESSWPNHTV-TGVSASCSHNGKSSFYRNLLWLTGKNGLYPNLSKSYVNNKEKEVLVLWGVHHPPNIGNQRALYHTENAYVSVVSSHYSRRFTPEIAKRPKVRDQEGRINYYWTLLEPGDTIIFEANGNLIAPWYAFALSRGFGSGIITSNAPMDECDAKCQTPQGAINSSLPFQNVHPVTIGECPKYVRSAKLRMVTGLRNIPSI"
new_caledonia_B = "GLFGAIAGFIEGGWTGMVDGWYGYHHQNEQGSGYAADQKSTQNAINGITNKVNSVIEKMNTQFTAVGKEFNKLERRMENLNKKVDDGFLDIWTYNAELLVLLENERTLDFHDSNVKNLYEKVKSQLKNNAKEIGNGCFEFYHKCNNECMESVKNGTYDYPKYSEESKLNREKIDGVKLE"
wisconsin_A = "STATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGGICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNDESFNWTGVTQNGTSSSCKRRSNNSFFSRLNWLTHLKFKYPALNVTMPNNEKFDKLYIWGVHHPVTDNDQIFLYAQASGRITVSTKRSQQTVIPNIGSRPRIRNIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCNSECITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEK"
wisconsin_B = "GIFGAIAGFIENGWEGMVDGWYGFRHQNSEGIGQAADLKSTQAAINQINGKLNRLIGKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFERTKKQLRENAEDMGNGCFKIYHKCDNACIGSIRNGTYDHDVYRDEALNNRFQIK"

one_letter = {'VAL': 'V', 'ILE': 'I', 'LEU': 'L', 'GLU': 'E', 'GLN': 'Q',
              'ASP': 'D', 'ASN': 'N', 'HIS': 'H', 'TRP': 'W', 'PHE': 'F', 'TYR': 'Y',
              'ARG': 'R', 'LYS': 'K', 'SER': 'S', 'THR': 'T', 'MET': 'M', 'ALA': 'A',
              'GLY': 'G', 'PRO': 'P', 'CYS': 'C'}
three_letter = dict([[v, k] for k, v in one_letter.items()])

def color_muts_H1(resn, resi, chain):
    if chain == "A":
        if one_letter[resn] != new_caledonia_A[int(resi) - 1]:
            print(one_letter[resn], resi, new_caledonia_A[int(resi) - 1], chain)
            cmd.select("to_color", f"resi {resi} and chain A")
            cmd.color("pink", "to_color")
    if chain == "B":
        if one_letter[resn] != new_caledonia_B[int(resi) - 1]:
            print(one_letter[resn], resi, new_caledonia_B[int(resi) - 1], chain)
            cmd.select("to_color", f"resi {resi} and chain B")
            cmd.color("pink", "to_color")


def color_muts_H3(resn, resi, chain):
    if chain == "A":
        if one_letter[resn] != wisconsin_A[int(resi) - 1]:
            print(one_letter[resn], resi,
                  wisconsin_A[int(resi) - 1], chain)
            cmd.select("to_color", f"resi {resi} and chain A")
            cmd.color("pink", "to_color")
    if chain == "B":
        if one_letter[resn] != wisconsin_B[int(resi) - 1]:
            print(one_letter[resn], resi,
                  wisconsin_B[int(resi) - 1], chain)
            cmd.select("to_color", f"resi {resi} and chain B")
            cmd.color("pink", "to_color")



cmd.select("HA", "chain A or chain B")
# cmd.iterate("HA and name CA", "color_muts_H1(resn, resi, chain)")
cmd.iterate("HA and name CA", "color_muts_H3(resn, resi, chain)")
