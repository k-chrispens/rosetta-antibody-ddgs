from pymol import cmd

pdbs = ["1DQJ", "1MHP", "1MLC", "1N8Z", "1VFB", "1YY9", "3GBN",
        "4FQY"]
relax_types = ["harm", "harm_025", "unconst", "crtsc", "clean"]
f = open("rmsds.csv", "w")
f.write("#PDB,type,RMSD\n")
for pdb in pdbs:
    cmd.load(f"./{pdb}.pdb")
    print(f"Loaded {pdb}")
    for type in relax_types:
        cmd.load(f"./{pdb}_{type}.pdb")
        print(f"Loaded {pdb}_{type}")
        rmsd = cmd.align(f"{pdb} and not hydro", f"{pdb}_{type} and not hydro", cycles = 0)
        print(rmsd)
        f.write(f"{pdb},{type},{rmsd[0]}\n")
        cmd.delete(f"{pdb}_{type}")
    cmd.delete(f"{pdb}")
f.close()
