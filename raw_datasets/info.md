**Filtering**
	1. AB Bind: took out homology modeled pdbs, took out mutations that were not on the antibody, took out mutations that were on both antibodies for PDB 1DVF that has two antibodies bound to each other.
	2. SKEMPI: SKEMPI search - fab, antibody, mab, fv for Protein 1 or 2 and webscraped the data to get all possible data, since the downloadable csv did not have all the data filled in. This is used in skempiscrape.py.
	3. When combining data, I favored the data from AB_bind first, then SiPDAB, then SKEMPI, then the rest without overlap.
	4. Database was *not* filtered on PDB resolution/clashscore.
	5. Sampled all the point mutations from Phillips, as well as combinations of interface residues up until 4 mutations.
		1. removed duplicate ddGs that were likely at limit values for the experimental analysis.
		2. PDBs in inputs are homology modeled with the hemagglutinin strain that is in Phillips et al.
	6. Interface residues are filtered as well.
