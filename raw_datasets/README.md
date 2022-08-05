## Datasets

This is a collection of all the datasets used in the process of producing the final dataset, which is contained in **use_this_data.csv**. The other datasets are either original datasets from SKEMPI, AB Bind, SiPMAB, Phillips et al., or Mason et al. 

To reproduce the dataset from the base data: 
(make sure you have installed and have activated the conda environment in the base directory)
1. Run **skempiscrape.py**.
2. Run **cleaning.py**.
3. Run **interface.py**.
4. Run the cells in **add_attributes.ipynb**
5. Remove duplicate ddG values (I did this in Excel, but it is easily done programmatically as well).

## Filtering

1. AB Bind: took out homology modeled PDBs, took out mutations that were not on the antibody, took out mutations that were on both antibodies for PDB 1DVF that has two antibodies bound to each other.
	
2. SKEMPI: SKEMPI search - fab, antibody, mab, fv for Protein 1 or 2 and webscraped the data to get all possible data, since the downloadable csv did not have all the data filled in. The link for the page is available and used in **skempiscrape.py**.
	
3. When combining data, I favored the data from AB_bind first, then SiPDAB, then SKEMPI, then the rest without overlap.
	
4. Database was *not* filtered on PDB resolution/clashscore.

5. Sampled all the point mutations from Phillips, as well as combinations of interface residues up until 4 mutations.
		
	1. removed duplicate ddGs that were likely at limit values for the experimental analysis.
	
	2. PDBs in inputs are homology modeled with the hemagglutinin strain that is in Phillips et al.
	
6. Filtered for interface residues.
