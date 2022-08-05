## Relaxing PDBs
These are the PDBs in the dataset, I relaxed these using various methods as documented below. The best seemed to be harmonic constraints within 9 Å and a 0.25 sd for the harmonic function.

Naming convention:
\_clean: constrain_relax_to_start_coords on side chains only.
\_crtsc: constrain_relax_to_start_coords on side chains and backbone.
\_harm: harmonic constraint (using atom_pair_constraints) on all pairs CA atoms within 9 of each other, set sd to be 0.5 Å, weight of 1. 
\_harm_025: same as above, but with sd set to 0.25 Å.
\_all: same as above, just renamed for my run on all PDBs so I could then do ddG predictions.
\_noHmg: non-homology modeled hemagglutinin for 4FQY and 3GBN.

## Additional notes
Removed duplicate chains in PDBs: 1BJ1, 1MHP, 1MLC, 3NGB, 3BN9, 2NZ9, 2B2X, and 1CZ8 required this. I re-relaxed these structures after removing these chains.

These PDBs must still be cleaned to get the fold trees set up properly. Run `pose.clean()` or `get_pose_and_cleanup(path_to_pdb)` to achieve this. The inputs folder in the main dir has all of the \_all files cleaned.
