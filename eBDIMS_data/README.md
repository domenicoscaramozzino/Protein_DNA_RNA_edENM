This folder contains the Fortran code to run eBDIMS transitions between two end states, as well as some pre-computed examples.

To use the code, first compile it with gfortran as:
gfortran eBDIMS_par.f90 -o eBDIMS_par -fopenmp

We use OpenMP for parallelization, using 16 threads by default. You can change this number modifying the parameter "num_threads" directly in the source code. Compilation without "-fopenmp" creates an executable that operates without core threading.

To run the simulation, simply import the two PDBs of the end states in your folder, e.g. 7SWQ.pdb and 7SWF.pdb, and type:
./eBDIMS_par 7SWQ ABD 7SWF ABD 500 99.9 0.8

The first PDB code in the input list will be the reference/starting structure, while the second is the target.
The chain lables ("ABD") are used to select which chains are used from both structures to select bead particles and drive the transition.
Note that bead correspondence for the transition will be computed by associating the corresponding chain not by chain label, but by chain ordering, e.g. typing "ABD" and "ABD" will create A-A, B-B, D-D correspondences, while typing "ABD" and "BAD" will create A-B, B-A, and D-D correspondences. This is done because often deposited PDB files have mismatches in chain labelling.
The first numeric input ("500") corresponds to the number of improvement steps for frame saving. Often, 500 works just fine. If you want more frequent frames, decrease this number.
The second numeric input ("99.9") is the % of transition convergence to end the simulation. Typically, 99.9% is enough to get close to the target.
The last numeric input ("0.8") is the RMSD convergence to end the simulation. Note that the eBDIMS code does not perform structural alignment, so adding this parameters makes sense only if you provide aligned input files - if not, RMSD values are inflated by rigid translations and rotations and should be corrected a posteriori.

More information can be found in [domenico](https://github.com/domenicoscaramozzino/eBDIMS2)
