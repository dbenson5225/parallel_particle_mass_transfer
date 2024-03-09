# parallel_particle_mass_transfer
Shared-memory mass-transfer particle-tracking

These codes are used in the paper "Parallelization of particle-mass-transfer algorithms on shared-memory, multi-core CPUs", submitted to Advances in Water Resources.  All of the codes solve an inter-particle mass transfer routine that is similar to traditional SPH but uses a physics-based kernel.  The many particle-particle dependencies require parallelization.

There are two directories with fortran code and one with matlab.  

The first fortran folder contains a sparse-matrix mass-transfer routine that uses classical domain decomposition to induce parallelism.  Each geometric subdomain also relys on "ghost particles" from neighboring domains within the mixing distance.  The codes are written to use shared-memory parallelism via OpenMP directives.  Therefore, one might use gfortran to compile the fortran in the sparse-matrix folder via the commands at the command prompt:

gfortran -g -fbacktrace -fcheck=all -fopenmp -O3 -c sparse_particle_module_2.f90 sparse_main_omp.f90 qtree_module.f90 

gfortran -g -fbacktrace -fcheck=all -fopenmp -O3 -o sparse.exe sparse_particle_module_2.o sp_main_omp.o qtree_module.o

Note that these commands have extra debugging directives that may be omitted.  Also, because the compilation looks for already-compiled modules, you can just run the first command twice and ignore error warnings.  Also, the code allows 1-, 2-, or 3-dimensional problems, so there are a few error warnings for accessing extra dimensions of a few arrays that don't exist - those are embedded within if statements and cannot be reached to give a real run-time error.  

After compilation, the code prompts for the name of an input file with problem parameters (e.g., domain size, number of particles, etc.) One is provided named "input.dat".  The code also prompts for a few output filenames.  A file is provided that may be edited called "filenames", so that after compilation, the user may just use re-directs and run it by entering (for Linix-based syntax):

./sparse.exe<filenames

The second fortran folder has code for a new non-geometrix domain decomposition based solely on quad-tree particle location and row-normalization.  This allows each row of the mass-transfer matrix to be calculated on its own on any processing core.  These codes would be similarly compiled by:

gfortran -g -fbacktrace -fcheck=all  -fopenmp -O3 -c list_particle_module_omp.f90 list_main_omp.f90 qtree_module.f90

gfortran -g -fbacktrace -fcheck=all  -fopenmp -O3 -o list.exe list_particle_module_omp.o list_main_omp.o qtree_module.o

Similar to the sparse-matrix instructions above, run the code by entering

./list.exe<filenames

MATLAB CODE
The matlab folder contains three .m files: One each for full-matrix, sparse-matrix, and list-based computations.  Note that these codes test the problem decomposition schemes in serial, i.e., no explicit multi-threading or multi-core matlab directives are used. Of course matlab does some multi-threading on its own.  The user may specify different levels of domain decomposition using e.g., the vector Sx = [1 3 5], which would test geometric division into 1 (no subdivision), 3, or 5 in each dimension.  

If you use these codes, please cite their description in "Benson, D. A., I. Pribec, N.B. Engdahl, S. Pankavich, and L. Schauer, Parallelization of particle-mass-transfer algorithms on shared-memory, multi-core CPUs, submitted, Advances in Water Resources, 2024. 
