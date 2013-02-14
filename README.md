POW
===


Parallel Optimization Workbench (POW) aims at simplifying the development of efficient optimization tools
able to deal with complex data structures, we developed.
POW allows the user to perform any kind of data manipulation required for the assessment of the fitness function.
The exploration of the search space is performed by an enhanced version of Particle Swarm Optimization (PSO
Kick and Reseed, PSO-KaR), working in a parallel fashion using MPI libraries.

Requirements
------------

The following python (>=2.5) packages are required:

- numpy
- scipy
- MDAnalysis (for flexible docking modules)
- mpi4py

mpi4py will also require the installation of OpenMPI


How to Launch
-------------

mpiexec -n 4 POW.py module input.dat

where "module" is the type of needed optimization (module), and "input.dat" your setup file


Available modules
-----------------

- DockDimer:       dock two proteins in a heterodimer, given experimental constraints
- DockSymmCircle : rigid/flexible assembly of n monomers according to a circular symmetry
                   and geometric constraints, possibly in the presence of a substrate
- Function:        generic function optimization
