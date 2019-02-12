# Dynamics around the Site Percolation Threshold on Hypercubic Lattices

Calculating the percolation threshold and dynamical exponents of the site percolation. The model is as the following: n-dimensional (hyper)cubic lattice, with each site is randomly assigned as occupied with the probability, or "cover fraction" p. An ant on a occupied site randomly walk through the neighboring occupied sites. At each time step, the ant randomly chooses a neighbor site (in d-dimensional hypercube, a site has 2*d neighbors), and moves onto that site if it is occupied, otherwise remains in situ. The neighboring occupied sites is named as a "cluster" in this model. Hence the ant one one cluster can only reach the sites on this cluster.

## Summary

- Language: C++
- Tested compiler Intel C++ compiler icpc version 19.0.0.117 (gcc version 4.9.0 compatibility)
- Platform: Linux or Mac
- Publication: Giulio Biroli, Patrick Charbonneau and Yi Hu. \"Dynamics around the Site Percolation Threshold on High-Dimensional Hypercubic Lattices\". Physical Review E 99, 022118 (2019). 
- URL: http://doi.org/10.1103/PhysRevE.99.022118
## Compiling options

Makefile has multiple labels which correspond to different programs and implementations:

- ant
  Run dynamical simulations of an ant on the lettice. There are different implementations based on different background data structure, including:
  - ant_ord
     The visited sites are stored in std::set, based on red-black tree, with O(log n) query time. It is adapted to the case of p>p_c where large clusters might be possible.
  - ant_und
     The visited sites are stored in std::unordered_set, based on hash table, with O(1) query time with more memory . It is adapted to the case of p <= p_c where large clusters is unlikely, otherwise the program may run out of the memory when a large cluster is encountered.
  - ant_hbd
     Solving dynamics on small clusters by transfer matrix and do simulations only on large clusters. Uses ant_und implementation as the parent class, which means the simulation part is the same as ant_und. One can easily modify the code to call ant_ord implementation if necessary.

- msdlim
  Calculate the infinite time limit of MSD on a lattice

- leath
  Calculate the infinite time limit of MSD by Leath algorithm, i.e., growing a cluster from the origin

- invade
  Invade percolation algorithm to compute the percolation threshold. An implementation of the method of S Mertens, and C Moore. Percolation Thresholds and Fisher Exponents in Hypercubic Lattices. Phys. Rev. E 98, 022120 (2018)

## File list

- AntLatticeWalk
  - ant.hpp
  - cluster_counter.hpp
  - leath_counter.hpp
  - ant_global.cpp
  - coordinate.h
  - leathsitenode.hpp
  - ant_global.hpp
  - grid_field.h
  - limited_queue.hpp
  - ant_main.cpp
  - invade.cpp
  - msdlimit.cpp
  - ant_ord.cpp
  - invade.hpp
  - read_input.cpp
  - ant_und.cpp
  - leath.cpp
  - read_input.h
  - cluster_counter.cpp
  - leath_counter.cpp
- DIM
  - *D/dim.h Headers define the dimension of systems when compiling
- HybridAnt
  - ant_hbd_main.cpp
  - ant_hbd.hpp
  - ant_mat.cpp
  - ant_hbd.cpp
- MSDDumpReader  Small data processing programs to compute the local $\mu_-$ 
  (MSDDumpReader in d<6; MSDDumpReaderB for the prefactor in d >= 6)

### Example input files
- inputfile_ant.dat
  for ant_ord, ant_und, ant_hbd
- inputfile_ivd.dat
  for invade
- inputfile_leath.dat
  for leath
- inputfile.dat
  for msdlim

## To compile

One may need to manually configure the environment variables INTEL_ROOT, MKL_ROOT to target the path of intel compiler and MKL librarys. Only ant_hbd need MKL_librarys; other programs should be compatible on g++, clang++ and other C++ compilers.

To compiler one of the target, for example, leath, for 3D systems, run command
make leath DIM=3D

and so on.

