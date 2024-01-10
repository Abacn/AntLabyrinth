# Random Lorentz gas codes

Contents

1. Void percolation threshold calculation in high dimensional periodic box
2. (cross-validation) Void percolation threshold calculation in 3-dimensional periodic box
3. Cage size calculation via static cavity reconstruction
4. Ballistic tracer dynamics in periodic box or in a spherical shell (dynamical cavity reconstruction)
5. Ballistic tracer dynamics in an inhomogeneous RLG in a spherical shell
6. Newtonian dynamics of MK/MKK model under flat/gaussdian shift in Dn periodic box

## Publications

The program led to the following publications. The study involved other programs in this repository are noted in front.

1. G Biroli, P Charbonneau, EI Corwin, Y Hu, H Ikeda, G Szamel, F Zamponi. Interplay between percolation and glassiness in the random Lorentz gas. Phys. Rev. E, 2021, 103, L030104. [DOI:10.1103/PhysRevE.103.L030104](https://doi.org/10.1103/physreve.103.l030104)
1. G Biroli, P Charbonneau, Y Hu, H Ikeda, G Szamel, F Zamponi. Mean-field caging in a random Lorentz gas. J. Phys. Chem. B, 2021, 125, 6244. [DOI:10.1021/acs.jpcb.1c02067](https://doi.org/10.1021/acs.jpcb.1c02067)
1. B Charbonneau, P Charbonneau, Y Hu, Z Yang. High-dimensional percolation criticality and hints of mean-field-like caging of the random Lorentz gas. Phys. Rev. E, 2021, 104, 024137. [DOI:10.1103/PhysRevE.104.024137](https://doi.org/10.1103/physreve.104.024137)
1. G Biroli, P Charbonneau, G Folena, Y Hu, F Zamponi. Local dynamical heterogeneity in simple glass formers. Phys. Rev. Letters, 2022, 128, 175501. [DOI:10.1103/PhysRevLett.128.175501](https://doi.org/10.1103/physrevlett.128.175501)
1. G Folena, G Biroli, P Charbonneau, Y Hu. Equilibrium fluctuations in mean-field disordered models.  Phys. Rev. E, 2022, 106, 024605. [DOI:10.1103/PhysRevE.106.024605](https://doi.org/10.1103/PhysRevE.106.024605)
1. (RLG & MKK) P Charbonneau, Y Hu, PK Morse. Dynamics and fluctuations of minimally-structured glass formers. 2023, [arXiv:2312.07643](https://arxiv.org/abs/2312.07643).

## Void percolation threshold in high dimension (quickhull)

- Source folder: hCombPercolation/

  Compilation prerequisite: 
  - boost library
- Makefile : hcombperco.make
  For general dimension simulation codes, Makefile takes an option DIM=dD (e.g. `make -f hcombperco.make DIM=4D`) to generate binary that simulates a specific dimension
- Example input file : inputfiles/input_RLG_hcombhull.dat

quickhull library modified from https://github.com/tomilov/quickhull


## Void percolation threshold in 3-dimension (CGAL periodic Delaunay package)

- Source folder: RLG_Delaunay_3D/

  Compilation prerequisite: 
  - boost library
  - Eigen3
  - CGAL (working with v4.14)
- Makefile : phip3d.make
- Example input file : inputfiles/input_RLG_Delaunay_3D.dat

## Static cavity reconstruction

- Source folder: RLG_CageShape/

  Compilation prerequisite: 
  - boost library
  - Eigen3
  - CGAL (working with v4.14)
- Makefile : rlgcshape.make
  For general dimension simulation codes, Makefile takes an option DIM=dD (e.g. `make -f rlgcshape.make DIM=4D`) to generate binary that simulates a specific dimension
- Example input file : inputfiles/input_CageShape.dat

## Tracer dynamics

- Source folder: dyshared/ and RLGDynamics/

  Compilation prerequisite: 
  - boost library
- Makefile : rlgdynamics.make
- Example input file (.ini file. Use of this is preferred)
  - inputfiles/input_RLG_dynamics_hcomb.ini : Dynamics in D_n periodic box
- Example input file (legacy format. Kept for rerun old simulations)
  - inputfiles/input_RLG_dynamics_box.dat : Dynamics in hypercubic periodic box
  - inputfiles/input_RLG_dynamics_hcomb.dat : Dynamics in D_n periodic box (and Leech lattice in 24D)
  - inputfiles/input_RLG_dynamics_sphere.dat : Dynamics in a shell

### Run the code

1. Compile. To compile 3D simulation program, for example, run the command
```
make -f rlgdynamics.make DIM=3D
```
This will generate a executable in the folder "../bin" named as "RLG_dy_3D"

2. Run. To run 3D simulation program, in "../bin" folfer run the command
```
./RLG_dy_3D [input file name]
```
The program takes the input filename as an optional parameter. If no input file is specified, it looks for input_RLG_dynamics.dat by default.

3. Get the results
When the program is running, at lease two files are generated, 'msd.dat' and 'dbgsummary.dat', see the next section for details. 'msd.dat' gets updated when every specified number of runs (repeatrun/50) are finished; and 'dbgsummary.dat' gets incremented when every run is finished.

It is useful to write your own scripts for batch run in computing clusters.


### Output file description for dynamics in a shell

- msd.dat: mean square displacements evolving with time. Columns (all columns are unscaled, tracer has unit velocity):
1   t - time
2   MSD - mean square displacement
3   [r^4] - mean quadruple displacement
4   [ <Delta>^2 ] - term in (old definition of) chi_het (<Delta> averaged with both initial and final points of the tracer)
5   [r_i^2 * r_j^2] - term in covariance of square displacement of orthogonal direction
6   [count] - number of cavities

- dbgsummary.dat: summary record for debug and additional information. Every row is a sample
Columns:
1   NCollision - Number of collision
2   Nescape - Number of tracer escaping from neighbor list
3   Finalt - Final time (either set final time or the time of tracer escaped from shell)
4   FinalDelta - Final square displacement at the final time
5   FinalMSD - MSD with only final points averaged
6   FinalMQD - Mean quadruple displacement with only final points averaged
7   InitFinalMSD - MSD with both initial and final points averaged
8   InitFinalMQD - Mean quadruple displacement with both initial and final points averaged
9   Nobstacle - Number of obstacles generated
10   statusCode - 0 - plateaued sufficiently; -1 - the tracer has escaped from the cage (r > r_max); 1 - MSD curve goes up; 2 - MSD curve goes down; 3 - not enough data to calculate MSD

For the cage size sample by sample, use column 5.

To calculate the chi_th of a single tracer, one can use <FinalMQD> - <FinalMSD>^2.

To calculate the chi_var of a single tracer, one can use <InitFinalMQD> - <InitFinalMSD>^2

## Tracer dynamics inhomogeneous RLG

- Source folder: dyshared/ and SqWellDynamics/
- Compilation prerequisites the same as [Tracer dynamics]
- Makefile : rlgdynamics.make , use the "SQRLG_dy" target
  For general dimension simulation codes, Makefile takes an option DIM=dD (e.g. make -f rlgdynamics.make DIM=4D) to generate binary that simulates a specific dimension
- Example input file : inputfiles/input_RLG_dynamics_sqwell.dat
