## Synopsis

**CaNS (Canonical Navier-Stokes)** is a code for massively-parallel numerical simulations of fluid flows. It aims at solving any fluid flow of an incompressible, Newtonian fluid that can benefit from a FFT-based solver for the second-order finite-difference Poisson equation in a 3D Cartesian grid. In two directions the grid is regular and the solver supports the following combination of (homogeneous) boundary conditions:

 * Neumann-Neumann
 * Dirichlet-Dirichlet
 * Neumann-Dirichlet
 * Periodic

In the third domain direction, the code is more flexible as it uses Gauss elimination. There the grid can also be non-uniform.

**Update**

CaNS now allows for choosing an implicit temporal discretization of the diffusion term of the N-S equations. This results in solving a Helmholtz equation for each velocity component. Since FFT-based solvers are also used, the same options described above for pressure BCs apply to the velocity boundary conditions.

**Reference**

P. Costa. *A FFT-based finite-difference solver for massively-parallel direct numerical simulations of turbulent flows.* [arXiv:1802.10323](https://arxiv.org/pdf/1802.10323.pdf) (2018).

## Features

Some features are:

 * Hybrid MPI/OpenMP parallelization
 * FFTW guru interface used for computing multi-dimensional vectors of 1D transforms
 * The right type of transformation (Fourier, Cosine, Sine, etc) is automatically determined form the input file
 * 2DECOMP&FFT routines used for performing global data transpositions and data i/o
 * A different canonical flow can be simulated just by changing the input files

Some examples of flows that this code can solve are:

 * periodic or developping channel
 * periodic or developping square duct
 * tri-periodic domain
 * lid-driven cavity

## Motivation

This project aimed first at being a modern alternative to the well-known FISHPACK routines (Paul Swarztrauber & Roland Sweet, NCAR) for solving a three-dimensional Helmholtz equation. After noticing some works simulating canonical flows with iterative solvers -- when faster direct solvers could have been used instead -- it seemed natural to create a versatile tool and make it available. This code can be used as a first base code for which solvers for more complex flows can be developed (e.g. extensions with ficticious domain methods).

## Method

The fluid flow is solved with a second-order finite-volume pressure correction scheme, discretized in a MAC grid arrangement. Time is advanced with a three-step low storage Runge-Kutta scheme. Optionally, for increased stability at low Reynolds numbers, at the price of higher computational demand, the diffusion term can be treated implicitly. See the arXiv reference above for details.

## Usage

### Input files

The input (header) files inside the `src/` folder, `setup.h90` and `bc.h90` setup a case. `setup.h90` sets most of the physical and computational parameters, and `bc.h90` the boundary conditions for the Pressure and velocity fields, together with some options for forcing the flow and inflow/outflow conditions. The comments in these files make them self-explanatory.

In the `examples/` folder are examples of these files for several canonical flows.

The files `out1d.h90`, `out2d.h90` and `out3d.h90` set which data are written in 1-, 2- and 3-dimensional output files, respectively. The corresponding output frequency is set in `setup.h90`.

### Compilation

After modifying files `setup.h90`, `bc.h90` and `out*d.h90`, the code should be compiled. 

The prerequisites are the following:

 * MPI
 * FFTW3
 * LAPACK
 * BLAS
 * OpenMP (optional)

The Makefile should be modified in agreement to the installation paths of each library. Also, the following preprocessor options are available:

 * `-DDEBUG`   : performs some basic checks for debugging purposes
 * `-DTIMING`  : wall-clock time per timestep is computed
 * `-DIMPDIFF` : diffusion term of the N-S equations is treated implicitly (thereby improving the stability of the numerical algorithm for viscous-dominated flows)

### Running the code

Run the executable with `mpirun` with a number of tasks and shared threads complying to what has been set in the input file `setup.h90`. Data will be written by defailt in a folder named `data`, which must be located where the executable is run.

## Notes

I appreciate any feedback that can improve the code. Also, feel free to send case files pertaining to flows not listed in the examples folder.

Please read the `ACKNOWLEDGEMENTS` and `LICENSE` files.

## Contributors

Pedro Costa (p.simoes.costa@gmail.com)
