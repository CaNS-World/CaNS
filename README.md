## Synopsis

**CaNS (Canonical Navier-Stokes)** is a code for massively-parallel numerical simulations of fluid flows. It aims at solving any fluid flow of an incompressible, Newtonian fluid that can benefit from a FFT-based solver for the second-order finite-difference Poisson equation in a 3D Cartesian grid. In two directions the grid is regular and the solver supports the following combination of (homogeneous) boundary conditions:

 * Neumann-Neumann
 * Dirichlet-Dirichlet
 * Neumann-Dirichlet
 * Periodic

In the third domain direction, the solver is more flexible as it uses Gauss elimination. There the grid can also be non-uniform (e.g. fine at the boundary and coarser in the center).

**Update**

CaNS now allows for choosing an implicit temporal discretization of the diffusion term of the N-S equations. This results in solving a Helmholtz equation for each velocity component. Since FFT-based solvers are also used, the same options described above for pressure boundary conditions apply to the velocity, in case implicit diffusion is used.

**Reference**

P. Costa. *A FFT-based finite-difference solver for massively-parallel direct numerical simulations of turbulent flows.* *Computers & Mathematics with Applications* 76: 1853--1862 (2018). [doi:10.1016/j.camwa.2018.07.034](https://doi.org/10.1016/j.camwa.2018.07.034) [[arXiv preprint]](https://arxiv.org/abs/1802.10323)

## News

The recent **many-GPU** version of CaNS can be found [**here**](https://github.com/maxcuda/CaNS).

08/02/2019 -- Input files corresponding to the simulations presented in the manuscript above have been added to `examples/`.

16/05/2019 -- Now a single input file, `dns.in`, can be used to run the executable without recompiling the source. The `examples/` folder has been updated accordingly. The implementation with the former input files (not maintained) can be found in branch `old_input_files`.

20/11/2019 -- Different simulation termination criteria, and a constant body force density have been added as option in the input file `dns.in`, see [`src/INFO_INPUT.md`](src/INFO_INPUT.md) for more details.

04/05/2020 -- Visualization of binary field data has been made easier, see [`src/INFO_VISU.md`](src/INFO_VISU.md) for more details.

13/05/2020 -- An implementation that allows for an arbritary default orientation of the pencil decomposition (i.e. x-, y- or z-aligned pencils) has been developed and kept for now in branch `x_pencil`. Using x-aligned pencils by default increases slightly the performance of CaNS by avoiding two *all-to-all* operations in the Poisson solver. To test this new feature, the code in branch `x_pencil` should be compiled with `-DDECOMP_X` (best performing) or `-DDECOMP_Y`.

## Features

Some features are:

 * Hybrid MPI/OpenMP parallelization
 * FFTW guru interface used for computing multi-dimensional vectors of 1D transforms
 * The right type of transformation (Fourier, Cosine, Sine, etc) is automatically determined form the input file
 * 2DECOMP&FFT routines used for performing global data transpositions and data I/O
 * A different canonical flow can be simulated just by changing the input files

Some examples of flows that this code can solve are:

 * periodic or developing channel
 * periodic or developing square duct
 * tri-periodic domain
 * lid-driven cavity

## Motivation

This project aimed first at being a modern alternative to the well-known FISHPACK routines (Paul Swarztrauber & Roland Sweet, NCAR) for solving a three-dimensional Helmholtz equation. After noticing some works simulating canonical flows with iterative solvers -- when faster direct solvers could have been used instead -- it seemed natural to create a versatile tool and make it available. This code can be used as a first base code for which solvers for more complex flows can be developed (e.g. extensions with fictitious domain methods).

## Method

The fluid flow is solved with a second-order finite-volume pressure correction scheme, discretized in a MAC grid arrangement. Time is advanced with a three-step low storage Runge-Kutta scheme. Optionally, for increased stability at low Reynolds numbers, at the price of higher computational demand, the diffusion term can be treated implicitly. See the reference above for details.

## Usage

### Input file

The input file `dns.in` sets the physical and computational parameters. In the `examples/` folder are examples of input files for several canonical flows. See `src/INFO_INPUT.md` for a detailed description of the input file.

Files `out1d.h90`, `out2d.h90` and `out3d.h90` in `src/` set which data are written in 1-, 2- and 3-dimensional output files, respectively. *The code should be recompiled after editing out?d.h90 files*.

### Compilation

The code should be compiled in `src/`. The prerequisites are the following:

 * MPI
 * FFTW3
 * OpenMP (optional)
 * LAPACK & BLAS (optional)

The Makefile in `src/` should be modified in agreement to the installation paths of each library. Also, the following preprocessor options are available:

 * `-DDEBUG`            : performs some basic checks for debugging purposes
 * `-DTIMING`           : wall-clock time per timestep is computed
 * `-DIMPDIFF`          : diffusion term of the N-S equations is integrated in time with an implicit discretization (thereby improving the stability of the numerical algorithm for viscous-dominated flows)
 * `-DSINGLE_PRECISION` : calculation will be carried out in single precision (the default precision is double)

Typing `make run` will compile the code and copy the executable `cans` and input file `dns.in` to a `run/` folder.

### Running the code

Run the executable with `mpirun` with a number of tasks and shared threads complying to what has been set in the input file `dns.bin`. Data will be written by default in a folder named `data/`, which must be located where the executable is run.

### Visualizing field data

See `src/INFO_VISU.md`.

## Notes

I appreciate any feedback that can improve the code. Also, feel free to send case files pertaining to flows not listed in the examples folder.

Please read the `ACKNOWLEDGEMENTS` and `LICENSE` files.

## Contributors

Pedro Costa (p.simoes.costa@gmail.com)
