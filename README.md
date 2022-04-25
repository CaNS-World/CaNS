## Synopsis

**CaNS (Canonical Navier-Stokes)** is a code for massively-parallel numerical simulations of fluid flows. It aims at solving any fluid flow of an incompressible, Newtonian fluid that can benefit from a FFT-based solver for the second-order finite-difference Poisson equation in a 3D Cartesian grid. In two directions the grid is regular and the solver supports the following combination of (homogeneous) boundary conditions:

 * Neumann-Neumann
 * Dirichlet-Dirichlet
 * Neumann-Dirichlet
 * Periodic

In the third domain direction, the solver is more flexible as it uses Gauss elimination. There the grid can also be non-uniform (e.g. fine at the boundary and coarser in the center).

CaNS also allows for choosing an implicit temporal discretization of the diffusion term of the N-S equations. This results in solving a Helmholtz equation for each velocity component. Since FFT-based solvers are also used, the same options described above for pressure boundary conditions apply to the velocity, in case implicit diffusion is active.

**Reference**

P. Costa. *A FFT-based finite-difference solver for massively-parallel direct numerical simulations of turbulent flows.* *Computers & Mathematics with Applications* 76: 1853--1862 (2018). [doi:10.1016/j.camwa.2018.07.034](https://doi.org/10.1016/j.camwa.2018.07.034) [[arXiv preprint]](https://arxiv.org/abs/1802.10323)

## News

[26/12/2021] Implicit temporal discretization of the diffusion term along only one direction (z) is now also supported.

[29/10/2021] **Major update** -- a few neat features have been incorporated in the the most recent version of *CaNS*:

* **x-aligned pencils are now used by default**, which results in improved speed and scalability. This behavior can be changed using the flags `-D_DECOMP_Y`/`-D_DECOMP_Z` for y- or z-aligned pencils;
* **support uneven partitioning of the computational subdomains**: the total number of grid points along one direction does not have to be divisible by the number of tasks;
* simplified `rk.f90` and `mom.f90` routines and the option of implicit diffusion (based on [*SNaC*](https://github.com/p-costa/SNaC));
* improved the routines for imposing boundary conditions, and the MPI I/O checkpointing  (based on [*SNaC*](https://github.com/p-costa/SNaC));
* ... and lots of minor improvements and polishing.

The **many-GPU** version of CaNS can be found [**here**](https://github.com/maxcuda/CaNS).

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

The input file `dns.in` sets the physical and computational parameters. In the `examples/` folder are examples of input files for several canonical flows. See [`src/INFO_INPUT.md`](src/INFO_INPUT.md) for a detailed description of the input file.

Files `out1d.h90`, `out2d.h90` and `out3d.h90` in `src/` set which data are written in 1-, 2- and 3-dimensional output files, respectively. *The code should be recompiled after editing out?d.h90 files*.

### Compilation

#### Prerequisites
The prerequisites for compiling *CaNS* are the following:

 * MPI
 * FFTW3
 * OpenMP (optional)
 * `awk` (to generate dependencies)

#### In short
For most systems, *CaNS* can be compiled from the root directory with the following commands `make library && make -j`.

#### Detailed instructions
The Makefile in root directory is used to compiled the code, and is expected to work out-of-the-box for most systems. The `build.conf` file in the root directory can be used to choose the Fortran compiler (MPI wrapper), a few pre-defined profiles depending on the nature of the run (e.g., production vs debugging), and pre-processing options:

```
#
# compiler and compiling profile
#
FCOMP=GNU           # other options: NVIDIA, INTEL
FFLAGS_OPT=1        # for production runs
FFLAGS_OPT_MAX=0    # for production runs (more aggressive optimization)
FFLAGS_DEBUG=0      # for debugging
FFLAGS_DEBUG_MAX=0  # for thorough debugging
#
# defines
#
DEBUG=1            # best = 1 (no performance penalty)
TIMING=1           # best = 1
IMPDIFF=0          #
IMPDIFF_1D=0       # = 1 requires IMPDIFF=1 too
DECOMP_X=1         # best = 1 if IMPDIFF_1D=0, else 1
DECOMP_Z=0         # best = 0 if IMPDIFF_1D=1, else 0
SINGLE_PRECISION=0 # not a good idea to change
```

In this file, `FCOMP` can be one of `GNU` (`gfortran`), `INTEL` (`ifort`), or `NVIDIA` (`nvfortran`); the predefined profiles for compiler options can be selected by choosing one of the `FFLAGS_*` option; finer control of the compiler flags may be achieved by building with, e.g., `make FFLAGS+=[OTHER_FLAGS]`, or by tweaking the profiles directly under `src/configs/flags.mk`. Finally, the following pre-processing options are available:

 * `DEBUG`            : performs some basic checks for debugging purposes
 * `TIMING`           : wall-clock time per timestep is computed
 * `IMPDIFF`          : diffusion term of the N-S equations is integrated in time with an implicit discretization (thereby improving the stability of the numerical algorithm for viscous-dominated flows)
 * `IMPDIFF_1D`       : same as above, but with implicit diffusion *only* along Z; this option needs to be combined with `IMPDIFF` (required) and `DECOMP_Z` (optional, but recommended for best performance)
 * `SINGLE_PRECISION` : calculation will be carried out in single precision (the default precision is double)

Typing `make library` will build the 2DECOMP&FFT library; then typing `make` will compile the code and copy the executable `cans` and input file `dns.in` to a `run/` folder.

Finally, the choice of compiler `FCOMP` (see `src/configs/flags.mk`), and profile flags `FFLAGS_*` (see `src/configs/flags.mk`) can easily be overloaded, for instance, as: `make FC=ftn FFLAGS=-O2`.

### Running the code

Run the executable with `mpirun` with a number of tasks and (optionally) OpenMP threads complying to what has been set in the input file `dns.in`. Data will be written by default in a folder named `data/`, which must be located where the executable is run.

### Visualizing field data

See [`src/INFO_VISU.md`](src/INFO_VISU.md).

## Notes

I appreciate any feedback that can improve the code. Also, feel free to send case files pertaining to flows not listed in the examples folder.

Please read the `ACKNOWLEDGEMENTS` and `LICENSE` files.

## Contributors

Pedro Costa (p.simoes.costa@gmail.com)
