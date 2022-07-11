# Compiling CaNS

For most systems, CaNS can be compiled from the root directory with the following commands `make library && make`, which will compile the 2DECOMP&FFT library and CaNS. `make clean` clears the CaNS build files, `make libclean` clears the 2DECOMP build, and `make allclean` clears both.

The `Makefile` in root directory is used to compiled the code, and is expected to work out-of-the-box for most systems. The `build.conf` file in the root directory can be used to choose the Fortran compiler (MPI wrapper), a few pre-defined profiles depending on the nature of the run (e.g., production vs debugging), and pre-processing options:

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
DEBUG=1                    # best = 1 (no performance penalty)
TIMING=1                   # best = 1
IMPDIFF=0                  #
IMPDIFF_1D=0               # = 1 requires IMPDIFF=1 too
DECOMP_X=1                 # best = 1 if IMPDIFF_1D=0, else 0
DECOMP_Z=0                 # best = 0 if IMPDIFF_1D=0, else 1
SINGLE_PRECISION=0         # not a good idea to change
SINGLE_PRECISION_POISSON=0 # downcast/upcast correction pressure precision to solve Poisson equation in single precision
```

In this file, `FCOMP` can be one of `GNU` (`gfortran`), `INTEL` (`ifort`), `NVIDIA` (`nvfortran`), or `CRAY` (`ftn`); the predefined profiles for compiler options can be selected by choosing one of the `FFLAGS_*` option; finer control of the compiler flags may be achieved by building with, e.g., `make FFLAGS+=[OTHER_FLAGS]`, or by tweaking the profiles directly under `src/configs/flags.mk`. Similarly, the library paths (e.g., for *FFTW*) may need to be adapted in the `Makefile` (`LIBS` variable) or by building with `make LIBS+='-L[PATH_TO_LIB] -l[NAME_OF_LIB]'`. Finally, the following pre-processing options are available:

 * `DEBUG`                    : performs some basic checks for debugging purposes
 * `TIMING`                   : wall-clock time per timestep is computed
 * `IMPDIFF`                  : diffusion term of the N-S equations is integrated in time with an implicit discretization (thereby improving the stability of the numerical algorithm for viscous-dominated flows)
 * `IMPDIFF_1D`               : same as above, but with implicit diffusion *only* along Z; this option needs to be combined with `IMPDIFF` (required) and `DECOMP_Z` (optional, but recommended for best performance)
 * `SINGLE_PRECISION`         : calculation will be carried out in single precision (the default precision is double)
 * `SINGLE_PRECISION_POISSON` : downcast/upcast correction pressure precision to solve Poisson equation in single precision; this option is compatible with explicit diffusion, or with `IMPDIFF=1 IMPDIFF_1D=1 DECOMP_Z=1`); note that `2DECOMP` needs to be built for single-precision transposes

Typing `make library` will build the 2DECOMP&FFT library; then typing `make` will compile the code and copy the executable `cans` and input file `dns.in` to a `run/` folder.

Finally, the choice of compiler `FCOMP` (see `src/configs/flags.mk`), and profile flags `FFLAGS_*` (see `src/configs/flags.mk`) can easily be overloaded, for instance, as: `make FC=ftn FFLAGS=-O2`.

Note: the older Makefile with explicit dependencies which was used in previous versions to compile CaNS is still present under `src/` (`makefile`). The pre-processing options above can be added there by appending `-D_[FEATURE]` to the variable `OTH` in the `makefile`.
