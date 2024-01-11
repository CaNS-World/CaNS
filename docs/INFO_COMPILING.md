# Compiling CaNS

For most systems, CaNS can be compiled from the root directory with the following commands `make libs && make`, which will compile the 2DECOMP&FFT/cuDecomp libraries and CaNS. `make clean` clears the CaNS build files, `make libsclean` clears the 2DECOMP/cuDecomp builds, and `make allclean` clears both.

The `Makefile` in root directory is used to compiled the code, and is expected to work out-of-the-box for most systems. The `build.conf` file in the root directory can be used to choose the Fortran compiler (MPI wrapper), a few pre-defined profiles depending on the nature of the run (e.g., production vs debugging), and pre-processing options:

```
#
# compiler and compiling profile
#
FCOMP=GNU           # options: GNU, NVIDIA, INTEL
FFLAGS_OPT=1        # for production runs
FFLAGS_OPT_MAX=0    # for production runs (more aggressive optimization)
FFLAGS_DEBUG=0      # for debugging
FFLAGS_DEBUG_MAX=0  # for thorough debugging
#
# defines
#
DEBUG=1                    # best = 1 (no significant performance penalty)
TIMING=1                   # best = 1
IMPDIFF=0                  #
IMPDIFF_1D=0               #
PENCIL_AXIS=1              # = 1/2/3 for X/Y/Z-aligned pencils
SINGLE_PRECISION=0         # perform the whole calculation in single precision
#
# GPU-related
#
GPU=0
USE_NVTX=0
```

In this file, `FCOMP` can be one of `GNU` (`gfortran`), `INTEL` (`ifort`), `NVIDIA` (`nvfortran`), or `CRAY` (`ftn`); the predefined profiles for compiler options can be selected by choosing one of the `FFLAGS_*` option; finer control of the compiler flags may be achieved by building with, e.g., `make FFLAGS+=[OTHER_FLAGS]`, or by tweaking the profiles directly under `configs/flags.mk`. Similarly, the library paths (e.g., for *FFTW*) may need to be adapted in the `Makefile` (`LIBS` variable) or by building with `make LIBS+='-L[PATH_TO_LIB] -l[NAME_OF_LIB]'`. Finally, the following pre-processing options are available:

 * `DEBUG`                    : performs some basic checks for debugging purposes
 * `TIMING`                   : wall-clock time per time step is computed
 * `IMPDIFF`                  : diffusion term of the N-S equations is integrated in time with an implicit discretization (thereby improving the stability of the numerical algorithm for viscous-dominated flows)
 * `IMPDIFF_1D`               : same as above, but with implicit diffusion *only* along Z; *for optimal parallel performance this option should be combined with* `PENCIL_AXIS=3`
 * `PENCIL_AXIS`              : sets the default pencil direction, one of [1,2,3] for [X,Y,Z]-aligned pencils; X-aligned is the default and should be optimal for all cases except for Z implicit diffusion, where using Z-pencils is recommended
 * `SINGLE_PRECISION`         : calculation will be carried out in single precision (the default precision is double)
 * `GPU`                      : enable GPU accelerated runs (requires the `FCOMP=NVIDIA`)
 * `USE_NVTX`                 : enable [NVTX](https://docs.nvidia.com/nsight-visual-studio-edition/nvtx) markers to tag certain code regions and assist with profiling

Typing `make libs` will build the 2DECOMP&FFT/cuDecomp libraries; then typing `make` will compile the code and copy the executable `cans` to a `run/` folder; `make run` will also copy the default input files `*.in` under `src/` to the same `run/` folder.

Note that cuDecomp needs to be dynamically linked before performing a GPU run. To do this, one should update the `LD_LIBRARY_PATH` environment variable as follows (from the root directory):
```
export LD_LIBRARY_PATH=$PWD/dependencies/cuDecomp/build/lib:$LD_LIBRARY_PATH
```

Finally, the choice of compiler `FCOMP` (see `configs/flags.mk`), and profile flags `FFLAGS_*` (see `configs/flags.mk`) can easily be overloaded, for instance, as: `make FC=ftn FFLAGS=-O2`. Linking and Include options can be changed in `configs/libs.mk`.
