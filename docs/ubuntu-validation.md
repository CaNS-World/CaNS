# Ubuntu validation

Validation target: `tud`

## Build status

With the current `configs/defaults/libs-default.mk`, the following strict-debug builds succeed on Ubuntu:

- default MPI-IO
- HDF5
- ADIOS2

Build flags used for the checks:

`FFLAGS_DEBUG_MAX=1 FFLAGS_DEBUG=0 FFLAGS_OPT=0 FFLAGS_OPT_MAX=0`

## Runtime status

The following `tests/start_stop` smoke runs completed successfully on Ubuntu:

- default MPI-IO, `mpirun -n 1`
- HDF5, `mpirun -n 1`
- ADIOS2, `mpirun -n 1`

The HDF5 runtime bounds failure was fixed in `src/load.f90` by changing the HDF5 metadata coordinate dummies in `io_field_hdf5()` to assumed-shape arrays. The previous failure came from passing compact subset coordinate vectors into dummies whose lower bounds were tied to field halo extents.

## Strict-debug warnings seen

The only write-time warnings observed in the successful Ubuntu smoke runs were:

- MPI-IO: `Fortran runtime warning: An array temporary was created for argument 'buf' of procedure 'mpi_file_write_all'`
- HDF5: `Fortran runtime warning: An array temporary was created for argument 'buf' of procedure 'h5dwrite_rkind_8_rank_3'`

No corresponding runtime error was observed after the HDF5 metadata fix.

## ADIOS2 matrix on Ubuntu

An ADIOS2 restart/I-O validation matrix was run on Ubuntu with:

- backends: `adios2`, `adios2_blosc`
- decompositions: `np1_p1`, `np6_3x2_p1`, `np6_2x3_p2`
- restart pairs:
  - `np1_p1 -> np6_3x2_p1`
  - `np6_3x2_p1 -> np1_p1`
  - `np6_2x3_p2 -> np6_3x2_p1`

The harness exercised:

- checkpoint/restart writes
- full 3D field output
- 2D slice output

Operationally, all ADIOS2 runs completed cleanly:

- no `ERROR:`
- no `Aborting...`
- no `NaN`
- no `SIGSEGV`

Each run directory produced the expected `.bp` outputs.

## ADIOS2 numeric comparison on Ubuntu

After installing Ubuntu ADIOS2 reader packages, the ADIOS2 restart matrix was rerun with numeric post-read comparison enabled.

An Ubuntu-specific issue was identified in the plain restart path: uncompressed ADIOS2 restart I/O used `SetMemorySelection` directly on haloed arrays, while the compressed path used a packed contiguous buffer. On Ubuntu system ADIOS2 2.9.2, the direct-memory path produced incorrect restart data after decomposition changes.

The fix implemented in `src/load.f90` is to force the packed path for ADIOS2 restart fields in both compressed and uncompressed checkpoint reads/writes.

Final results on Ubuntu after the fix:

- `adios2`: all three restart pairs passed with exact agreement for both `Velocity_Y` and `Pressure_P`
- `adios2_blosc`: all three restart pairs passed with exact agreement for both `Velocity_Y` and `Pressure_P`

Validated restart pairs:

- `np1_p1 -> np6_3x2_p1`
- `np6_3x2_p1 -> np1_p1`
- `np6_2x3_p2 -> np6_3x2_p1`

For all of the above, the observed maximum absolute difference was `0` for both compared fields.

## ADIOS2 packaging notes on Ubuntu

Ubuntu exposes the binary reader tool as `bpls.mpi` instead of `bpls`, so the validation command environment needed a small shim.

Ubuntu's Python package layout for ADIOS2 is also MPI-only in this install:

- `/usr/lib/python3/dist-packages/adios2/adios2_mpi*.so` is present
- `adios2_serial` is not present

The top-level `adios2` Python package falls back to `adios2_serial` unless MPI is explicitly considered active, so `import adios2` from plain `/usr/bin/python3` is currently not reliable on this machine without additional environment control.

## ADIOS2 compression note

For this small Ubuntu validation case, the total output sizes for `adios2` and `adios2_blosc` were essentially identical, so this matrix is not informative about compression benefit. It does, however, confirm that both compressed and uncompressed ADIOS2 restart paths are now correct on Ubuntu after forcing packed restart I/O.
