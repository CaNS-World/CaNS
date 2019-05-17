module mod_types
  use mpi
#ifdef SINGLE_PRECISION
  integer, parameter, public :: rp = KIND(1.0)
  integer, parameter, public :: MPI_REAL_RP = MPI_REAL
#else
  integer, parameter, public :: rp = KIND(0.0D0)
  integer, parameter, public :: MPI_REAL_RP = MPI_DOUBLE_PRECISION
#endif
end module mod_types
