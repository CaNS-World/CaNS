module mod_types
  use mpi, only: MPI_REAL,MPI_DOUBLE_PRECISION
  integer, parameter, public :: sp = selected_real_kind(6 , 37)
  integer, parameter, public :: dp = selected_real_kind(15,307)
#if defined(_SINGLE_PRECISION)
  integer, parameter, public :: rp = sp
  integer, parameter, public :: MPI_REAL_RP = MPI_REAL
#else
  integer, parameter, public :: rp = dp
  integer, parameter, public :: MPI_REAL_RP = MPI_DOUBLE_PRECISION
#endif
#if defined(_SINGLE_PRECISION_POISSON)
  integer, parameter, public :: gp = sp
#else
  integer, parameter, public :: gp = rp
#endif
end module mod_types
