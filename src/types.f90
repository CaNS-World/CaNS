module mod_types
  use mpi, only: MPI_REAL,MPI_DOUBLE_PRECISION
#ifdef SINGLE_PRECISION
  integer, parameter, public :: rp = selected_real_kind(6 , 37)
  integer, parameter, public :: MPI_REAL_RP = MPI_REAL
#else
  integer, parameter, public :: rp = selected_real_kind(15,307)
  integer, parameter, public :: MPI_REAL_RP = MPI_DOUBLE_PRECISION
#endif
end module mod_types
