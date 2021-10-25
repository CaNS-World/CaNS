module mod_common_mpi
  use mpi
  use mod_param, only: dims
  implicit none
  integer :: myid
  integer :: comm_cart,ierr
  integer :: halo(3)
end module mod_common_mpi
