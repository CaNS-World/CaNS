module mod_common_mpi
  use mod_param, only: dims
  implicit none
  integer :: myid
  integer :: comm_cart,ierr
  integer :: halo(3)
  integer :: ipencil
end module mod_common_mpi
