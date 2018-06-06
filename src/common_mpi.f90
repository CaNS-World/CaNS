module mod_common_mpi
  use mpi
  use mod_param, only: dims
  implicit none
  integer :: myid
  integer :: left,right,front,back
  integer, dimension(2) :: coord
  integer :: comm_cart,ierr
  integer :: xhalo,yhalo
  integer :: status(MPI_STATUS_SIZE)
end module mod_common_mpi
