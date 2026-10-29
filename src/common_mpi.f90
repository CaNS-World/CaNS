! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_common_mpi
  use decomp_2d, only: decomp_info
  implicit none
  public
  integer :: myid,ierr
  integer, allocatable, dimension(:,:) :: halo ! MPI halo datatypes, halo(idir,nh) for halo width nh
  integer :: ipencil_axis
  type(decomp_info) :: dinfo_dtdma
end module mod_common_mpi
