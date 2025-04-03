! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_common_mpi
#if defined(_POISSON_PCR_TDMA)
  use decomp_2d, only: decomp_info
#endif
  implicit none
  public
  integer :: myid,ierr
  integer :: halo(3)
  integer :: ipencil_axis
#if defined(_POISSON_PCR_TDMA)
  type(decomp_info) :: dinfo_ptdma
#endif
end module mod_common_mpi
