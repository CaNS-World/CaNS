! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_common_cudecomp
#if defined(_OPENACC)
  use mod_types
  !@cuf use cudafor
  use cudecomp
  use openacc
  use mod_param, only: cudecomp_is_t_in_place
  implicit none
  public
  integer :: cudecomp_real_rp
  type(cudecompHandle)     :: handle
  type(cudecompGridDesc)   :: gd_halo,gd_poi,gd_poi_io
  type(cudecompPencilInfo) :: ap_x,ap_y,ap_z,ap_x_poi,ap_y_poi,ap_z_poi
  type(cudecompGridDesc)   :: gd_ptdma
  type(cudecompPencilInfo) :: ap_y_ptdma,ap_z_ptdma
  !
  ! workspace stuff
  !
  integer(i8) :: wsize_fft
  real(rp), pointer, contiguous, dimension(:) :: work     ,work_cuda
  real(rp), pointer, contiguous, dimension(:) :: work_halo,work_halo_cuda
  !@cuf attributes(device) :: work_cuda,work_halo_cuda
  real(rp), pointer, contiguous, dimension(:) :: work_ptdma,work_ptdma_cuda
  !@cuf attributes(device) :: work_halo_cuda,work_ptdma_cuda
  real(rp), target, allocatable, dimension(:) :: solver_buf_0,solver_buf_1
  real(rp), allocatable, dimension(:,:,:) :: pz_aux_1,pz_aux_2
  integer(acc_handle_kind) :: istream_acc_queue_1
#endif
end module mod_common_cudecomp
