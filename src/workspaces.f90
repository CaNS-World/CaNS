! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_workspaces
#if defined(_OPENACC)
  use mod_types
  use mod_common_cudecomp, only: work,work_cuda,work_halo,work_halo_cuda,work_ptdma,work_ptdma_cuda
  use mod_utils, only: f_sizeof
  implicit none
  private
  public init_wspace_arrays,set_cufft_wspace
contains
  subroutine init_wspace_arrays
    use mod_common_cudecomp, only: handle,gd_halo,gd_poi,gd_ptdma, &
                                   ap_x_poi,ap_y_poi,ap_z_poi, &
                                   solver_buf_0,solver_buf_1, &
                                   ap_z,pz_aux_1, &
                                   istream_acc_queue_1
    use mod_fft            , only: wsize_fft
    use mod_param          , only: ng,cudecomp_is_t_in_place,cbcpre,ipencil => ipencil_axis,is_poisson_pcr_tdma
    use cudecomp
    use openacc
    implicit none
    integer :: istat
    integer(i8) :: wsize,max_wsize
    integer :: nh(3)
    !
    ! allocate cuDecomp workspace buffer for transposes (reused for FFTs)
    !
    max_wsize = -1
    wsize     = wsize_fft
    max_wsize = max(max_wsize,wsize)
    istat = cudecompGetTransposeWorkspaceSize(handle,gd_poi,wsize)
    max_wsize = max(max_wsize,wsize)
    !
    ! one can append more checks here (e.g., for spectra calculation)
    !
    allocate(work(max_wsize))
    istat = cudecompMalloc(handle,gd_poi,work_cuda,max_wsize)
    call acc_map_data(work,work_cuda,max_wsize*f_sizeof(work(1)))
    !
    ! allocate cuDecomp workspace buffer for halos
    ! (needs to be different due to the possible need of an NVSHMEM allocator
    !  in one of the descriptors, rather than a simple cudaMaloc)
    !
    nh(:) = 1
    istat = cudecompGetHaloWorkspaceSize(handle,gd_halo,ipencil,nh,max_wsize)
    allocate(work_halo(max_wsize))
    !
    istat = cudecompMalloc(handle,gd_halo,work_halo_cuda,max_wsize)
    call acc_map_data(work_halo,work_halo_cuda,max_wsize*f_sizeof(work_halo(1)))
    !
    ! allocate transpose buffers
    !
    wsize = max(ap_x_poi%size,ap_y_poi%size,ap_z_poi%size)
    allocate(solver_buf_0(wsize),solver_buf_1(wsize))
    !$acc enter data create(solver_buf_0,solver_buf_1)
    if(cbcpre(0,3)//cbcpre(1,3) == 'PP') then
      allocate(pz_aux_1(ap_z%shape(1),ap_z%shape(2),ap_z%shape(3)))
      !$acc enter data create(pz_aux_1)
    end if
    if(is_poisson_pcr_tdma) then
      if(.not.allocated(pz_aux_1)) then
        allocate(pz_aux_1(ap_z%shape(1),ap_z%shape(2),ap_z%shape(3)))
        !$acc enter data create(pz_aux_1)
      end if
      !
      ! allocate pcr-tdma transpose workspaces: separate buffer is needed because work is used along with work_ptdma
      !
      istat = cudecompGetTransposeWorkspaceSize(handle,gd_ptdma,wsize)
      wsize = max(wsize,(3*(ng(3)+1))) ! work_ptdma also use as a buffer with this size in `gaussel_ptdma_gpu_fast_1d`
      allocate(work_ptdma(wsize))
      istat = cudecompMalloc(handle,gd_ptdma,work_ptdma_cuda,wsize)
      call acc_map_data(work_ptdma,work_ptdma_cuda,wsize*f_sizeof(work_ptdma(1)))
    end if
    !
    istream_acc_queue_1 = acc_get_cuda_stream(1) ! fetch CUDA stream of OpenACC queue 1
  end subroutine init_wspace_arrays
  subroutine set_cufft_wspace(arrplan,istream)
    use cufft
    use openacc
    !
    ! to be done after initializing all FFTs and allocating work
    !
    implicit none
    integer                 , intent(in)           :: arrplan(:)
    integer(acc_handle_kind), intent(in), optional :: istream
    integer :: istat,i
    do i=1,size(arrplan)
      !$acc host_data use_device(work)
      istat = cufftSetWorkArea(arrplan(i),work)
      !$acc end host_data
      if(present(istream)) then
        istat = cufftSetStream(arrplan(i),istream)
      end if
    end do
  end subroutine set_cufft_wspace
#endif
end module mod_workspaces
