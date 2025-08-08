! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2025 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_workspaces
#if defined(_OPENACC) || defined(_OPENMP)
  use mod_types
  use mod_common_cudecomp, only: work,work_cuda,work_halo,work_halo_cuda,work_ptdma,work_ptdma_cuda
  use mod_utils          , only: f_sizeof
  implicit none
  private
  public init_wspace_arrays,set_cufft_wspace,cudecomp_finalize
contains
  subroutine init_wspace_arrays
    use mod_common_cudecomp, only: handle,gd_halo,gd_poi,gd_poi_io,gd_ptdma, &
                                   ap_x_poi,ap_y_poi,ap_z_poi, &
                                   solver_buf_0,solver_buf_1, &
                                   ap_z,pz_aux_1, &
                                   istream_acc_queue_1,istream_acc_queue_1_comm_lib
    use mod_fft            , only: wsize_fft
    use mod_param          , only: ng,cudecomp_is_t_in_place,cbcpre,ipencil => ipencil_axis,is_poisson_pcr_tdma, &
                                   is_use_diezdecomp
#if !defined(_USE_DIEZDECOMP)
    use cudecomp
#else
    use diezdecomp
#endif
#if   defined(_OPENACC)
    use openacc
#elif defined(_OPENMP)
    use omp_lib
#endif
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
    istat = cudecompGetTransposeWorkspaceSize(handle,gd_poi_io,wsize)
    max_wsize = max(max_wsize,wsize)
    !
    ! one can append more checks here (e.g., for spectra calculation)
    !
    allocate(work(max_wsize))
    !$acc        enter data create(   work) if(is_use_diezdecomp)
    !$omp target enter data map(alloc:work) if(is_use_diezdecomp)
#if !defined(_USE_DIEZDECOMP)
    istat = cudecompMalloc(handle,gd_poi,work_cuda,max_wsize)
#if   defined(_OPENACC)
    call acc_map_data(work,work_cuda,max_wsize*f_sizeof(work(1)))
#elif defined(_OPENMP)
    istat = omp_target_associate_ptr(c_loc(work(1)),c_loc(work_cuda(1)),max_wsize*f_sizeof(work(1)), &
                                     0_i8,omp_get_default_device())
#endif
#endif
    !
    ! allocate cuDecomp workspace buffer for halos
    ! (needs to be different due to the possible need of an NVSHMEM allocator
    !  in one of the descriptors, rather than a simple cudaMaloc)
    !
    nh(:) = 1
    istat = cudecompGetHaloWorkspaceSize(handle,gd_halo,ipencil,nh,max_wsize)
    allocate(work_halo(max_wsize))
    !
    !$acc        enter data create(   work_halo) if(is_use_diezdecomp)
    !$omp target enter data map(alloc:work_halo) if(is_use_diezdecomp)
#if !defined(_USE_DIEZDECOMP)
    istat = cudecompMalloc(handle,gd_halo,work_halo_cuda,max_wsize)
#if   defined(_OPENACC)
    call acc_map_data(work_halo,work_halo_cuda,max_wsize*f_sizeof(work_halo(1)))
#elif defined(_OPENMP)
    istat = omp_target_associate_ptr(c_loc(work_halo(1)),c_loc(work_halo_cuda(1)),max_wsize*f_sizeof(work_halo(1)), &
                                     0_i8,omp_get_default_device())
#endif
#endif
    !
    ! allocate transpose buffers
    !
    wsize = max(ap_x_poi%size,ap_y_poi%size,ap_z_poi%size)
    allocate(solver_buf_0(wsize),solver_buf_1(wsize))
    !$acc        enter data create(   solver_buf_0,solver_buf_1)
    !$omp target enter data map(alloc:solver_buf_0,solver_buf_1)
    if(cbcpre(0,3)//cbcpre(1,3) == 'PP') then
      allocate(pz_aux_1(ap_z%shape(1),ap_z%shape(2),ap_z%shape(3)))
      !$acc        enter data create(   pz_aux_1)
      !$omp target enter data map(alloc:pz_aux_1)
    end if
    if(is_poisson_pcr_tdma) then
      if(.not.allocated(pz_aux_1)) then
        allocate(pz_aux_1(ap_z%shape(1),ap_z%shape(2),ap_z%shape(3)))
        !$acc        enter data create(   pz_aux_1)
        !$omp target enter data map(alloc:pz_aux_1)
      end if
      !
      ! allocate pcr-tdma transpose workspaces: separate buffer is needed because work is used along with work_ptdma
      !
      istat = cudecompGetTransposeWorkspaceSize(handle,gd_ptdma,wsize)
      wsize = max(wsize,(3*(ng(3)+1))) ! work_ptdma also use as a buffer with this size in `gaussel_ptdma_gpu_fast_1d`
      allocate(work_ptdma(wsize))
      !$acc        enter data create(   work_ptdma) if(is_use_diezdecomp)
      !$omp target enter data map(alloc:work_ptdma) if(is_use_diezdecomp)
#if !defined(_USE_DIEZDECOMP)
      istat = cudecompMalloc(handle,gd_ptdma,work_ptdma_cuda,wsize)
#if   defined(_OPENACC)
      call acc_map_data(work_ptdma,work_ptdma_cuda,wsize*f_sizeof(work_ptdma(1)))
#elif defined(_OPENMP)
      istat = omp_target_associate_ptr(c_loc(work_ptdma(1)),c_loc(work_ptdma_cuda(1)),wsize*f_sizeof(work_ptdma(1)), &
                                       0_i8,omp_get_default_device())
#endif
#endif
    end if
    !
#if   defined(_OPENACC)
#if !defined(_USE_HIP)
    istream_acc_queue_1 = acc_get_cuda_stream(1) ! fetch CUDA stream of OpenACC queue 1
#else
    istream_acc_queue_1 = acc_get_hip_stream(1)
#endif
#elif defined(_OPENMP)
    istream_acc_queue_1 = 0
#endif
    istream_acc_queue_1_comm_lib = istream_acc_queue_1
#if defined(_USE_DIEZDECOMP)
    istream_acc_queue_1_comm_lib = 1
#endif
  end subroutine init_wspace_arrays
  subroutine set_cufft_wspace(arrplan,istream)
#if !defined(_USE_HIP)
    use cufft
#else
    use, intrinsic :: iso_c_binding, only: C_PTR,c_loc
    use hipfort_hipfft, only: cufftSetWorkArea => hipfftSetWorkArea_, &
                              cufftSetStream   => hipfftSetStream_
#endif
    use mod_common_cudecomp, only: cuda_stream_kind
    !
    ! to be done after initializing all FFTs and allocating work
    !
    implicit none
#if !defined(_USE_HIP)
    integer    , intent(in) :: arrplan(:)
#else
    type(C_PTR), intent(in) :: arrplan(:)
#endif
    integer(cuda_stream_kind), target, intent(in), optional :: istream
    integer :: istat,i
    do i=1,size(arrplan)
      !$acc   host_data use_device(     work)
      !$omp target data use_device_addr(work)
#if !defined(_USE_HIP)
      istat = cufftSetWorkArea(arrplan(i),work)
      if(present(istream)) then
        istat = cufftSetStream(arrplan(i),istream)
      end if
#else
      istat = cufftSetWorkArea(arrplan(i),c_loc(work))
      !if(present(istream)) then
      !  istat = cufftSetStream(arrplan(i),c_loc(istream))
      !end if
#endif
      !$omp end target data
      !$acc end   host_data
    end do
  end subroutine set_cufft_wspace
  subroutine cudecomp_finalize
#if !defined(_USE_DIEZDECOMP)
    use cudecomp
    use mod_common_cudecomp, only: handle,gd_halo,gd_poi,gd_poi_io,gd_ptdma
    use mod_param          , only: is_poisson_pcr_tdma
    !
    implicit none
    integer :: istat
    istat = cudecompGridDescDestroy(handle,gd_halo)
    istat = cudecompGridDescDestroy(handle,gd_poi)
    istat = cudecompGridDescDestroy(handle,gd_poi_io)
    if(is_poisson_pcr_tdma) istat = cudecompGridDescDestroy(handle,gd_ptdma)
    istat = cudecompFinalize(handle)
#endif
  end subroutine cudecomp_finalize
#endif
end module mod_workspaces
