! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_param
use mod_types
#if defined(_OPENACC)
use cudecomp
#endif
implicit none
public
!
! parameters
!
real(rp), parameter :: pi = acos(-1._rp)
#if !defined(_EPS_EXACT_ZERO) /* recommended */
real(rp), parameter :: eps = epsilon(1._rp)
#else
real(rp), parameter :: eps = 0._rp
#endif
real(rp), parameter :: small = epsilon(1._rp)*10**(precision(1._rp)/2)
character(len=100), parameter :: datadir = 'data/'
real(rp), parameter, dimension(2,3) :: rkcoeff = reshape([32._rp/60._rp,  0._rp        , &
                                                          25._rp/60._rp, -17._rp/60._rp, &
                                                          45._rp/60._rp, -25._rp/60._rp], shape(rkcoeff))
real(rp), parameter, dimension(3)   :: rkcoeff12 = rkcoeff(1,:)+rkcoeff(2,:)
!
! variables to be determined from the input file
!
integer , protected, dimension(3) :: ng
real(rp), protected, dimension(3) :: l
integer , protected :: gtype
real(rp), protected :: gr
real(rp), protected :: cfl,dtmax,dt_f
real(rp), protected :: visci
!
character(len=100), protected :: inivel
logical, protected :: is_wallturb
!
integer , protected :: nstep
real(rp), protected :: time_max,tw_max
logical , protected, dimension(3) :: stop_type
logical , protected :: restart,is_overwrite_save
integer , protected :: nsaves_max
integer , protected :: icheck,iout0d,iout1d,iout2d,iout3d,isave
!
integer , dimension(2) :: dims
!
integer, dimension(0:1,3) :: nb
logical, dimension(0:1,3) :: is_bound
character(len=1), protected, dimension(0:1,3,3) ::  cbcvel
real(rp)        , protected, dimension(0:1,3,3) ::   bcvel
character(len=1), protected, dimension(0:1,3)   ::  cbcpre
real(rp)        , protected, dimension(0:1,3)   ::   bcpre
!
real(rp), protected, dimension(3) :: bforce
logical , protected, dimension(3) :: is_forced
real(rp), protected, dimension(3) :: velf
!
real(rp), protected, dimension(3) :: dl,dli
real(rp), protected :: visc
!
#if defined(_OPENACC)
!
! cuDecomp input parameters
!
integer, protected :: cudecomp_t_comm_backend,cudecomp_h_comm_backend
logical, protected :: cudecomp_is_t_comm_autotune ,cudecomp_is_h_comm_autotune , &
                      cudecomp_is_t_enable_nccl   ,cudecomp_is_h_enable_nccl   , &
                      cudecomp_is_t_enable_nvshmem,cudecomp_is_h_enable_nvshmem, &
                      cudecomp_is_t_in_place
#endif
contains
  subroutine read_input(myid)
    use mpi
    implicit none
    character(len=*), parameter :: input_file = 'input.nml'
    integer, intent(in) :: myid
    integer :: iunit,ierr
    character(len=1024) :: c_iomsg
    namelist /dns/ &
                  ng, &
                  l, &
                  gtype,gr, &
                  cfl,dtmax,dt_f, &
                  visci, &
                  inivel, &
                  is_wallturb, &
                  nstep,time_max,tw_max, &
                  stop_type, &
                  restart,is_overwrite_save,nsaves_max, &
                  icheck,iout0d,iout1d,iout2d,iout3d,isave, &
                  cbcvel,cbcpre,bcvel,bcpre, &
                  bforce, &
                  is_forced, &
                  velf, &
                  dims
#if defined(_OPENACC)
    namelist /cudecomp/ &
                       cudecomp_t_comm_backend,cudecomp_is_t_enable_nccl,cudecomp_is_t_enable_nvshmem, &
                       cudecomp_h_comm_backend,cudecomp_is_h_enable_nccl,cudecomp_is_h_enable_nvshmem
#endif
    !
    ! defaults
    !
    dt_f = -1.
    open(newunit=iunit,file=input_file,status='old',action='read',iostat=ierr,iomsg=c_iomsg)
      if(ierr /= 0) then
        if(myid == 0) print*, 'Error reading the input file: ', trim(c_iomsg)
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        close(iunit)
        error stop
      end if
      read(iunit,nml=dns,iostat=ierr,iomsg=c_iomsg)
      if(ierr /= 0) then
        if(myid == 0) print*, 'Error reading dns namelist: ', trim(c_iomsg)
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        close(iunit)
        error stop
      end if
      !
      dl(:) = l(:)/(1.*ng(:))
      dli(:) = dl(:)**(-1)
      visc = visci**(-1)
#if defined(_OPENACC)
      !
      ! reading cuDecomp parameters, if these are set
      !
      ! defaults
      !
      cudecomp_is_t_comm_autotune  = .true.
      cudecomp_is_h_comm_autotune  = .true.
      cudecomp_is_t_enable_nccl    = .true.
      cudecomp_is_h_enable_nccl    = .true.
      cudecomp_is_t_enable_nvshmem = .true.
      cudecomp_is_h_enable_nvshmem = .true.
      rewind(iunit)
      read(iunit,nml=cudecomp,iostat=ierr,iomsg=c_iomsg)
      if(ierr /= 0) then
        if(myid == 0) print*, 'Error reading cudecomp namelist: ', trim(c_iomsg)
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        close(iunit)
        error stop
      end if
      !
      if(cudecomp_t_comm_backend >= 1 .and. cudecomp_t_comm_backend <= 7) then
        cudecomp_is_t_comm_autotune = .false. ! do not autotune if backend is prescribed
        select case(cudecomp_t_comm_backend)
        case(1)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P
        case(2)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P_PL
        case(3)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_A2A
        case(4)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_NCCL
        case(5)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_NCCL_PL
        case(6)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_NVSHMEM
        case(7)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_NVSHMEM_PL
        case default
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P
        end select
      end if
      if(cudecomp_h_comm_backend >= 1 .and. cudecomp_h_comm_backend <= 4) then
        cudecomp_is_h_comm_autotune = .false. ! do not autotune if backend is prescribed
        select case(cudecomp_h_comm_backend)
        case(1)
          cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_MPI
        case(2)
          cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_MPI_BLOCKING
        case(3)
          cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_NCCL
        case(4)
          cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_NVSHMEM
        case(5)
          cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_NVSHMEM_BLOCKING
        case default
          cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_MPI
        end select
      end if
      !
      ! manually set cuDecomp out-of-place transposes by default
      !
      cudecomp_is_t_in_place = .false.
#endif
    close(iunit)
  end subroutine read_input
end module mod_param
