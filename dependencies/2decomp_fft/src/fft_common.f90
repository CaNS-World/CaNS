!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2021 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This file contains common code shared by all FFT engines

  integer, parameter, public :: DECOMP_2D_FFT_FORWARD = -1
  integer, parameter, public :: DECOMP_2D_FFT_BACKWARD = 1
  
  ! Physical space data can be stored in either X-pencil or Z-pencil
  integer, parameter, public :: PHYSICAL_IN_X = 1
  integer, parameter, public :: PHYSICAL_IN_Z = 3 

  integer, save :: format                 ! input X-pencil or Z-pencil
  
  ! The libary can only be initialised once
  logical, save :: initialised = .false. 

  ! Global size of the FFT
  integer, save :: nx_fft, ny_fft, nz_fft

  ! 2D processor grid
  integer, save, dimension(2) :: dims

  ! Decomposition objects
  TYPE(DECOMP_INFO), save :: ph  ! physical space
  TYPE(DECOMP_INFO), save :: sp  ! spectral space

  ! Workspace to store the intermediate Y-pencil data
  ! *** TODO: investigate how to use only one workspace array
  complex(mytype), allocatable, dimension(:,:,:) :: wk2_c2c, wk2_r2c
  complex(mytype), allocatable, dimension(:,:,:) :: wk13

  public :: decomp_2d_fft_init, decomp_2d_fft_3d, &
       decomp_2d_fft_finalize, decomp_2d_fft_get_size
  
  ! Declare generic interfaces to handle different inputs
  
  interface decomp_2d_fft_init
     module procedure fft_init_noarg
     module procedure fft_init_arg
     module procedure fft_init_general
  end interface
  
  interface decomp_2d_fft_3d
     module procedure fft_3d_c2c
     module procedure fft_3d_r2c
     module procedure fft_3d_c2r
  end interface

  
contains
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialise the FFT module
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fft_init_noarg
    
    implicit none
    
    call fft_init_arg(PHYSICAL_IN_X)  ! default input is X-pencil data
    
    return
  end subroutine fft_init_noarg

  subroutine fft_init_arg(pencil)     ! allow to handle Z-pencil input

    implicit none

    integer, intent(IN) :: pencil

    call fft_init_general(pencil, nx_global, ny_global, nz_global)

    return
  end subroutine fft_init_arg

  ! Initialise the FFT library to perform arbitrary size transforms
  subroutine fft_init_general(pencil, nx, ny, nz)

    implicit none

    integer, intent(IN) :: pencil
    integer, intent(IN) :: nx, ny, nz

    logical, dimension(2) :: dummy_periods
    integer, dimension(2) :: dummy_coords
    integer :: status, errorcode, ierror

    if (initialised) then
       errorcode = 4
       call decomp_2d_abort(errorcode, &
            'FFT library should only be initialised once')
    end if
    
    format = pencil
    nx_fft = nx
    ny_fft = ny
    nz_fft = nz

    ! determine the processor grid in use
    call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, &
         dims, dummy_periods, dummy_coords, ierror)

    ! for c2r/r2c interface:
    ! if in physical space, a real array is of size: nx*ny*nz
    ! in spectral space, the complex array is of size:
    !         (nx/2+1)*ny*nz, if PHYSICAL_IN_X
    !      or nx*ny*(nz/2+1), if PHYSICAL_IN_Z

    call decomp_info_init(nx, ny, nz, ph)
    if (format==PHYSICAL_IN_X) then
       call decomp_info_init(nx/2+1, ny, nz, sp)
    else if (format==PHYSICAL_IN_Z) then
       call decomp_info_init(nx, ny, nz/2+1, sp)
    end if

    allocate(wk2_c2c(ph%ysz(1),ph%ysz(2),ph%ysz(3)), STAT=status)
    allocate(wk2_r2c(sp%ysz(1),sp%ysz(2),sp%ysz(3)), STAT=status)
    if (format==PHYSICAL_IN_X) then
       allocate(wk13(sp%xsz(1),sp%xsz(2),sp%xsz(3)), STAT=status)
    else if (format==PHYSICAL_IN_Z) then
       allocate(wk13(sp%zsz(1),sp%zsz(2),sp%zsz(3)), STAT=status)
    end if
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising FFT')
    end if

    call init_fft_engine
    
    initialised = .true.
    
    return
  end subroutine fft_init_general

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Final clean up
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_fft_finalize
    
    implicit none

    call decomp_info_finalize(ph)
    call decomp_info_finalize(sp)

    deallocate(wk2_c2c, wk2_r2c, wk13)

    call finalize_fft_engine

    initialised = .false.

    return
  end subroutine decomp_2d_fft_finalize


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Return the size, starting/ending index of the distributed array 
  !  whose global size is (nx/2+1)*ny*nz, for defining data structures
  !  in r2c and c2r interfaces
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_fft_get_size(istart, iend, isize)
    
    implicit none
    integer, dimension(3), intent(OUT) :: istart, iend, isize
    
    if (format==PHYSICAL_IN_X) then
       istart = sp%zst
       iend   = sp%zen
       isize  = sp%zsz
    else if (format==PHYSICAL_IN_Z) then
       istart = sp%xst
       iend   = sp%xen
       isize  = sp%xsz
    end if
    
    return
  end subroutine decomp_2d_fft_get_size
