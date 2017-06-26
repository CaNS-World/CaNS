!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2012 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This file contains the routines that transpose data from X to Z pencil

  subroutine transpose_x_to_z_real(src, dst, opt_decomp)!complex(src, dst, opt_decomp)

    implicit none
    
    !complex(mytype), dimension(:,:,:), intent(IN) :: src
    !complex(mytype), dimension(:,:,:), intent(OUT) :: dst
    real(mytype), dimension(:,:,:), intent(IN) :: src
    real(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

#if !defined(MPI3)
    call MPI_Alltoallw(src,decomp%xcnts_xz,decomp%xdispls_xz,decomp%xtypes_xz, &
      dst,decomp%zcnts_xz,decomp%zdispls_xz,decomp%ztypes_xz,MPI_COMM_WORLD,ierror)
#endif

#ifdef MPI3
    call MPI_Neighbor_alltoallw( &
      src,decomp%xcnts_xz(decomp%xranks),decomp%xdispls_xz(decomp%xranks),decomp%xtypes_xz(decomp%xranks), &
      dst,decomp%zcnts_xz(decomp%zranks),decomp%zdispls_xz(decomp%zranks),decomp%ztypes_xz(decomp%zranks), &
      decomp%xtozNeighborComm,ierror)
#endif

    return
  end subroutine transpose_x_to_z_real!complex

