!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2011 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This module provides parallel IO facilities for applications based on
! 2D decomposition.

module decomp_2d_io

  use decomp_2d
  use MPI
#ifdef T3PIO
  use t3pio
#endif

  implicit none

  private        ! Make everything private unless declared public

  public :: decomp_2d_write_one, decomp_2d_read_one, &
       decomp_2d_write_var, decomp_2d_read_var, &
       decomp_2d_write_scalar, decomp_2d_read_scalar, &
       decomp_2d_write_plane, decomp_2d_write_every

  ! Generic interface to handle multiple data types and decompositions

  interface decomp_2d_write_one
     module procedure write_one_real
     module procedure write_one_complex
  end interface decomp_2d_write_one

  interface decomp_2d_read_one
     module procedure read_one_real
     module procedure read_one_complex
  end interface decomp_2d_read_one

  interface decomp_2d_write_var
     module procedure write_var_real
     module procedure write_var_complex
  end interface decomp_2d_write_var

  interface decomp_2d_read_var
     module procedure read_var_real
     module procedure read_var_complex
  end interface decomp_2d_read_var

  interface decomp_2d_write_scalar
     module procedure write_scalar_real
     module procedure write_scalar_complex
     module procedure write_scalar_integer
  end interface decomp_2d_write_scalar

  interface decomp_2d_read_scalar
     module procedure read_scalar_real
     module procedure read_scalar_complex
     module procedure read_scalar_integer
  end interface decomp_2d_read_scalar

  interface decomp_2d_write_plane
     module procedure write_plane_3d_real
     module procedure write_plane_3d_complex
!     module procedure write_plane_2d
  end interface decomp_2d_write_plane

  interface decomp_2d_write_every
     module procedure write_every_real
     module procedure write_every_complex
  end interface decomp_2d_write_every
  
contains
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Using MPI-IO library to write a single 3D array to a file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_one_real(ipencil,var,filename,opt_decomp)
    
    implicit none
    
    integer, intent(IN) :: ipencil
    real(mytype), dimension(:,:,:), intent(IN) :: var
    character(len=*), intent(IN) :: filename
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, fh, data_type, info, gs

    data_type = real_type

#include "io_write_one.f90"

    return
  end subroutine write_one_real


  subroutine write_one_complex(ipencil,var,filename,opt_decomp)
    
    implicit none
    
    integer, intent(IN) :: ipencil
    complex(mytype), dimension(:,:,:), intent(IN) :: var
    character(len=*), intent(IN) :: filename
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, fh, data_type, info, gs

    data_type = complex_type

#include "io_write_one.f90"
    
    return
  end subroutine write_one_complex


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Using MPI-IO library to read from a file a single 3D array
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_one_real(ipencil,var,filename,opt_decomp)
    
    implicit none
    
    integer, intent(IN) :: ipencil
    real(mytype), dimension(:,:,:), intent(INOUT) :: var
    character(len=*), intent(IN) :: filename
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer(kind=MPI_OFFSET_KIND) :: disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, fh, data_type

    data_type = real_type
    
#include "io_read_one.f90"
    
    return
  end subroutine read_one_real


  subroutine read_one_complex(ipencil,var,filename,opt_decomp)
    
    implicit none
    
    integer, intent(IN) :: ipencil
    complex(mytype), dimension(:,:,:), intent(INOUT) :: var
    character(len=*), intent(IN) :: filename
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer(kind=MPI_OFFSET_KIND) :: disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, fh, data_type
    
    data_type = complex_type

#include "io_read_one.f90"
    
    return
  end subroutine read_one_complex


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write a 3D array as part of a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated after the writing
  !  operation to prepare the writing of next chunk of data.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_var_real(fh,disp,ipencil,var,opt_decomp)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: ipencil
    real(mytype), dimension(:,:,:), intent(IN) :: var
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, data_type

    data_type = real_type

#include "io_write_var.f90"

    return
  end subroutine write_var_real


  subroutine write_var_complex(fh,disp,ipencil,var,opt_decomp)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: ipencil
    complex(mytype), dimension(:,:,:), intent(IN) :: var
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, data_type

    data_type = complex_type

#include "io_write_var.f90"

    return
  end subroutine write_var_complex


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Read a 3D array as part of a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated after the reading
  !  operation to prepare the reading of next chunk of data.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_var_real(fh,disp,ipencil,var,opt_decomp)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: ipencil
    real(mytype), dimension(:,:,:), intent(INOUT) :: var
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, data_type

    data_type = real_type

#include "io_read_var.f90"

    return
  end subroutine read_var_real


  subroutine read_var_complex(fh,disp,ipencil,var,opt_decomp)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: ipencil
    complex(mytype), dimension(:,:,:), intent(INOUT) :: var
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, data_type

    data_type = complex_type

#include "io_read_var.f90"

    return
  end subroutine read_var_complex

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write scalar variables as part of a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated after the reading
  !  operation to prepare the reading of next chunk of data.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_scalar_real(fh,disp,n,var)

    implicit none

    integer, intent(IN) :: fh             ! file handle
    integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement
    integer, intent(IN) :: n              ! number of scalars
    real(mytype), dimension(n), &
         intent(IN) :: var                ! array of scalars

    integer :: m, ierror

    call MPI_FILE_SET_VIEW(fh,disp,real_type, &
         real_type,'native',MPI_INFO_NULL,ierror)
    if (nrank==0) then
       m = n ! only one rank needs to write
    else
       m = 0
    end if
    call MPI_FILE_WRITE_ALL(fh, var, m, real_type, &
         MPI_STATUS_IGNORE, ierror)
    disp = disp + n*mytype_bytes

    return
  end subroutine write_scalar_real


  subroutine write_scalar_complex(fh,disp,n,var)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: n
    complex(mytype), dimension(n), intent(IN) :: var

    integer :: m, ierror

    call MPI_FILE_SET_VIEW(fh,disp,complex_type, &
         complex_type,'native',MPI_INFO_NULL,ierror)
    if (nrank==0) then
       m = n
    else
       m = 0
    end if
    call MPI_FILE_WRITE_ALL(fh, var, m, complex_type, &
         MPI_STATUS_IGNORE, ierror)
    disp = disp + n*mytype_bytes*2

    return
  end subroutine write_scalar_complex


  subroutine write_scalar_integer(fh,disp,n,var)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: n
    integer, dimension(n), intent(IN) :: var

    integer :: m, ierror

    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         MPI_INTEGER,'native',MPI_INFO_NULL,ierror)
    if (nrank==0) then
       m = n
    else
       m = 0
    end if
    call MPI_FILE_WRITE_ALL(fh, var, m, MPI_INTEGER, &
         MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_SIZE(MPI_INTEGER,m,ierror)
    disp = disp + n*m

    return
  end subroutine write_scalar_integer


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Read scalar variables as part of a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated after the reading
  !  operation to prepare the reading of next chunk of data.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_scalar_real(fh,disp,n,var)

    implicit none

    integer, intent(IN) :: fh             ! file handle
    integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement
    integer, intent(IN) :: n              ! number of scalars
    real(mytype), dimension(n), &
         intent(INOUT) :: var             ! array of scalars

    integer :: ierror

    call MPI_FILE_SET_VIEW(fh,disp,real_type, &
         real_type,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh, var, n, real_type, &
         MPI_STATUS_IGNORE, ierror)
    disp = disp + n*mytype_bytes

    return
  end subroutine read_scalar_real
  

  subroutine read_scalar_complex(fh,disp,n,var)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: n
    complex(mytype), dimension(n), intent(INOUT) :: var

    integer :: ierror

    call MPI_FILE_SET_VIEW(fh,disp,complex_type, &
         complex_type,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh, var, n, complex_type, &
         MPI_STATUS_IGNORE, ierror)
    disp = disp + n*mytype_bytes*2

    return
  end subroutine read_scalar_complex


  subroutine read_scalar_integer(fh,disp,n,var)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: n
    integer, dimension(n), intent(INOUT) :: var

    integer :: m, ierror

    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         MPI_INTEGER,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh, var, n, MPI_INTEGER, &
         MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_SIZE(MPI_INTEGER,m,ierror)
    disp = disp + n*m

    return
  end subroutine read_scalar_integer


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write a 2D slice of the 3D data to a file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_plane_3d_real(ipencil,var,iplane,n,filename, &
       opt_decomp)

    implicit none

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    real(mytype), dimension(:,:,:), intent(IN) :: var
    integer, intent(IN) :: iplane !(x-plane=1; y-plane=2; z-plane=3)
    integer, intent(IN) :: n ! which plane to write (global coordinate)
    character(len=*), intent(IN) :: filename
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    real(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    real(mytype), allocatable, dimension(:,:,:) :: wk2d
    TYPE(DECOMP_INFO) :: decomp
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh, data_type

    data_type = real_type

#include "io_write_plane.f90"

    return
  end subroutine write_plane_3d_real


  subroutine write_plane_3d_complex(ipencil,var,iplane,n, &
       filename,opt_decomp)

    implicit none

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    complex(mytype), dimension(:,:,:), intent(IN) :: var
    integer, intent(IN) :: iplane !(x-plane=1; y-plane=2; z-plane=3)
    integer, intent(IN) :: n ! which plane to write (global coordinate)
    character(len=*), intent(IN) :: filename
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    complex(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    complex(mytype), allocatable, dimension(:,:,:) :: wk2d
    TYPE(DECOMP_INFO) :: decomp
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh, data_type

    data_type = complex_type

#include "io_write_plane.f90"

    return
  end subroutine write_plane_3d_complex


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write a 2D array to a file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !************** TO DO ***************
!  subroutine write_plane_2d(ipencil,var,filename)
!    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
!    real(mytype), dimension(:,:), intent(IN) :: var ! 2D array
!    character(len=*), intent(IN) :: filename
!
!    if (ipencil==1) then
!       ! var should be defined as var(xsize(2)
!
!    else if (ipencil==2) then
!
!    else if (ipencil==3) then
!
!    end if
!
!    return
!  end subroutine write_plane_2d


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write 3D array data for every specified mesh point
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_every_real(ipencil,var,iskip,jskip,kskip, &
       filename, from1)

    implicit none

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    real(mytype), dimension(:,:,:), intent(IN) :: var
    integer, intent(IN) :: iskip,jskip,kskip 
    character(len=*), intent(IN) :: filename
    logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
                                  ! .false. - save n,2n,3n...

    real(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh, key,color,newcomm, data_type
    integer, dimension(3) :: xsz,ysz,zsz,xst,yst,zst,xen,yen,zen,skip

    data_type = real_type

#include "io_write_every.f90"
    
    return
  end subroutine write_every_real


  subroutine write_every_complex(ipencil,var,iskip,jskip,kskip, &
       filename, from1)

    implicit none

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    complex(mytype), dimension(:,:,:), intent(IN) :: var
    integer, intent(IN) :: iskip,jskip,kskip 
    character(len=*), intent(IN) :: filename
    logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
                                  ! .false. - save n,2n,3n...

    complex(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh, key,color,newcomm, data_type
    integer, dimension(3) :: xsz,ysz,zsz,xst,yst,zst,xen,yen,zen,skip

    data_type = complex_type

#include "io_write_every.f90"
    
    return
  end subroutine write_every_complex


end module decomp_2d_io
