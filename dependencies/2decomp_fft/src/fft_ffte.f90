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

! This is the FFTE implementation of the FFT library
! Note that FFTE only support transform size with prime factor 2/3/5.
! Other sizes are computed by the generic engine.

module decomp_2d_fft
  
  use decomp_2d  ! 2D decomposition module
  use glassman
  
  implicit none
  
  private        ! Make everything private unless declared public

  ! engine-specific global variables
  complex(mytype), allocatable, dimension(:) :: buf, scratch

  logical, save :: x_goodsize, y_goodsize, z_goodsize

  ! common code used for all engines, including global variables, 
  ! generic interface definitions and several subroutines
#include "fft_common.f90"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time initialisations for the FFT engine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_fft_engine

    implicit none

    integer :: cbuf_size, errorcode

    if (KIND(0.0) == mytype) then
       errorcode = 5
       call decomp_2d_abort(errorcode, &
            'FFTE engine only suppport double-precision calculations')
    end if

    if (nrank==0) then
       write(*,*) ' '
       write(*,*) '***** Using the FFTE engine *****'
       write(*,*) ' '
    end if

    cbuf_size = max(ph%xsz(1), ph%ysz(2))
    cbuf_size = max(cbuf_size, ph%zsz(3))
    allocate(buf(cbuf_size))
    allocate(scratch(2*cbuf_size))  ! twice space for vffte library

    ! check the transform sizes
    x_goodsize = check_size(ph%xsz(1))
    y_goodsize = check_size(ph%ysz(2))
    z_goodsize = check_size(ph%zsz(3))

    return
  end subroutine init_fft_engine

  ! factorisation routine
#include "factor.f90"

  logical function check_size(nsize)

    integer, intent(IN) :: nsize

    integer, allocatable, dimension(:) :: factors
    integer :: nfact, i

    i = ceiling(log(real(nsize))*1.442695)  ! log_2_N
    allocate(factors(i))

    call primefactors(nsize,factors,nfact)

    check_size = .true.

    do i=1,nfact
       if (factors(i) > 5) then
          check_size = .false.
          exit
       end if
    end do

    return
  end function check_size


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time finalisations for the FFT engine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine finalize_fft_engine

    implicit none

    deallocate(buf,scratch)

    return
  end subroutine finalize_fft_engine


  ! Following routines calculate multiple one-dimensional FFTs to form 
  ! the basis of three-dimensional FFTs.

  ! c2c transform, multiple 1D FFTs in x direction
  subroutine c2c_1m_x(inout, isign, decomp)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: i,j,k
    
    if (x_goodsize) then
       call zfft1d(buf,decomp%xsz(1),0,scratch)
    end if

    do k=1,decomp%xsz(3)
       do j=1,decomp%xsz(2)
          do i=1,decomp%xsz(1)
             buf(i) = inout(i,j,k)
          end do
          if (x_goodsize) then
             call zfft1d(buf,decomp%xsz(1),isign,scratch)
          else
             call spcfft(buf,decomp%xsz(1),isign,scratch)
          end if
          do i=1,decomp%xsz(1)
             inout(i,j,k) = buf(i)
          end do
       end do
    end do

    ! Note FFTE backward transform is scaled, unscale here
    if (x_goodsize .AND. isign==DECOMP_2D_FFT_BACKWARD) then
       inout = inout * real(decomp%xsz(1), kind=mytype)
    end if

    return

  end subroutine c2c_1m_x

  ! c2c transform, multiple 1D FFTs in y direction
  subroutine c2c_1m_y(inout, isign, decomp)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: i,j,k

    if (y_goodsize) then
       call zfft1d(buf,decomp%ysz(2),0,scratch)
    end if

    do k=1,decomp%ysz(3)
       do i=1,decomp%ysz(1)
          do j=1,decomp%ysz(2)
             buf(j) = inout(i,j,k)
          end do
          if (y_goodsize) then
             call zfft1d(buf,decomp%ysz(2),isign,scratch)
          else
             call spcfft(buf,decomp%ysz(2),isign,scratch)
          end if
          do j=1,decomp%ysz(2)
             inout(i,j,k) = buf(j)
          end do
       end do
    end do

    if (y_goodsize .AND. isign==DECOMP_2D_FFT_BACKWARD) then
       inout = inout * real(decomp%ysz(2), kind=mytype)
    end if

    return

  end subroutine c2c_1m_y

  ! c2c transform, multiple 1D FFTs in z direction
  subroutine c2c_1m_z(inout, isign, decomp)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: i,j,k

    if (z_goodsize) then
       call zfft1d(buf,decomp%zsz(3),0,scratch)
    end if

    do j=1,decomp%zsz(2)
       do i=1,decomp%zsz(1)
          do k=1,decomp%zsz(3)
             buf(k) = inout(i,j,k)
          end do
          if (z_goodsize) then
             call zfft1d(buf,decomp%zsz(3),isign,scratch)
          else
             call spcfft(buf,decomp%zsz(3),isign,scratch)
          end if
          do k=1,decomp%zsz(3)
             inout(i,j,k) = buf(k)
          end do
       end do
    end do

    if (z_goodsize .AND. isign==DECOMP_2D_FFT_BACKWARD) then
       inout = inout * real(decomp%zsz(3), kind=mytype)
    end if

    return

  end subroutine c2c_1m_z

  ! r2c transform, multiple 1D FFTs in x direction
  subroutine r2c_1m_x(input, output)

    implicit none

    real(mytype), dimension(:,:,:), intent(IN)  ::  input
    complex(mytype), dimension(:,:,:), intent(OUT) :: output

    integer :: i,j,k, s1,s2,s3, d1

    s1 = size(input,1)
    s2 = size(input,2)
    s3 = size(input,3)
    d1 = size(output,1)

    if (x_goodsize) then
       call zfft1d(buf,s1,0,scratch)
    end if

    do k=1,s3
       do j=1,s2
          ! The base FFT engine is c2c only, 
          ! needing some pre- and post-processing for r2c
          ! pack real input in complex storage
          do i=1,s1
             buf(i) = cmplx(input(i,j,k),0._mytype, kind=mytype)
          end do
          if (x_goodsize) then
             call zfft1d(buf,s1,-1,scratch)
          else
             call spcfft(buf,s1,-1,scratch)
          end if
          ! note d1 ~ s1/2+1
          ! simply drop the redundant part of the complex output
          do i=1,d1
             output(i,j,k) = buf(i)
          end do
       end do
    end do

    return

  end subroutine r2c_1m_x

  ! r2c transform, multiple 1D FFTs in z direction
  subroutine r2c_1m_z(input, output)

    implicit none

    real(mytype), dimension(:,:,:), intent(IN)  ::  input
    complex(mytype), dimension(:,:,:), intent(OUT) :: output

    integer :: i,j,k, s1,s2,s3, d3

    s1 = size(input,1)
    s2 = size(input,2)
    s3 = size(input,3)
    d3 = size(output,3)

    if (z_goodsize) then
       call zfft1d(buf,s3,0,scratch)
    end if

    do j=1,s2
       do i=1,s1
          do k=1,s3
             buf(k) = cmplx(input(i,j,k),0._mytype, kind=mytype)
          end do
          if (z_goodsize) then
             call zfft1d(buf,s3,-1,scratch)
          else
             call spcfft(buf,s3,-1,scratch)
          end if
          do k=1,d3
             output(i,j,k) = buf(k)
          end do
       end do
    end do

    return

  end subroutine r2c_1m_z

  ! c2r transform, multiple 1D FFTs in x direction
  subroutine c2r_1m_x(input, output)

    implicit none

    complex(mytype), dimension(:,:,:), intent(IN)  ::  input
    real(mytype), dimension(:,:,:), intent(OUT) :: output

    integer :: i,j,k, d1,d2,d3

    d1 = size(output,1)
    d2 = size(output,2)
    d3 = size(output,3)

    if (x_goodsize) then
       call zfft1d(buf,d1,0,scratch)
    end if

    do k=1,d3
       do j=1,d2
          ! The base FFT engine is c2c only, 
          ! needing some pre- and post-processing for c2r
          do i=1,d1/2+1
             buf(i) = input(i,j,k)
          end do
          ! expanding to a full-size complex array
          ! For odd N, the storage is:
          !  1, 2, ...... N/2+1   integer division rounded down
          !     N, ...... N/2+2   => a(i) is conjugate of a(N+2-i)
          ! For even N, the storage is:
          !  1, 2, ...... N/2  , N/2+1
          !     N, ...... N/2+2  again a(i) conjugate of a(N+2-i)
          do i=d1/2+2,d1
             buf(i) =  conjg(buf(d1+2-i))
          end do
          if (x_goodsize) then
             call zfft1d(buf,d1,1,scratch)
          else
             call spcfft(buf,d1,1,scratch)
          end if
          do i=1,d1
             ! simply drop imaginary part
             output(i,j,k) = real(buf(i), kind=mytype)
          end do
       end do
    end do

    if (x_goodsize) then
       output = output * real(d1, kind=mytype)
    end if

    return

  end subroutine c2r_1m_x

  ! c2r transform, multiple 1D FFTs in z direction
  subroutine c2r_1m_z(input, output)

    implicit none

    complex(mytype), dimension(:,:,:), intent(IN)  ::  input
    real(mytype), dimension(:,:,:), intent(OUT) :: output

    integer :: i,j,k, d1,d2,d3

    d1 = size(output,1)
    d2 = size(output,2)
    d3 = size(output,3)

    if (z_goodsize) then
       call zfft1d(buf,d3,0,scratch)
    end if

    do j=1,d2
       do i=1,d1
          do k=1,d3/2+1
             buf(k) = input(i,j,k)
          end do
          do k=d3/2+2,d3
             buf(k) =  conjg(buf(d3+2-k))
          end do
          if (z_goodsize) then
             call zfft1d(buf,d3,1,scratch)
          else
             call spcfft(buf,d3,1,scratch)
          end if
          do k=1,d3
             output(i,j,k) = real(buf(k), kind=mytype)
          end do
       end do
    end do

    if (z_goodsize) then
       output = output * real(d3, kind=mytype)
    end if

    return

  end subroutine c2r_1m_z


#include "fft_common_3d.f90"

  
end module decomp_2d_fft
