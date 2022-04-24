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

! This is the FFTPACK5 implementation of the FFT library

! Note this engine code has been modified to suit the packaging of
! FFTPACK 5.1 where single- and double-precision routines are pack
! in separate files but using the same symbols. The preprocessing
! branches for precisions remain, in case the packaging is changed
! to the old way, that rename double-precision routines (c->z; r->d).
!
! Do not use the 5.0 version as it contains serious bugs and gives
! incorrect results for transform sizes with large prime numbers.
! http://people.sc.fsu.edu/~jburkardt/f_src/fftpack5.1/fftpack5.1.html

module decomp_2d_fft
  
  use decomp_2d  ! 2D decomposition module
  
  implicit none
  
  private        ! Make everything private unless declared public

  ! engine-specific global variables

  ! FFTPACK work arrays
  complex(mytype), allocatable, dimension(:) :: wsave, work
  integer, save :: lensav, lenwrk

  ! work arrays to save intermediate 1D data sets
  complex(mytype), allocatable, dimension(:) :: buf_c
  real   (mytype), allocatable, dimension(:) :: buf_r

  ! common code used for all engines, including global variables, 
  ! generic interface definitions and several subroutines
#include "fft_common.f90"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time initialisations for the FFT engine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_fft_engine

    implicit none

    integer :: maxdim

    if (nrank==0) then
       write(*,*) ' '
       write(*,*) '***** Using the FFTPACK5 engine *****'
       write(*,*) ' '
    end if

    maxdim = max(ph%xsz(1),max(ph%ysz(2),ph%zsz(3)))
    lensav = 2*maxdim + int(log(real(maxdim))/log(2.)) + 4
    lenwrk = 2*maxdim

    allocate(wsave(lensav))
    allocate(work(lenwrk))
    allocate(buf_c(maxdim))

    if (format==PHYSICAL_IN_X) then
       allocate(buf_r(ph%xsz(1)))
    else if(format==PHYSICAL_IN_Z) then
       allocate(buf_r(ph%zsz(3)))
    end if

    return
  end subroutine init_fft_engine


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time finalisations for the FFT engine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine finalize_fft_engine

    implicit none

    deallocate(wsave, work, buf_c, buf_r)

    return
  end subroutine finalize_fft_engine


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 3D FFT - complex to complex
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fft_3d_c2c(in, out, isign)
    
    implicit none
    
    complex(mytype), dimension(:,:,:), intent(IN) :: in
    complex(mytype), dimension(:,:,:), intent(OUT) :: out
    integer, intent(IN) :: isign   

    complex(mytype), allocatable, dimension(:,:,:) :: wk1
    integer :: i,j,k, ier
 
    if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_FORWARD .OR.  &
         format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_BACKWARD) then
       
       ! ===== 1D FFTs in X =====
       allocate (wk1(ph%xsz(1),ph%xsz(2),ph%xsz(3)))
#ifdef DOUBLE_PREC
       call cfft1i(ph%xsz(1), wsave, lensav, ier)
#else
       call cfft1i(ph%xsz(1), wsave, lensav, ier)
#endif
       do k=1,ph%xsz(3)
          do j=1,ph%xsz(2)
             do i=1,ph%xsz(1)
                buf_c(i) = in(i,j,k)
             enddo
             if (isign==DECOMP_2D_FFT_FORWARD) then
#ifdef DOUBLE_PREC
                call cfft1f(ph%xsz(1),1,buf_c,ph%xsz(1), &
                     wsave,lensav,work,lenwrk,ier)
#else
                call cfft1f(ph%xsz(1),1,buf_c,ph%xsz(1), &
                     wsave,lensav,work,lenwrk,ier)
#endif
             else
#ifdef DOUBLE_PREC
                call cfft1b(ph%xsz(1),1,buf_c,ph%xsz(1), &
                     wsave,lensav,work,lenwrk,ier)
#else
                call cfft1b(ph%xsz(1),1,buf_c,ph%xsz(1), &
                     wsave,lensav,work,lenwrk,ier)
#endif
             end if
             do i=1,ph%xsz(1)
                wk1(i,j,k) = buf_c(i)
             end do
          end do
       end do

       ! ===== Swap X --> Y =====
       call transpose_x_to_y(wk1,wk2_c2c,ph)
       
       ! ===== 1D FFTs in Y =====
#ifdef DOUBLE_PREC
       call cfft1i(ph%ysz(2), wsave, lensav, ier)
#else
       call cfft1i(ph%ysz(2), wsave, lensav, ier)
#endif
       do k=1,ph%ysz(3)
          do i=1,ph%ysz(1)
             do j=1,ph%ysz(2)
                buf_c(j) = wk2_c2c(i,j,k)
             end do
             if (isign==DECOMP_2D_FFT_FORWARD) then
#ifdef DOUBLE_PREC
                call cfft1f(ph%ysz(2),1,buf_c,ph%ysz(2), &
                     wsave,lensav,work,lenwrk,ier)
#else
                call cfft1f(ph%ysz(2),1,buf_c,ph%ysz(2), &
                     wsave,lensav,work,lenwrk,ier)
#endif
             else
#ifdef DOUBLE_PREC
                call cfft1b(ph%ysz(2),1,buf_c,ph%ysz(2), &
                     wsave,lensav,work,lenwrk,ier)
#else
                call cfft1b(ph%ysz(2),1,buf_c,ph%ysz(2), &
                     wsave,lensav,work,lenwrk,ier)
#endif
             end if
             do j=1,ph%ysz(2)
                wk2_c2c(i,j,k) = buf_c(j)
             end do
          end do
       end do

       ! ===== Swap Y --> Z =====
       call transpose_y_to_z(wk2_c2c,out,ph)

       ! ===== 1D FFTs in Z =====
#ifdef DOUBLE_PREC       
       call cfft1i(ph%zsz(3), wsave, lensav, ier)
#else
       call cfft1i(ph%zsz(3), wsave, lensav, ier)
#endif
       do j=1,ph%zsz(2)
          do i=1,ph%zsz(1)
             do k=1,ph%zsz(3)
                buf_c(k) = out(i,j,k)
             end do
             if (isign==DECOMP_2D_FFT_FORWARD) then
#ifdef DOUBLE_PREC
                call cfft1f(ph%zsz(3),1,buf_c,ph%zsz(3), &
                     wsave,lensav,work,lenwrk,ier)
#else
                call cfft1f(ph%zsz(3),1,buf_c,ph%zsz(3), &
                     wsave,lensav,work,lenwrk,ier)
#endif
             else
#ifdef DOUBLE_PREC
                call cfft1b(ph%zsz(3),1,buf_c,ph%zsz(3), &
                     wsave,lensav,work,lenwrk,ier)
#else
                call cfft1b(ph%zsz(3),1,buf_c,ph%zsz(3), &
                     wsave,lensav,work,lenwrk,ier)
#endif
             end if
             do k=1,ph%zsz(3)
                out(i,j,k) = buf_c(k)
             end do
          end do
       end do

    else if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_BACKWARD &
         .OR. & 
         format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_FORWARD) then

       ! ===== 1D FFTs in Z =====
       allocate (wk1(ph%zsz(1),ph%zsz(2),ph%zsz(3)))
#ifdef DOUBLE_PREC
       call cfft1i(ph%zsz(3), wsave, lensav, ier)
#else
       call cfft1i(ph%zsz(3), wsave, lensav, ier)
#endif
       do j=1,ph%zsz(2)
          do i=1,ph%zsz(1)
             do k=1,ph%zsz(3)
                buf_c(k) = in(i,j,k)
             end do
             if (isign==DECOMP_2D_FFT_FORWARD) then
#ifdef DOUBLE_PREC
                call cfft1f(ph%zsz(3),1,buf_c,ph%zsz(3), &
                     wsave,lensav,work,lenwrk,ier)
#else
                call cfft1f(ph%zsz(3),1,buf_c,ph%zsz(3), &
                     wsave,lensav,work,lenwrk,ier)
#endif
             else
#ifdef DOUBLE_PREC
                call cfft1b(ph%zsz(3),1,buf_c,ph%zsz(3), &
                     wsave,lensav,work,lenwrk,ier)
#else
                call cfft1b(ph%zsz(3),1,buf_c,ph%zsz(3), &
                     wsave,lensav,work,lenwrk,ier)
#endif
             endif
             do k=1,ph%zsz(3)
                wk1(i,j,k) = buf_c(k)
             end do
          end do
       end do
       
       ! ===== Swap Z --> Y =====
       call transpose_z_to_y(wk1,wk2_c2c,ph)
       
       ! ===== 1D FFTs in Y =====
#ifdef DOUBLE_PREC
       call cfft1i(ph%ysz(2), wsave, lensav, ier)
#else
       call cfft1i(ph%ysz(2), wsave, lensav, ier)
#endif
       do k=1,ph%ysz(3)
          do i=1,ph%ysz(1)
             do j=1,ph%ysz(2)
                buf_c(j) = wk2_c2c(i,j,k)
             end do
             if (isign==DECOMP_2D_FFT_FORWARD) then
#ifdef DOUBLE_PREC
                call cfft1f(ph%ysz(2),1,buf_c,ph%ysz(2), &
                     wsave,lensav,work,lenwrk,ier)
#else
                call cfft1f(ph%ysz(2),1,buf_c,ph%ysz(2), &
                     wsave,lensav,work,lenwrk,ier)
#endif
             else
#ifdef DOUBLE_PREC
                call cfft1b(ph%ysz(2),1,buf_c,ph%ysz(2), &
                     wsave,lensav,work,lenwrk,ier)
#else
                call cfft1b(ph%ysz(2),1,buf_c,ph%ysz(2), &
                     wsave,lensav,work,lenwrk,ier)
#endif
             endif
             do j=1,ph%ysz(2)
                wk2_c2c(i,j,k) = buf_c(j)
             end do
          end do
       end do
       
       ! ===== Swap Y --> X =====
       call transpose_y_to_x(wk2_c2c,out,ph)
       
       ! ===== 1D FFTs in X =====
#ifdef DOUBLE_PREC
       call cfft1i(ph%xsz(1), wsave, lensav, ier)
#else
       call cfft1i(ph%xsz(1), wsave, lensav, ier)
#endif
       do k=1,ph%xsz(3)
          do j=1,ph%xsz(2)
             do i=1,ph%xsz(1)
                buf_c(i) = out(i,j,k)
             enddo
             if (isign==DECOMP_2D_FFT_FORWARD) then
#ifdef DOUBLE_PREC
                call cfft1f(ph%xsz(1),1,buf_c,ph%xsz(1), &
                     wsave,lensav,work,lenwrk,ier)
#else
                call cfft1f(ph%xsz(1),1,buf_c,ph%xsz(1), &
                     wsave,lensav,work,lenwrk,ier)
#endif
             else
#ifdef DOUBLE_PREC
                call cfft1b(ph%xsz(1),1,buf_c,ph%xsz(1), &
                     wsave,lensav,work,lenwrk,ier)
#else
                call cfft1b(ph%xsz(1),1,buf_c,ph%xsz(1), &
                     wsave,lensav,work,lenwrk,ier)
#endif
             endif
             do i=1,ph%xsz(1)
                out(i,j,k) = buf_c(i)
             end do
          end do
       end do
       
    end if

    ! Undo the scaling done in the forward FFTs by the FFT engine
    if (isign==DECOMP_2D_FFT_FORWARD) then
       out = out * real(nx_global, kind=mytype) * &
            real(ny_global, kind=mytype) * &
            real(nz_global, kind=mytype)
    end if

    return
  end subroutine fft_3d_c2c

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 3D forward FFT - real to complex
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fft_3d_r2c(in_r, out_c)
    
    implicit none
    
    real(mytype), dimension(:,:,:), intent(IN) :: in_r
    complex(mytype), dimension(:,:,:), intent(OUT) :: out_c
    
    integer :: i,j,k, ier

    if (format==PHYSICAL_IN_X) then

       ! ===== 1D FFTs in X =====
#ifdef DOUBLE_PREC
       call rfft1i(ph%xsz(1),wsave,lensav,ier)
#else
       call rfft1i(ph%xsz(1),wsave,lensav,ier)
#endif
       do k=1,ph%xsz(3)
          do j=1,ph%xsz(2)
             do i=1,ph%xsz(1)
                buf_r(i) = in_r(i,j,k)
             end do
#ifdef DOUBLE_PREC
             call rfft1f(ph%xsz(1),1,buf_r,ph%xsz(1), &
                  wsave,lensav,work,lenwrk,ier)
#else
             call rfft1f(ph%xsz(1),1,buf_r,ph%xsz(1), &
                  wsave,lensav,work,lenwrk,ier)
#endif
             ! convert to complex array of size N/2+1
             call to_complex(buf_r, ph%xsz(1), wk13(:,j,k))
          end do
       end do

       ! ===== Swap X --> Y =====
       call transpose_x_to_y(wk13,wk2_r2c,sp)

       ! ===== 1D FFTs in Y =====
#ifdef DOUBLE_PREC
       call cfft1i(sp%ysz(2), wsave, lensav, ier)
#else
       call cfft1i(sp%ysz(2), wsave, lensav, ier)
#endif
       do k=1,sp%ysz(3)
          do i=1,sp%ysz(1)
             do j=1,sp%ysz(2)
                buf_c(j) = wk2_r2c(i,j,k)
             end do
#ifdef DOUBLE_PREC
             call cfft1f(sp%ysz(2),1,buf_c,sp%ysz(2), &
                  wsave,lensav,work,lenwrk,ier)
#else
             call cfft1f(sp%ysz(2),1,buf_c,sp%ysz(2), &
                  wsave,lensav,work,lenwrk,ier)
#endif
             do j=1,sp%ysz(2)
                ! also remove the scaling here
                wk2_r2c(i,j,k) = buf_c(j) * real(sp%ysz(2), kind=mytype)
             end do
          end do
       end do

       ! ===== Swap Y --> Z =====
       call transpose_y_to_z(wk2_r2c,out_c,sp)

       ! ===== 1D FFTs in Z =====
#ifdef DOUBLE_PREC
       call cfft1i(sp%zsz(3), wsave, lensav, ier)
#else
       call cfft1i(sp%zsz(3), wsave, lensav, ier)
#endif
       do j=1,sp%zsz(2)
          do i=1,sp%zsz(1)
             do k=1,sp%zsz(3)
                buf_c(k) = out_c(i,j,k)
             end do
#ifdef DOUBLE_PREC
             call cfft1f(sp%zsz(3),1,buf_c,sp%zsz(3), &
                  wsave,lensav,work,lenwrk,ier)
#else
             call cfft1f(sp%zsz(3),1,buf_c,sp%zsz(3), &
                  wsave,lensav,work,lenwrk,ier)
#endif
             do k=1,sp%zsz(3)
                ! also remove the scaling here
                out_c(i,j,k) = buf_c(k) * real(sp%zsz(3), kind=mytype)
             end do
          end do
       end do

    else if (format==PHYSICAL_IN_Z) then

       ! ===== 1D FFTs in Z =====
#ifdef DOUBLE_PREC
       call rfft1i(ph%zsz(3),wsave,lensav,ier)
#else
       call rfft1i(ph%zsz(3),wsave,lensav,ier)
#endif
       do j=1,ph%zsz(2)
          do i=1,ph%zsz(1)
             do k=1,ph%zsz(3)
                buf_r(k) = in_r(i,j,k)
             end do
#ifdef DOUBLE_PREC
             call rfft1f(ph%zsz(3),1,buf_r,ph%zsz(3), &
                  wsave,lensav,work,lenwrk,ier)
#else
             call rfft1f(ph%zsz(3),1,buf_r,ph%zsz(3), &
                  wsave,lensav,work,lenwrk,ier)
#endif
             ! convert to complex array of size N/2+1
             ! note the use of temp array here, 
             ! cannot write directly into wk13 because of stride
             call to_complex(buf_r, ph%zsz(3), buf_c)
             do k=1,sp%zsz(3) 
                wk13(i,j,k) = buf_c(k)
             end do
          end do
       end do

       ! ===== Swap Z --> Y =====
       call transpose_z_to_y(wk13,wk2_r2c,sp)

       ! ===== 1D FFTs in Y =====
#ifdef DOUBLE_PREC
       call cfft1i(sp%ysz(2), wsave, lensav, ier)
#else
       call cfft1i(sp%ysz(2), wsave, lensav, ier)
#endif
       do k=1,sp%ysz(3)
          do i=1,sp%ysz(1)
             do j=1,sp%ysz(2)
                buf_c(j) = wk2_r2c(i,j,k)
             end do
#ifdef DOUBLE_PREC
             call cfft1f(sp%ysz(2),1,buf_c,sp%ysz(2), &
                  wsave,lensav,work,lenwrk,ier)
#else
             call cfft1f(sp%ysz(2),1,buf_c,sp%ysz(2), &
                  wsave,lensav,work,lenwrk,ier)
#endif
             do j=1,sp%ysz(2)
                ! also remove the scaling here
                wk2_r2c(i,j,k) = buf_c(j) * real(sp%ysz(2), kind=mytype)
             end do
          end do
       end do

       ! ===== Swap Y --> X =====
       call transpose_y_to_x(wk2_r2c,out_c,sp)

       ! ===== 1D FFTs in X =====
#ifdef DOUBLE_PREC
       call cfft1i(sp%xsz(1), wsave, lensav, ier)
#else
       call cfft1i(sp%xsz(1), wsave, lensav, ier)
#endif
       do k=1,sp%xsz(3)
          do j=1,sp%xsz(2)
             do i=1,sp%xsz(1)
                buf_c(i) = out_c(i,j,k)
             end do
#ifdef DOUBLE_PREC
             call cfft1f(sp%xsz(1),1,buf_c,sp%xsz(1), &
                  wsave,lensav,work,lenwrk,ier)
#else
             call cfft1f(sp%xsz(1),1,buf_c,sp%xsz(1), &
                  wsave,lensav,work,lenwrk,ier)
#endif
             do i=1,sp%xsz(1)
                ! also remove the scaling here
                out_c(i,j,k) = buf_c(i) * real(sp%xsz(1), kind=mytype)
             end do
          end do
       end do

    end if

    return
  end subroutine fft_3d_r2c

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 3D inverse FFT - complex to real
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fft_3d_c2r(in_c, out_r)
    
    implicit none
    
    complex(mytype), dimension(:,:,:), intent(IN) :: in_c
    real(mytype), dimension(:,:,:), intent(OUT) :: out_r
    
    complex(mytype), allocatable, dimension(:,:,:) :: wk1
    integer :: i,j,k, ier

    if (format==PHYSICAL_IN_X) then

       ! ===== 1D FFTs in Z =====
       allocate (wk1(sp%zsz(1),sp%zsz(2),sp%zsz(3)))   
#ifdef DOUBLE_PREC
       call cfft1i(sp%zsz(3),wsave,lensav,ier)
#else
       call cfft1i(sp%zsz(3),wsave,lensav,ier)
#endif
       do j=1,sp%zsz(2)
          do i=1,sp%zsz(1)
             do k=1,sp%zsz(3)
                buf_c(k) = in_c(i,j,k)
             end do
#ifdef DOUBLE_PREC
             call cfft1b(sp%zsz(3),1,buf_c,sp%zsz(3), &
                  wsave,lensav,work,lenwrk,ier)
#else
             call cfft1b(sp%zsz(3),1,buf_c,sp%zsz(3), &
                  wsave,lensav,work,lenwrk,ier)
#endif
             do k=1,sp%zsz(3)
                wk1(i,j,k) = buf_c(k)
             end do
          end do
       end do

       ! ===== Swap Z --> Y =====
       call transpose_z_to_y(wk1,wk2_r2c,sp)

       ! ===== 1D FFTs in Y =====
#ifdef DOUBLE_PREC
       call cfft1i(sp%ysz(2),wsave,lensav,ier)
#else
       call cfft1i(sp%ysz(2),wsave,lensav,ier)
#endif
       do k=1,sp%ysz(3)
          do i=1,sp%ysz(1)
             do j=1,sp%ysz(2)
                buf_c(j) = wk2_r2c(i,j,k)
             end do
#ifdef DOUBLE_PREC
             call cfft1b(sp%ysz(2),1,buf_c,sp%ysz(2), &
                  wsave,lensav,work,lenwrk,ier)
#else
             call cfft1b(sp%ysz(2),1,buf_c,sp%ysz(2), &
                  wsave,lensav,work,lenwrk,ier)
#endif
             do j=1,sp%ysz(2)
                wk2_r2c(i,j,k) = buf_c(j)
             end do
          end do
       end do

       ! ===== Swap Y --> X =====
       call transpose_y_to_x(wk2_r2c,wk13,sp)

       ! ===== 1D FFTs in X =====
#ifdef DOUBLE_PREC
       call rfft1i(ph%xsz(1),wsave,lensav,ier)
#else
       call rfft1i(ph%xsz(1),wsave,lensav,ier)
#endif
       do k=1,ph%xsz(3)
          do j=1,ph%xsz(2)
             call from_complex(wk13(:,j,k), ph%xsz(1), buf_r)
#ifdef DOUBLE_PREC
             call rfft1b(ph%xsz(1),1,buf_r,ph%xsz(1), &
                  wsave,lensav,work,lenwrk,ier)
#else
             call rfft1b(ph%xsz(1),1,buf_r,ph%xsz(1), &
                  wsave,lensav,work,lenwrk,ier)
#endif
             do i=1,ph%xsz(1)
                out_r(i,j,k) = buf_r(i)
             end do
          end do
       end do

    else if (format==PHYSICAL_IN_Z) then

       ! ===== 1D FFTs in X =====
       allocate(wk1(sp%xsz(1),sp%xsz(2),sp%xsz(3)))
#ifdef DOUBLE_PREC
       call cfft1i(sp%xsz(1),wsave,lensav,ier)
#else
       call cfft1i(sp%xsz(1),wsave,lensav,ier)
#endif
       do k=1,sp%xsz(3)
          do j=1,sp%xsz(2)
             do i=1,sp%xsz(1)
                buf_c(i) = in_c(i,j,k)
             end do
#ifdef DOUBLE_PREC
             call cfft1b(sp%xsz(1),1,buf_c,sp%xsz(1), &
                  wsave,lensav,work,lenwrk,ier)
#else
             call cfft1b(sp%xsz(1),1,buf_c,sp%xsz(1), &
                  wsave,lensav,work,lenwrk,ier)
#endif
             do i=1,sp%xsz(1)
                wk1(i,j,k) = buf_c(i)
             end do
          end do
       end do

       ! ===== Swap X --> Y =====
       call transpose_x_to_y(wk1,wk2_r2c,sp)

       ! ===== 1D FFTs in Y =====
#ifdef DOUBLE_PREC
       call cfft1i(sp%ysz(2),wsave,lensav,ier)
#else
       call cfft1i(sp%ysz(2),wsave,lensav,ier)
#endif
       do k=1,sp%ysz(3)
          do i=1,sp%ysz(1)
             do j=1,sp%ysz(2)
                buf_c(j) = wk2_r2c(i,j,k)
             end do
#ifdef DOUBLE_PREC
             call cfft1b(sp%ysz(2),1,buf_c,sp%ysz(2), &
                  wsave,lensav,work,lenwrk,ier)
#else
             call cfft1b(sp%ysz(2),1,buf_c,sp%ysz(2), &
                  wsave,lensav,work,lenwrk,ier)
#endif
             do j=1,sp%ysz(2)
                wk2_r2c(i,j,k) = buf_c(j)
             end do
          end do
       end do

       ! ===== Swap Y --> Z =====
       call transpose_y_to_z(wk2_r2c,wk13,sp)

       ! ===== 1D FFTs in Z =====
#ifdef DOUBLE_PREC
       call rfft1i(ph%zsz(3),wsave,lensav,ier)
#else
       call rfft1i(ph%zsz(3),wsave,lensav,ier)
#endif
       do j=1,ph%zsz(2)
          do i=1,ph%zsz(1)
             do k=1,sp%zsz(3)
                buf_c(k) = wk13(i,j,k)
             end do
             call from_complex(buf_c, ph%zsz(3), buf_r)
#ifdef DOUBLE_PREC
             call rfft1b(ph%zsz(3),1,buf_r,ph%zsz(3), &
                  wsave,lensav,work,lenwrk,ier)
#else
             call rfft1b(ph%zsz(3),1,buf_r,ph%zsz(3), &
                  wsave,lensav,work,lenwrk,ier)
#endif
             do k=1,ph%zsz(3)
                out_r(i,j,k) = buf_r(k)
             end do
          end do
       end do

    end if

    return
  end subroutine fft_3d_c2r


  ! Two utility routines to convert internal storage format

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Rearrange data to half-plus-1-complex-storage format
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine to_complex(in,n,out)
    
    implicit none
    
    integer, intent(IN) :: n
    real(mytype), dimension(0:n-1), intent(IN) :: in
    complex(mytype), dimension(0:n/2), intent(OUT) :: out
    
    integer :: i
    
    ! From FFTPACK5's document, the 1D r2c result is returned in the
    ! same real input array and the storage format appears to be:
    !   -- a real number:  0th result (real by definition)
    !   -- real part and imaginary part of complex numbers
    !   -- another real number for the (N/2)th result, if N is even 
    ! Also all complex numbers appear to have a factor of 2 in them.
    
    out(0) = cmplx(in(0), 0._mytype, kind=mytype)

    if (n == n/2*2) then ! n is even
       do i=1,n/2-1
          out(i) = cmplx(in(2*i-1)*0.5_mytype, -in(2*i)*0.5_mytype, &
               kind=mytype)
       end do
       out(n/2) = cmplx(in(n-1), 0._mytype, kind=mytype)
    else                 ! n is odd
       do i=1,n/2
          out(i) = cmplx(in(2*i-1)*0.5_mytype, -in(2*i)*0.5_mytype, &
               kind=mytype)
       end do
    end if

    ! also remove the scaling as we need unscaled 1D FFT
    out = out * real(n, kind=mytype)
    
    return
  end subroutine to_complex

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! The opposite of the above
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine from_complex(in,n,out)
    
    implicit none
    
    integer, intent(IN) :: n
    complex(mytype), dimension(0:n/2), intent(IN) :: in
    real(mytype), dimension(0:n-1), intent(OUT) :: out
    
    integer :: i

    out(0) = real(in(0), kind=mytype)

    if (n == n/2*2) then ! n is even
       do i=1,n/2-1
          out(2*i-1) = real(in(i), kind=mytype)*2.0_mytype
          out(2*i)   = -aimag(in(i))*2.0_mytype
       end do
       out(n-1) = real(in(n/2), kind=mytype)
    else                 ! n is odd
       do i=1,n/2
          out(2*i-1) = real(in(i), kind=mytype)*2.0_mytype
          out(2*i)   = -aimag(in(i))*2.0_mytype
       end do
    end if

    return
  end subroutine from_complex

  
end module decomp_2d_fft
