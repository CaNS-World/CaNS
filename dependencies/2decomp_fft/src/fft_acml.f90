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

! This is the ACML implementation of the FFT library

module decomp_2d_fft
  
  use decomp_2d  ! 2D decomposition module
  
  implicit none
  
  private        ! Make everything private unless declared public

  ! engine-specific global variables
  integer, save :: plan_type = 0    ! *** TODO: debug mode 100

  real(mytype), parameter :: scale = 1.0_mytype ! compute unscaled FFT

  ! ACML plans
  complex(mytype), allocatable, dimension(:) :: &
       comm_c2c_x, comm_c2c_y, comm_c2c_z, &  ! for 3D c2c
       comm_c2c_y2, comm_c2c_z2, & ! for y/z sub-steps of 3D r2c
       comm_c2c_x2                 ! for y/x sub-steps of 3D r2c
  real(mytype), allocatable, dimension(:) :: comm_r2c_x, comm_c2r_x
  real(mytype), allocatable, dimension(:) :: comm_r2c_z, comm_c2r_z

  ! common code used for all engines, including global variables, 
  ! generic interface definitions and several subroutines
#include "fft_common.f90"

  ! subroutines that generate ACML plans
  ! for several types of 1D multiple FFTs.
#include "acml_plan.f90"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time initialisations for the FFT engine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_fft_engine

    implicit none

    if (nrank==0) then
       write(*,*) ' '
       write(*,*) '***** Using the ACML engine *****'
       write(*,*) ' '
    end if

    ! *** Read ACML doc carefully!
    ! Seems that c2c plans can be shared by forward/backward transforms
    ! r2c & c2r must use separate plans (otherwise problem with mode 100)

    ! For C2C transforms
    call c2c_1m_x_plan(comm_c2c_x, ph)
    call c2c_1m_y_plan(comm_c2c_y, ph)
    call c2c_1m_z_plan(comm_c2c_z, ph)

    ! For R2C/C2R tranforms
    if (format == PHYSICAL_IN_X) then
       call r2c_1m_x_plan(comm_r2c_x, ph)
       call c2r_1m_x_plan(comm_c2r_x, ph)
       call c2c_1m_y_plan(comm_c2c_y2, sp)
       call c2c_1m_z_plan(comm_c2c_z2, sp)
    else if (format == PHYSICAL_IN_Z) then
       call r2c_1m_z_plan(comm_r2c_z, ph)
       call c2r_1m_z_plan(comm_c2r_z, ph)
       call c2c_1m_y_plan(comm_c2c_y2, sp)
       call c2c_1m_x_plan(comm_c2c_x2, sp)
    end if

    return
  end subroutine init_fft_engine


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time finalisations for the FFT engine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine finalize_fft_engine

    implicit none

    deallocate(comm_c2c_x)
    deallocate(comm_c2c_y)
    deallocate(comm_c2c_z)

    if (format == PHYSICAL_IN_X) then
       deallocate(comm_r2c_x)
       deallocate(comm_c2r_x)
       deallocate(comm_c2c_y2)
       deallocate(comm_c2c_z2)
    else if (format==PHYSICAL_IN_Z) then
       deallocate(comm_r2c_z)
       deallocate(comm_c2r_z)
       deallocate(comm_c2c_y2)
       deallocate(comm_c2c_x2)
    end if

    return
  end subroutine finalize_fft_engine



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 3D FFT - complex to complex
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fft_3d_c2c(in, out, isign)
    
    implicit none
    
    complex(mytype), dimension(:,:,:), intent(INOUT) :: in
    complex(mytype), dimension(:,:,:), intent(OUT) :: out
    integer, intent(IN) :: isign   

#ifndef OVERWRITE
    complex(mytype), allocatable, dimension(:,:,:) :: wk1
#endif

    integer :: k

    integer :: info ! error flag 
    
    if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_FORWARD .OR.  &
         format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_BACKWARD) then
       
       ! ===== 1D FFTs in X =====
#ifndef OVERWRITE
       allocate (wk1(ph%xsz(1),ph%xsz(2),ph%xsz(3)))
#endif

#ifdef OVERWRITE

#ifdef DOUBLE_PREC
       call zfft1mx(isign,scale,.true.,ph%xsz(2)*ph%xsz(3),ph%xsz(1), &
            in,1,ph%xsz(1),in,1,ph%xsz(1),comm_c2c_x,info)
#else
       call cfft1mx(isign,scale,.true.,ph%xsz(2)*ph%xsz(3),ph%xsz(1), &
            in,1,ph%xsz(1),in,1,ph%xsz(1),comm_c2c_x,info)
#endif

#else

#ifdef DOUBLE_PREC
       call zfft1mx(isign,scale,.false.,ph%xsz(2)*ph%xsz(3),ph%xsz(1), &
            in,1,ph%xsz(1),wk1,1,ph%xsz(1),comm_c2c_x,info)
#else
       call cfft1mx(isign,scale,.false.,ph%xsz(2)*ph%xsz(3),ph%xsz(1), &
            in,1,ph%xsz(1),wk1,1,ph%xsz(1),comm_c2c_x,info)
#endif

#endif

       ! ===== Swap X --> Y =====
#ifdef OVERWRITE
       call transpose_x_to_y(in,wk2_c2c,ph)
#else
       call transpose_x_to_y(wk1,wk2_c2c,ph)
#endif
       
       ! ===== 1D FFTs in Y =====
       do k=1, ph%ysz(3)  ! loop through Z-planes
#ifdef DOUBLE_PREC
          call zfft1mx(isign,scale,.true.,ph%ysz(1),ph%ysz(2), &
               wk2_c2c(:,:,k),ph%ysz(1),1,wk2_c2c(:,:,k),ph%ysz(1), &
               1,comm_c2c_y,info)
#else
          call cfft1mx(isign,scale,.true.,ph%ysz(1),ph%ysz(2), &
               wk2_c2c(:,:,k),ph%ysz(1),1,wk2_c2c(:,:,k),ph%ysz(1), &
               1,comm_c2c_y,info)
#endif
       end do

       ! ===== Swap Y --> Z =====
       call transpose_y_to_z(wk2_c2c,out,ph)

       ! ===== 1D FFTs in Z =====
#ifdef DOUBLE_PREC
       call zfft1mx(isign,scale,.true.,ph%zsz(1)*ph%zsz(2),ph%zsz(3), &
            out,ph%zsz(1)*ph%zsz(2),1,out,ph%zsz(1)*ph%zsz(2),1, &
            comm_c2c_z,info)
#else
       call cfft1mx(isign,scale,.true.,ph%zsz(1)*ph%zsz(2),ph%zsz(3), &
            out,ph%zsz(1)*ph%zsz(2),1,out,ph%zsz(1)*ph%zsz(2),1, &
            comm_c2c_z,info)
#endif       

    else if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_BACKWARD &
         .OR. & 
         format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_FORWARD) then

       ! ===== 1D FFTs in Z =====
#ifndef OVERWRITE
       allocate (wk1(ph%zsz(1),ph%zsz(2),ph%zsz(3)))
#endif

#ifdef OVERWRITE

#ifdef DOUBLE_PREC
       call zfft1mx(isign,scale,.true.,ph%zsz(1)*ph%zsz(2),ph%zsz(3), &
            in,ph%zsz(1)*ph%zsz(2),1,in,ph%zsz(1)*ph%zsz(2),1, &
            comm_c2c_z,info)
#else
       call cfft1mx(isign,scale,.true.,ph%zsz(1)*ph%zsz(2),ph%zsz(3), &
            in,ph%zsz(1)*ph%zsz(2),1,in,ph%zsz(1)*ph%zsz(2),1, &
            comm_c2c_z,info)
#endif

#else

#ifdef DOUBLE_PREC
       call zfft1mx(isign,scale,.false.,ph%zsz(1)*ph%zsz(2),ph%zsz(3), &
            in,ph%zsz(1)*ph%zsz(2),1,wk1,ph%zsz(1)*ph%zsz(2),1, &
            comm_c2c_z,info)
#else
       call cfft1mx(isign,scale,.false.,ph%zsz(1)*ph%zsz(2),ph%zsz(3), &
            in,ph%zsz(1)*ph%zsz(2),1,wk1,ph%zsz(1)*ph%zsz(2),1, &
            comm_c2c_z,info)
#endif

#endif

       ! ===== Swap Z --> Y =====
#ifdef OVERWRITE
       call transpose_z_to_y(in,wk2_c2c,ph)
#else
       call transpose_z_to_y(wk1,wk2_c2c,ph)
#endif

       ! ===== 1D FFTs in Y =====
       do k=1, ph%ysz(3)  ! loop through Z-planes
#ifdef DOUBLE_PREC
          call zfft1mx(isign,scale,.true.,ph%ysz(1),ph%ysz(2), &
               wk2_c2c(:,:,k),ph%ysz(1),1,wk2_c2c(:,:,k),ph%ysz(1),1, &
               comm_c2c_y,info)
#else
          call cfft1mx(isign,scale,.true.,ph%ysz(1),ph%ysz(2), &
               wk2_c2c(:,:,k),ph%ysz(1),1,wk2_c2c(:,:,k),ph%ysz(1),1, &
               comm_c2c_y,info)
#endif
       end do

       ! ===== Swap Y --> X =====
       call transpose_y_to_x(wk2_c2c,out,ph)
       
       ! ===== 1D FFTs in X =====
#ifdef DOUBLE_PREC
       call zfft1mx(isign,scale,.true.,ph%xsz(2)*ph%xsz(3),ph%xsz(1), &
            out,1,ph%xsz(1),out,1,ph%xsz(1),comm_c2c_x,info)
#else
       call cfft1mx(isign,scale,.true.,ph%xsz(2)*ph%xsz(3),ph%xsz(1), &
            out,1,ph%xsz(1),out,1,ph%xsz(1),comm_c2c_x,info)
#endif      

    end if

    return
  end subroutine fft_3d_c2c
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 3D forward FFT - real to complex
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fft_3d_r2c(in_r, out_c)
    
    implicit none
    
    real(mytype), dimension(:,:,:), intent(INOUT) :: in_r
    complex(mytype), dimension(:,:,:), intent(OUT) :: out_c
    
#ifndef OVERWRITE
    real(mytype), allocatable, dimension(:,:,:) :: in_r_copy
#endif
    complex(mytype), allocatable, dimension(:,:,:) :: wk1
    integer :: i,j,k
    integer :: info 

    real(mytype), allocatable, dimension(:) :: tmp, zbuf_r
    complex(mytype), allocatable, dimension(:) :: zbuf_c

    if (format==PHYSICAL_IN_X) then

       ! ===== 1D FFTs in X =====
       allocate(wk1(sp%xsz(1),sp%xsz(2),sp%xsz(3)))
#ifndef OVERWRITE
       ! need a copy of input as dz/scfftm overwrite the input
       allocate(in_r_copy(ph%xsz(1),ph%xsz(2),ph%xsz(3)))
       in_r_copy = in_r
#endif

#ifdef OVERWRITE

#ifdef DOUBLE_PREC
       call dzfftm(ph%xsz(2)*ph%xsz(3),ph%xsz(1),in_r,comm_r2c_x,info)
#else
       call scfftm(ph%xsz(2)*ph%xsz(3),ph%xsz(1),in_r,comm_r2c_x,info)
#endif

#else 

#ifdef DOUBLE_PREC
       call dzfftm(ph%xsz(2)*ph%xsz(3),ph%xsz(1),in_r_copy, &
            comm_r2c_x,info)
#else
       call scfftm(ph%xsz(2)*ph%xsz(3),ph%xsz(1),in_r_copy, &
            comm_r2c_x,info)
#endif

#endif

       ! translate Hermitian storage into half-plus-1-complex-storage
       allocate(tmp(ph%xsz(1)))
       do k=1,ph%xsz(3)
          do j=1,ph%xsz(2)
#ifdef OVERWRITE
             tmp(1:ph%xsz(1)) = in_r(1:ph%xsz(1),j,k)
#else
             tmp(1:ph%xsz(1)) = in_r_copy(1:ph%xsz(1),j,k)
#endif
             call to_complex(tmp, ph%xsz(1), wk1(:,j,k))
          end do
       end do
#ifdef OVERWRITE
       deallocate(tmp)
#else
       deallocate(tmp, in_r_copy)
#endif

       ! remove the scaling 
       wk1 = wk1 * sqrt(real(ph%xsz(1), kind=mytype))

       ! ===== Swap X --> Y =====
       call transpose_x_to_y(wk1,wk2_r2c,sp)

       ! ===== 1D FFTs in Y =====
       do k=1,sp%ysz(3) ! compute one z-plane at a time
#ifdef DOUBLE_PREC
          call zfft1mx(-1,scale,.true.,sp%ysz(1),sp%ysz(2),  &
               wk2_r2c(:,:,k),sp%ysz(1),1,wk2_r2c(:,:,k),sp%ysz(1),1, &
               comm_c2c_y2,info)
#else
          call cfft1mx(-1,scale,.true.,sp%ysz(1),sp%ysz(2),  &
               wk2_r2c(:,:,k),sp%ysz(1),1,wk2_r2c(:,:,k),sp%ysz(1),1, &
               comm_c2c_y2,info)
#endif
       end do
       
       ! ===== Swap Y --> Z =====
       call transpose_y_to_z(wk2_r2c,out_c,sp)

       ! ===== 1D FFTs in Z =====
#ifdef DOUBLE_PREC
       call zfft1mx(-1,scale,.true.,sp%zsz(1)*sp%zsz(2),sp%zsz(3), &
            out_c,sp%zsz(1)*sp%zsz(2),1,out_c,sp%zsz(1)*sp%zsz(2),1, &
            comm_c2c_z2, info)
#else
       call cfft1mx(-1,scale,.true.,sp%zsz(1)*sp%zsz(2),sp%zsz(3), &
            out_c,sp%zsz(1)*sp%zsz(2),1,out_c,sp%zsz(1)*sp%zsz(2),1, &
            comm_c2c_z2, info)
#endif


    else if (format==PHYSICAL_IN_Z) then

       ! ===== 1D FFTs in Z =====
       ! note that there is no expert driver to perform multiple 1D r2c
       ! FFTs in ACML --> need to copy numbers into buffers and then
       ! make multiple calls to simple driver.
       allocate(zbuf_r(ph%zsz(3)))
       allocate(zbuf_c(sp%zsz(3)))
       allocate(wk1(sp%zsz(1),sp%zsz(2),sp%zsz(3)))
       do j=1,ph%zsz(2)
          do i=1,ph%zsz(1)
             ! copy data set along Z into a buffer
             do k=1,ph%zsz(3)
                zbuf_r(k) = in_r(i,j,k)
             end do
             ! 1D r2c FFT 
#ifdef DOUBLE_PREC
             call dzfft(1,ph%zsz(3),zbuf_r,comm_r2c_z,info)
#else
             call scfft(1,ph%zsz(3),zbuf_r,comm_r2c_z,info)
#endif
             ! The result is in Hermitian storage, convert to complex
             call to_complex(zbuf_r,ph%zsz(3),zbuf_c)
             ! copy back into 3D data structure
             do k=1,sp%zsz(3) 
                wk1(i,j,k) = zbuf_c(k)
             end do
          end do
       end do
       ! remove ACML scaling, the parallel library never scales anything
       wk1 = wk1 * sqrt(real(ph%zsz(3), kind=mytype))

       ! ===== Swap Z --> Y =====
       call transpose_z_to_y(wk1,wk2_r2c,sp)

       ! ===== 1D FFTs in Y =====
       do k=1,sp%ysz(3) ! compute one z-plane at a time
#ifdef DOUBLE_PREC
          call zfft1mx(-1,scale,.true.,sp%ysz(1),sp%ysz(2),  &
               wk2_r2c(:,:,k),sp%ysz(1),1,wk2_r2c(:,:,k),sp%ysz(1),1, &
               comm_c2c_y2,info)
#else
          call cfft1mx(-1,scale,.true.,sp%ysz(1),sp%ysz(2),  &
               wk2_r2c(:,:,k),sp%ysz(1),1,wk2_r2c(:,:,k),sp%ysz(1),1, &
               comm_c2c_y2,info)
#endif
       end do      

       ! ===== Swap Y --> X =====
       call transpose_y_to_x(wk2_r2c,out_c,sp)

       ! ===== 1D FFTs in X =====
#ifdef DOUBLE_PREC
       call zfft1mx(-1,scale,.true.,sp%xsz(2)*sp%xsz(3),sp%xsz(1), &
            out_c,1,sp%xsz(1),out_c, 1, sp%xsz(1),comm_c2c_x2, info)
#else
       call cfft1mx(-1,scale,.true.,sp%xsz(2)*sp%xsz(3),sp%xsz(1), &
            out_c,1,sp%xsz(1),out_c, 1, sp%xsz(1),comm_c2c_x2, info)
#endif

    end if
    
    return
  end subroutine fft_3d_r2c
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 3D inverse FFT - complex to real
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fft_3d_c2r(in_c, out_r)
    
    implicit none
    
    complex(mytype), dimension(:,:,:), intent(INOUT) :: in_c
    real(mytype), dimension(:,:,:), intent(OUT) :: out_r
    
#ifndef OVERWRITE
    complex(mytype), allocatable, dimension(:,:,:) :: wk1
#endif
    complex(mytype), allocatable, dimension(:,:,:) :: wk3
    integer :: i,j,k
    integer :: info 

    real(mytype), allocatable, dimension(:) :: tmp, zbuf_r 
    complex(mytype), allocatable, dimension(:) :: zbuf_c

    if (format==PHYSICAL_IN_X) then

       ! ===== 1D FFTs in Z =====
#ifndef OVERWRITE
       allocate(wk1(sp%zsz(1),sp%zsz(2),sp%zsz(3)))
#endif

#ifdef OVERWRITE

#ifdef DOUBLE_PREC
       call zfft1mx(1,scale,.true.,sp%zsz(1)*sp%zsz(2),sp%zsz(3),  &
            in_c,sp%zsz(1)*sp%zsz(2),1,in_c, sp%zsz(1)*sp%zsz(2), 1, &
            comm_c2c_z2, info)
#else
       call cfft1mx(1,scale,.true.,sp%zsz(1)*sp%zsz(2),sp%zsz(3),  &
            in_c,sp%zsz(1)*sp%zsz(2),1,in_c, sp%zsz(1)*sp%zsz(2), 1, &
            comm_c2c_z2, info)
#endif 

#else

#ifdef DOUBLE_PREC
       call zfft1mx(1,scale,.false.,sp%zsz(1)*sp%zsz(2),sp%zsz(3),  &
            in_c,sp%zsz(1)*sp%zsz(2),1,wk1, sp%zsz(1)*sp%zsz(2), 1, &
            comm_c2c_z2, info)
#else
       call cfft1mx(1,scale,.false.,sp%zsz(1)*sp%zsz(2),sp%zsz(3),  &
            in_c,sp%zsz(1)*sp%zsz(2),1,wk1, sp%zsz(1)*sp%zsz(2), 1, &
            comm_c2c_z2, info)
#endif 

#endif

       ! ===== Swap Z --> Y =====
#ifdef OVERWRITE
       call transpose_z_to_y(in_c,wk2_r2c,sp)
#else
       call transpose_z_to_y(wk1,wk2_r2c,sp)
#endif

       ! ===== 1D FFTs in Y =====
       do k=1,sp%ysz(3) ! compute one z-plane at a time
#ifdef DOUBLE_PREC
          call zfft1mx(1,scale,.true.,sp%ysz(1),sp%ysz(2),  &
               wk2_r2c(:,:,k),sp%ysz(1),1,wk2_r2c(:,:,k),sp%ysz(1),1, &
               comm_c2c_y2,info)
#else
          call cfft1mx(1,scale,.true.,sp%ysz(1),sp%ysz(2),  &
               wk2_r2c(:,:,k),sp%ysz(1),1,wk2_r2c(:,:,k),sp%ysz(1),1, &
               comm_c2c_y2,info)
#endif
       end do

       ! ===== Swap Y --> X =====
       allocate (wk3(sp%xsz(1),sp%xsz(2),sp%xsz(3)))
       call transpose_y_to_x(wk2_r2c,wk3,sp)
       
       ! translate half-plus-1-complex-storage into Hermitian storage
       allocate(tmp(ph%xsz(1)))
       do k=1,ph%xsz(3)
          do j=1,ph%xsz(2)
             call to_hermitian(wk3(:,j,k), ph%xsz(1), tmp)
             out_r(1:ph%xsz(1),j,k) = tmp(1:ph%xsz(1))
          end do
       end do
       deallocate(tmp)
       out_r = out_r * sqrt(real(ph%xsz(1),kind=mytype))

       ! ===== 1D FFTs in X =====
#ifdef DOUBLE_PREC
       call zdfftm(ph%xsz(2)*ph%xsz(3),ph%xsz(1),out_r,comm_c2r_x,info)
#else
       call csfftm(ph%xsz(2)*ph%xsz(3),ph%xsz(1),out_r,comm_c2r_x,info)
#endif

    else if (format==PHYSICAL_IN_Z) then

       ! ===== 1D FFTs in X =====
#ifndef OVERWRITE
       allocate(wk1(sp%xsz(1),sp%xsz(2),sp%xsz(3)))
#endif

#ifdef OVERWRITE

#ifdef DOUBLE_PREC
       call zfft1mx(1,scale,.true.,sp%xsz(2)*sp%xsz(3),sp%xsz(1), &
            in_c,1,sp%xsz(1),in_c,1,sp%xsz(1),comm_c2c_x2,info)
#else
       call cfft1mx(1,scale,.true.,sp%xsz(2)*sp%xsz(3),sp%xsz(1), &
            in_c,1,sp%xsz(1),in_c,1,sp%xsz(1),comm_c2c_x2,info)
#endif

#else

#ifdef DOUBLE_PREC
       call zfft1mx(1,scale,.false.,sp%xsz(2)*sp%xsz(3),sp%xsz(1), &
            in_c,1,sp%xsz(1),wk1,1,sp%xsz(1),comm_c2c_x2,info)
#else
       call cfft1mx(1,scale,.false.,sp%xsz(2)*sp%xsz(3),sp%xsz(1), &
            in_c,1,sp%xsz(1),wk1,1,sp%xsz(1),comm_c2c_x2,info)
#endif

#endif

       ! ===== Swap X --> Y =====
#ifdef OVERWRITE
       call transpose_x_to_y(in_c,wk2_r2c,sp)
#else
       call transpose_x_to_y(wk1,wk2_r2c,sp)
#endif

       ! ===== 1D FFTs in Y =====

       do k=1,sp%ysz(3) ! compute one z-plane at a time
#ifdef DOUBLE_PREC
          call zfft1mx(1,scale,.true.,sp%ysz(1),sp%ysz(2),  &
               wk2_r2c(:,:,k),sp%ysz(1),1,wk2_r2c(:,:,k),sp%ysz(1),1, &
               comm_c2c_y2,info)
#else
          call cfft1mx(1,scale,.true.,sp%ysz(1),sp%ysz(2),  &
               wk2_r2c(:,:,k),sp%ysz(1),1,wk2_r2c(:,:,k),sp%ysz(1),1, &
               comm_c2c_y2,info)
#endif
       end do

       ! ===== Swap Y --> Z =====
       allocate (wk3(sp%zsz(1),sp%zsz(2),sp%zsz(3)))
       call transpose_y_to_z(wk2_r2c,wk3,sp)

       ! ===== 1D FFTs in Z =====
       ! note that there is no expert driver to perform multiple 1D c2r
       ! FFTs in ACML --> need to copy numbers into buffers and then
       ! make multiple calls to simple driver.
       allocate(zbuf_r(ph%zsz(3)))
       allocate(zbuf_c(sp%zsz(3)))
       do j=1,ph%zsz(2)
          do i=1,ph%zsz(1)
             ! copy data set along Z into a buffer
             do k=1,sp%zsz(3) 
                zbuf_c(k) = wk3(i,j,k)
             end do
             ! convert complex storage to Hermitian storage
             call to_hermitian(zbuf_c,ph%zsz(3),zbuf_r)
             ! 1D c2r FFT
#ifdef DOUBLE_PREC
             call zdfft(1,ph%zsz(3),zbuf_r,comm_c2r_z,info)
#else
             call csfft(1,ph%zsz(3),zbuf_r,comm_c2r_z,info)
#endif
             ! copy back into 3D data structure
             do k=1,ph%zsz(3)
                out_r(i,j,k) = zbuf_r(k)
             end do
          end do
       end do
       ! remove ACML scaling, the parallel library never scales anything
       out_r = out_r * sqrt(real(ph%zsz(3),kind=mytype))

    end if

    return
  end subroutine fft_3d_c2r


  ! Two utility routines to convert internal storage format

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Rearrange Hermitian storage data as half-plus-1-complex-storage data
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine to_complex(in,n,out)
    
    implicit none
    
    integer, intent(IN) :: n
    real(mytype), dimension(0:n-1), intent(IN) :: in
    complex(mytype), dimension(0:n/2), intent(OUT) :: out
    
    integer :: i
    
    ! Let out(i) be the full-complex spectral representation of real input
    ! Let 'in' be an Hermitian storage array of length 'n' (0-based)
    !  - in(i) contains the real part of out(i) for i=0,...,n/2
    !  - in(n−i) contains the imaginary part of out(i) for i=1,...,(n−1)/2

    out(0) = cmplx(in(0),0._mytype,kind=mytype)
    
    if (n == n/2*2) then   ! n is even, n/2 is larger than (n-1)/2 by 1
       do i=1,(n-1)/2
          out(i) = cmplx(in(i),in(n-i),kind=mytype)
       end do
       out(n/2) = cmplx(in(n/2),0._mytype,kind=mytype)
    else                   ! n is odd, then n/2 == (n-1)/2 
       do i=1,n/2
          out(i) = cmplx(in(i),in(n-i),kind=mytype)
       end do
    end if
    
    return
  end subroutine to_complex

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! The opposite of the above
  !   Also need to be conjudated before inverse transform
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine to_hermitian(in,n,out)
    
    implicit none
    
    integer, intent(IN) :: n
    complex(mytype), dimension(0:n/2), intent(IN) :: in
    real(mytype), dimension(0:n-1), intent(OUT) :: out
    
    integer :: i

    out(0) = real(in(0),kind=mytype)
    
    if (n == n/2*2) then   ! n is even, n/2 is larger than (n-1)/2 by 1
       do i=1,(n-1)/2
          out(i) = real(in(i),kind=mytype)
          out(n-i) = -aimag(in(i))
       end do
       out(n/2) = real(in(n/2), kind=mytype)
    else                   ! n is odd, then n/2 == (n-1)/2 
       do i=1,n/2
          out(i) = real(in(i),kind=mytype)
          out(n-i) = -aimag(in(i))
       end do
    end if

    return
  end subroutine to_hermitian

  
end module decomp_2d_fft
