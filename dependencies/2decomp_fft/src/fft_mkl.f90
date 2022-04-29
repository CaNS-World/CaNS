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

! This is the Intel MKL implementation of the FFT library

module decomp_2d_fft
  
  use decomp_2d  ! 2D decomposition module
  use MKL_DFTI   ! MKL FFT module
  
  implicit none
  
  private        ! Make everything private unless declared public

  ! engine-specific global variables
  
  ! Descriptors for MKL FFT, one for each set of 1D FFTs
  !  for c2c transforms
  type(DFTI_DESCRIPTOR), pointer :: c2c_x, c2c_y, c2c_z     
  !  for r2c/c2r transforms, PHYSICAL_IN_X
  type(DFTI_DESCRIPTOR), pointer :: r2c_x, c2c_y2, c2c_z2, c2r_x
  !  for r2c/c2r transforms, PHYSICAL_IN_Z
  type(DFTI_DESCRIPTOR), pointer :: r2c_z, c2c_x2, c2r_z

  ! common code used for all engines, including global variables, 
  ! generic interface definitions and several subroutines
#include "fft_common.f90"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time initialisations for the FFT engine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_fft_engine

    implicit none

    if (nrank==0) then
       write(*,*) ' '
       write(*,*) '***** Using the MKL engine *****'
       write(*,*) ' '
    end if

    ! For C2C transforms
    call c2c_1m_x_plan(c2c_x, ph)
    call c2c_1m_y_plan(c2c_y, ph)
    call c2c_1m_z_plan(c2c_z, ph)

    ! For R2C/C2R tranfroms with physical space in X-pencil
    if (format == PHYSICAL_IN_X) then
       call r2c_1m_x_plan(r2c_x, ph, sp, -1)
       call c2c_1m_y_plan(c2c_y2, sp)
       call c2c_1m_z_plan(c2c_z2, sp)
       call r2c_1m_x_plan(c2r_x, ph, sp,  1)

    ! For R2C/C2R tranfroms with physical space in Z-pencil
    else if (format == PHYSICAL_IN_Z) then
       call r2c_1m_z_plan(r2c_z, ph, sp, -1)
       call c2c_1m_y_plan(c2c_y2, sp)
       call c2c_1m_x_plan(c2c_x2, sp)
       call r2c_1m_z_plan(c2r_z, ph, sp,  1)
    end if

    return
  end subroutine init_fft_engine

  
  ! Return an MKL plan for multiple 1D c2c FFTs in X direction
  subroutine c2c_1m_x_plan(desc, decomp)

    implicit none

    type(DFTI_DESCRIPTOR), pointer :: desc
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: status

#ifdef DOUBLE_PREC
    status = DftiCreateDescriptor(desc, DFTI_DOUBLE, &
         DFTI_COMPLEX, 1, decomp%xsz(1))
#else       
    status = DftiCreateDescriptor(desc, DFTI_SINGLE, &
         DFTI_COMPLEX, 1, decomp%xsz(1))
#endif
    status = DftiSetValue(desc, DFTI_NUMBER_OF_TRANSFORMS, &
         decomp%xsz(2)*decomp%xsz(3))
    status = DftiSetValue(desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    status = DftiSetValue(desc, DFTI_INPUT_DISTANCE, decomp%xsz(1))
    status = DftiSetValue(desc, DFTI_OUTPUT_DISTANCE, decomp%xsz(1))
    status = DftiCommitDescriptor(desc)

    return
  end subroutine c2c_1m_x_plan

  ! Return an MKL plan for multiple 1D c2c FFTs in Y direction
  subroutine c2c_1m_y_plan(desc, decomp)

    implicit none

    type(DFTI_DESCRIPTOR), pointer :: desc
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: status, strides(2)

    ! Due to memory pattern of 3D arrays, 1D FFTs along Y have to be
    ! done one Z-plane at a time. So plan for 2D data sets here.

#ifdef DOUBLE_PREC
    status = DftiCreateDescriptor(desc, DFTI_DOUBLE, &
         DFTI_COMPLEX, 1, decomp%ysz(2))
#else
    status = DftiCreateDescriptor(desc, DFTI_SINGLE, &
         DFTI_COMPLEX, 1, decomp%ysz(2))
#endif
    status = DftiSetValue(desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    status = DftiSetValue(desc, DFTI_NUMBER_OF_TRANSFORMS, decomp%ysz(1))
    status = DftiSetValue(desc, DFTI_INPUT_DISTANCE, 1)
    status = DftiSetValue(desc, DFTI_OUTPUT_DISTANCE, 1)
    strides(1) = 0
    strides(2) = decomp%ysz(1)
    status = DftiSetValue(desc, DFTI_INPUT_STRIDES, strides)
    status = DftiSetValue(desc, DFTI_OUTPUT_STRIDES, strides)
    status = DftiCommitDescriptor(desc)

    return
  end subroutine c2c_1m_y_plan

  ! Return an MKL plan for multiple 1D c2c FFTs in Z direction
  subroutine c2c_1m_z_plan(desc, decomp)

    implicit none

    type(DFTI_DESCRIPTOR), pointer :: desc
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: status, strides(2)

#ifdef DOUBLE_PREC
    status = DftiCreateDescriptor(desc, DFTI_DOUBLE, &
         DFTI_COMPLEX, 1, decomp%zsz(3))
#else
    status = DftiCreateDescriptor(desc, DFTI_SINGLE, &
         DFTI_COMPLEX, 1, decomp%zsz(3))
#endif
    status = DftiSetValue(desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    status = DftiSetValue(desc, DFTI_NUMBER_OF_TRANSFORMS, &
         decomp%zsz(1)*decomp%zsz(2))
    status = DftiSetValue(desc, DFTI_INPUT_DISTANCE, 1)
    status = DftiSetValue(desc, DFTI_OUTPUT_DISTANCE, 1)
    strides(1) = 0
    strides(2) = decomp%zsz(1)*decomp%zsz(2)
    status = DftiSetValue(desc, DFTI_INPUT_STRIDES, strides)
    status = DftiSetValue(desc, DFTI_OUTPUT_STRIDES, strides)
    status = DftiCommitDescriptor(desc)    

    return
  end subroutine c2c_1m_z_plan

  ! Return an MKL plan for multiple 1D r2c FFTs in X direction
  subroutine r2c_1m_x_plan(desc, decomp_ph, decomp_sp, direction)

    implicit none

    type(DFTI_DESCRIPTOR), pointer :: desc
    TYPE(DECOMP_INFO), intent(IN) :: decomp_ph, decomp_sp
    integer, intent(IN) :: direction ! (-1=r2c; 1=c2r)

    integer :: status

    ! c2r and r2c plans are almost the same, just swap input/output

#ifdef DOUBLE_PREC
    status = DftiCreateDescriptor(desc, DFTI_DOUBLE, &
         DFTI_REAL, 1, decomp_ph%xsz(1))
#else
    status = DftiCreateDescriptor(desc, DFTI_SINGLE, &
         DFTI_REAL, 1, decomp_ph%xsz(1))
#endif
    status = DftiSetValue(desc, DFTI_NUMBER_OF_TRANSFORMS, &
         decomp_ph%xsz(2)*decomp_ph%xsz(3))
    status = DftiSetValue(desc, DFTI_CONJUGATE_EVEN_STORAGE, &
         DFTI_COMPLEX_COMPLEX)
    status = DftiSetValue(desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    if (direction == -1) then  ! r2c
       status = DftiSetValue(desc, DFTI_INPUT_DISTANCE, &
            decomp_ph%xsz(1))
       status = DftiSetValue(desc, DFTI_OUTPUT_DISTANCE, &
            decomp_sp%xsz(1))
    else if (direction == 1) then  ! c2r
       status = DftiSetValue(desc, DFTI_INPUT_DISTANCE, &
            decomp_sp%xsz(1))
       status = DftiSetValue(desc, DFTI_OUTPUT_DISTANCE, &
            decomp_ph%xsz(1))
    end if
    status = DftiCommitDescriptor(desc)
    
    return
  end subroutine r2c_1m_x_plan

  ! Return an MKL plan for multiple 1D r2c FFTs in Z direction
  subroutine r2c_1m_z_plan(desc, decomp_ph, decomp_sp, direction)

    implicit none

    type(DFTI_DESCRIPTOR), pointer :: desc
    TYPE(DECOMP_INFO), intent(IN) :: decomp_ph, decomp_sp
    integer, intent(IN) :: direction ! (-1=r2c; 1=c2r)

    integer :: status, strides(2)

    ! c2r and r2c plans are almost the same, just swap input/output
    
#ifdef DOUBLE_PREC
    status = DftiCreateDescriptor(desc, DFTI_DOUBLE, &
         DFTI_REAL, 1, decomp_ph%zsz(3))
#else
    status = DftiCreateDescriptor(desc, DFTI_SINGLE, &
         DFTI_REAL, 1, decomp_ph%zsz(3))
#endif
    status = DftiSetValue(desc, DFTI_NUMBER_OF_TRANSFORMS, &
         decomp_ph%zsz(1)*decomp_ph%zsz(2))
    status = DftiSetValue(desc, DFTI_CONJUGATE_EVEN_STORAGE, &
         DFTI_COMPLEX_COMPLEX)
    status = DftiSetValue(desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    status = DftiSetValue(desc, DFTI_INPUT_DISTANCE, 1)
    status = DftiSetValue(desc, DFTI_OUTPUT_DISTANCE, 1)
    strides(1) = 0
    strides(2) = decomp_ph%zsz(1)*decomp_ph%zsz(2)
    if (direction == -1) then
       status = DftiSetValue(desc, DFTI_INPUT_STRIDES, strides)
    else if (direction == 1) then
       status = DftiSetValue(desc, DFTI_OUTPUT_STRIDES, strides)
    end if
    strides(2) = decomp_sp%zsz(1)*decomp_sp%zsz(2)
    if (direction == -1) then
       status = DftiSetValue(desc, DFTI_OUTPUT_STRIDES, strides)
    else if (direction == 1) then
       status = DftiSetValue(desc, DFTI_INPUT_STRIDES, strides)
    end if
    status = DftiCommitDescriptor(desc)
    
    return
  end subroutine r2c_1m_z_plan


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time finalisations for the FFT engine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine finalize_fft_engine

    implicit none

    integer :: status
    
    status = DftiFreeDescriptor(c2c_x)
    status = DftiFreeDescriptor(c2c_y)
    status = DftiFreeDescriptor(c2c_z)
    if (format==PHYSICAL_IN_X) then
       status = DftiFreeDescriptor(r2c_x)
       status = DftiFreeDescriptor(c2c_z2)
       status = DftiFreeDescriptor(c2r_x)
    else if (format==PHYSICAL_IN_Z) then
       status = DftiFreeDescriptor(r2c_z)
       status = DftiFreeDescriptor(c2c_x2)
       status = DftiFreeDescriptor(c2r_z)
    end if
    status = DftiFreeDescriptor(c2c_y2)

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

    complex(mytype), allocatable, dimension(:,:,:) :: wk1,wk2,wk2b,wk3
    integer :: k, status

    if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_FORWARD .OR.  &
         format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_BACKWARD) then
       
       ! ===== 1D FFTs in X =====
       allocate (wk1(ph%xsz(1),ph%xsz(2),ph%xsz(3)))
!       if (isign==DECOMP_2D_FFT_FORWARD) then
!          status = DftiComputeForward(c2c_x, in(:,1,1), wk1(:,1,1))
!       else if (isign==DECOMP_2D_FFT_BACKWARD) then
!          status = DftiComputeBackward(c2c_x, in(:,1,1), wk1(:,1,1))
!       end if
       status = wrapper_c2c(c2c_x, in, wk1, isign)

       ! ===== Swap X --> Y =====
       allocate (wk2(ph%ysz(1),ph%ysz(2),ph%ysz(3)))
       call transpose_x_to_y(wk1,wk2,ph)
       
       ! ===== 1D FFTs in Y =====
       allocate (wk2b(ph%ysz(1),ph%ysz(2),ph%ysz(3))) 
       do k=1,ph%xsz(3) ! one Z-plane at a time
!          if (isign==DECOMP_2D_FFT_FORWARD) then
!             status = DftiComputeForward(c2c_y, wk2(:,1,k), wk2b(:,1,k))
!          else if (isign==DECOMP_2D_FFT_BACKWARD) then
!             status = DftiComputeBackward(c2c_y, wk2(:,1,k), wk2b(:,1,k))
!          end if
          status = wrapper_c2c(c2c_y, wk2(1,1,k), wk2b(1,1,k), isign)
       end do

       ! ===== Swap Y --> Z =====
       allocate (wk3(ph%zsz(1),ph%zsz(2),ph%zsz(3)))
       call transpose_y_to_z(wk2b,wk3,ph)

       ! ===== 1D FFTs in Z =====
!       if (isign==DECOMP_2D_FFT_FORWARD) then
!          status = DftiComputeForward(c2c_z, wk3(:,1,1), out(:,1,1))
!       else if (isign==DECOMP_2D_FFT_BACKWARD) then
!          status = DftiComputeBackward(c2c_z, wk3(:,1,1), out(:,1,1))
!       end if
       status = wrapper_c2c(c2c_z, wk3, out, isign)

    else if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_BACKWARD &
         .OR. & 
         format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_FORWARD) then

       ! ===== 1D FFTs in Z =====
       allocate (wk1(ph%zsz(1),ph%zsz(2),ph%zsz(3)))
!       if (isign==DECOMP_2D_FFT_FORWARD) then
!          status = DftiComputeForward(c2c_z, in(:,1,1), wk1(:,1,1))
!       else if (isign==DECOMP_2D_FFT_BACKWARD) then
!          status = DftiComputeBackward(c2c_z, in(:,1,1), wk1(:,1,1))
!       end if
       status = wrapper_c2c(c2c_z, in, wk1, isign)
       
       ! ===== Swap Z --> Y =====
       allocate (wk2(ph%ysz(1),ph%ysz(2),ph%ysz(3)))
       call transpose_z_to_y(wk1,wk2,ph)
       
       ! ===== 1D FFTs in Y =====
       allocate (wk2b(ph%ysz(1),ph%ysz(2),ph%ysz(3)))
       do k=1,ph%xsz(3) ! one Z-plane at a time
!          if (isign==DECOMP_2D_FFT_FORWARD) then
!             status = DftiComputeForward(c2c_y, wk2(:,1,k), wk2b(:,1,k))
!          else if (isign==DECOMP_2D_FFT_BACKWARD) then
!             status = DftiComputeBackward(c2c_y, wk2(:,1,k), wk2b(:,1,k))
!          end if
          status = wrapper_c2c(c2c_y, wk2(1,1,k), wk2b(1,1,k), isign)
       end do
       
       ! ===== Swap Y --> X =====
       allocate (wk3(ph%xsz(1),ph%xsz(2),ph%xsz(3)))
       call transpose_y_to_x(wk2b,wk3,ph)
       
       ! ===== 1D FFTs in X =====
!       if (isign==DECOMP_2D_FFT_FORWARD) then
!          status = DftiComputeForward(c2c_x, wk3(:,1,1), out(:,1,1))
!       else if (isign==DECOMP_2D_FFT_BACKWARD) then
!          status = DftiComputeBackward(c2c_x, wk3(:,1,1), out(:,1,1))
!       end if
       status = wrapper_c2c(c2c_x, wk3, out, isign)
       
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
    
    complex(mytype), allocatable, dimension(:,:,:) :: wk1,wk2,wk2b,wk3
    integer :: k, status, isign

    isign = DECOMP_2D_FFT_FORWARD

    if (format==PHYSICAL_IN_X) then

       ! ===== 1D FFTs in X =====
       allocate(wk1(sp%xsz(1),sp%xsz(2),sp%xsz(3)))
!       status = DftiComputeForward(r2c_x, in_r(:,1,1), wk1(:,1,1))       
       status = wrapper_r2c(r2c_x, in_r, wk1)

       ! ===== Swap X --> Y =====
       allocate (wk2(sp%ysz(1),sp%ysz(2),sp%ysz(3)))
       call transpose_x_to_y(wk1,wk2,sp)

       ! ===== 1D FFTs in Y =====
       allocate (wk2b(sp%ysz(1),sp%ysz(2),sp%ysz(3)))
       do k=1,sp%ysz(3)
!          status = DftiComputeForward(c2c_y2, wk2(:,1,k), wk2b(:,1,k))
          status = wrapper_c2c(c2c_y2, wk2(1,1,k), wk2b(1,1,k), isign)
       end do

       ! ===== Swap Y --> Z =====
       allocate (wk3(sp%zsz(1),sp%zsz(2),sp%zsz(3)))
       call transpose_y_to_z(wk2b,wk3,sp)

       ! ===== 1D FFTs in Z =====
!       status = DftiComputeForward(c2c_z2, wk3(:,1,1), out_c(:,1,1))
       status = wrapper_c2c(c2c_z2, wk3, out_c, isign)
                
    else if (format==PHYSICAL_IN_Z) then

       ! ===== 1D FFTs in Z =====
       allocate(wk1(sp%zsz(1),sp%zsz(2),sp%zsz(3)))
!       status = DftiComputeForward(r2c_z, in_r(:,1,1), wk1(:,1,1))
       status = wrapper_r2c(r2c_z, in_r, wk1)

       ! ===== Swap Z --> Y =====
       allocate (wk2(sp%ysz(1),sp%ysz(2),sp%ysz(3)))
       call transpose_z_to_y(wk1,wk2,sp)

       ! ===== 1D FFTs in Y =====
       allocate (wk2b(sp%ysz(1),sp%ysz(2),sp%ysz(3)))
       do k=1,sp%ysz(3)
!          status = DftiComputeForward(c2c_y2, wk2(:,1,k), wk2b(:,1,k))
          status = wrapper_c2c(c2c_y2, wk2(1,1,k), wk2b(1,1,k), isign)
       end do

       ! ===== Swap Y --> X =====
       allocate (wk3(sp%xsz(1),sp%xsz(2),sp%xsz(3)))
       call transpose_y_to_x(wk2b,wk3,sp)

       ! ===== 1D FFTs in X =====
!       status = DftiComputeForward(c2c_x2, wk3(:,1,1), out_c(:,1,1))
       status = wrapper_c2c(c2c_x2, wk3, out_c, isign)

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
    
    complex(mytype), allocatable, dimension(:,:,:) :: wk1,wk2,wk2b,wk3
    integer :: k, status, isign

    isign = DECOMP_2D_FFT_BACKWARD

    if (format==PHYSICAL_IN_X) then

       ! ===== 1D FFTs in Z =====
       allocate (wk1(sp%zsz(1),sp%zsz(2),sp%zsz(3)))
!       status = DftiComputeBackward(c2c_z2, in_c(:,1,1), wk1(:,1,1))
       status = wrapper_c2c(c2c_z2, in_c, wk1, isign)

       ! ===== Swap Z --> Y =====
       allocate (wk2(sp%ysz(1),sp%ysz(2),sp%ysz(3)))
       call transpose_z_to_y(wk1,wk2,sp)

       ! ===== 1D FFTs in Y =====
       allocate (wk2b(sp%ysz(1),sp%ysz(2),sp%ysz(3)))
       do k=1,sp%ysz(3)
!          status = DftiComputeBackward(c2c_y2, wk2(:,1,k), wk2b(:,1,k))
          status = wrapper_c2c(c2c_y2, wk2(1,1,k), wk2b(1,1,k), isign)
       end do

       ! ===== Swap Y --> X =====
       allocate (wk3(sp%xsz(1),sp%xsz(2),sp%xsz(3)))
       call transpose_y_to_x(wk2b,wk3,sp)

       ! ===== 1D FFTs in X =====
!       status = DftiComputeBackward(c2r_x, wk3(:,1,1), out_r(:,1,1))
       status =  wrapper_c2r(c2r_x, wk3, out_r)

    else if (format==PHYSICAL_IN_Z) then

       ! ===== 1D FFTs in X =====
       allocate(wk1(sp%xsz(1),sp%xsz(2),sp%xsz(3)))
!       status = DftiComputeBackward(c2c_x2, in_c(:,1,1), wk1(:,1,1))
       status = wrapper_c2c(c2c_x2, in_c, wk1, isign)

       ! ===== Swap X --> Y =====
       allocate (wk2(sp%ysz(1),sp%ysz(2),sp%ysz(3)))
       call transpose_x_to_y(wk1,wk2,sp)

       ! ===== 1D FFTs in Y =====
       allocate (wk2b(sp%ysz(1),sp%ysz(2),sp%ysz(3)))
       do k=1,sp%ysz(3)
!          status = DftiComputeBackward(c2c_y2, wk2(:,1,k), wk2b(:,1,k))
          status = wrapper_c2c(c2c_y2, wk2(1,1,k), wk2b(1,1,k), isign)
       end do

       ! ===== Swap Y --> Z =====
       allocate (wk3(sp%zsz(1),sp%zsz(2),sp%zsz(3)))
       call transpose_y_to_z(wk2b,wk3,sp)

       ! ===== 1D FFTs in Z =====
!       status = DftiComputeBackward(c2r_z, wk3(:,1,1), out_r(:,1,1))
       status =  wrapper_c2r(c2r_z, wk3, out_r)

    end if

    return
  end subroutine fft_3d_c2r


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Wrapper functions so that one can pass 3D arrays to DftiCompute
  !  -- MKL accepts only 1D arrays as input/output for its multi-
  !     dimensional FFTs.
  !  -- Using EQUIVALENCE as suggested by MKL documents is impossible 
  !     for allocated arrays, not to mention bad coding style
  !  -- All code commented out above may well work but not safe. There 
  !     is no guarantee that compiler wouldn't make copies of 1D arrays
  !     (which would contain only one slice of the original 3D data) 
  !     rather than referring to the same memory address, i.e. 3D array
  !     A and 1D array A(:,1,1) may refer to different memory location.
  !  -- Using the following wrappers is safe and standard conforming.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function wrapper_c2c(desc, in, out, isign)

    implicit none

    type(DFTI_DESCRIPTOR), pointer :: desc
    complex(mytype), dimension(*) :: in, out
    integer :: isign, status

    if (isign == DECOMP_2D_FFT_FORWARD) then
       status = DftiComputeForward(desc, in, out)
    else if (isign == DECOMP_2D_FFT_BACKWARD) then
       status = DftiComputeBackward(desc, in, out)
    end if

    wrapper_c2c = status

    return
  end function wrapper_c2c

  integer function wrapper_r2c(desc, in, out)
    
    implicit none
    
    type(DFTI_DESCRIPTOR), pointer :: desc
    real(mytype), dimension(*) :: in 
    complex(mytype), dimension(*) :: out

    wrapper_r2c = DftiComputeForward(desc, in, out)

    return
  end function wrapper_r2c

  integer function wrapper_c2r(desc, in, out)

    implicit none

    type(DFTI_DESCRIPTOR), pointer :: desc
    complex(mytype), dimension(*) :: in 
    real(mytype), dimension(*) :: out

    wrapper_c2r = DftiComputeBackward(desc, in, out)

    return
  end function wrapper_c2r
  
end module decomp_2d_fft
