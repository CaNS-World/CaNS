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

! This file contains 3D c2c/r2c/c2r transform subroutines which are
! identical for several FFT engines 

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

    if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_FORWARD .OR.  &
         format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_BACKWARD) then
       
       ! ===== 1D FFTs in X =====
#ifdef OVERWRITE
       call c2c_1m_x(in,isign,ph)
#else
       allocate (wk1(ph%xsz(1),ph%xsz(2),ph%xsz(3)))
       wk1 = in
       call c2c_1m_x(wk1,isign,ph)
#endif

       ! ===== Swap X --> Y; 1D FFTs in Y =====

       if (dims(1)>1) then
#ifdef OVERWRITE
          call transpose_x_to_y(in,wk2_c2c,ph)
#else
          call transpose_x_to_y(wk1,wk2_c2c,ph)
#endif
          call c2c_1m_y(wk2_c2c,isign,ph)
       else
#ifdef OVERWRITE
          call c2c_1m_y(in,isign,ph)
#else
          call c2c_1m_y(wk1,isign,ph)
#endif
       end if

       ! ===== Swap Y --> Z; 1D FFTs in Z =====
       if (dims(1)>1) then
          call transpose_y_to_z(wk2_c2c,out,ph)
       else
#ifdef OVERWRITE
          call transpose_y_to_z(in,out,ph)
#else
          call transpose_y_to_z(wk1,out,ph)
#endif
       end if
       call c2c_1m_z(out,isign,ph)

    else if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_BACKWARD &
         .OR. & 
         format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_FORWARD) then

       ! ===== 1D FFTs in Z =====
#ifdef OVERWRITE
       call c2c_1m_z(in,isign,ph)
#else
       allocate (wk1(ph%zsz(1),ph%zsz(2),ph%zsz(3)))
       wk1 = in
       call c2c_1m_z(wk1,isign,ph)
#endif

       ! ===== Swap Z --> Y; 1D FFTs in Y =====
       if (dims(1)>1) then
#ifdef OVERWRITE
          call transpose_z_to_y(in,wk2_c2c,ph)
#else
          call transpose_z_to_y(wk1,wk2_c2c,ph)
#endif
          call c2c_1m_y(wk2_c2c,isign,ph)
       else  ! out==wk2_c2c if 1D decomposition
#ifdef OVERWRITE
          call transpose_z_to_y(in,out,ph)
#else
          call transpose_z_to_y(wk1,out,ph)
#endif
          call c2c_1m_y(out,isign,ph)
       end if

       ! ===== Swap Y --> X; 1D FFTs in X =====
       if (dims(1)>1) then
          call transpose_y_to_x(wk2_c2c,out,ph)
       end if
       call c2c_1m_x(out,isign,ph)
       
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

    if (format==PHYSICAL_IN_X) then

       ! ===== 1D FFTs in X =====
       call r2c_1m_x(in_r,wk13)

       ! ===== Swap X --> Y; 1D FFTs in Y =====
       if (dims(1)>1) then
          call transpose_x_to_y(wk13,wk2_r2c,sp)
          call c2c_1m_y(wk2_r2c,-1,sp)
       else
          call c2c_1m_y(wk13,-1,sp)
       end if

       ! ===== Swap Y --> Z; 1D FFTs in Z =====
       if (dims(1)>1) then
          call transpose_y_to_z(wk2_r2c,out_c,sp)
       else
          call transpose_y_to_z(wk13,out_c,sp)
       end if
       call c2c_1m_z(out_c,-1,sp)
                
    else if (format==PHYSICAL_IN_Z) then

       ! ===== 1D FFTs in Z =====
       call r2c_1m_z(in_r,wk13)

       ! ===== Swap Z --> Y; 1D FFTs in Y =====
       if (dims(1)>1) then
          call transpose_z_to_y(wk13,wk2_r2c,sp)
          call c2c_1m_y(wk2_r2c,-1,sp)
       else  ! out_c==wk2_r2c if 1D decomposition
          call transpose_z_to_y(wk13,out_c,sp)
          call c2c_1m_y(out_c,-1,sp)
       end if

       ! ===== Swap Y --> X; 1D FFTs in X =====
       if (dims(1)>1) then
          call transpose_y_to_x(wk2_r2c,out_c,sp)
       end if
       call c2c_1m_x(out_c,-1,sp)

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

    if (format==PHYSICAL_IN_X) then

       ! ===== 1D FFTs in Z =====
#ifdef OVERWRITE
       call c2c_1m_z(in_c,1,sp)       
#else
       allocate(wk1(sp%zsz(1),sp%zsz(2),sp%zsz(3)))
       wk1 = in_c
       call c2c_1m_z(wk1,1,sp)
#endif

       ! ===== Swap Z --> Y; 1D FFTs in Y =====
#ifdef OVERWRITE
       call transpose_z_to_y(in_c,wk2_r2c,sp)
#else
       call transpose_z_to_y(wk1,wk2_r2c,sp)
#endif
       call c2c_1m_y(wk2_r2c,1,sp)

       ! ===== Swap Y --> X; 1D FFTs in X =====
       if (dims(1)>1) then
          call transpose_y_to_x(wk2_r2c,wk13,sp)
          call c2r_1m_x(wk13,out_r)
       else
          call c2r_1m_x(wk2_r2c,out_r)
       end if

    else if (format==PHYSICAL_IN_Z) then

       ! ===== 1D FFTs in X =====
#ifdef OVERWRITE
       call c2c_1m_x(in_c,1,sp)
#else
       allocate(wk1(sp%xsz(1),sp%xsz(2),sp%xsz(3)))
       wk1 = in_c
       call c2c_1m_x(wk1,1,sp)
#endif

       ! ===== Swap X --> Y; 1D FFTs in Y =====
       if (dims(1)>1) then
#ifdef OVERWRITE
          call transpose_x_to_y(in_c,wk2_r2c,sp)
#else
          call transpose_x_to_y(wk1,wk2_r2c,sp)
#endif
          call c2c_1m_y(wk2_r2c,1,sp)
       else  ! in_c==wk2_r2c if 1D decomposition
#ifdef OVERWRITE
          call c2c_1m_y(in_c,1,sp)
#else
          call c2c_1m_y(wk1,1,sp)
#endif
       end if

       ! ===== Swap Y --> Z; 1D FFTs in Z =====
       if (dims(1)>1) then
          call transpose_y_to_z(wk2_r2c,wk13,sp)
       else
#ifdef OVERWRITE
          call transpose_y_to_z(in_c,wk13,sp)
#else
          call transpose_y_to_z(wk1,wk13,sp)
#endif
       end if
       call c2r_1m_z(wk13,out_r)

    end if

    return
  end subroutine fft_3d_c2r
