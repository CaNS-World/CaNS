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

! This is an experimental CUFFT implementation of the FFT library.
! The CPUs are still responsible for the data transpositions while
! only the 1D FFTs computations are passed to the NVidia GPUs.

module cuda_alloc	

  use iso_c_binding

  interface
     integer (C_INT) function cudaMallocHost(buffer, size)  &
          bind(C,name="cudaMallocHost")
       use iso_c_binding
       implicit none
       type (C_PTR) :: buffer
       integer (C_SIZE_T), value :: size
     end function cudaMallocHost

     integer (C_INT) function cudaFreeHost(buffer)  &
          bind(C,name="cudaFreeHost")
       use iso_c_binding
       implicit none
       type (C_PTR), value :: buffer
     end function cudaFreeHost
  end interface

end module cuda_alloc


module decomp_2d_fft
  
  use decomp_2d  ! 2D decomposition module
  use glassman
  
  implicit none
  
  private        ! Make everything private unless declared public

  ! engine-specific global variables
  complex(mytype), allocatable, dimension(:) :: buf, scratch

  ! common code used for all engines, including global variables, 
  ! generic interface definitions and several subroutines
#include "fft_common.f90"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time initialisations for the FFT engine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_fft_engine

    implicit none

    integer :: cbuf_size

    if (nrank==0) then
       write(*,*) ' '
       write(*,*) '***** Using the CUFFT engine *****'
       write(*,*) ' '
    end if

    cbuf_size = max(ph%xsz(1), ph%ysz(2))
    cbuf_size = max(cbuf_size, ph%zsz(3))
    allocate(buf(cbuf_size))
    allocate(scratch(cbuf_size))

    return
  end subroutine init_fft_engine


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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Multiple C2C 1D FFTs
  !   - assuming input/output are normal 3D arrays
  !   - CUFFT only supports multiple 1D FFTs with stride 1.
  !   - The 'type' parameter: 1=X FFT, 2=Y FFT, 3=Z FFT.
  !   - For type 2/3, this routine rearrange storage to ensure stride 1
  ! *** TO DO *** the memory-copying operations for stride-1 data is
  !   not necessary and the cost can be absorbed by the communication
  !   routines packing/unpacking the MPI ALLTOALLV buffers
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine c2c_1m(c1, c2, isign, type)
    
    use iso_c_binding
    use cuda_alloc
    
    implicit none

    complex(mytype), dimension(:,:,:), intent(IN) :: c1
    complex(mytype), dimension(:,:,:), intent(OUT) :: c2
    integer, intent(IN) :: isign, type
    
    integer :: Nx, Ny, Nz, res, i,j,k
    complex(mytype), allocatable, dimension(:,:,:) :: wk1, wk2
    complex(mytype), parameter :: size1 = dcmplx(0.0,0.0)
    type(C_PTR) :: cptr_c1d, cptr_c2d
    complex(mytype), dimension(:,:,:), pointer :: c1d, c2d

    Nx = size(c1,1)
    Ny = size(c1,2)
    Nz = size(c1,3)
    
    if (type == 1) then
       res = cudaMallocHost(cptr_c1d, Nx*Ny*Nz*c_sizeof(size1))
       call c_f_pointer(cptr_c1d, c1d, [Nx,Ny,Nz])
       res = cudaMallocHost(cptr_c2d, Nx*Ny*Nz*c_sizeof(size1))
       call c_f_pointer(cptr_c2d, c2d, [Nx,Ny,Nz])
       c1d = c1
       call fft_1m_c2c(Nx, Ny*Nz, c1d, c2d, isign)
       c2 = c2d
    else if (type == 2) then
       res = cudaMallocHost(cptr_c1d, Nx*Ny*Nz*c_sizeof(size1))
       call c_f_pointer(cptr_c1d, c1d, [Ny,Nx,Nz])
       res = cudaMallocHost(cptr_c2d, Nx*Ny*Nz*c_sizeof(size1))
       call c_f_pointer(cptr_c2d, c2d, [Ny,Nx,Nz])
       allocate(wk1(Ny, Nx, Nz))
       do k=1,Nz                                                               
          do j=1,Ny                                                            
             do i=1,Nx
                wk1(j,i,k) = c1(i,j,k)
             end do
          end do
       end do
       c1d = wk1
       call fft_1m_c2c(Ny, Nx*Nz, c1d, c2d, isign)
       wk1 = c2d
       do k=1,Nz
          do j=1,Ny
             do i=1,Nx
                c2(i,j,k) = wk1(j,i,k)
             end do
          end do
       end do
       deallocate(wk1)
    else if (type == 3) then
       res = cudaMallocHost(cptr_c1d, Nx*Ny*Nz*c_sizeof(size1))
       call c_f_pointer(cptr_c1d, c1d, [Nz,Nx,Ny])
       res = cudaMallocHost(cptr_c2d, Nx*Ny*Nz*c_sizeof(size1))
       call c_f_pointer(cptr_c2d, c2d, [Nz,Nx,Ny])
       allocate(wk2(Nz, Nx, Ny))
       do k=1,Nz
          do j=1,Ny
             do i=1,Nx
                wk2(k,i,j) = c1(i,j,k)
             end do
          end do
       end do
       c1d = wk2
       call fft_1m_c2c(Nz, Nx*Ny, c1d, c2d, isign)
       wk2 = c2d
       do k=1,Nz
          do j=1,Ny
             do i=1,Nx
                c2(i,j,k) = wk2(k,i,j)
             end do
          end do
       end do
       deallocate(wk2)
    end if

    res = cudaFreeHost(cptr_c1d)
    res = cudaFreeHost(cptr_c2d)
    
  end subroutine c2c_1m


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Multiple R2C 1D FFTs
  !   - assuming input/output are normal 3D arrays
  !   - and multiple 1D r2c applied on the first dimension
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine r2c_1m_x(r, c) 
    
    use iso_c_binding
    use cuda_alloc
    
    implicit none
    
    real(mytype), dimension(:,:,:), intent(IN) :: r 
    complex(mytype), dimension(:,:,:), intent(OUT) :: c
    
    integer :: Nx, Ny, Nz, res
    complex(mytype), parameter :: size1 = dcmplx(0.0,0.0)
    real(mytype), parameter :: size2 = dble(0.0)
    type(C_PTR) :: cptr_cd, cptr_rd
    complex(mytype), dimension(:,:,:), pointer :: cd
    real(mytype), dimension(:,:,:), pointer :: rd
    
    Nx = size(r,1)
    Ny = size(r,2)
    Nz = size(r,3)
    
    res = cudaMallocHost(cptr_cd,(Nx/2+1)*Ny*Nz*c_sizeof(size1))
    call c_f_pointer(cptr_cd,cd, [Nx/2+1,Ny,Nz])
    res = cudaMallocHost(cptr_rd, Nx*Ny*Nz*c_sizeof(size2))
    call c_f_pointer(cptr_rd,rd, [Nx,Ny,Nz])
    rd = r
    call fft_1m_r2c(Nx, Ny*Nz, rd, cd)
    c = cd
    res = cudaFreeHost(cptr_rd)
    res = cudaFreeHost(cptr_cd)
    
  end subroutine r2c_1m_x


  ! r2c transform, multiple 1D FFTs in z direction
  ! *** TO DO ***  doing this on host for the moment because the
  !                input/output is not stride-1
  subroutine r2c_1m_z(input, output)

    implicit none

    real(mytype), dimension(:,:,:), intent(IN)  ::  input
    complex(mytype), dimension(:,:,:), intent(OUT) :: output

    integer :: i,j,k, s1,s2,s3, d3

    s1 = size(input,1)
    s2 = size(input,2)
    s3 = size(input,3)
    d3 = size(output,3)

    do j=1,s2
       do i=1,s1
          ! Glassman's FFT is c2c only, 
          ! needing some pre- and post-processing for r2c
          ! pack real input in complex storage
          do k=1,s3
             buf(k) = cmplx(input(i,j,k),0._mytype, kind=mytype)
          end do
          call spcfft(buf,s3,-1,scratch)
          ! note d3 ~ s3/2+1
          ! simply drop the redundant part of the complex output
          do k=1,d3
             output(i,j,k) = buf(k)
          end do
       end do
    end do

    return

  end subroutine r2c_1m_z


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Multiple R2C 1D FFTs
  !   - assuming input/output are normal 3D arrays
  !   - and multiple 1D r2c applied on the first dimension
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine c2r_1m_x(c, r)
    
    use iso_c_binding
    use cuda_alloc
    
    implicit none

    complex(mytype), dimension(:,:,:), intent(IN) :: c
    real(mytype), dimension(:,:,:), intent(OUT) :: r

    integer :: Nx, Ny, Nz, res
    complex(mytype), parameter :: size1 = dcmplx(0.0,0.0)
    real(mytype), parameter :: size2 = dble(0.0)
    type(C_PTR) :: cptr_cd, cptr_rd
    complex(mytype), dimension(:,:,:), pointer :: cd
    real(mytype), dimension(:,:,:), pointer :: rd
    
    Nx = size(r,1)
    Ny = size(r,2)
    Nz = size(r,3)

    res = cudaMallocHost(cptr_cd,(Nx/2+1)*Ny*Nz*c_sizeof(size1))
    call c_f_pointer(cptr_cd,cd, [Nx/2+1,Ny,Nz])
    res = cudaMallocHost(cptr_rd, Nx*Ny*Nz*c_sizeof(size2))
    call c_f_pointer(cptr_rd,rd, [Nx,Ny,Nz])
    cd = c
    call fft_1m_c2r(Nx, Ny*Nz, cd, rd)
    r = rd
    res = cudaFreeHost(cptr_rd)
    res = cudaFreeHost(cptr_cd)
    
  end subroutine c2r_1m_x


  ! c2r transform, multiple 1D FFTs in z direction
  ! *** TO DO ***  doing this on host for the moment because the
  !                input/output is not stride-1
  subroutine c2r_1m_z(input, output)

    implicit none

    complex(mytype), dimension(:,:,:), intent(IN)  ::  input
    real(mytype), dimension(:,:,:), intent(OUT) :: output

    integer :: i,j,k, d1,d2,d3

    d1 = size(output,1)
    d2 = size(output,2)
    d3 = size(output,3)

    do j=1,d2
       do i=1,d1
          do k=1,d3/2+1
             buf(k) = input(i,j,k)
          end do
          do k=d3/2+2,d3
             buf(k) =  conjg(buf(d3+2-k))
          end do
          call spcfft(buf,d3,1,scratch)
          do k=1,d3
             output(i,j,k) = real(buf(k), kind=mytype)
          end do
       end do
    end do

    return

  end subroutine c2r_1m_z



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 3D FFT - complex to complex
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fft_3d_c2c(in, out, isign)
    
    implicit none
    
    complex(mytype), dimension(:,:,:), intent(IN) :: in
    complex(mytype), dimension(:,:,:), intent(OUT) :: out
    integer, intent(IN) :: isign   

    complex(mytype), allocatable, dimension(:,:,:) :: wk1,wk2,wk2b,wk3

    if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_FORWARD .OR.  &
         format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_BACKWARD) then
       
       ! ===== 1D FFTs in X =====
       allocate (wk1(ph%xsz(1),ph%xsz(2),ph%xsz(3)))
       call c2c_1m(in,wk1,isign,1)

       ! ===== Swap X --> Y =====
       allocate (wk2(ph%ysz(1),ph%ysz(2),ph%ysz(3)))
       call transpose_x_to_y(wk1,wk2,ph)
       
       ! ===== 1D FFTs in Y =====
       allocate (wk2b(ph%ysz(1),ph%ysz(2),ph%ysz(3))) 
       call c2c_1m(wk2,wk2b,isign,2)

       ! ===== Swap Y --> Z =====
       allocate (wk3(ph%zsz(1),ph%zsz(2),ph%zsz(3)))
       call transpose_y_to_z(wk2b,wk3,ph)

       ! ===== 1D FFTs in Z =====
       call c2c_1m(wk3,out,isign,3)

    else if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_BACKWARD &
         .OR. & 
         format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_FORWARD) then

       ! ===== 1D FFTs in Z =====
       allocate (wk1(ph%zsz(1),ph%zsz(2),ph%zsz(3)))
       call c2c_1m(in,wk1,isign,3)
       
       ! ===== Swap Z --> Y =====
       allocate (wk2(ph%ysz(1),ph%ysz(2),ph%ysz(3)))
       call transpose_z_to_y(wk1,wk2,ph)
       
       ! ===== 1D FFTs in Y =====
       allocate (wk2b(ph%ysz(1),ph%ysz(2),ph%ysz(3)))
       call c2c_1m(wk2,wk2b,isign,2)
       
       ! ===== Swap Y --> X =====
       allocate (wk3(ph%xsz(1),ph%xsz(2),ph%xsz(3)))
       call transpose_y_to_x(wk2b,wk3,ph)
       
       ! ===== 1D FFTs in X =====
       call c2c_1m(wk3,out,isign,1)
       
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

    if (format==PHYSICAL_IN_X) then

       ! ===== 1D FFTs in X =====
       allocate(wk1(sp%xsz(1),sp%xsz(2),sp%xsz(3))) 
       call r2c_1m_x(in_r,wk1)

       ! ===== Swap X --> Y =====
       allocate (wk2(sp%ysz(1),sp%ysz(2),sp%ysz(3)))
       call transpose_x_to_y(wk1,wk2,sp)

       ! ===== 1D FFTs in Y =====
       allocate (wk2b(sp%ysz(1),sp%ysz(2),sp%ysz(3)))
       call c2c_1m(wk2,wk2b,-1,2)

       ! ===== Swap Y --> Z =====
       allocate (wk3(sp%zsz(1),sp%zsz(2),sp%zsz(3)))
       call transpose_y_to_z(wk2b,wk3,sp)

       ! ===== 1D FFTs in Z =====
       call c2c_1m(wk3,out_c,-1,3)
                
    else if (format==PHYSICAL_IN_Z) then

       ! ===== 1D FFTs in Z =====
       allocate(wk1(sp%zsz(1),sp%zsz(2),sp%zsz(3)))
       call r2c_1m_z(in_r,wk1)

       ! ===== Swap Z --> Y =====
       allocate (wk2(sp%ysz(1),sp%ysz(2),sp%ysz(3)))
       call transpose_z_to_y(wk1,wk2,sp)

       ! ===== 1D FFTs in Y =====
       allocate (wk2b(sp%ysz(1),sp%ysz(2),sp%ysz(3)))
       call c2c_1m(wk2,wk2b,-1,2)

       ! ===== Swap Y --> X =====
       allocate (wk3(sp%xsz(1),sp%xsz(2),sp%xsz(3)))
       call transpose_y_to_x(wk2b,wk3,sp)

       ! ===== 1D FFTs in X =====
       call c2c_1m(wk3,out_c,-1,1)

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

    if (format==PHYSICAL_IN_X) then

       ! ===== 1D FFTs in Z =====
       allocate (wk1(sp%zsz(1),sp%zsz(2),sp%zsz(3)))
       call c2c_1m(in_c,wk1,1,3)

       ! ===== Swap Z --> Y =====
       allocate (wk2(sp%ysz(1),sp%ysz(2),sp%ysz(3)))
       call transpose_z_to_y(wk1,wk2,sp)

       ! ===== 1D FFTs in Y =====
       allocate (wk2b(sp%ysz(1),sp%ysz(2),sp%ysz(3)))
       call c2c_1m(wk2,wk2b,1,2)

       ! ===== Swap Y --> X =====
       allocate (wk3(sp%xsz(1),sp%xsz(2),sp%xsz(3)))
       call transpose_y_to_x(wk2b,wk3,sp)

       ! ===== 1D FFTs in X =====
       call c2r_1m_x(wk3,out_r)

    else if (format==PHYSICAL_IN_Z) then

       ! ===== 1D FFTs in X =====
       allocate(wk1(sp%xsz(1),sp%xsz(2),sp%xsz(3)))
       call c2c_1m(in_c,wk1,1,1)

       ! ===== Swap X --> Y =====
       allocate (wk2(sp%ysz(1),sp%ysz(2),sp%ysz(3)))
       call transpose_x_to_y(wk1,wk2,sp)

       ! ===== 1D FFTs in Y =====
       allocate (wk2b(sp%ysz(1),sp%ysz(2),sp%ysz(3)))
       call c2c_1m(wk2,wk2b,1,2)

       ! ===== Swap Y --> Z =====
       allocate (wk3(sp%zsz(1),sp%zsz(2),sp%zsz(3)))
       call transpose_y_to_z(wk2b,wk3,sp)

       ! ===== 1D FFTs in Z =====
       call c2r_1m_z(wk3,out_r)

    end if

    return
  end subroutine fft_3d_c2r

  
end module decomp_2d_fft
