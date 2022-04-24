!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main test program for the FFT r2c/c2r interface
!  also demonstrate the use of the IO library 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program fft_test_r2c

  use decomp_2d
  use decomp_2d_fft
  use glassman
  use decomp_2d_io
  
  use MPI
  
  implicit none
  !include "fftw3.f"
  
  integer, parameter :: nx=4, ny=2, nz=3
  integer, parameter :: p_row=2, p_col=2
  
  real(mytype), allocatable, dimension(:,:,:) :: in, in2
  complex(mytype), allocatable, dimension(:,:,:) :: out
  
  integer, dimension(3) :: fft_start, fft_end, fft_size
  
  real(mytype), dimension(nx,ny,nz) :: in_global, in_g2, in_g3
  complex(mytype), dimension(nx/2+1,ny,nz) :: out_global
  
  integer (kind=MPI_OFFSET_KIND) :: filesize, disp
  
  real(mytype) :: err
  !integer*8 :: plan
  integer :: fh, ierror, i,j,k, n,iol
  
  call MPI_INIT(ierror)
  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  call decomp_2d_fft_init
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Compute a small problem all on rank 0 as reference
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call random_number(in_global)
  
  if (nrank==0) then
     write(*,*) '*** Reference serial computation on rank 0 only'
     write(*,*) ' global real input'
     do i=1,nx
        write(*,20) ((in_global(i,j,k),j=1,ny),k=1,nz)
     end do
     
     ! Using a 3D FFT routine supplied by this library
     call glassman_3d_r2c(in_global,nx,ny,nz,out_global)
     
     ! If using FFTW library:
     !  - uncomment the FFTW include file & plan above
     !  - uncomment the follwing function calls
     !  - change names to dfftw... for double precision 
     !call sfftw_plan_dft_r2c_3d(plan,nx,ny,nz, &
     !     in_global,out_global,FFTW_ESTIMATE)
     !call sfftw_execute_dft_r2c(plan,in_global,out_global)
     
     write(*,*) ' global complex output'
     do i=1,nx/2+1
        write(*,10) ((out_global(i,j,k),j=1,ny),k=1,nz)
     end do
  end if
10 format(1x,6(:,'(',F5.2,',',F5.2,')'))
20 format(1x,6F5.2)
  
  ! File for testing IO
  call MPI_FILE_OPEN(MPI_COMM_WORLD, 'fftdata', &
       MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
       fh, ierror)
  filesize = 0_MPI_OFFSET_KIND
  call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
  disp = 0_MPI_OFFSET_KIND
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test the real-to-complex interface (r2c) 
  
  !  input is X-pencil real    data whose global size is nx*ny*nz
  ! output is Z-pencil complex data whose global size is (nx/2+1)*ny*nz
  allocate (in(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
  
  call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)
  allocate (out(fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))
  
  ! each processor gets its local portion of global data
  do k=xstart(3),xend(3)
     do j=xstart(2),xend(2)
        do i=xstart(1),xend(1)
           in(i,j,k) = in_global(i,j,k)
        end do
     end do
  end do
  
  ! write input to file
  call decomp_2d_write_var(fh,disp,1,in)
  
  if (nrank==0) then
     write(*,*) ' '
     write(*,*) '*** Distributed computation (X-pencil input)'
     write(*,*) ' real input held by rank 0:'
     write(*,20) in
  end if
  
  ! compute r2c transform 
  call decomp_2d_fft_3d(in,out)
  
  if (nrank==0) then
     write(*,*) ' - after forward transform'
     write(*,*) ' complex output held by rank 0:'
     write(*,10) out
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test the complex-to-real interface (c2r) 

  allocate (in2(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
  
  ! compute c2r transform
  call decomp_2d_fft_3d(out,in2)
  
  ! normalisation
  in2 = in2 / real(nx) / real(ny) / real(nz)
  
  ! write the data recovered by inverse FFT to file
  call decomp_2d_write_var(fh,disp,1,in2)
  
  if (nrank==0) then
     write(*,*) ' - after backward transform and normalisation'
     write(*,*) ' real output held by rank 0:'
     write(*,20) in2
  end if
  
  deallocate(in,in2,out)
  call decomp_2d_fft_finalize
  
  call MPI_FILE_CLOSE(fh,ierror)
  
  ! check on rank 0 if input data is properly recovered
  ! this also tests the IO routines
  if (nrank==0) then
     in_g2(1,1,1) = real(0., mytype)
     inquire(iolength=iol) in_g2(1,1,1)
     OPEN(10, FILE='fftdata', FORM='unformatted', &
          ACCESS='DIRECT', RECL=iol)
     n=1
     do k=1,nz
        do j=1,ny
           do i=1,nx
              read(10,rec=n) in_g2(i,j,k)
              n=n+1
           end do
        end do
     end do
     do k=1,nz
        do j=1,ny
           do i=1,nx
              read(10,rec=n) in_g3(i,j,k)
              n=n+1
           end do
        end do
     end do
     err = 0._mytype
     do k=1,nz
        do j=1,ny
           do i=1,nx
              err = err + (in_g2(i,j,k)-in_g3(i,j,k))**2
           end do
        end do
     end do
     err = err / real(nx,mytype) / real(ny,mytype) / real(nz,mytype)
     write(*,*) ' error / mesh point: ', sqrt(err)
  end if


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Repeat the above but using Z-pencil input
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call decomp_2d_fft_init(PHYSICAL_IN_Z)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test the real-to-complex interface (r2c) 
  
  allocate (in(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
  call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)
  allocate (out(fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))
  
  ! each processor gets its local portion of global data
  do k=zstart(3),zend(3)
     do j=zstart(2),zend(2)
        do i=zstart(1),zend(1)
           in(i,j,k) = in_global(i,j,k)
        end do
     end do
  end do
  if (nrank==0) then
     write(*,*) ' '
     write(*,*) '*** Distributed computation (Z-pencil input)'
     write(*,*) ' real input held by rank 0:'
     write(*,20) in
  end if
  
  ! compute r2c transform 
  call decomp_2d_fft_3d(in,out)
  
  if (nrank==0) then
     write(*,*) ' - after forward transform'
     write(*,*) ' complex output held by rank 0:'
     write(*,10) out
  end if
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test the complex-to-real interface (c2r) 
  
  allocate (in2(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
  
  ! compute c2r transform
  call decomp_2d_fft_3d(out,in2)
  
  ! normalisation
  in2 = in2 / real(nx) / real(ny) / real(nz)
  
  if (nrank==0) then
     write(*,*) ' - after backward transform and normalisation'
     write(*,*) ' real output held by rank 0:'
     write(*,20) in2
  end if
  
  deallocate(in,in2,out)
  
  call decomp_2d_fft_finalize
  call decomp_2d_finalize
  call MPI_FINALIZE(ierror)
  
end program fft_test_r2c
