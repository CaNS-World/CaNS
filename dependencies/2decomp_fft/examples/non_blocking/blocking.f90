!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This program computes multiple distributed 3D FFTs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program blocking

  use decomp_2d
  use MPI

  implicit none

  include "fftw3.f"

  integer, parameter :: nx=16, ny=16, nz=16
  integer, parameter :: p_row=2, p_col=2
  integer, parameter :: NFFT=20   ! number of independent FFTs

  integer :: i,j,k, m, nmax, ierror
  real(mytype) :: tmp1, tmp2

  double precision :: t1, t2

  ! FFTW plans for the 1D forward/backward transforms
  integer*8, save :: x_plan_f, x_plan_b
  integer*8, save :: y_plan_f, y_plan_b
  integer*8, save :: z_plan_f, z_plan_b

  ! dummy array used for planning
  complex(mytype), allocatable, dimension(:) :: buf1, buf2

  ! input/output of the FFT
  complex(mytype), allocatable, dimension(:,:,:) :: in, out, wk2


  call MPI_INIT(ierror)
  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  
  ! ===== planning 1D FFT in x =====
  allocate(buf1(xsize(1)), buf2(xsize(1)))

#ifdef DOUBLE_PREC
  call dfftw_plan_dft_1d(x_plan_f, xsize(1), buf1, buf2,  &
       FFTW_FORWARD,  FFTW_MEASURE)
  call dfftw_plan_dft_1d(x_plan_b, xsize(1), buf1, buf2,  &
       FFTW_BACKWARD, FFTW_MEASURE)
#else
  call sfftw_plan_dft_1d(x_plan_f, xsize(1), buf1, buf2,  &
       FFTW_FORWARD,  FFTW_MEASURE)
  call sfftw_plan_dft_1d(x_plan_b, xsize(1), buf1, buf2,  &
       FFTW_BACKWARD, FFTW_MEASURE)
#endif

  deallocate(buf1,buf2)

  ! ===== planning 1D FFT in Y =====
  allocate(buf1(ysize(2)), buf2(ysize(2)))
  
#ifdef DOUBLE_PREC
  call dfftw_plan_dft_1d(y_plan_f, ysize(2), buf1, buf2,  &
       FFTW_FORWARD,  FFTW_MEASURE)
  call dfftw_plan_dft_1d(y_plan_b, ysize(2), buf1, buf2,  &
       FFTW_BACKWARD, FFTW_MEASURE)
#else
  call sfftw_plan_dft_1d(y_plan_f, ysize(2), buf1, buf2,  &
       FFTW_FORWARD,  FFTW_MEASURE)
  call sfftw_plan_dft_1d(y_plan_b, ysize(2), buf1, buf2,  &
       FFTW_BACKWARD, FFTW_MEASURE)
#endif
  
  deallocate(buf1,buf2)
  
  ! ===== planning 1D FFT in Z =====
  allocate(buf1(zsize(3)), buf2(zsize(3)))
  
#ifdef DOUBLE_PREC
  call dfftw_plan_dft_1d(z_plan_f, zsize(3), buf1, buf2,  &
       FFTW_FORWARD,  FFTW_MEASURE)
  call dfftw_plan_dft_1d(z_plan_b, zsize(3), buf1, buf2,  &
       FFTW_BACKWARD, FFTW_MEASURE)
#else
  call sfftw_plan_dft_1d(z_plan_f, zsize(3), buf1, buf2,  &
       FFTW_FORWARD,  FFTW_MEASURE)
  call sfftw_plan_dft_1d(z_plan_b, zsize(3), buf1, buf2,  &
       FFTW_BACKWARD, FFTW_MEASURE)
#endif
  
  deallocate(buf1,buf2)


  allocate( in(xsize(1),xsize(2),xsize(3)))   ! x-pencil input
  allocate(out(zsize(1),zsize(2),zsize(3)))   ! z-pencil output
  allocate(wk2(ysize(1),ysize(2),ysize(3)))   ! y-pencil intermediate

  ! 1D temp buffer
  nmax = max(xsize(1),max(ysize(2),zsize(3)))
  allocate (buf1(nmax))
  allocate (buf2(nmax))

  t1 = MPI_WTIME()

  do m=1,NFFT

     do k=1,xsize(3)
        do j=1,xsize(2)
           do i=1,xsize(1)
              tmp1 = real(xstart(1)+i-1, mytype) / real(nx, mytype) &
                   * real(xstart(2)+j-1, mytype) / real(ny, mytype) &
                   * real(xstart(3)+k-1, mytype) / real(nz, mytype) &
                   * real(m, mytype) / real(NFFT, mytype)
              in(i,j,k) = cmplx(tmp1, 0._mytype, mytype)
           end do
        end do
     end do

     ! This shows how to perform 3D FFT by using the FFTW basic interface.
     ! Copy data to/from 1D buffers and loop through all 1D FFTs.

     ! 1D FFT in X
     do k=1,xsize(3)
        do j=1,xsize(2)
           do i=1,xsize(1)
              buf1(i) = in(i,j,k)
           end do
#ifdef DOUBLE_PREC
           call dfftw_execute_dft(x_plan_f, buf1, buf2)
#else
           call sfftw_execute_dft(x_plan_f, buf1, buf2)
#endif
           do i=1,xsize(1)
              in(i,j,k) = buf2(i)
           end do
        end do
     end do

     ! ===== Swap X --> Y =====
     call transpose_x_to_y(in,wk2)

     ! ===== 1D FFTs in Y =====
     do k=1,ysize(3)
        do i=1,ysize(1)
           do j=1,ysize(2)
              buf1(j) = wk2(i,j,k)
           end do
#ifdef DOUBLE_PREC
           call dfftw_execute_dft(y_plan_f, buf1, buf2)
#else
           call sfftw_execute_dft(y_plan_f, buf1, buf2)
#endif
           do j=1,ysize(2)
              wk2(i,j,k) = buf2(j)
           end do
        end do
     end do

     ! ===== Swap Y --> Z =====
     call transpose_y_to_z(wk2,out)

     ! ===== 1D FFTs in Z =====
     do j=1,zsize(2)
        do i=1,zsize(1)
           do k=1,zsize(3)
              buf1(k) = out(i,j,k)
           end do
#ifdef DOUBLE_PREC
           call dfftw_execute_dft(z_plan_f, buf1, buf2)
#else
           call sfftw_execute_dft(z_plan_f, buf1, buf2)
#endif
           do k=1,zsize(3)
              out(i,j,k) = buf2(k)
           end do
        end do
     end do

     if (nrank==0) write(*,*) 'TEST ', m, out(1:2,1:2,1:2)

  end do ! NFFT

  t2 = MPI_WTIME() - t1
  call MPI_ALLREDUCE(t2,t1,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
       MPI_COMM_WORLD,ierror)
  t1 = t1 / real(nproc, mytype)

  if (nrank==0) then
     write(*,*) 'Average Forward FFT Time(sec): ', t1
  end if

  ! clean up
#ifdef DOUBLE_PREC
  call dfftw_destroy_plan(x_plan_f)
  call dfftw_destroy_plan(x_plan_b)
  call dfftw_destroy_plan(y_plan_f)
  call dfftw_destroy_plan(y_plan_b)
  call dfftw_destroy_plan(z_plan_f)
  call dfftw_destroy_plan(z_plan_b)
#else
  call sfftw_destroy_plan(x_plan_f)
  call sfftw_destroy_plan(x_plan_b)
  call sfftw_destroy_plan(y_plan_f)
  call sfftw_destroy_plan(y_plan_b)
  call sfftw_destroy_plan(z_plan_f)
  call sfftw_destroy_plan(z_plan_b)
#endif

  call decomp_2d_finalize 
  call MPI_FINALIZE(ierror)
  deallocate(in,out,wk2,buf1,buf2)
  

end program blocking
