!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This program computes multiple distributed 3D FFTs using the
! non-blocking communication API provided by 2DECOMP&FFT to demonstrate
! the overlap of computation and communication
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program non_blocking

  use decomp_2d
  use MPI

  implicit none

  include "fftw3.f"

  integer, parameter :: nx=16, ny=16, nz=16
  integer, parameter :: p_row=2, p_col=2
  integer, parameter :: NFFT=20   ! number of independent FFTs
  integer, parameter :: NTEST=5   ! number of MPI_TEST per 1D FFT

  integer :: i,j,k, m,n, nmax, ierror
  real(mytype) :: tmp1, tmp2

  double precision :: t1, t2

  integer, dimension(NFFT) :: handles

  ! FFTW plans for the 1D forward/backward transforms
  integer*8, save :: x_plan_f, x_plan_b
  integer*8, save :: y_plan_f, y_plan_b
  integer*8, save :: z_plan_f, z_plan_b

  ! dummy array used for planning
  complex(mytype), allocatable, dimension(:) :: buf1, buf2

  ! input/output of the FFT
  complex(mytype), allocatable, dimension(:,:,:,:) :: in, out
  ! intermediate Y-pencil storage
  complex(mytype), allocatable, dimension(:,:,:,:) :: wk2
  ! send/recv buffers for non-blocking ALLTOALLV
  complex(mytype), allocatable, dimension(:,:,:,:) :: sbuf, rbuf


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


  allocate( in(xsize(1),xsize(2),xsize(3),NFFT))   ! x-pencil input
  allocate(out(zsize(1),zsize(2),zsize(3),NFFT))   ! z-pencil output
  allocate(wk2(ysize(1),ysize(2),ysize(3),NFFT))   ! y-pencil intermediate

  ! send/recv buffers
  allocate(sbuf(ysize(1),ysize(2),ysize(3),0:1))   ! y-pencil
  allocate(rbuf(zsize(1),zsize(2),zsize(3),0:1))   ! z-pencil

  ! 1D temp buffer
  nmax = max(xsize(1),max(ysize(2),zsize(3)))
  allocate (buf1(nmax))
  allocate (buf2(nmax))

  t1 = MPI_WTIME()

  ! first data set: V_1
  do k=1,xsize(3)
     do j=1,xsize(2)
        do i=1,xsize(1)
           tmp1 = real(xstart(1)+i-1, mytype) / real(nx, mytype) &
                * real(xstart(2)+j-1, mytype) / real(ny, mytype) &
                * real(xstart(3)+k-1, mytype) / real(nz, mytype) &
                * 1._mytype / real(NFFT, mytype)
           in(i,j,k,1) = cmplx(tmp1, 0._mytype, mytype)
        end do
     end do
  end do
  
  ! 1D FFT in X for V_1
  do k=1,xsize(3)
     do j=1,xsize(2)
        do i=1,xsize(1)
           buf1(i) = in(i,j,k,1)
        end do
#ifdef DOUBLE_PREC
        call dfftw_execute_dft(x_plan_f, buf1, buf2)
#else
        call sfftw_execute_dft(x_plan_f, buf1, buf2)
#endif
        do i=1,xsize(1)
           in(i,j,k,1) = buf2(i)
        end do
     end do
  end do

  ! transpose X->Y for V_1 (blocking)
  call transpose_x_to_y(in(:,:,:,1),wk2(:,:,:,1))

  ! 1D FFT in Y for V_1
  do k=1,ysize(3)
     do i=1,ysize(1)
        do j=1,ysize(2)
           buf1(j) = wk2(i,j,k,1)
        end do
#ifdef DOUBLE_PREC
        call dfftw_execute_dft(y_plan_f, buf1, buf2)
#else
        call sfftw_execute_dft(y_plan_f, buf1, buf2)
#endif
        do j=1,ysize(2)
           wk2(i,j,k,1) = buf2(j)
        end do
     end do
  end do

  ! initiate transpose Y->Z for V_1
  call transpose_y_to_z_start(handles(1),wk2(:,:,:,1),out(:,:,:,1), &
       sbuf(:,:,:,1),rbuf(:,:,:,1))
  
  ! Loop through remaining data sets
  do m=2,NFFT

     ! new data set: V_m
     do k=1,xsize(3)
        do j=1,xsize(2)
           do i=1,xsize(1)
              tmp1 = real(xstart(1)+i-1, mytype) / real(nx, mytype) &
                   * real(xstart(2)+j-1, mytype) / real(ny, mytype) &
                   * real(xstart(3)+k-1, mytype) / real(nz, mytype) &
                   * real(m, mytype) / real(NFFT, mytype)
              in(i,j,k,m) = cmplx(tmp1, 0._mytype, mytype)
           end do
        end do
     end do

     ! 1D FFT in X for V_m
     do k=1,xsize(3)
        do j=1,xsize(2)
           do i=1,xsize(1)
              buf1(i) = in(i,j,k,m)
           end do
#ifdef DOUBLE_PREC
           call dfftw_execute_dft(x_plan_f, buf1, buf2)
#else
           call sfftw_execute_dft(x_plan_f, buf1, buf2)
#endif
           do n=1,NTEST
              call transpose_test(handles(m-1))
           end do
           do i=1,xsize(1)
              in(i,j,k,m) = buf2(i)
           end do
        end do
     end do

     ! transpose X->Y for V_m (blocking)
     call transpose_x_to_y(in(:,:,:,m),wk2(:,:,:,m))

     ! 1D FFT in Y for V_m
     do k=1,ysize(3)
        do i=1,ysize(1)
           do j=1,ysize(2)
              buf1(j) = wk2(i,j,k,m)
           end do
#ifdef DOUBLE_PREC
           call dfftw_execute_dft(y_plan_f, buf1, buf2)
#else
           call sfftw_execute_dft(y_plan_f, buf1, buf2)
#endif
           do n=1,NTEST
              call transpose_test(handles(m-1))
           end do
           do j=1,ysize(2)
              wk2(i,j,k,m) = buf2(j)
           end do
        end do
     end do

     ! NOTE: there are two transposes ongoing at any time
     ! The mod operation is to use different sets of send/recv buffers

     ! initiate transpose Y->Z for V_m
     call transpose_y_to_z_start(handles(m),wk2(:,:,:,m),out(:,:,:,m), &
          sbuf(:,:,:,mod(m,2)),rbuf(:,:,:,mod(m,2)))

     ! wait for transpose complete for V_m-1
     call transpose_y_to_z_wait(handles(m-1),wk2(:,:,:,m-1),out(:,:,:,m-1), &
          sbuf(:,:,:,mod(m-1,2)),rbuf(:,:,:,mod(m-1,2)))

     ! 1D FFT in Z for V_m-1
     do j=1,zsize(2)
        do i=1,zsize(1)
           do k=1,zsize(3)
              buf1(k) = out(i,j,k,m-1)
           end do
#ifdef DOUBLE_PREC
           call dfftw_execute_dft(z_plan_f, buf1, buf2)
#else
           call sfftw_execute_dft(z_plan_f, buf1, buf2)
#endif
           do k=1,zsize(3)
              out(i,j,k,m-1) = buf2(k)
           end do
        end do
     end do

     if (nrank==0) write(*,*) 'TEST ', m-1, out(1:2,1:2,1:2,m-1)

  end do ! NTEST Loop

  ! wait for transpose complete for V_NFFT
  call transpose_y_to_z_wait(handles(NFFT),wk2(:,:,:,NFFT),out(:,:,:,NFFT), &
       sbuf(:,:,:,mod(NFFT,2)),rbuf(:,:,:,mod(NFFT,2)))
  
  ! 1D FFT in Z for V_NFFT
  do j=1,zsize(2)
     do i=1,zsize(1)
        do k=1,zsize(3)
           buf1(k) = out(i,j,k,NFFT)
        end do
#ifdef DOUBLE_PREC
        call dfftw_execute_dft(z_plan_f, buf1, buf2)
#else
        call sfftw_execute_dft(z_plan_f, buf1, buf2)
#endif
        do k=1,zsize(3)
           out(i,j,k,NFFT) = buf2(k)
        end do
     end do
  end do

  if (nrank==0) write(*,*) 'TEST ', NFFT, out(1:2,1:2,1:2,NFFT)
  
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
  deallocate(in,out,wk2,buf1,buf2,sbuf,rbuf)
  

end program non_blocking
