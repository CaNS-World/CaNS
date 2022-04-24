program fft_timing

  use decomp_2d
  use decomp_2d_fft
  use MPI
  
  implicit none
  
  integer, parameter :: nx=17, ny=13, nz=11
  integer, parameter :: p_row=0, p_col=0
  
  integer, parameter :: NTEST = 10  ! repeat test this times
  
  complex(mytype), allocatable, dimension(:,:,:) :: in, out
  real(mytype), allocatable, dimension(:,:,:) :: in_r
  
  integer, dimension(3) :: fft_start, fft_end, fft_size
  
  real(mytype) :: dr,di, err, err_all, n1,flops
  integer :: ierror, i,j,k,m
  double precision :: t1, t2, t3 ,t4
  
  call MPI_INIT(ierror)
  call decomp_2d_init(nx,ny,nz,p_row,p_col)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test the c2c interface
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call decomp_2d_fft_init(PHYSICAL_IN_Z) ! non-default Z-pencil input

  !  input is Z-pencil data
  ! output is X-pencil data
  allocate (in(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
  allocate (out(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
  ! initilise input
  do k=zstart(3),zend(3)
     do j=zstart(2),zend(2)
        do i=zstart(1),zend(1)
           dr = real(i,mytype)/real(nx,mytype)*real(j,mytype) &
                /real(ny,mytype)*real(k,mytype)/real(nz,mytype)
           di = dr
           in(i,j,k) = cmplx(dr,di,mytype)
        end do
     end do
  end do
  
  t2 = 0.0D0
  t4 = 0.0D0
  do m=1,NTEST
     
     ! forward FFT
     t1 = MPI_WTIME()
     call decomp_2d_fft_3d(in, out, DECOMP_2D_FFT_FORWARD)
     t2 = t2 + MPI_WTIME() - t1

     ! inverse FFT
     t3 = MPI_WTIME()
     call decomp_2d_fft_3d(out, in, DECOMP_2D_FFT_BACKWARD)
     t4 = t4 + MPI_WTIME() - t3
  
     ! normalisation - note 2DECOMP&FFT doesn't normalise
     in = in / real(nx,mytype) / real(ny,mytype) /real(nz,mytype)

  end do
  
  call MPI_ALLREDUCE(t2,t1,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
       MPI_COMM_WORLD,ierror)
  t1 = t1 / real(nproc,mytype)
  call MPI_ALLREDUCE(t4,t3,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
       MPI_COMM_WORLD,ierror)
  t3 = t3 / real(nproc,mytype)
  
  ! checking accuracy
  err = 0.
  do k=zstart(3),zend(3)
     do j=zstart(2),zend(2)
        do i=zstart(1),zend(1)
           dr = real(i,mytype)/real(nx,mytype)*real(j,mytype) &
                /real(ny,mytype)*real(k,mytype)/real(nz,mytype)
           di = dr
           dr = dr - real(in(i,j,k),mytype)
           di = di - aimag(in(i,j,k))
           err = err + sqrt(dr*dr + di*di)
        end do
     end do
  end do
  call MPI_ALLREDUCE(err,err_all,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
  err_all = err_all / real(nx,mytype) / real(ny,mytype) / real(nz,mytype)

  if (nrank==0) then
     write(*,*) '===== c2c interface ====='
     write(*,*) 'error / mesh point: ', err_all 
     write(*,*) 'time (sec): ', t1,t3
     n1 = real(nx,mytype) * real(ny,mytype) * real(nz,mytype)
     n1 = n1 ** (1._mytype/3._mytype)
     ! 5n*log(n) flops per 1D FFT of size n using Cooley-Tukey algorithm
     flops = 5._mytype * n1 * log(n1) / log(2.0_mytype)
     ! 3 sets of 1D FFTs for 3 directions, each having n^2 1D FFTs
     flops = flops * 3._mytype * n1**2  
     flops = 2._mytype * flops / ((t1+t3)/real(NTEST,mytype))
     write(*,*) 'GFLOPS : ', flops / 1000._mytype**3
  end if
  
  deallocate(in,out)
  call decomp_2d_fft_finalize


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test the r2c/c2r interface
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call decomp_2d_fft_init
  
  allocate (in_r(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
  call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)
  allocate (out(fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))
  
  ! initilise input
  do k=xstart(3),xend(3)
     do j=xstart(2),xend(2)
        do i=xstart(1),xend(1)
           in_r(i,j,k) = real(i,mytype)/real(nx,mytype)*real(j,mytype) &
                /real(ny,mytype)*real(k,mytype)/real(nz,mytype)
        end do
     end do
  end do
  
  t2 = 0.0D0
  t4 = 0.0D0
  do m=1,NTEST
  
     ! 3D r2c FFT
     t1 = MPI_WTIME()
     call decomp_2d_fft_3d(in_r, out)
     t2 = t2 + MPI_WTIME() - t1
  
     ! 3D inverse FFT
     t3 = MPI_WTIME()
     call decomp_2d_fft_3d(out, in_r)
     t4 = t4 + MPI_WTIME() - t3
  
     ! normalisation - note 2DECOMP&FFT doesn't normalise
     do k=xstart(3),xend(3)
        do j=xstart(2),xend(2)
           do i=xstart(1),xend(1)
              in_r(i,j,k) = in_r(i,j,k) &
                   / (real(nx,mytype)*real(ny,mytype)*real(nz,mytype))
           end do
        end do
     end do

  end do
  
  call MPI_ALLREDUCE(t2,t1,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
       MPI_COMM_WORLD,ierror)
  t1 = t1 / real(nproc,mytype)
  call MPI_ALLREDUCE(t4,t3,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
       MPI_COMM_WORLD,ierror)
  t3 = t3 / real(nproc,mytype)
  
  ! checking accuracy
  err = 0.
  do k=xstart(3),xend(3)
     do j=xstart(2),xend(2)
        do i=xstart(1),xend(1)
           dr = real(i,mytype)/real(nx,mytype)*real(j,mytype) &
                /real(ny,mytype)*real(k,mytype)/real(nz,mytype)
           err = err + abs(in_r(i,j,k)-dr)
        end do
     end do
  end do
  call MPI_ALLREDUCE(err,err_all,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
  err_all = err_all / real(nx,mytype) / real(ny,mytype) / real(nz,mytype)
  
  if (nrank==0) then
     write(*,*) '===== r2c/c2r interface ====='
     write(*,*) 'error / mesh point: ', err_all
     write(*,*) 'time (sec): ', t1,t3
  end if
  
  deallocate(in_r,out)
  call decomp_2d_fft_finalize
  call decomp_2d_finalize
  call MPI_FINALIZE(ierror)

end program fft_timing

