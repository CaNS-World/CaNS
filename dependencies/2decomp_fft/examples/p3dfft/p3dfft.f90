!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This test computes a 3D FFT using 2DECOMP&FFT and P3DFFT separately,
! crosschecks the result and compares the performance of the libraries.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program p3dfft_test

  use p3dfft
  use decomp_2d
  use decomp_2d_fft
  use MPI
  
  implicit none
  
  integer, parameter :: wp=KIND(0.0d0)
  integer, parameter :: nx=128, ny=128, nz=128
  ! use 0*0 for auto-tuning, particularly good for large core count
  integer, parameter :: p_row=0, p_col=0 
  integer, parameter :: ntest = 10
  
  real(wp), dimension(:,:,:), allocatable :: BEG,FIN
  complex(wp), dimension(:,:,:), allocatable :: AEND,AEND2
  
  real(wp) :: cdiff, ccdiff
  double precision :: rtime1, rtime2
  integer :: ierr, ndims(2), nsize, proc_id
  integer :: dims(2), dummy_coords(2)
  logical :: dummy_periods(2)
  integer :: istart(3), iend(3), isize(3)
  integer :: fstart(3), fend(3), fsize(3)
  integer :: i, m, x, y, z, imin
  
  call MPI_INIT (ierr)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 2DECOMP&FFT
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, &
       dims, dummy_periods, dummy_coords, ierr)
  
  if (nrank == 0) then
     print *,' '
     print *,'----------------------------------------------'
     print *,'2DECOMP&FFT result'
  endif
  
  call decomp_2d_fft_init
  
  allocate (BEG(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
  allocate (FIN(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
  call decomp_2d_fft_get_size(istart,iend,isize)
  allocate(AEND(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3)))
  allocate(AEND2(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3)))
  
  do z=xstart(3),xend(3)
     do y=xstart(2),xend(2)
        do x=xstart(1),xend(1)
           BEG(x,y,z) = real(x,wp)/real(nx,wp)*real(y,wp) &
                /real(ny,wp)*real(z,wp)/real(nz,wp)
           FIN(x,y,z) = BEG(x,y,z) ! keep a copy
        end do
     end do
  end do
  
  rtime1 = 0.0D0
  do m = 1, ntest
     
     BEG = FIN ! Note if OVERWRITE flag used BEG may be overwritten.
     ! So here restore exact input from the copy.
     
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     rtime1 = rtime1 - MPI_WTIME()
     call decomp_2d_fft_3d(BEG,AEND)
     rtime1 = rtime1 + MPI_WTIME()
     
     ! keep a copy of the result for comparison with P3DFFT
     if (m==1) AEND2=AEND
     
     if (nrank==0 .AND. m==1) then
        imin = min(isize(1),min(isize(2),isize(3)))
        write(*,*) 'Result of forward transform:'
        write(*,*) (AEND(i,i,i),i=1,min(10,imin))
     end if
     
     AEND = AEND / real(nx,wp) / real(ny,wp) / real(nz,wp)
     
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     rtime1 = rtime1 - MPI_WTIME()
     call decomp_2d_fft_3d(AEND,BEG)
     rtime1 = rtime1 + MPI_WTIME()
     
  end do
  
  cdiff = 0.0_wp
  ccdiff = 0.0_wp
  do z=xstart(3),xend(3)
     do y=xstart(2),xend(2)
        do x=xstart(1),xend(1)
           if (cdiff < abs(BEG(x,y,z)-FIN(x,y,z))) then
              cdiff = abs(BEG(x,y,z)-FIN(x,y,z))
           end if
        end do
     end do
  end do
  call MPI_Reduce(cdiff,ccdiff,1,real_type,MPI_MAX,0, &
       MPI_COMM_WORLD,ierr)
  if (nrank == 0) then
     write(*,*) 'after backward transform'
     write (*,*) ' max diff =',ccdiff
  end if
  
  cdiff = 0.0_wp
  ccdiff = 0.0_wp
  do z=xstart(3),xend(3)
     do y=xstart(2),xend(2)
        do x=xstart(1),xend(1)
           cdiff = cdiff + (BEG(x,y,z)-FIN(x,y,z))**2
        end do
     end do
  end do
  call MPI_Reduce(cdiff,ccdiff,1,real_type,MPI_SUM,0, &
       MPI_COMM_WORLD,ierr)
  ccdiff = sqrt(ccdiff/real(nx,wp)/real(ny,wp)/real(nz,wp))
  if (nrank == 0) then
     write (*,*) 'aver diff =',ccdiff
  end if
  
  call MPI_Reduce(rtime1,rtime2,1,MPI_REAL8,MPI_MAX,0, &
       MPI_COMM_WORLD,ierr)
  if (nrank==0) then
     write(*,*)'cpu time per forward+backward transforms', &
          rtime2/dble(ntest)
  end if
  
  ! save some information for P3DFFT
  proc_id = nrank
  ndims(1) = dims(1)
  ndims(2) = dims(2)
  
  call decomp_2d_fft_finalize
  call decomp_2d_finalize
  
  deallocate(BEG, FIN, AEND)
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! P3DFFT
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call MPI_COMM_SIZE(MPI_COMM_WORLD,nsize,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,proc_id,ierr)
  
  if (ndims(1)*ndims(2) /= nsize) then
     call MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
  end if
  
  if (proc_id == 0) then
     print *,' '
     print *,'----------------------------------------------'
     print *,'P3DFFT result'
  endif
  
  call p3dfft_setup(ndims,nx,ny,nz,MPI_COMM_WORLD,nx,ny,nz,.true.)
  call p3dfft_get_dims(istart,iend,isize,1)
  call p3dfft_get_dims(fstart,fend,fsize,2)
  
  allocate (BEG(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3)))
  allocate (AEND(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)))
  allocate (FIN(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3)))
  
  do z=istart(3),iend(3)
     do y=istart(2),iend(2)
        do x=istart(1),iend(1)
           BEG(x,y,z) = real(x,wp)/real(nx,wp)*real(y,wp) &
                /real(ny,wp)*real(z,wp)/real(nz,wp)
        end do
     end do
  end do
  !
  ! transform from physical space to wavenumber space
  ! (XgYiZj to XiYjZg)
  ! then transform back to physical space
  ! (XiYjZg to XgYiZj)
  !
  rtime1 = 0.0D0               
  do m=1, ntest
     
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     
     ! Forward transform
     rtime1 = rtime1 - MPI_wtime()
     ! The third parameter introduce in version 2.5.1 of P3DFFT
     ! API not compatible with earlier versions
     call p3dfft_ftran_r2c(BEG,AEND,'fft')
     rtime1 = rtime1 + MPI_wtime()
     
     ! check against 2DECOMP&FFT
     if (m==1) then
        cdiff = 0.0_wp
        do z=fstart(3),fend(3)
           do y=fstart(2),fend(2)
              do x=fstart(1),fend(1)
                 cdiff = cdiff + &
                      (real(AEND(x,y,z),wp)-real(AEND2(x,y,z),wp))**2 &
                      +(aimag(AEND(x,y,z))-aimag(AEND2(x,y,z)))**2
              end do
           end do
        end do
        call MPI_Reduce(cdiff,ccdiff,1,real_type,MPI_SUM,0, &
             MPI_COMM_WORLD,ierr)
        ccdiff = sqrt(ccdiff/real(nx/2+1,wp)/real(ny,wp)/real(nz,wp))
        if (proc_id==0) then
           write(*,*) ' '
           write(*,*) '==== Crosscheck (2DECOMP&FFT vs. P3DFFT) ===='
           write(*,*) 'aver diff = ', ccdiff
           write(*,*) ' '
        end if
     end if
     
     if (proc_id==0 .AND. m==1) then
        imin = min(fsize(1),min(fsize(2),fsize(3)))
        write(*,*) 'Result of forward transform:'
        write(*,*) (AEND(i,i,i),i=1,min(10,imin))
     end if
     
     AEND = AEND / real(nx,wp) / real(ny,wp) / real(nz,wp)
     
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     
     ! Backward transform
     rtime1 = rtime1 - MPI_wtime()
     call p3dfft_btran_c2r(AEND,FIN,'tff')       
     rtime1 = rtime1 + MPI_wtime()
     
  end do
  
  cdiff = 0.0_wp
  ccdiff = 0.0_wp
  do z=istart(3),iend(3)
     do y=istart(2),iend(2)
        do x=istart(1),iend(1)
           if (cdiff < abs(BEG(x,y,z)-FIN(x,y,z))) then
              cdiff = abs(BEG(x,y,z)-FIN(x,y,z))
           end if
        end do
     end do
  end do
  call MPI_Reduce(cdiff,ccdiff,1,real_type,MPI_MAX,0, &
       MPI_COMM_WORLD,ierr)
  if (proc_id == 0) then
     write(*,*) 'after backward transform'
     write (*,*) ' max diff =',ccdiff
  end if
  
  cdiff = 0.0_wp
  ccdiff = 0.0_wp
  do z=istart(3),iend(3)
     do y=istart(2),iend(2)
        do x=istart(1),iend(1)
           cdiff = cdiff + (BEG(x,y,z)-FIN(x,y,z))**2
        end do
     end do
  end do
  call MPI_Reduce(cdiff,ccdiff,1,real_type,MPI_SUM,0, &
       MPI_COMM_WORLD,ierr)
  ccdiff = sqrt(ccdiff/real(nx,wp)/real(ny,wp)/real(nz,wp))
  if (proc_id == 0) then
     write (*,*) 'aver diff =',ccdiff
  end if
  
  call MPI_Reduce(rtime1,rtime2,1,mpi_real8,MPI_MAX,0, &
       MPI_COMM_WORLD,ierr)
  
  if (proc_id == 0) then
     write(*,*)'cpu time per forward+backward transforms', &
          rtime2/dble(ntest)
  end if
  
  call p3dfft_clean
  deallocate(BEG, FIN, AEND, AEND2)
  
  
  call MPI_FINALIZE (ierr)
  
end program p3dfft_test
