program tecplot_view
  
  use decomp_2d
  use decomp_2d_io
  use MPI

  implicit none

  integer, parameter :: nx=17, ny=13, nz=11
  integer, parameter :: p_row=4, p_col=3

  real(mytype), dimension(nx,ny,nz) :: data1
  
  integer, dimension(3) :: lstart, lend, lsize

  integer :: i,j,k, m, ierror
  
  call MPI_INIT(ierror)
  call decomp_2d_init(nx,ny,nz,p_row,p_col)

  ! a copy of global data saved on every rank
  m = 1
  do k=1,nz
     do j=1,ny
        do i=1,nx
           data1(i,j,k) = real(m, mytype)
           m = m+1
        end do
     end do
  end do

  ! master rank generated a Tecplot view of the global data 
  if (nrank==0) then
     open(10,file='data0.dat',form='formatted')
     write(10,*) 'TITLE="Tecplot Output"'
     write(10,*) 'VARIABLES= "X" "Y" "Z" "VAR"'
     write(10,*) 'ZONE F=POINT I=',nx,' J=',ny,' ,K=',nz
     do k=1,nz
        do j=1,ny
           do i=1,nx
              write(10,*) i,j,k, data1(i,j,k)
           end do
        end do
     end do
     close(10)
  end if
  
  ! Generate Tecplot views of the decompositions
  ! -------------------------------------------- 
  ! For each pencil orientation there are two ways decomposing: 
  ! p_row*p_col or p_col*p_row. One set is used in 2DECOMP and is 
  ! described by the library's global variables. The other set is
  ! generated here for visualisation.
  
  ! (/ 1,2,3 /)
  call tecplot(nx, ny, nz, data1, xstart, xend, xsize, 'data1.dat')
  ! (/ 2,1,3 /)
  call tecplot(nx, ny, nz, data1, ystart, yend, ysize, 'data2.dat')
  ! (/ 2,3,1 /)
  call tecplot(nx, ny, nz, data1, zstart, zend, zsize, 'data3b.dat')

  call partition(nx, ny, nz, (/ 1,3,2 /), lstart, lend, lsize)
  call tecplot(nx, ny, nz, data1, lstart, lend, lsize, 'data1b.dat')
  call partition(nx, ny, nz, (/ 3,1,2 /), lstart, lend, lsize)
  call tecplot(nx, ny, nz, data1, lstart, lend, lsize, 'data2b.dat')
  call partition(nx, ny, nz, (/ 3,2,1 /), lstart, lend, lsize)
  call tecplot(nx, ny, nz, data1, lstart, lend, lsize, 'data3.dat')

  call decomp_2d_finalize
  call MPI_FINALIZE(ierror)

end program tecplot_view


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generate Tecplot files to visualise the 2D decompositions
!   - each rank corresponds to a Tecplot 'zone'
!   - rank 0 handles all I/O
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine tecplot(nx, ny, nz, data1, lstart, lend, lsize, filename)

  use decomp_2d
  use MPI

  implicit none

  integer, intent(IN) :: nx ,ny ,nz
  real(mytype), dimension(nx,ny,nz), intent(IN) :: data1
  integer, dimension(3), intent(IN) :: lstart, lend, lsize
  character(len=*), intent(IN) :: filename

  real(mytype), dimension(nx,ny,nz) :: data1b
  real(mytype), allocatable, dimension(:,:,:) :: local_data, rbuf
  integer, dimension(9) :: sbuf1, rbuf1
  integer, dimension(MPI_STATUS_SIZE) :: status
  character(len=7) :: tempstr
  integer :: i,j,k, m, ierror

  ! data1 holds the first copy of global data, generated locally
  ! data1b holds a second copy of global data, collected via communication

  ! each rank holds its local data
  allocate (local_data(lstart(1):lend(1),lstart(2):lend(2), &
       lstart(3):lend(3)))
  do k=lstart(3),lend(3)
     do j=lstart(2),lend(2)
        do i=lstart(1),lend(1)
           local_data(i,j,k) = data1(i,j,k)
        end do
     end do
  end do

  if (nrank==0) then

     ! master writes file header, collect data from each slave process,
     ! and write as seraprate Tecplot zones
     open(10,file=filename,form='formatted')
     write(10,*) 'TITLE="Tecplot Output"'
     write(10,*) 'VARIABLES= "X" "Y" "Z" "VAR"'
     write(10,*) 'ZONE F=POINT T="Rank 00" I=',lsize(1),' J=',lsize(2), &
          ' ,K=',lsize(3)
     do k=lstart(3),lend(3)
        do j=lstart(2),lend(2)
           do i=lstart(1),lend(1)
              write(10,*) i,j,k, local_data(i,j,k)
              ! master copies its local data to the second global array
              data1b(i,j,k)=local_data(i,j,k)
           end do
        end do
     end do

     ! loop through all other ranks to receive data and write
     do m=1,nproc-1
        CALL MPI_RECV(rbuf1,9,MPI_INTEGER,m,m,MPI_COMM_WORLD, &
             status,ierror)
        write(tempstr,100)'Rank ',m
100     format(A,I2.2)
        write(10,*) 'ZONE F=POINT T="', tempstr, '" I=',rbuf1(3), &
             ' J=',rbuf1(6), ' ,K=',rbuf1(9)
        allocate (rbuf(rbuf1(1):rbuf1(2),rbuf1(4):rbuf1(5), &
             rbuf1(7):rbuf1(8)))
        CALL MPI_RECV(rbuf,rbuf1(3)*rbuf1(6)*rbuf1(9),real_type,m, &
             m+nproc,MPI_COMM_WORLD,status,ierror)
        do k=rbuf1(7),rbuf1(8)
           do j=rbuf1(4),rbuf1(5)
              do i=rbuf1(1),rbuf1(2)
                 write(10,*) i,j,k, rbuf(i,j,k)
                 ! data received copied to global array
                 data1b(i,j,k)=rbuf(i,j,k)
              end do
           end do
        end do
        deallocate(rbuf)
     end do

     close (10)

     ! check if data set collected via communication is correct
     do k=1,nz
        do j=1,ny
           do i=1,nx
              if (abs(data1b(i,j,k)-data1(i,j,k)) > 1.0e-5) then
                 stop "error"
              end if
           end do
        end do
     end do
     
  else

     ! slaves send data to the master
     sbuf1(1) = lstart(1)
     sbuf1(2) = lend(1)
     sbuf1(3) = lsize(1)
     sbuf1(4) = lstart(2)
     sbuf1(5) = lend(2)
     sbuf1(6) = lsize(2)
     sbuf1(7) = lstart(3)
     sbuf1(8) = lend(3)
     sbuf1(9) = lsize(3)
     CALL MPI_SEND(sbuf1,9,MPI_INTEGER,0,nrank,MPI_COMM_WORLD,ierror)
     CALL MPI_SEND(local_data,lsize(1)*lsize(2)*lsize(3),real_type,0, &
          nrank+nproc,MPI_COMM_WORLD,ierror)

  endif

  deallocate(local_data)
  return

end subroutine tecplot
