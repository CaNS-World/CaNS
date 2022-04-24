program io_plane_test

  use decomp_2d
  use decomp_2d_io

  implicit none

  integer, parameter :: nx=17, ny=13, nz=11
  integer, parameter :: p_row=4, p_col=3

  real(mytype), dimension(nx,ny,nz) :: data1
  real(mytype), allocatable, dimension(:,:,:) :: u1, u2, u3

  real(mytype), allocatable, dimension(:,:,:) :: work
  
  integer :: i,j,k, m, ierror, iol

  call MPI_INIT(ierror)
  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  
  ! ***** global data *****
  m = 1
  do k=1,nz
     do j=1,ny
        do i=1,nx
           data1(i,j,k) = real(m,mytype)
           m = m+1
        end do
     end do
  end do

  allocate(u1(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
  allocate(u2(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
  allocate(u3(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))

  ! original X-pensil based data 
  do k=xstart(3),xend(3)
    do j=xstart(2),xend(2)
      do i=xstart(1),xend(1)
        u1(i,j,k) = data1(i,j,k)
      end do
    end do
  end do
  call decomp_2d_write_plane(1,u1,1,nx/2,'x_pencil-x_plane.dat')
  call decomp_2d_write_plane(1,u1,2,ny/2,'x_pencil-y_plane.dat')
  call decomp_2d_write_plane(1,u1,3,nz/2,'x_pencil-z_plane.dat')

  ! Y-pencil data
  call transpose_x_to_y(u1,u2)
  call decomp_2d_write_plane(2,u2,1,nx/2,'y_pencil-x_plane.dat')
  call decomp_2d_write_plane(2,u2,2,ny/2,'y_pencil-y_plane.dat')
  call decomp_2d_write_plane(2,u2,3,nz/2,'y_pencil-z_plane.dat')

  ! Z-pencil data
  call transpose_y_to_z(u2,u3)
  call decomp_2d_write_plane(3,u3,1,nx/2,'z_pencil-x_plane.dat')
  call decomp_2d_write_plane(3,u3,2,ny/2,'z_pencil-y_plane.dat')
  call decomp_2d_write_plane(3,u3,3,nz/2,'z_pencil-z_plane.dat')

  ! Attemp to read the files
  if (nrank==0) then
     inquire(iolength=iol) data1(1,1,1)

     ! X-plane
     allocate(work(1,ny,nz))
     open(10, FILE='x_pencil-x_plane.dat', FORM='unformatted', &
          ACCESS='DIRECT', RECL=iol)
     m=1
     do k=1,nz
        do j=1,ny
           read(10,rec=m) work(1,j,k)
           m=m+1
        end do
     end do
     write(*,*) ' '
     write(*,'(15I5)') int(work)
     close(10)
     deallocate(work)

     ! Y-plane
     allocate(work(nx,1,nz))
     open(10, FILE='x_pencil-y_plane.dat', FORM='unformatted', &
          ACCESS='DIRECT', RECL=iol)
     m=1
     do k=1,nz
        do i=1,nx
           read(10,rec=m) work(i,1,k)
           m=m+1
        end do
     end do
     write(*,*) ' '
     write(*,'(15I5)') int(work)
     close(10)
     deallocate(work)

     ! Z-plane
     allocate(work(nx,ny,1))
     open(10, FILE='x_pencil-z_plane.dat', FORM='unformatted', &
          ACCESS='DIRECT', RECL=iol)
     m=1
     do j=1,ny
        do i=1,nx
           read(10,rec=m) work(i,j,1)
           m=m+1
        end do
     end do
     write(*,*) ' '
     write(*,'(15I5)') int(work)
     close(10)
     deallocate(work)

  end if

  call decomp_2d_finalize 
  call MPI_FINALIZE(ierror)
  deallocate(u1,u2,u3)
  
end program io_plane_test
