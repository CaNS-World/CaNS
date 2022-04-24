program io_read

  use decomp_2d
  use decomp_2d_io

  implicit none

  integer, parameter :: nx=17, ny=13, nz=11
  ! use different number of processes
  integer, parameter :: p_row=3, p_col=2

#ifdef COMPLEX_TEST
  complex(mytype), dimension(nx,ny,nz) :: data1

  complex(mytype), allocatable, dimension(:,:,:) :: u1b, u2b, u3b
#else
  real(mytype), dimension(nx,ny,nz) :: data1

  real(mytype), allocatable, dimension(:,:,:) :: u1b, u2b, u3b
#endif

  real(mytype), parameter :: eps = 1.0E-7_mytype
 
  integer :: i,j,k, m, ierror
  
  call MPI_INIT(ierror)
  call decomp_2d_init(nx,ny,nz,p_row,p_col)

  ! ***** global data *****
  m = 1
  do k=1,nz
     do j=1,ny
        do i=1,nx
#ifdef COMPLEX_TEST
           data1(i,j,k) = cmplx(real(m,mytype), real(nx*ny*nz-m,mytype)) 
#else
           data1(i,j,k) = real(m,mytype)
#endif
           m = m+1
        end do
     end do
  end do

  allocate(u1b(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
  allocate(u2b(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
  allocate(u3b(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))

  ! read back to different arrays
  call decomp_2d_read_one(1,u1b,'u1.dat')
  call decomp_2d_read_one(2,u2b,'u2.dat')
  call decomp_2d_read_one(3,u3b,'u3.dat')

  ! Check against the global data array
  do k=xstart(3),xend(3)
    do j=xstart(2),xend(2)
      do i=xstart(1),xend(1)
        if (abs((data1(i,j,k)-u1b(i,j,k))) > eps) stop 4
      end do
    end do
  end do

  do k=ystart(3),yend(3)
    do j=ystart(2),yend(2)
      do i=ystart(1),yend(1)
        if (abs((data1(i,j,k)-u2b(i,j,k))) > eps) stop 5
      end do
    end do
  end do
  
  do k=zstart(3),zend(3)
    do j=zstart(2),zend(2)
      do i=zstart(1),zend(1)
        if (abs((data1(i,j,k)-u3b(i,j,k))) > eps) stop 6
      end do
    end do
  end do

  call decomp_2d_finalize 
  call MPI_FINALIZE(ierror)
  deallocate(u1b,u2b,u3b)

end program io_read
