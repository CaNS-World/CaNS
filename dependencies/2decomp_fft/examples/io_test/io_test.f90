program io_test

  use decomp_2d
  use decomp_2d_io

  implicit none

  integer, parameter :: nx=17, ny=13, nz=11
  integer, parameter :: p_row=4, p_col=3

#ifdef COMPLEX_TEST
  complex(mytype), dimension(nx,ny,nz) :: data1

  complex(mytype), allocatable, dimension(:,:,:) :: u1, u2, u3
  complex(mytype), allocatable, dimension(:,:,:) :: u1b, u2b, u3b
#else
  real(mytype), dimension(nx,ny,nz) :: data1

  real(mytype), allocatable, dimension(:,:,:) :: u1, u2, u3
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

  allocate(u1(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
  allocate(u2(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
  allocate(u3(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))

  allocate(u1b(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
  allocate(u2b(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
  allocate(u3b(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))

  ! original x-pensil based data 
  do k=xstart(3),xend(3)
    do j=xstart(2),xend(2)
      do i=xstart(1),xend(1)
        u1(i,j,k) = data1(i,j,k)
      end do
    end do
  end do

  ! transpose
  call transpose_x_to_y(u1,u2)
  call transpose_y_to_z(u2,u3)

  ! write to disk
  call decomp_2d_write_one(1,u1,'u1.dat')
  call decomp_2d_write_one(2,u2,'u2.dat')
  call decomp_2d_write_one(3,u3,'u3.dat')

  ! read back to different arrays
  call decomp_2d_read_one(1,u1b,'u1.dat')
  call decomp_2d_read_one(2,u2b,'u2.dat')
  call decomp_2d_read_one(3,u3b,'u3.dat')

  ! compare  
  do k=xstart(3),xend(3)
    do j=xstart(2),xend(2)
      do i=xstart(1),xend(1)
        if (abs((u1(i,j,k)-u1b(i,j,k))) > eps) stop 1
      end do
    end do
  end do

  do k=ystart(3),yend(3)
    do j=ystart(2),yend(2)
      do i=ystart(1),yend(1)
        if (abs((u2(i,j,k)-u2b(i,j,k))) > eps) stop 2
      end do
    end do
  end do

  do k=zstart(3),zend(3)
    do j=zstart(2),zend(2)
      do i=zstart(1),zend(1)
        if (abs((u3(i,j,k)-u3b(i,j,k))) > eps) stop 3
      end do
    end do
  end do

  ! Also check against the global data array
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
  deallocate(u1,u2,u3)
  deallocate(u1b,u2b,u3b)

end program io_test
