!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This example calculates the divergency of a random field using
!   (1) global transposition
!   (2) halo-cell exchange
! The two method should give identical results
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program halo_test

  use decomp_2d
  use decomp_2d_io

  implicit none

  integer, parameter :: nx=171, ny=132, nz=113
  integer, parameter :: p_row=4, p_col=3

  real(mytype), allocatable, dimension(:,:,:) :: u1, u2, u3
  real(mytype), allocatable, dimension(:,:,:) :: v1, v2, v3
  real(mytype), allocatable, dimension(:,:,:) :: w1, w2, w3
  real(mytype), allocatable, dimension(:,:,:) :: wk2, wk3
  real(mytype), allocatable, dimension(:,:,:) :: uh, vh, wh
  real(mytype), allocatable, dimension(:,:,:) :: div1, div2, div3, div4

  integer :: i,j,k, ierror, n

  integer, allocatable, dimension(:) :: seed

  call MPI_INIT(ierror)
  call decomp_2d_init(nx,ny,nz,p_row,p_col)

  ! initialise u,v,w with random numbers in X-pencil
#ifdef HALO_GLOBAL
  allocate(u1(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
  allocate(v1(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
  allocate(w1(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
#else
  allocate(u1(xsize(1), xsize(2), xsize(3)))
  allocate(v1(xsize(1), xsize(2), xsize(3)))
  allocate(w1(xsize(1), xsize(2), xsize(3)))
#endif

  call random_seed(size = n)
  allocate(seed(n))
  seed = nrank
  call random_seed(put=seed)
  call random_number(u1)
  call random_number(v1)
  call random_number(w1)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calculate divergency using global transposition
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! du/dx calculated on X-pencil
#ifdef HALO_GLOBAL
  allocate(div1(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
#else
  allocate(div1(xsize(1), xsize(2), xsize(3)))
#endif

  div1 = 0.

#ifdef HALO_GLOBAL
  do k=xstart(3),xend(3)
     do j=xstart(2),xend(2)
#else
  do k=1,xsize(3)
     do j=1,xsize(2)
#endif     
        do i=2,xsize(1)-1
           div1(i,j,k) = u1(i+1,j,k)-u1(i-1,j,k)
        end do
     end do
  end do

  ! dv/dy calculated on Y-pencil
#ifdef HALO_GLOBAL
  allocate(v2(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
  allocate(wk2(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
#else
  allocate(v2(ysize(1), ysize(2), ysize(3)))
  allocate(wk2(ysize(1), ysize(2), ysize(3)))
#endif
  call transpose_x_to_y(v1,v2)
  call transpose_x_to_y(div1,wk2)
#ifdef HALO_GLOBAL
  do k=ystart(3),yend(3)
     do i=ystart(1),yend(1)
#else
  do k=1,ysize(3)
     do i=1,ysize(1)
#endif
        do j=2,ysize(2)-1
           wk2(i,j,k) = wk2(i,j,k) + v2(i,j+1,k)-v2(i,j-1,k)
        end do
     end do
  end do

  ! dw/dz calculated on Z-pencil
#ifdef HALO_GLOBAL
  allocate(w2(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
  allocate(w3(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
  allocate(wk3(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
#else
  allocate(w2(ysize(1), ysize(2), ysize(3)))
  allocate(w3(zsize(1), zsize(2), zsize(3)))
  allocate(wk3(zsize(1), zsize(2), zsize(3)))
#endif
  call transpose_x_to_y(w1,w2)
  call transpose_y_to_z(w2,w3)
  call transpose_y_to_z(wk2,wk3)
#ifdef HALO_GLOBAL
  do j=zstart(2),zend(2)
     do i=zstart(1),zend(1)
#else
  do j=1,zsize(2)
     do i=1,zsize(1)
#endif
        do k=2,zsize(3)-1
           wk3(i,j,k) = wk3(i,j,k) + w3(i,j,k+1)-w3(i,j,k-1)
        end do
     end do
  end do
  
  ! result in X-pencil
  call transpose_z_to_y(wk3,wk2)
  call transpose_y_to_x(wk2,div1)

  if (nrank==0) then
     write(*,*) 'Calculated via global transposition'
     write(*,*) (div1(i,i,i), i=2,13)
  end if
#ifdef HALO_GLOBAL
  call decomp_2d_write_one(1,div1,'div1g.dat')
#else
  call decomp_2d_write_one(1,div1,'div1.dat')
#endif

  deallocate(v2,w2,w3,wk2,wk3)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calculate divergency using halo-cell exchange (data in X-pencil)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef HALO_GLOBAL
  allocate(div2(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
#else
  allocate(div2(xsize(1), xsize(2), xsize(3)))
#endif
  ! du/dx
#ifdef HALO_GLOBAL
  do k=xstart(3),xend(3)
     do j=xstart(2),xend(2)
#else
  do k=1,xsize(3)
     do j=1,xsize(2)
#endif
        do i=2,xsize(1)-1
           div2(i,j,k) = u1(i+1,j,k)-u1(i-1,j,k) 
        end do
     end do
  end do

  ! dv/dy
#ifdef HALO_GLOBAL
  call update_halo(v1,vh,1,opt_global=.true.)
#else
  call update_halo(v1,vh,1)
#endif
#ifdef HALO_GLOBAL
  do k=xstart(3),xend(3)
     do j=xstart(2),xend(2)
#else
  do k=1,xsize(3)
     do j=1,xsize(2)
#endif
        do i=2,xsize(1)-1
           div2(i,j,k) = div2(i,j,k) + vh(i,j+1,k)-vh(i,j-1,k)
        end do
     end do
  end do

  ! dw/dz
#ifdef HALO_GLOBAL
  call update_halo(w1,wh,1,opt_global=.true.)
#else
  call update_halo(w1,wh,1)
#endif
#ifdef HALO_GLOBAL
  do k=xstart(3),xend(3)
     do j=xstart(2),xend(2)
#else
  do k=1,xsize(3)
     do j=1,xsize(2)
#endif
        do i=2,xsize(1)-1
           div2(i,j,k) = div2(i,j,k) + wh(i,j,k+1)-wh(i,j,k-1)
        end do
     end do
  end do

  if (nrank==0) then
     write(*,*) '-----------------------------------------------'
     write(*,*) 'Calculated via halo exchange (data in X-pencil)'
     write(*,*) (div2(i,i,i), i=2,13)
  end if
#ifdef HALO_GLOBAL
  call decomp_2d_write_one(1,div1,'div2g.dat')
#else
  call decomp_2d_write_one(1,div1,'div2.dat')
#endif
           
  deallocate(vh,wh)

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calculate divergency using halo-cell exchange (data in Y-pencil)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
#ifdef HALO_GLOBAL
  allocate(u2(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
  allocate(v2(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
  allocate(w2(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
#else
  allocate(u2(ysize(1), ysize(2), ysize(3)))
  allocate(v2(ysize(1), ysize(2), ysize(3)))
  allocate(w2(ysize(1), ysize(2), ysize(3)))
#endif
  call transpose_x_to_y(u1,u2)
  call transpose_x_to_y(v1,v2)
  call transpose_x_to_y(w1,w2)

#ifdef HALO_GLOBAL
  allocate(div3(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
  allocate(wk2(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
#else
  allocate(div3(xsize(1), xsize(2), xsize(3)))
  allocate(wk2(ysize(1), ysize(2), ysize(3)))
#endif

  ! du/dx
#ifdef HALO_GLOBAL
  call update_halo(u2,uh,1,opt_global=.true.)
#else
  call update_halo(u2,uh,1)
#endif
#ifdef HALO_GLOBAL
  do k=ystart(3),yend(3)
     do i=ystart(1),yend(1)
#else
  do k=1,ysize(3)
     do i=1,ysize(1)
#endif
        do j=2,ysize(2)-1
           wk2(i,j,k) = uh(i+1,j,k)-uh(i-1,j,k)
        end do
     end do
  end do

  ! dv/dy
#ifdef HALO_GLOBAL
  do k=ystart(3),yend(3)
     do i=ystart(1),yend(1)
#else
  do k=1,ysize(3)
     do i=1,ysize(1)
#endif
        do j=2,ysize(2)-1
           wk2(i,j,k) = wk2(i,j,k) + v2(i,j+1,k)-v2(i,j-1,k)
        end do
     end do
  end do

  ! dw/dz
#ifdef HALO_GLOBAL
  call update_halo(w2,wh,1,opt_global=.true.)
#else
  call update_halo(w2,wh,1)
#endif
#ifdef HALO_GLOBAL
  do k=ystart(3),yend(3)
     do i=ystart(1),yend(1)
#else
  do k=1,ysize(3)
     do i=1,ysize(1)
#endif
        do j=2,ysize(2)-1
           wk2(i,j,k) = wk2(i,j,k) + wh(i,j,k+1)-wh(i,j,k-1)
        end do
     end do
  end do

  call transpose_y_to_x(wk2,div3)

  if (nrank==0) then
     write(*,*) '-----------------------------------------------'
     write(*,*) 'Calculated via halo exchange (data in Y-pencil)'
     write(*,*) (div3(i,i,i), i=2,13)
  end if
#ifdef HALO_GLOBAL
  call decomp_2d_write_one(1,div1,'div3g.dat')
#else
  call decomp_2d_write_one(1,div1,'div3.dat')
#endif

  deallocate(uh,wh,wk2)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calculate divergency using halo-cell exchange (data in Z-pencil)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef HALO_GLOBAL
  allocate(u3(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
  allocate(v3(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
  allocate(w3(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
#else
  allocate(u3(zsize(1), zsize(2), zsize(3)))
  allocate(v3(zsize(1), zsize(2), zsize(3)))
  allocate(w3(zsize(1), zsize(2), zsize(3)))
#endif
  call transpose_y_to_z(u2,u3)
  call transpose_y_to_z(v2,v3)
  call transpose_y_to_z(w2,w3)

#ifdef HALO_GLOBAL
  allocate(div4(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
  allocate(wk2(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
  allocate(wk3(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
#else
  allocate(div4(xsize(1), xsize(2), xsize(3)))
  allocate(wk2(ysize(1), ysize(2), ysize(3)))
  allocate(wk3(zsize(1), zsize(2), zsize(3)))
#endif

  ! du/dx
#ifdef HALO_GLOBAL
  call update_halo(u3,uh,1,opt_global=.true.)
#else
  call update_halo(u3,uh,1)
#endif
#ifdef HALO_GLOBAL
  do j=zstart(2),zend(2)
     do i=zstart(1),zend(1)
#else
  do j=1,zsize(2)
     do i=1,zsize(1)
#endif
        do k=2,zsize(3)-1
           wk3(i,j,k) = uh(i+1,j,k)-uh(i-1,j,k)
        end do
     end do
  end do

  ! dv/dy
#ifdef HALO_GLOBAL
  call update_halo(v3,vh,1,opt_global=.true.)
#else
  call update_halo(v3,vh,1)
#endif
#ifdef HALO_GLOBAL
  do j=zstart(2),zend(2)
     do i=zstart(1),zend(1)
#else
  do j=1,zsize(2)
     do i=1,zsize(1)
#endif
        do k=2,zsize(3)-1
           wk3(i,j,k) = wk3(i,j,k) + vh(i,j+1,k)-vh(i,j-1,k)
        end do
     end do
  end do

  ! dw/dz
#ifdef HALO_GLOBAL
  do j=zstart(2),zend(2)
     do i=zstart(1),zend(1)
#else
  do j=1,zsize(2)
     do i=1,zsize(1)
#endif
        do k=2,zsize(3)-1
           wk3(i,j,k) = wk3(i,j,k) + w3(i,j,k+1)-w3(i,j,k-1)
        end do
     end do
  end do

  call transpose_z_to_y(wk3,wk2)
  call transpose_y_to_x(wk2,div4)

  if (nrank==0) then
     write(*,*) '-----------------------------------------------'
     write(*,*) 'Calculated via halo exchange (data in Z-pencil)'
     write(*,*) (div4(i,i,i), i=2,13)
  end if
#ifdef HALO_GLOBAL
  call decomp_2d_write_one(1,div1,'div4g.dat')
#else
  call decomp_2d_write_one(1,div1,'div4.dat')
#endif

  deallocate(uh,vh,wk2,wk3,u1,v1,w1,u2,v2,w2,u3,v3,w3)

  deallocate(div1,div2,div3,div4)


  call decomp_2d_finalize 
  call MPI_FINALIZE(ierror)

end program halo_test
