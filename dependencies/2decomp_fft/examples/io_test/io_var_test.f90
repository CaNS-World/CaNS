! Sample application to test the read/write_var sets of routines
! in the IO library

program io_var_test

  use decomp_2d
  use decomp_2d_io
  use MPI

  implicit none

  integer, parameter :: nx=17, ny=13, nz=11
  integer :: p_row, p_col

  real(mytype), parameter :: eps = 1.0E-7

  ! for global data
  real(mytype), dimension(nx,ny,nz) :: data1
  real(mytype), dimension(nx*2,ny*2,nz*2) :: data1_large
  complex(mytype), dimension(nx,ny,nz) :: cdata1

  ! for distributed data
  real(mytype), allocatable, dimension(:,:,:) :: u1, u2, u3
  real(mytype), allocatable, dimension(:,:,:) :: u1l, u2l, u3l
  complex(mytype), allocatable, dimension(:,:,:) :: cu1, cu2, cu3

  ! another copy
  real(mytype), allocatable, dimension(:,:,:) :: u1_b, u2_b, u3_b
  real(mytype), allocatable, dimension(:,:,:) :: u1l_b, u2l_b, u3l_b
  complex(mytype), allocatable, dimension(:,:,:) :: cu1_b, cu2_b, cu3_b

  real(mytype), allocatable, dimension(:) :: tmp
  complex(mytype), allocatable, dimension(:) :: ctmp
  integer, allocatable, dimension(:) :: itmp

  TYPE(DECOMP_INFO) :: large

  integer :: i,j,k, m, ierror, fh
  character(len=15) :: filename, arg
  integer (kind=MPI_OFFSET_KIND) :: filesize, disp

  call MPI_INIT(ierror)

  i = command_argument_count()
  if (i/=2) then
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierror)
  else
     call get_command_argument(1, arg)
     read(arg, '(I10)') i
     p_row = i
     call get_command_argument(2, arg)
     read(arg, '(I10)') i
     p_col = i
  end if
  
  call decomp_2d_init(nx,ny,nz,p_row,p_col)

  ! also create a data set over a large domain
  call decomp_info_init(nx*2, ny*2, nz*2, large)

  ! initialise global data
  m = 1
  do k=1,nz
     do j=1,ny
        do i=1,nx
           data1(i,j,k) = real(m,mytype)
           cdata1(i,j,k) = cmplx(real(m,mytype),real(m,mytype))
           m = m+1
        end do
     end do
  end do

  m = 1
  do k=1,nz*2
     do j=1,ny*2
        do i=1,nx*2
           data1_large(i,j,k) = real(m,mytype)
           m = m+1
        end do
     end do
  end do

  ! allocate memory
  allocate(u1(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
  allocate(u2(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
  allocate(u3(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
  allocate(u1_b(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
  allocate(u2_b(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
  allocate(u3_b(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))

  allocate(cu1(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
  allocate(cu2(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
  allocate(cu3(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
  allocate(cu1_b(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
  allocate(cu2_b(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
  allocate(cu3_b(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))

  allocate(u1l(large%xst(1):large%xen(1), large%xst(2):large%xen(2), &
       large%xst(3):large%xen(3)))
  allocate(u2l(large%yst(1):large%yen(1), large%yst(2):large%yen(2), &
       large%yst(3):large%yen(3)))
  allocate(u3l(large%zst(1):large%zen(1), large%zst(2):large%zen(2), &
       large%zst(3):large%zen(3)))
  allocate(u1l_b(large%xst(1):large%xen(1), large%xst(2):large%xen(2), &
       large%xst(3):large%xen(3)))
  allocate(u2l_b(large%yst(1):large%yen(1), large%yst(2):large%yen(2), &
       large%yst(3):large%yen(3)))
  allocate(u3l_b(large%zst(1):large%zen(1), large%zst(2):large%zen(2), &
       large%zst(3):large%zen(3)))

  ! distribute the data
  do k=xstart(3),xend(3)
     do j=xstart(2),xend(2)
        do i=xstart(1),xend(1)
           u1(i,j,k) = data1(i,j,k)
           cu1(i,j,k) = cdata1(i,j,k)
        end do
     end do
  end do
  do k=large%xst(3),large%xen(3)
     do j=large%xst(2),large%xen(2)
        do i=large%xst(1),large%xen(1)
           u1l(i,j,k) = data1_large(i,j,k)
        end do
     end do
  end do

  ! transpose
  call transpose_x_to_y(u1,u2)
  call transpose_y_to_z(u2,u3)
  call transpose_x_to_y(u1l,u2l,large)
  call transpose_y_to_z(u2l,u3l,large)
  call transpose_x_to_y(cu1,cu2)
  call transpose_y_to_z(cu2,cu3)

  ! open file for IO
  write(filename,'(A,I3.3)') 'io_var_data.', nproc
  call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
       MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
       fh, ierror)
  filesize = 0_MPI_OFFSET_KIND
  call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
  disp = 0_MPI_OFFSET_KIND

  ! test writing scalar data
  allocate(tmp(2))
  tmp(1) = 1._mytype
  tmp(2) = 2._mytype
  allocate(ctmp(3))
  ctmp(1) = cmplx(1.0,1.0,mytype)
  ctmp(2) = cmplx(2.0,2.0,mytype)
  ctmp(3) = cmplx(3.0,3.0,mytype)
  allocate(itmp(3))
  call decomp_2d_write_scalar(fh,disp,2,tmp)
  call decomp_2d_write_scalar(fh,disp,3,ctmp)
  call decomp_2d_write_scalar(fh,disp,3,(/nx,ny,nz/))
  
  ! test the IO routines by writing all data to disk
  call decomp_2d_write_var(fh,disp,1,u1)
  call decomp_2d_write_var(fh,disp,2,u2)
  call decomp_2d_write_var(fh,disp,3,u3)
  call decomp_2d_write_var(fh,disp,1,u1l,large)
  call decomp_2d_write_var(fh,disp,2,u2l,large)
  call decomp_2d_write_var(fh,disp,3,u3l,large)
  call decomp_2d_write_var(fh,disp,1,cu1)
  call decomp_2d_write_var(fh,disp,2,cu2)
  call decomp_2d_write_var(fh,disp,3,cu3)

  call MPI_FILE_CLOSE(fh,ierror)

  if (nrank==0) write(*,*) 'disp=',disp

  ! read data back in from file
  call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
       MPI_MODE_RDONLY, MPI_INFO_NULL, &
       fh, ierror)
  ! skip the scalars (2 real, 3 cmplx, 3 int)
#ifdef DOUBLE_PREC
  ! if double precision: 2*8+3*8*2+3*4
  disp = 76_MPI_OFFSET_KIND
#else
  ! if single precision: 2*4+3*4*2+3*4
  disp = 44_MPI_OFFSET_KIND
#endif
  
  call decomp_2d_read_var(fh,disp,1,u1_b)
  call decomp_2d_read_var(fh,disp,2,u2_b)
  call decomp_2d_read_var(fh,disp,3,u3_b)
  call decomp_2d_read_var(fh,disp,1,u1l_b,large)
  call decomp_2d_read_var(fh,disp,2,u2l_b,large)
  call decomp_2d_read_var(fh,disp,3,u3l_b,large)
  call decomp_2d_read_var(fh,disp,1,cu1_b)
  call decomp_2d_read_var(fh,disp,2,cu2_b)
  call decomp_2d_read_var(fh,disp,3,cu3_b)

  disp = 0_MPI_OFFSET_KIND
  call decomp_2d_read_scalar(fh,disp,2,tmp)
  call decomp_2d_read_scalar(fh,disp,3,ctmp)
  call decomp_2d_read_scalar(fh,disp,3,itmp)
  if (nrank==0) then
     write(*,'(2F8.3)') tmp
     write(*,20) ctmp
20   format(3(:,'(',F5.2,',',F5.2,')'))
     write(*,'(A,3I5)') 'nx,ny,nz', itmp
  end if

  call MPI_FILE_CLOSE(fh,ierror)
  deallocate(tmp, ctmp, itmp)

  ! validate the data
  do k=xstart(3),xend(3)
     do j=xstart(2),xend(2)
        do i=xstart(1),xend(1)
           if (abs(u1(i,j,k)-u1_b(i,j,k)) > eps) stop 1
           if (abs(cu1(i,j,k)-cu1_b(i,j,k)) > eps) stop 2
        end do
     end do
  end do

  do k=ystart(3),yend(3)
     do j=ystart(2),yend(2)
        do i=ystart(1),yend(1)
           if (abs(u2(i,j,k)-u2_b(i,j,k)) > eps) stop 3
           if (abs(cu2(i,j,k)-cu2_b(i,j,k)) > eps) stop 4
        end do
     end do
  end do

  do k=zstart(3),zend(3)
     do j=zstart(2),zend(2)
        do i=zstart(1),zend(1)
           if (abs(u3(i,j,k)-u3_b(i,j,k)) > eps) stop 5
           if (abs(cu3(i,j,k)-cu3_b(i,j,k)) > eps) stop 6
        end do
     end do
  end do

  do k=large%xst(3),large%xen(3)
     do j=large%xst(2),large%xen(2)
        do i=large%xst(1),large%xen(1)
           if (abs(u1l(i,j,k)-u1l_b(i,j,k)) > eps) stop 7
        end do
     end do
  end do

  do k=large%yst(3),large%yen(3)
     do j=large%yst(2),large%yen(2)
        do i=large%yst(1),large%yen(1)
           if (abs(u2l(i,j,k)-u2l_b(i,j,k)) > eps) stop 8
        end do
     end do
  end do

  do k=large%zst(3),large%zen(3)
     do j=large%zst(2),large%zen(2)
        do i=large%zst(1),large%zen(1)
           if (abs(u3l(i,j,k)-u3l_b(i,j,k)) > eps) stop 9
        end do
     end do
  end do

  if (nrank==0) write(*,*) 'passed self test'

  ! clean up
  call decomp_info_finalize(large)
  call decomp_2d_finalize
  call MPI_FINALIZE(ierror)
  deallocate(u1,u2,u3)
  deallocate(u1l,u2l,u3l)
  deallocate(cu1,cu2,cu3)
  deallocate(u1_b,u2_b,u3_b)
  deallocate(u1l_b,u2l_b,u3l_b)
  deallocate(cu1_b,cu2_b,cu3_b)
  
end program io_var_test
