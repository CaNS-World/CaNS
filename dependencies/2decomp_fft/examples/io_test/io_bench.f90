program io_bench

  use decomp_2d
  use decomp_2d_io
  use MPI

  implicit none
  
  integer, parameter :: nx=100, ny=100, nz=100
  integer, parameter :: p_row=4, p_col=4

  real(mytype), allocatable, dimension(:,:,:) :: u1
  
  double precision :: t1, t2
  integer :: ierror

  call MPI_INIT(ierror)
  call decomp_2d_init(nx,ny,nz,p_row,p_col)

  allocate(u1(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
  call random_number(u1)

  t1 = MPI_WTIME()
  call decomp_2d_write_one(1,u1,'io.dat')
  t2 = MPI_WTIME()

  if (nrank==0) write(*,*) 'I/O time: ', t2-t1

  call decomp_2d_finalize 
  call MPI_FINALIZE(ierror)
  deallocate(u1)

end program io_bench
  
