module mod_initmpi
  use mpi
  use decomp_2d
  use mod_param     , only: dims
  use mod_common_mpi, only: myid,ijk_min,dims_xyz,ipencil,comm_cart,is_bound,halo,nb,ierr
  use mod_types
  implicit none
  private
  public initmpi
  contains
  subroutine initmpi(n,bc)
    implicit none
    integer, intent(in), dimension(3) :: n
    character(len=1), intent(in), dimension(0:1,3) :: bc
    integer :: ntx,nty,ntz
    logical, dimension(3) :: periods
    integer, dimension(3) :: coord
    !
    periods(:) = .false.
    if( bc(0,1)//bc(1,1).eq.'PP' ) periods(1) = .true.
    if( bc(0,2)//bc(1,2).eq.'PP' ) periods(2) = .true.
    if( bc(0,3)//bc(1,3).eq.'PP' ) periods(3) = .true.
    call decomp_2d_init(n(1),n(2),n(3),dims(1),dims(2),periods)
    myid = nrank
    !
    dims_xyz(1,1) = 1
    dims_xyz(2,1) = dims(1)
    dims_xyz(3,1) = dims(2)
    dims_xyz(1,2) = dims(1)
    dims_xyz(2,2) = 1
    dims_xyz(3,2) = dims(2)
    dims_xyz(1,3) = dims(1)
    dims_xyz(2,3) = dims(2)
    dims_xyz(3,3) = 1
#ifdef DECOMP_X
    dims(1:3) = dims_xyz(1:3,1)
    comm_cart = DECOMP_2D_COMM_CART_X
    coord(1) = (xstart(1)-1)*dims(1)/n(1)
    coord(2) = (xstart(2)-1)*dims(2)/n(2)
    coord(3) = (xstart(3)-1)*dims(3)/n(3)
    nb(0,1) = MPI_PROC_NULL
    nb(1,1) = MPI_PROC_NULL
    call MPI_CART_SHIFT(comm_cart,0,1,nb(0,2),nb(1,2),ierr)
    call MPI_CART_SHIFT(comm_cart,1,1,nb(0,3),nb(1,3),ierr)
    ipencil = 1
#elif DECOMP_Y
    dims(1:3) = dims_xyz(1:3,2)
    comm_cart = DECOMP_2D_COMM_CART_Y
    coord(1) = (ystart(1)-1)*dims(1)/n(1)
    coord(2) = (ystart(2)-1)*dims(2)/n(2)
    coord(3) = (ystart(3)-1)*dims(3)/n(3)
    call MPI_CART_SHIFT(comm_cart,0,1,nb(0,1),nb(1,1),ierr)
    nb(0,2) = MPI_PROC_NULL
    nb(1,2) = MPI_PROC_NULL
    call MPI_CART_SHIFT(comm_cart,1,1,nb(0,3),nb(1,3),ierr)
    ipencil = 2
#else
!#elif DECOMP_Z
    dims(1:3) = dims_xyz(1:3,3)
    comm_cart = DECOMP_2D_COMM_CART_Z
    coord(1) = (zstart(1)-1)*dims(1)/n(1)
    coord(2) = (zstart(2)-1)*dims(2)/n(2)
    coord(3) = (zstart(3)-1)*dims(3)/n(3)
    call MPI_CART_SHIFT(comm_cart,0,1,nb(0,1),nb(1,1),ierr)
    call MPI_CART_SHIFT(comm_cart,1,1,nb(0,2),nb(1,2),ierr)
    nb(0,3) = MPI_PROC_NULL
    nb(1,3) = MPI_PROC_NULL
    ipencil = 3
#endif
    ijk_min(:) = coord(:)*n(:)/dims(:)
    is_bound(:,:) = .false.
    !
    ! TODO: Define ng, n_x, n_y, n_z, n, ijk_start_x, ijk_start_y, ijk_start_z, ijk_start
    !
    !
    ! best make a loop below instead
    !
    if(nb(0,1).eq. MPI_PROC_NULL) is_bound(0,1) = .true.
    if(nb(1,1).eq. MPI_PROC_NULL) is_bound(1,1) = .true.
    if(nb(0,2).eq. MPI_PROC_NULL) is_bound(0,2) = .true. 
    if(nb(1,2).eq. MPI_PROC_NULL) is_bound(1,2) = .true.
    if(nb(0,3).eq. MPI_PROC_NULL) is_bound(0,3) = .true. 
    if(nb(1,3).eq. MPI_PROC_NULL) is_bound(1,3) = .true.
    !
    ntx = n(1)/dims(1)+2
    nty = n(2)/dims(2)+2
    ntz = n(3)/dims(3)+2
    !
    ! definitions of datatypes for velocity and pressure b.c.'s
    ! note: array(i,j,k) is basically a 1-dim array;
    !       k is the outer and i is the inner loop counter =>
    !         * for fixed i, (nty)*(ntz) blocks of 1 element,
    !           with (ntx) elements between start and end
    !         * for fixed j, (ntz) blocks of (ntx) elements,
    !           with (ntx)*(nty) elements between start and end
    !         * for fixed k, (1) blocks of (ntx*nty) elements,
    !           with (ntx)*(nty)*(ntz) elements between start and end
    !
    call MPI_TYPE_VECTOR(nty*ntz,1      ,ntx        ,MPI_REAL_RP,halo(1),ierr)
    call MPI_TYPE_VECTOR(ntz    ,ntx    ,ntx*nty    ,MPI_REAL_RP,halo(2),ierr)
    call MPI_TYPE_VECTOR(1      ,ntx*nty,ntx*nty*ntz,MPI_REAL_RP,halo(3),ierr)
    call MPI_TYPE_COMMIT(halo(1),ierr)
    call MPI_TYPE_COMMIT(halo(2),ierr)
    call MPI_TYPE_COMMIT(halo(3),ierr)
    return
  end subroutine initmpi
end module mod_initmpi
