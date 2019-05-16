module mod_initmpi
  use mpi
  use decomp_2d
  use mod_param     , only: dims
  use mod_common_mpi, only: myid,coord,comm_cart,left,right,front,back,xhalo,yhalo,ierr
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
    !
    periods(:) = .false.
    if( bc(0,1)//bc(1,1).eq.'PP' ) periods(1) = .true.
    if( bc(0,2)//bc(1,2).eq.'PP' ) periods(2) = .true.
    if( bc(0,3)//bc(1,3).eq.'PP' ) periods(3) = .true.
    call decomp_2d_init(n(1),n(2),n(3),dims(1),dims(2),periods)
    myid = nrank
    comm_cart = DECOMP_2D_COMM_CART_Z
    coord(1) = (zstart(1)-1)*dims(1)/n(1)
    coord(2) = (zstart(2)-1)*dims(2)/n(2)
    !
    call MPI_CART_SHIFT(comm_cart,0,1,left,right,ierr)
    call MPI_CART_SHIFT(comm_cart,1,1,front,back,ierr)
    !
    ntx = n(1)/dims(1)+2
    nty = n(2)/dims(2)+2
    ntz = n(3)+2
    !
    ! Definitions of datatypes for velocity and pressure b.c.'s
    ! Note: array(i,j,k) is basically a 1-dim array;
    !       k is the outer and i is the inner loop counter =>
    !         * for fixed i, (j1+1)*(k1+1) blocks of 1 element,
    !           with (i1+1) elements between start and end
    !         * for fixed j, (k1+1) blocks of (i1+1) elements,
    !           with (i1+1)*(j1+1) elements between start and end
    !
    call MPI_TYPE_VECTOR(nty*ntz,1  ,ntx    ,MPI_REAL_RP,xhalo,ierr)
    call MPI_TYPE_VECTOR(ntz    ,ntx,ntx*nty,MPI_REAL_RP,yhalo,ierr)
    call MPI_TYPE_COMMIT(xhalo,ierr)
    call MPI_TYPE_COMMIT(yhalo,ierr)
    return
  end subroutine initmpi
end module mod_initmpi
