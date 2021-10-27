module mod_initmpi
  use mpi
  use decomp_2d
  use mod_param     , only: dims
  use mod_common_mpi, only: comm_cart,myid,ierr,halo,ipencil
  use mod_types
  implicit none
  private
  public initmpi
  contains
  subroutine initmpi(n,bc,n_z,lo,hi,nb,is_bound)
    implicit none
    integer, intent(in), dimension(3) :: n
    character(len=1), intent(in ), dimension(0:1,3) :: bc
    integer, intent(out), dimension(3    ) :: n_z,lo,hi
    integer, intent(out), dimension(0:1,3) :: nb
    logical, intent(out), dimension(0:1,3) :: is_bound
    logical, dimension(3) :: periods
    integer :: l1,l2,l
    !
    periods(:) = .false.
    if( bc(0,1)//bc(1,1).eq.'PP' ) periods(1) = .true.
    if( bc(0,2)//bc(1,2).eq.'PP' ) periods(2) = .true.
    if( bc(0,3)//bc(1,3).eq.'PP' ) periods(3) = .true.
    call decomp_2d_init(n(1),n(2),n(3),dims(1),dims(2),periods)
    myid = nrank
    n_z(:) = zsize(:)
#ifdef DECOMP_X
    ipencil=1
    l1 = 2
    l2 = 3
    comm_cart = DECOMP_2D_COMM_CART_X
    lo(:) = xstart(:)
    hi(:) = xend(:)
#elif  DECOMP_Y
    ipencil=2
    l1 = 1
    l2 = 3
    comm_cart = DECOMP_2D_COMM_CART_Y
    lo(:) = ystart(:)
    hi(:) = yend(:)
#else /*DECOMP_Z*/
    ipencil=3
    l1 = 1
    l2 = 2
    comm_cart = DECOMP_2D_COMM_CART_Z
    lo(:) = zstart(:)
    hi(:) = zend(:)
#endif
    nb(:,:) = MPI_PROC_NULL
    call MPI_CART_SHIFT(comm_cart,0,1,nb(0,l1),nb(1,l1),ierr)
    call MPI_CART_SHIFT(comm_cart,1,1,nb(0,l2),nb(1,l2),ierr)
    is_bound(:,:) = .false.
    where(nb(:,:).eq.MPI_PROC_NULL) is_bound(:,:) = .true.
    do l=1,3
      call makehalo(l,1,hi(:)-lo(:)+1,halo(l))
    enddo
  end subroutine initmpi
  subroutine makehalo(idir,nh,n,halo)
    implicit none
    integer, intent(in ) :: idir,nh
    integer, intent(in ), dimension(3) :: n
    integer, intent(out) :: halo
    integer, dimension(3) :: nn
    nn(:) = n(:) + 2*nh
    select case(idir)
    case(1)
      call MPI_TYPE_VECTOR(nn(2)*nn(3),nh            ,nn(1)            ,MPI_REAL_RP,halo,ierr)
    case(2)
      call MPI_TYPE_VECTOR(      nn(3),nh*nn(1)      ,nn(1)*nn(2)      ,MPI_REAL_RP,halo,ierr)
    case(3)
      call MPI_TYPE_VECTOR(          1,nh*nn(1)*nn(2),nn(1)*nn(2)*nn(3),MPI_REAL_RP,halo,ierr)
    end select
    call MPI_TYPE_COMMIT(halo,ierr)
  end subroutine makehalo
end module mod_initmpi
