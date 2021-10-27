module mod_chkdiv
  use mpi
  use mod_common_mpi, only: myid,ierr
  use mod_types
  implicit none
  private
  public chkdiv
  contains
  subroutine chkdiv(lo,hi,dli,dzfi,u,v,w,divtot,divmax)
    !
    ! checks the divergence of the velocity field
    !
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(3) :: dli
    real(rp), intent(in), dimension(lo(3)-1:) :: dzfi
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: u,v,w
    real(rp), intent(out) :: divtot,divmax
    real(rp) :: dxi,dyi,div!,dzi
    integer :: i,j,k
    !
    dxi = dli(1)
    dyi = dli(2)
    !dzi = dli(3)
    divtot = 0.
    divmax = 0.
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(lo,hi,u,v,w,dxi,dyi,dzfi) &
    !$OMP PRIVATE(i,j,k,div) &
    !$OMP REDUCTION(+:divtot) &
    !$OMP REDUCTION(max:divmax)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          div = (w(i,j,k)-w(i,j,k-1))*dzfi(k) + &
                (v(i,j,k)-v(i,j-1,k))*dyi     + &
                (u(i,j,k)-u(i-1,j,k))*dxi
          divmax = max(divmax,abs(div))
          divtot = divtot + div
          !if(abs(div).ge.1.e-12) print*,div,'Large divergence at grid cell: ',i,j,k,div
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    call mpi_allreduce(MPI_IN_PLACE,divtot,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,divmax,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
  end subroutine chkdiv
end module mod_chkdiv
