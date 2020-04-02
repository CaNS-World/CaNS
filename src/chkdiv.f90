module mod_chkdiv
  use mpi
  use mod_common_mpi, only: myid,ierr,ijk_start
  use mod_types
  implicit none
  private
  public chkdiv
  contains
  subroutine chkdiv(n,dli,dzfi,u,v,w,divtot,divmax)
    !
    ! checks the divergence of the velocity field
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3) :: dli
    real(rp), intent(in), dimension(0:) :: dzfi
    real(rp), intent(in), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out) :: divtot,divmax
    real(rp) :: dxi,dyi,div!,dzi,div
    integer :: i,j,k,im,jm,km
    integer :: ii,jj,kk
    !
    dxi = dli(1)
    dyi = dli(2)
    !dzi = dli(3)
    divtot = 0.
    divmax = 0.
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,u,v,w,dxi,dyi,dzfi) &
    !$OMP PRIVATE(i,j,k,im,jm,km,div) &
    !$OMP REDUCTION(+:divtot) &
    !$OMP REDUCTION(max:divmax)
    do k=1,n(3)
       km = k-1
       do j=1,n(2)
          jm = j-1
          do i=1,n(1)
             im = i-1
             div = (w(i,j,k)-w(i,j,km))*dzfi(k) + &
                   (v(i,j,k)-v(i,jm,k))*dyi     + &
                   (u(i,j,k)-u(im,j,k))*dxi
             divmax = max(divmax,abs(div))
             divtot = divtot + div
             ii = ijk_start(1)+i
             jj = ijk_start(2)+j
             kk = ijk_start(3)+k
             if(abs(div).ge.1.e-12) print*,div,'Large divergence at grid cell: ',ii,jj,kk,div
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    call mpi_allreduce(MPI_IN_PLACE,divtot,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,divmax,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(myid.eq.0) print*, 'Total divergence = ', divtot, '| Maximum divergence = ', divmax
    return
  end subroutine chkdiv
end module mod_chkdiv
