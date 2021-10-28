module mod_fillps
  use mod_types
  implicit none
  private
  public fillps
  contains
  subroutine fillps(n,dli,dzfi,dti,up,vp,wp,p)
    !
    !  fill the right-hand side of the Poisson equation for the correction pressure.
    !
    !  the discrete divergence is:
    !
    !  w(i,j,k)-w(i,j,k-1)   v(i,j,k)-v(i,j-1,k)   u(i,j,k)-u(i-1,j,k)
    !  ------------------- + ------------------- + -------------------  = div
    !          dz                    dy                    dx
    !
    implicit none
    integer , intent(in ), dimension(3) :: n
    real(rp), intent(in ), dimension(3) :: dli
    real(rp), intent(in ), dimension(0:) :: dzfi
    real(rp), intent(in ) :: dti
    real(rp), intent(in ), dimension(0:,0:,0:) :: up,vp,wp
    real(rp), intent(out), dimension(0:,0:,0:) :: p
    real(rp) :: dtidxi,dtidyi!,dtidzi
    real(rp), dimension(0:n(3)+1) :: dtidzfi
    integer :: i,j,k
    !
    dtidxi = dti*dli(1)
    dtidyi = dti*dli(2)
    !dtidzi = dti*dli(3)
    dtidzfi(:) = dti*dzfi(:)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,p,up,vp,wp,dtidzfi,dtidyi,dtidxi) &
    !$OMP PRIVATE(i,j,k)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          p(i,j,k) = ( &
                      (wp(i,j,k)-wp(i,j,k-1))*dtidzfi(k)+ &
                      (vp(i,j,k)-vp(i,j-1,k))*dtidyi    + &
                      (up(i,j,k)-up(i-1,j,k))*dtidxi    )
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    !
  end subroutine fillps
end module mod_fillps
