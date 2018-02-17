module mod_correc
  implicit none
  private
  public correc
  contains
  !
  ! corrects the velocity so that it is divergence free
  !
  subroutine correc(n,dli,dzci,dt,p,up,vp,wp,u,v,w)
    implicit none
    integer, intent(in), dimension(3) :: n
    real(8), intent(in), dimension(3) :: dli
    real(8), intent(in), dimension(0:) :: dzci
    real(8), intent(in) :: dt
    real(8), intent(in) , dimension(0:,0:,0:) :: p,up,vp,wp
    real(8), intent(out), dimension(0:,0:,0:) :: u,v,w
    real(8) :: factori,factorj
    real(8), dimension(0:n(3)+1) :: factork
    integer :: i,j,k,ip,jp,kp
    !
    !factor = rkcoeffab(rkiter)*dt
    !
    factori = dt*dli(1)
    factorj = dt*dli(2)
    factork = dt*dzci!dli(3)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,factori,factorj,factork,u,v,w,up,vp,wp,p) &
    !$OMP PRIVATE(i,j,k,ip,jp,kp)
    do k=1,n(3)
      kp = k+1
      do j=1,n(2)
        jp = j+1
        do i=1,n(1)
          ip = i+1
          u(i,j,k) = up(i,j,k) - factori*(   p(ip,j,k)-p(i,j,k))
          v(i,j,k) = vp(i,j,k) - factorj*(   p(i,jp,k)-p(i,j,k))
          w(i,j,k) = wp(i,j,k) - factork(k)*(p(i,j,kp)-p(i,j,k))
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    return
  end subroutine correc
end module mod_correc
