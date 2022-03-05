module mod_correc
  use mod_types
  implicit none
  private
  public correc
  contains
  subroutine correc(n,dli,dzci,dt,p,u,v,w)
    !
    ! corrects the velocity so that it is divergence free
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3 ) :: dli
    real(rp), intent(in), dimension(0:) :: dzci
    real(rp), intent(in) :: dt
    real(rp), intent(in   ), dimension(0:,0:,0:) :: p
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    real(rp) :: factori,factorj
    real(rp), dimension(0:n(3)+1) :: factork
    integer :: i,j,k
    !
    !factor = rkcoeffab(rkiter)*dt
    !
    factori = dt*dli(1)
    factorj = dt*dli(2)
    factork = dt*dzci!dli(3)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,factori,u,p) &
    !$OMP PRIVATE(i,j,k)
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)
          u(i,j,k) = u(i,j,k) - factori*(   p(i+1,j,k)-p(i,j,k))
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,factorj,v,p) &
    !$OMP PRIVATE(i,j,k)
    do k=0,n(3)+1
      do j=0,n(2)
        do i=0,n(1)+1
          v(i,j,k) = v(i,j,k) - factorj*(   p(i,j+1,k)-p(i,j,k))
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,factork,w,p) &
    !$OMP PRIVATE(i,j,k)
    do k=0,n(3)
      do j=0,n(2)+1
        do i=0,n(1)+1
          w(i,j,k) = w(i,j,k) - factork(k)*(p(i,j,k+1)-p(i,j,k))
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine correc
end module mod_correc
