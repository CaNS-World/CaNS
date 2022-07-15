!  subroutine dgtsv_homebrewed(n,a,b,c,p)
    implicit none
    integer , intent(in) :: n
    MYREAL, intent(in   ), dimension(:) :: a,b,c
    MYREAL, intent(inout), dimension(:) :: p
    MYREAL, dimension(n) :: d
    MYREAL :: z
    integer :: l
    !
    ! Gauss elimination
    !
    z = 1./(b(1)+eps)
    d(1) = c(1)*z
    p(1) = p(1)*z
    do l=2,n
      z    = 1./(b(l)-a(l)*d(l-1)+eps)
      d(l) = c(l)*z
      p(l) = (p(l)-a(l)*p(l-1))*z
    end do
    !
    ! backward substitution
    !
    do l=n-1,1,-1
      p(l) = p(l) - d(l)*p(l+1)
    end do
!  end subroutine dgtsv_homebrewed
