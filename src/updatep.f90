module mod_updatep
  use mod_types
  implicit none
  private
  public updatep
  contains
  subroutine updatep(n,dli,dzci,dzfi,alpha,u,v,w,pp,p)
    !
    ! updates the final pressure
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    real(rp), intent(in   ), dimension(3 ) :: dli
    real(rp), intent(in   ), dimension(0:) :: dzci,dzfi
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u,v,w,pp
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    real(rp), intent(in   ) :: alpha
    real(rp) :: dxi,dyi,alphai
    integer :: i,j,k
    !
#if defined(_IMPDIFF)
      alphai = alpha**(-1)
      dxi = dli(1); dyi = dli(2)
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP SHARED(n,p,pp,dxi,dyi,dzfi,dzci,alphai)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            p(i,j,k) = p(i,j,k) + pp(i,j,k) + alphai*( &
#if !defined(_IMPDIFF_1D)
                        (pp(i+1,j,k)-2.*pp(i,j,k)+pp(i-1,j,k))*(dxi**2) + &
                        (pp(i,j+1,k)-2.*pp(i,j,k)+pp(i,j-1,k))*(dyi**2) + &
#endif
                        ((pp(i,j,k+1)-pp(i,j,k  ))*dzci(k  ) - &
                         (pp(i,j,k  )-pp(i,j,k-1))*dzci(k-1))*dzfi(k) )
          end do
        end do
      end do
#else
      !$OMP WORKSHARE
      p(1:n(1),1:n(2),1:n(3)) = p(1:n(1),1:n(2),1:n(3)) + pp(1:n(1),1:n(2),1:n(3))
      !$OMP END WORKSHARE
#endif
  end subroutine updatep
end module mod_updatep
