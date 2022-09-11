! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_updatep
  use mod_types
  implicit none
  private
  public updatep
  contains
  subroutine updatep(n,dli,dzci,dzfi,alpha,pp,p)
    !
    ! updates the final pressure
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    real(rp), intent(in   ), dimension(3 ) :: dli
    real(rp), intent(in   ), dimension(0:) :: dzci,dzfi
    real(rp), intent(in   ), dimension(0:,0:,0:) :: pp
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    real(rp), intent(in   ) :: alpha
    real(rp) :: dxi,dyi
    integer :: i,j,k
    !
#if defined(_IMPDIFF)
    dxi = dli(1); dyi = dli(2)
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,p,pp,dxi,dyi,dzfi,dzci,alpha)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          p(i,j,k) = p(i,j,k) + pp(i,j,k) + alpha*( &
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
    !$acc kernels default(present) async(1)
    !$OMP PARALLEL WORKSHARE
    p(1:n(1),1:n(2),1:n(3)) = p(1:n(1),1:n(2),1:n(3)) + pp(1:n(1),1:n(2),1:n(3))
    !$OMP END PARALLEL WORKSHARE
    !$acc end kernels
#endif
  end subroutine updatep
end module mod_updatep
