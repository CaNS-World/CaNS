! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
!  subroutine dgtsv_homebrewed(n,a,b,c,p)
    implicit none
    integer , intent(in) :: n
    real(wp), intent(in   ), dimension(:) :: a,b,c
    real(wp), intent(inout), dimension(:) :: p
    real(wp), dimension(n) :: d
    real(wp) :: z
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
