! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
!  subroutine gaussel_periodic(nx,ny,n,nh,a,b,c,p,lambdaxy)
    implicit none
    integer , intent(in) :: nx,ny,n,nh
    real(wp), intent(in), dimension(:) :: a,b,c
    real(wp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    real(wp), intent(in), dimension(nx,ny), optional :: lambdaxy
    real(wp), dimension(n) :: bb,p1,p2
    integer :: i,j
    !
    ! solve tridiagonal system
    !
    if(present(lambdaxy)) then
      !$OMP PARALLEL DEFAULT(none) &
      !$OMP PRIVATE(bb,p1,p2) &
      !$OMP SHARED(nx,ny,n,a,b,c,lambdaxy,p)
      !$OMP DO COLLAPSE(2)
      do j=1,ny
        do i=1,nx
          bb(:)  = b(:) + lambdaxy(i,j)
          p1(1:n-1) = p(i,j,1:n-1)
          call dgtsv_homebrewed(n-1,a,bb,c,p1)
          p2(:) = 0.
          p2(1  ) = -a(1  )
          p2(n-1) = -c(n-1)
          call dgtsv_homebrewed(n-1,a,bb,c,p2)
          p(i,j,n) = (p(i,j,n) - c(n)*p1(1) - a(n)*p1(n-1)) / &
                     (bb(   n) + c(n)*p2(1) + a(n)*p2(n-1)+eps)
          p(i,j,1:n-1) = p1(1:n-1) + p2(1:n-1)*p(i,j,n)
        end do
      end do
      !$OMP END PARALLEL
    else
      !$OMP PARALLEL DEFAULT(none) &
      !$OMP PRIVATE(p1,p2) &
      !$OMP SHARED(nx,ny,n,a,b,c,p)
      !$OMP DO COLLAPSE(2)
      do j=1,ny
        do i=1,nx
          p1(1:n-1) = p(i,j,1:n-1)
          call dgtsv_homebrewed(n-1,a,b,c,p1)
          p2(:) = 0.
          p2(1  ) = -a(1  )
          p2(n-1) = -c(n-1)
          call dgtsv_homebrewed(n-1,a,b,c,p2)
          p(i,j,n) = (p(i,j,n) - c(n)*p1(1) - a(n)*p1(n-1)) / &
                     (b(    n) + c(n)*p2(1) + a(n)*p2(n-1)+eps)
          p(i,j,1:n-1) = p1(1:n-1) + p2(1:n-1)*p(i,j,n)
        end do
      end do
      !$OMP END PARALLEL
    end if
!  end subroutine gaussel_periodic
