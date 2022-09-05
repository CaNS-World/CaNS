! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
!  subroutine gaussel(nx,ny,n,nh,a,b,c,p,lambdaxy)
    implicit none
    integer , intent(in) :: nx,ny,n,nh
    real(wp), intent(in), dimension(:) :: a,b,c
    real(wp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    real(wp), intent(in), dimension(nx,ny), optional :: lambdaxy
    real(wp), dimension(n) :: bb
    integer :: i,j
    !
    ! solve tridiagonal system
    !
    if(present(lambdaxy)) then
      !$OMP PARALLEL DEFAULT(none) &
      !$OMP PRIVATE(bb) &
      !$OMP SHARED(nx,ny,n,a,b,c,lambdaxy,p)
      !$OMP DO COLLAPSE(2)
      do j=1,ny
        do i=1,nx
          bb(:) = b(1:n) + lambdaxy(i,j)
          call dgtsv_homebrewed(n,a,bb,c,p(i,j,1:n))
        end do
      end do
      !$OMP END PARALLEL
    else
      !$OMP PARALLEL DEFAULT(none) &
      !$OMP SHARED(nx,ny,n,a,b,c,p)
      !$OMP DO COLLAPSE(2)
      do j=1,ny
        do i=1,nx
          call dgtsv_homebrewed(n,a,b,c,p(i,j,1:n))
        end do
      end do
      !$OMP END PARALLEL
    end if
!  end subroutine gaussel
