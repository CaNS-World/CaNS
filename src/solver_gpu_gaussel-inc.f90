! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
!  subroutine gaussel_gpu(nx,ny,n,nh,a,b,c,p,d,lambdaxy)
    implicit none
    integer , intent(in) :: nx,ny,n,nh
    real(wp), intent(in), dimension(:) :: a,b,c
    real(wp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    !dir$ ignore_tkr d
    real(wp),                dimension(nx,ny,n) :: d
    real(wp), intent(in), dimension(:,:), optional :: lambdaxy
    real(wp) :: z,lxy
    integer :: i,j,k
    !
    !solve tridiagonal system
    !
    if(present(lambdaxy)) then
      !$acc parallel loop gang vector collapse(2) default(present) private(lxy,z) async(1)
      do j=1,ny
        do i=1,nx
          lxy = lambdaxy(i,j)
          !
          ! call dgtsv_homebrewed(n,a,bb,c,p)
          !
          z = 1./(b(1)+lxy+eps)
          d(i,j,1) = c(1)*z
          p(i,j,1) = p(i,j,1)*z
          !$acc loop seq
          do k=2,n
            z = 1./(b(k)+lxy-a(k)*d(i,j,k-1)+eps)
            d(i,j,k) = c(k)*z
            p(i,j,k) = (p(i,j,k)-a(k)*p(i,j,k-1))*z
          end do
          !
          ! backward substitution
          !
          !$acc loop seq
          do k=n-1,1,-1
            p(i,j,k) = p(i,j,k) - d(i,j,k)*p(i,j,k+1)
          end do
        end do
      end do
    else
      !$acc parallel loop gang vector collapse(2) default(present) private(z) async(1)
      do j=1,ny
        do i=1,nx
          !
          ! call dgtsv_homebrewed(n,a,b,c,p)
          !
          z = 1./(b(1)+eps)
          d(i,j,1) = c(1)*z
          p(i,j,1) = p(i,j,1)*z
          !$acc loop seq
          do k=2,n
            z = 1./(b(k)-a(k)*d(i,j,k-1)+eps)
            d(i,j,k) = c(k)*z
            p(i,j,k) = (p(i,j,k)-a(k)*p(i,j,k-1))*z
          end do
          !
          ! backward substitution
          !
          !$acc loop seq
          do k=n-1,1,-1
            p(i,j,k) = p(i,j,k) - d(i,j,k)*p(i,j,k+1)
          end do
        end do
      end do
    end if
!  end subroutine gaussel_gpu
