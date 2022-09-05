! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
!  subroutine gaussel_periodic_gpu(nx,ny,n,nh,a,b,c,p,d,p1,p2,lambdaxy)
    implicit none
    integer , intent(in) :: nx,ny,n,nh
    real(wp), intent(in), dimension(:) :: a,b,c
    real(wp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    !dir$ ignore_tkr d,p1,p2
    real(wp),                dimension(nx,ny,n) :: d,p1,p2
    real(wp), intent(in), dimension(:,:), optional :: lambdaxy
    real(wp) :: z,lxy
    integer :: i,j,k
    !
    ! solve tridiagonal system
    !
    if(present(lambdaxy)) then
      !$acc parallel loop gang vector collapse(2) default(present) private(lxy,z) async(1)
      do j=1,ny
        do i=1,nx
          lxy = lambdaxy(i,j)
          !$acc loop seq
          do k=1,n-1
            p1(i,j,k) = p(i,j,k)
          end do
          !
          ! call dgtsv_homebrewed(n-1,a,bb,c,p1)
          !
          z = 1./(b(1)+lxy+eps)
          d(i,j,1)  = c(1)*z
          p1(i,j,1) = p1(i,j,1)*z
          !$acc loop seq
          do k=2,n-1
            z         = 1./(b(k)+lxy-a(k)*d(i,j,k-1)+eps)
            d(i,j,k)  = c(k)*z
            p1(i,j,k) = (p1(i,j,k)-a(k)*p1(i,j,k-1))*z
          end do
          !
          ! backward substitution
          !
          !$acc loop seq
          do k=n-2,1,-1
            p1(i,j,k) = p1(i,j,k) - d(i,j,k)*p1(i,j,k+1)
          end do
          !
          !$acc loop seq
          do k=1,n
            p2(i,j,k) = 0.
          end do
          p2(i,j,1  ) = -a(1  )
          p2(i,j,n-1) = -c(n-1)
          !
          ! call dgtsv_homebrewed(n-1,a,bb,c,p2)
          !
          z = 1./(b(1)+lxy+eps)
          d( i,j,1) = c(1)*z
          p2(i,j,1) = p2(i,j,1)*z
          !$acc loop seq
          do k=2,n-1
            z        = 1./(b(k)+lxy-a(k)*d(i,j,k-1)+eps)
            d(i,j,k) = c(k)*z
            p2(i,j,k) = (p2(i,j,k)-a(k)*p2(i,j,k-1))*z
          end do
          !
          ! backward substitution
          !
          !$acc loop seq
          do k=n-2,1,-1
            p2(i,j,k) = p2(i,j,k) - d(i,j,k)*p2(i,j,k+1)
          end do
          p(i,j,n) = (p(i,j,n)       - c(n)*p1(i,j,1) - a(n)*p1(i,j,n-1)) / &
                     (b(    n) + lxy + c(n)*p2(i,j,1) + a(n)*p2(i,j,n-1)+eps)
          !$acc loop seq
          do k=1,n-1
            p(i,j,k) = p1(i,j,k) + p2(i,j,k)*p(i,j,n)
          end do
        end do
      end do
    else
      !$acc parallel loop gang vector collapse(2) default(present) private(z) async(1)
      do j=1,ny
        do i=1,nx
          !$acc loop seq
          do k=1,n-1
            p1(i,j,k) = p(i,j,k)
          end do
          !
          ! call dgtsv_homebrewed(n-1,a,b,c,p1)
          !
          z = 1./(b(1)+eps)
          d(i,j,1)  = c(1)*z
          p1(i,j,1) = p1(i,j,1)*z
          !$acc loop seq
          do k=2,n-1
            z         = 1./(b(k)-a(k)*d(i,j,k-1)+eps)
            d(i,j,k)  = c(k)*z
            p1(i,j,k) = (p1(i,j,k)-a(k)*p1(i,j,k-1))*z
          end do
          !
          ! backward substitution
          !
          !$acc loop seq
          do k=n-2,1,-1
            p1(i,j,k) = p1(i,j,k) - d(i,j,k)*p1(i,j,k+1)
          end do
          !
          !!$acc loop seq
          !do k=1,n
          !  p2(i,j,k) = 0.
          !end do
          p2(i,j,1  ) = -a(1  )
          p2(i,j,n-1) = -c(n-1)
          !
          ! call dgtsv_homebrewed(n-1,a,b,c,p2)
          !
          z = 1./(b(1)+eps)
          d( i,j,1) = c(1)*z
          p2(i,j,1) = p2(i,j,1)*z
          !$acc loop seq
          do k=2,n-1
            z        = 1./(b(k)-a(k)*d(i,j,k-1)+eps)
            d(i,j,k) = c(k)*z
           !p2(i,j,k) = (p2(i,j,k)-a(k)*p2(i,j,k-1))*z
            p2(i,j,k) = (         -a(k)*p2(i,j,k-1))*z
          end do
          !
          ! backward substitution
          !
          !$acc loop seq
          do k=n-2,1,-1
            p2(i,j,k) = p2(i,j,k) - d(i,j,k)*p2(i,j,k+1)
          end do
          p(i,j,n) = (p(i,j,n) - c(n)*p1(i,j,1) - a(n)*p1(i,j,n-1)) / &
                     (b(    n) + c(n)*p2(i,j,1) + a(n)*p2(i,j,n-1)+eps)
          !$acc loop seq
          do k=1,n-1
            p(i,j,k) = p1(i,j,k) + p2(i,j,k)*p(i,j,n)
          end do
        end do
      end do
    end if
!  end subroutine gaussel_periodic_gpu
