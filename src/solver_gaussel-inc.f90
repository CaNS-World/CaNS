!  subroutine gaussel(nx,ny,n,nh,a,b,c,p,lambdaxy)
    implicit none
    integer , intent(in) :: nx,ny,n,nh
    MYREAL, intent(in), dimension(:) :: a,b,c
    MYREAL, intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    MYREAL, intent(in), dimension(nx,ny), optional :: lambdaxy
    MYREAL, dimension(n) :: bb
    integer :: i,j
    !
    ! solve tridiagonal system
    !
    if(present(lambdaxy)) then
      !$OMP PARALLEL DEFAULT(none) &
      !$OMP PRIVATE(i,j,bb) &
      !$OMP SHARED(nx,ny,n,a,b,c,lambdaxy,p)
      !$OMP DO COLLAPSE(2)
      do j=1,ny
        do i=1,nx
          bb(:) = b(1:n) + lambdaxy(i,j)
          call dgtsv_homebrewed(n,a,bb,c,p(i,j,1:n))
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL
    else
      !$OMP PARALLEL DEFAULT(none) &
      !$OMP PRIVATE(i,j) &
      !$OMP SHARED(nx,ny,n,a,b,c,p)
      !$OMP DO COLLAPSE(2)
      do j=1,ny
        do i=1,nx
          call dgtsv_homebrewed(n,a,b,c,p(i,j,1:n))
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL
    endif
!  end subroutine gaussel
