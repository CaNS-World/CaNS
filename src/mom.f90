module mod_mom
  use mpi
  use mod_common_mpi, only: ierr
  use mod_types
  implicit none
  private
  public momx_a,momy_a,momz_a, &
         momx_d,momy_d,momz_d, &
         momx_p,momy_p,momz_p, cmpt_wallshear
  contains
  !
  subroutine momx_a(nx,ny,nz,dxi,dyi,dzci,dzfi,u,v,w,dudt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in   ) :: u,v,w
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    integer :: i,j,k
    real(rp) :: uuip,uuim,uvjp,uvjm,uwkp,uwkm
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(uuip,uuim,uvjp,uvjm,uwkp,uwkm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,u,v,w,dudt,dzci,dzfi)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          uuip  = 0.25*( u(i+1,j,k)+u(i,j,k) )*( u(i+1,j  ,k  )+u(i,j  ,k  ) )
          uuim  = 0.25*( u(i-1,j,k)+u(i,j,k) )*( u(i-1,j  ,k  )+u(i,j  ,k  ) )
          uvjp  = 0.25*( u(i,j+1,k)+u(i,j,k) )*( v(i+1,j  ,k  )+v(i,j  ,k  ) )
          uvjm  = 0.25*( u(i,j-1,k)+u(i,j,k) )*( v(i+1,j-1,k  )+v(i,j-1,k  ) )
          uwkp  = 0.25*( u(i,j,k+1)+u(i,j,k) )*( w(i+1,j  ,k  )+w(i,j  ,k  ) )
          uwkm  = 0.25*( u(i,j,k-1)+u(i,j,k) )*( w(i+1,j  ,k-1)+w(i,j  ,k-1) )
          !
          ! Momentum balance
          !
          dudt(i,j,k) = dudt(i,j,k) + &
                        dxi*(     -uuip + uuim ) + &
                        dyi*(     -uvjp + uvjm ) + &
                        dzfi(k)*( -uwkp + uwkm )
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine momx_a
  !
  subroutine momy_a(nx,ny,nz,dxi,dyi,dzci,dzfi,u,v,w,dvdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in   ) :: u,v,w
    real(rp), dimension( :, :, :), intent(inout) :: dvdt
    integer :: i,j,k
    real(rp) :: uvip,uvim,vvjp,vvjm,wvkp,wvkm
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(uvip,uvim,vvjp,vvjm,wvkp,wvkm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzci,dzfi,u,v,w,dvdt)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          uvip  = 0.25*( u(i  ,j,k  )+u(i  ,j+1,k  ) )*( v(i,j,k)+v(i+1,j,k) )
          uvim  = 0.25*( u(i-1,j,k  )+u(i-1,j+1,k  ) )*( v(i,j,k)+v(i-1,j,k) )
          vvjp  = 0.25*( v(i  ,j,k  )+v(i  ,j+1,k  ) )*( v(i,j,k)+v(i,j+1,k) )
          vvjm  = 0.25*( v(i  ,j,k  )+v(i  ,j-1,k  ) )*( v(i,j,k)+v(i,j-1,k) )
          wvkp  = 0.25*( w(i  ,j,k  )+w(i  ,j+1,k  ) )*( v(i,j,k+1)+v(i,j,k) )
          wvkm  = 0.25*( w(i  ,j,k-1)+w(i  ,j+1,k-1) )*( v(i,j,k-1)+v(i,j,k) )
          !
          ! Momentum balance
          !
          dvdt(i,j,k) = dvdt(i,j,k) + &
                        dxi*(     -uvip + uvim ) + &
                        dyi*(     -vvjp + vvjm ) + &
                        dzfi(k)*( -wvkp + wvkm )
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine momy_a
  !
  subroutine momz_a(nx,ny,nz,dxi,dyi,dzci,dzfi,u,v,w,dwdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in   ) :: u,v,w
    real(rp), dimension( :, :, :), intent(inout) :: dwdt
    integer :: i,j,k
    real(rp) :: uwip,uwim,vwjp,vwjm,wwkp,wwkm
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(uwip,uwim,vwjp,vwjm,wwkp,wwkm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzci,dzfi,u,v,w,dwdt)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          uwip  = 0.25*( w(i,j,k)+w(i+1,j,k) )*( u(i  ,j  ,k)+u(i  ,j  ,k+1) )
          uwim  = 0.25*( w(i,j,k)+w(i-1,j,k) )*( u(i-1,j  ,k)+u(i-1,j  ,k+1) )
          vwjp  = 0.25*( w(i,j,k)+w(i,j+1,k) )*( v(i  ,j  ,k)+v(i  ,j  ,k+1) )
          vwjm  = 0.25*( w(i,j,k)+w(i,j-1,k) )*( v(i  ,j-1,k)+v(i  ,j-1,k+1) )
          wwkp  = 0.25*( w(i,j,k)+w(i,j,k+1) )*( w(i  ,j  ,k)+w(i  ,j  ,k+1) )
          wwkm  = 0.25*( w(i,j,k)+w(i,j,k-1) )*( w(i  ,j  ,k)+w(i  ,j  ,k-1) )
          !
          ! Momentum balance
          !
          dwdt(i,j,k) = dwdt(i,j,k) + &
                        dxi*(     -uwip + uwim ) + &
                        dyi*(     -vwjp + vwjm ) + &
                        dzci(k)*( -wwkp + wwkm )
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine momz_a
  !
  subroutine momx_d(nx,ny,nz,dxi,dyi,dzci,dzfi,visc,u,dudt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,visc
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in   ) :: u
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    real(rp) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm
    integer :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(dudxp,dudxm,dudyp,dudym,dudzp,dudzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzci,dzfi,u,dudt,visc)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dudxp = (u(i+1,j,k)-u(i,j,k))*dxi
          dudxm = (u(i,j,k)-u(i-1,j,k))*dxi
          dudyp = (u(i,j+1,k)-u(i,j,k))*dyi
          dudym = (u(i,j,k)-u(i,j-1,k))*dyi
          dudzp = (u(i,j,k+1)-u(i,j,k))*dzci(k  )
          dudzm = (u(i,j,k)-u(i,j,k-1))*dzci(k-1)
          dudt(i,j,k) = dudt(i,j,k) + &
                        (dudxp-dudxm)*visc*dxi + &
                        (dudyp-dudym)*visc*dyi + &
                        (dudzp-dudzm)*visc*dzfi(k)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine momx_d
  !
  subroutine momy_d(nx,ny,nz,dxi,dyi,dzci,dzfi,visc,v,dvdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,visc
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in   ) :: v
    real(rp), dimension( :, :, :), intent(inout) :: dvdt
    real(rp) :: dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm
    integer :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzci,dzfi,v,dvdt,visc)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dvdxp = (v(i+1,j,k)-v(i,j,k))*dxi
          dvdxm = (v(i,j,k)-v(i-1,j,k))*dxi
          dvdyp = (v(i,j+1,k)-v(i,j,k))*dyi
          dvdym = (v(i,j,k)-v(i,j-1,k))*dyi
          dvdzp = (v(i,j,k+1)-v(i,j,k))*dzci(k  )
          dvdzm = (v(i,j,k)-v(i,j,k-1))*dzci(k-1)
          dvdt(i,j,k) = dvdt(i,j,k) + &
                        (dvdxp-dvdxm)*visc*dxi+ &
                        (dvdyp-dvdym)*visc*dyi+ &
                        (dvdzp-dvdzm)*visc*dzfi(k)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine momy_d
  !
  subroutine momz_d(nx,ny,nz,dxi,dyi,dzci,dzfi,visc,w,dwdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,visc
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in   ) :: w
    real(rp), dimension( :, :, :), intent(inout) :: dwdt
    integer :: i,j,k
    real(rp) :: dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzci,dzfi,w,dwdt,visc)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dwdxp = (w(i+1,j,k)-w(i,j,k))*dxi
          dwdxm = (w(i,j,k)-w(i-1,j,k))*dxi
          dwdyp = (w(i,j+1,k)-w(i,j,k))*dyi
          dwdym = (w(i,j,k)-w(i,j-1,k))*dyi
          dwdzp = (w(i,j,k+1)-w(i,j,k))*dzfi(k+1)
          dwdzm = (w(i,j,k)-w(i,j,k-1))*dzfi(k  )
          dwdt(i,j,k) = dwdt(i,j,k) + &
                        (dwdxp-dwdxm)*visc*dxi+ &
                        (dwdyp-dwdym)*visc*dyi+ &
                        (dwdzp-dwdzm)*visc*dzci(k)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine momz_d
  !
  subroutine momx_p(nx,ny,nz,dxi,bforce,p,dudt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi
    real(rp), intent(in) :: bforce
    real(rp), dimension(0:,0:,0:), intent(in) :: p
    real(rp), dimension(:,:,:), intent(out) :: dudt
    integer :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(nx,ny,nz,dxi,bforce,p,dudt)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dudt(i,j,k) = - dxi*( p(i+1,j,k)-p(i,j,k) ) + bforce
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine momx_p
  !
  subroutine momy_p(nx,ny,nz,dyi,bforce,p,dvdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dyi
    real(rp), intent(in) :: bforce
    real(rp), dimension(0:,0:,0:), intent(in) :: p
    real(rp), dimension(:,:,:), intent(out) :: dvdt
    integer :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(nx,ny,nz,dyi,bforce,p,dvdt)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dvdt(i,j,k) = - dyi*( p(i,j+1,k)-p(i,j,k) ) + bforce
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine momy_p
  !
  subroutine momz_p(nx,ny,nz,dzci,bforce,p,dwdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in), dimension(0:) :: dzci
    real(rp), intent(in) :: bforce
    real(rp), dimension(0:,0:,0:), intent(in) :: p
    real(rp), dimension(:,:,:), intent(out) :: dwdt
    integer :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(nx,ny,nz,dzci,bforce,p,dwdt)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dwdt(i,j,k) = - dzci(k)*( p(i,j,k+1)-p(i,j,k) ) + bforce
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine momz_p
  subroutine cmpt_wallshear(n,is_bound,l,dli,dzci,dzfi,visc,u,v,w,taux,tauy,tauz)
    implicit none
    integer , intent(in ), dimension(3) :: n
    logical , intent(in ), dimension(0:1,3) :: is_bound
    real(rp), intent(in ), dimension(3)     :: l,dli
    real(rp), intent(in ), dimension(0:)    :: dzci,dzfi
    real(rp), intent(in )                   :: visc
    real(rp), intent(in ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out), dimension(3) :: taux,tauy,tauz
    real(rp) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm, &
                dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm, &
                dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
    integer :: i,j,k,nx,ny,nz
    !
    nx = n(1)
    ny = n(2)
    nz = n(3)
    taux(:) = 0.
    if(is_bound(0,2)) then
      do k=1,nz
        do i=1,nx
          dudyp = (u(i,1 ,k)-u(i,0   ,k))*dli(2)*visc
          taux(2) = taux(2) + dudyp/(dli(1)*dzfi(k)*l(1)*l(3))
        end do
      end do
    end if
    if(is_bound(1,2)) then
      do k=1,nz
        do i=1,nx
          dudym = (u(i,ny,k)-u(i,ny+1,k))*dli(2)*visc
          taux(2) = taux(2) + dudyp/(dli(1)*dzfi(k)*l(1)*l(3))
        end do
      end do
    end if
    if(is_bound(0,3)) then
      do j=1,ny
        do i=1,nx
          dudzp = (u(i,j,1 )-u(i,j,0   ))*dzci(0)*visc
          taux(3) = taux(3) + dudzp/(dli(1)*dli(2)*l(1)*l(2))
        end do
      end do
    end if
    if(is_bound(1,3)) then
      do j=1,ny
        do i=1,nx
          dudzm = (u(i,j,nz)-u(i,j,nz+1))*dzci(nz)*visc
          taux(3) = taux(3) + dudzm/(dli(1)*dli(2)*l(1)*l(2))
        end do
      end do
    end if
    call mpi_allreduce(MPI_IN_PLACE,taux(1),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    !
    tauy(:) = 0.
    if(is_bound(0,1)) then
      do k=1,nz
        do j=1,ny
          dvdxp = (v(1  ,j,k)-v(0  ,j,k))*dli(1)*visc
          tauy(1) = tauy(1) + dvdxp/(dli(2)*dzfi(1)*l(2)*l(3))
        end do
      end do
    end if
    if(is_bound(1,1)) then
      do k=1,nz
        do j=1,ny
          dvdxm = (v(nx,j,k)-v(nx+1,j,k))*dli(1)*visc
          tauy(1) = tauy(1) + dvdxm/(dli(2)*dzfi(1)*l(2)*l(3))
        end do
      end do
    end if
    if(is_bound(0,3)) then
      do j=1,ny
        do i=1,nx
          dvdzp = (v(i,j,1 )-v(i,j,0   ))*dzci(0)*visc
          tauy(3) = tauy(3) + dvdzp/(dli(1)*dli(2)*l(1)*l(2))
        end do
      end do
    end if
    if(is_bound(1,3)) then
      do j=1,ny
        do i=1,nx
          dvdzm = (v(i,j,nz)-v(i,j,nz+1))*dzci(nz)*visc
          tauy(3) = tauy(3) + dvdzm/(dli(1)*dli(2)*l(1)*l(2))
        end do
      end do
    end if
    call mpi_allreduce(MPI_IN_PLACE,tauy(1),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    !
    tauz(:) = 0.
    if(is_bound(0,1)) then
      do k=1,nz
        do j=1,ny
          dwdxp = (w(1 ,j,k)-w(0   ,j,k))*dli(1)*visc
          tauz(1) = tauz(1) + dwdxp/(dli(2)*dzfi(k)*l(2)*l(3))
        end do
      end do
    end if
    if(is_bound(1,1)) then
      do k=1,nz
        do j=1,ny
          dwdxm = (w(nx,j,k)-w(nx+1,j,k))*dli(1)*visc
          tauz(1) = tauz(1) + dwdxm/(dli(2)*dzfi(k)*l(2)*l(3))
        end do
      end do
    end if
    if(is_bound(0,2)) then
      do k=1,nz
        do i=1,nx
          dwdyp = (w(i,1,k )-w(i,0   ,k))*dli(2)*visc
          tauz(2) = tauz(2) + dwdyp/(dli(1)*dzfi(k)*l(1)*l(3))
        end do
      end do
    end if
    if(is_bound(1,2)) then
      do k=1,nz
        do i=1,nx
          dwdym = (w(i,ny,k)-w(i,ny+1,k))*dli(2)*visc
          tauz(2) = tauz(2) + dwdym/(dli(1)*dzfi(k)*l(1)*l(3))
        end do
      end do
    end if
    call mpi_allreduce(MPI_IN_PLACE,tauz(1),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  end subroutine cmpt_wallshear
end module mod_mom
