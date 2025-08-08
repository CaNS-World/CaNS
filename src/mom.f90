! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2025 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_mom
  use mpi
  use mod_common_mpi, only: ierr
  use mod_types
  implicit none
  private
  public momx_a,momy_a,momz_a, &
         momx_d,momy_d,momz_d, &
         momx_p,momy_p,momz_p, &
         cmpt_wallshear,bulk_forcing, &
         momx_d_xy,momy_d_xy,momz_d_xy, &
         momx_d_z ,momy_d_z ,momz_d_z, &
         mom_xyz_ad
  contains
  !
  subroutine momx_a(nx,ny,nz,dxi,dyi,dzfi,u,v,w,dudt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(0:) :: dzfi
    real(rp), dimension(0:,0:,0:), intent(in   ) :: u,v,w
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    integer :: i,j,k
    real(rp) :: uuip,uuim,vujp,vujm,wukp,wukm
    !
    !$acc parallel     loop collapse(3) default(present) private(uuip,uuim,vujp,vujm,wukp,wukm) async(1)
    !$omp target teams loop collapse(3)                  private(uuip,uuim,vujp,vujm,wukp,wukm)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          uuip  = 0.25*( u(i+1,j  ,k  )+u(i,j  ,k  ) )*( u(i,j,k)+u(i+1,j,k) )
          uuim  = 0.25*( u(i-1,j  ,k  )+u(i,j  ,k  ) )*( u(i,j,k)+u(i-1,j,k) )
          vujp  = 0.25*( v(i+1,j  ,k  )+v(i,j  ,k  ) )*( u(i,j,k)+u(i,j+1,k) )
          vujm  = 0.25*( v(i+1,j-1,k  )+v(i,j-1,k  ) )*( u(i,j,k)+u(i,j-1,k) )
          wukp  = 0.25*( w(i+1,j  ,k  )+w(i,j  ,k  ) )*( u(i,j,k)+u(i,j,k+1) )
          wukm  = 0.25*( w(i+1,j  ,k-1)+w(i,j  ,k-1) )*( u(i,j,k)+u(i,j,k-1) )
          !
          ! Momentum balance
          !
          dudt(i,j,k) = dudt(i,j,k) + &
                        dxi*(     -uuip + uuim ) + &
                        dyi*(     -vujp + vujm ) + &
                        dzfi(k)*( -wukp + wukm )
        end do
      end do
    end do
  end subroutine momx_a
  !
  subroutine momy_a(nx,ny,nz,dxi,dyi,dzfi,u,v,w,dvdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(0:) :: dzfi
    real(rp), dimension(0:,0:,0:), intent(in   ) :: u,v,w
    real(rp), dimension( :, :, :), intent(inout) :: dvdt
    integer :: i,j,k
    real(rp) :: uvip,uvim,vvjp,vvjm,wvkp,wvkm
    !
    !$acc parallel     loop collapse(3) default(present) private(uvip,uvim,vvjp,vvjm,wvkp,wvkm) async(1)
    !$omp target teams loop collapse(3)                  private(uvip,uvim,vvjp,vvjm,wvkp,wvkm)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          uvip  = 0.25*( u(i  ,j,k  )+u(i  ,j+1,k  ) )*( v(i,j,k)+v(i+1,j,k) )
          uvim  = 0.25*( u(i-1,j,k  )+u(i-1,j+1,k  ) )*( v(i,j,k)+v(i-1,j,k) )
          vvjp  = 0.25*( v(i  ,j,k  )+v(i  ,j+1,k  ) )*( v(i,j,k)+v(i,j+1,k) )
          vvjm  = 0.25*( v(i  ,j,k  )+v(i  ,j-1,k  ) )*( v(i,j,k)+v(i,j-1,k) )
          wvkp  = 0.25*( w(i  ,j,k  )+w(i  ,j+1,k  ) )*( v(i,j,k)+v(i,j,k+1) )
          wvkm  = 0.25*( w(i  ,j,k-1)+w(i  ,j+1,k-1) )*( v(i,j,k)+v(i,j,k-1) )
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
  end subroutine momy_a
  !
  subroutine momz_a(nx,ny,nz,dxi,dyi,dzci,u,v,w,dwdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(0:) :: dzci
    real(rp), dimension(0:,0:,0:), intent(in   ) :: u,v,w
    real(rp), dimension( :, :, :), intent(inout) :: dwdt
    integer :: i,j,k
    real(rp) :: uwip,uwim,vwjp,vwjm,wwkp,wwkm
    !
    !$acc parallel     loop collapse(3) default(present) private(uwip,uwim,vwjp,vwjm,wwkp,wwkm) async(1)
    !$omp target teams loop collapse(3)                  private(uwip,uwim,vwjp,vwjm,wwkp,wwkm)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          uwip  = 0.25*( u(i  ,j  ,k)+u(i  ,j  ,k+1) )*( w(i,j,k)+w(i+1,j,k) )
          uwim  = 0.25*( u(i-1,j  ,k)+u(i-1,j  ,k+1) )*( w(i,j,k)+w(i-1,j,k) )
          vwjp  = 0.25*( v(i  ,j  ,k)+v(i  ,j  ,k+1) )*( w(i,j,k)+w(i,j+1,k) )
          vwjm  = 0.25*( v(i  ,j-1,k)+v(i  ,j-1,k+1) )*( w(i,j,k)+w(i,j-1,k) )
          wwkp  = 0.25*( w(i  ,j  ,k)+w(i  ,j  ,k+1) )*( w(i,j,k)+w(i,j,k+1) )
          wwkm  = 0.25*( w(i  ,j  ,k)+w(i  ,j  ,k-1) )*( w(i,j,k)+w(i,j,k-1) )
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
    !$acc parallel     loop collapse(3) default(present) private(dudxp,dudxm,dudyp,dudym,dudzp,dudzm) async(1)
    !$omp target teams loop collapse(3)                  private(dudxp,dudxm,dudyp,dudym,dudzp,dudzm)
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
    !$acc parallel     loop collapse(3) default(present) private(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm) async(1)
    !$omp target teams loop collapse(3)                  private(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm)
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
                        (dvdxp-dvdxm)*visc*dxi + &
                        (dvdyp-dvdym)*visc*dyi + &
                        (dvdzp-dvdzm)*visc*dzfi(k)
        end do
      end do
    end do
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
    !$acc parallel     loop collapse(3) default(present) private(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm) async(1)
    !$omp target teams loop collapse(3)                  private(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm)
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
                        (dwdxp-dwdxm)*visc*dxi + &
                        (dwdyp-dwdym)*visc*dyi + &
                        (dwdzp-dwdzm)*visc*dzci(k)
        end do
      end do
    end do
  end subroutine momz_d
  !
  subroutine momx_p(nx,ny,nz,dxi,bforce,p,dudt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi
    real(rp), intent(in) :: bforce
    real(rp), dimension(0:,0:,0:), intent(in ) :: p
    real(rp), dimension( :, :, :), intent(out) :: dudt
    integer :: i,j,k
    !
    !$acc parallel     loop collapse(3) default(present) async(1)
    !$omp target teams loop collapse(3)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dudt(i,j,k) = - dxi*( p(i+1,j,k)-p(i,j,k) ) + bforce
        end do
      end do
    end do
  end subroutine momx_p
  !
  subroutine momy_p(nx,ny,nz,dyi,bforce,p,dvdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dyi
    real(rp), intent(in) :: bforce
    real(rp), dimension(0:,0:,0:), intent(in ) :: p
    real(rp), dimension( :, :, :), intent(out) :: dvdt
    integer :: i,j,k
    !
    !$acc parallel     loop collapse(3) default(present) async(1)
    !$omp target teams loop collapse(3)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dvdt(i,j,k) = - dyi*( p(i,j+1,k)-p(i,j,k) ) + bforce
        end do
      end do
    end do
  end subroutine momy_p
  !
  subroutine momz_p(nx,ny,nz,dzci,bforce,p,dwdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in), dimension(0:) :: dzci
    real(rp), intent(in) :: bforce
    real(rp), dimension(0:,0:,0:), intent(in ) :: p
    real(rp), dimension( :, :, :), intent(out) :: dwdt
    integer :: i,j,k
    !
    !$acc parallel     loop collapse(3) default(present) async(1)
    !$omp target teams loop collapse(3)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dwdt(i,j,k) = - dzci(k)*( p(i,j,k+1)-p(i,j,k) ) + bforce
        end do
      end do
    end do
  end subroutine momz_p
  !
  subroutine momx_d_z(nx,ny,nz,dzci,dzfi,visc,u,dudt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: visc
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in   ) :: u
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    real(rp) :: dudzp,dudzm
    integer :: i,j,k
    !
    !$acc parallel     loop collapse(3) default(present) private(dudzp,dudzm) async(1)
    !$omp target teams loop collapse(3)                  private(dudzp,dudzm)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dudzp = (u(i,j,k+1)-u(i,j,k))*dzci(k  )
          dudzm = (u(i,j,k)-u(i,j,k-1))*dzci(k-1)
          dudt(i,j,k) = dudt(i,j,k) + &
                        (dudzp-dudzm)*visc*dzfi(k)
        end do
      end do
    end do
  end subroutine momx_d_z
  !
  subroutine momy_d_z(nx,ny,nz,dzci,dzfi,visc,v,dvdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: visc
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in   ) :: v
    real(rp), dimension( :, :, :), intent(inout) :: dvdt
    real(rp) :: dvdzp,dvdzm
    integer :: i,j,k
    !
    !$acc parallel     loop collapse(3) default(present) private(dvdzp,dvdzm) async(1)
    !$omp target teams loop collapse(3)                  private(dvdzp,dvdzm)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dvdzp = (v(i,j,k+1)-v(i,j,k))*dzci(k  )
          dvdzm = (v(i,j,k)-v(i,j,k-1))*dzci(k-1)
          dvdt(i,j,k) = dvdt(i,j,k) + &
                        (dvdzp-dvdzm)*visc*dzfi(k)
        end do
      end do
    end do
  end subroutine momy_d_z
  !
  subroutine momz_d_z(nx,ny,nz,dzci,dzfi,visc,w,dwdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: visc
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in   ) :: w
    real(rp), dimension( :, :, :), intent(inout) :: dwdt
    integer :: i,j,k
    real(rp) :: dwdzp,dwdzm
    !
    !$acc parallel     loop collapse(3) default(present) private(dwdzp,dwdzm) async(1)
    !$omp target teams loop collapse(3)                  private(dwdzp,dwdzm)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dwdzp = (w(i,j,k+1)-w(i,j,k))*dzfi(k+1)
          dwdzm = (w(i,j,k)-w(i,j,k-1))*dzfi(k  )
          dwdt(i,j,k) = dwdt(i,j,k) + &
                        (dwdzp-dwdzm)*visc*dzci(k)
        end do
      end do
    end do
  end subroutine momz_d_z
  !
  subroutine momx_d_xy(nx,ny,nz,dxi,dyi,visc,u,dudt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,visc
    real(rp), dimension(0:,0:,0:), intent(in   ) :: u
    real(rp), dimension( :, :, :), intent(inout) :: dudt
    real(rp) :: dudxp,dudxm,dudyp,dudym
    integer :: i,j,k
    !
    !$acc parallel     loop collapse(3) default(present) private(dudxp,dudxm,dudyp,dudym) async(1)
    !$omp target teams loop collapse(3)                  private(dudxp,dudxm,dudyp,dudym)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dudxp = (u(i+1,j,k)-u(i,j,k))*dxi
          dudxm = (u(i,j,k)-u(i-1,j,k))*dxi
          dudyp = (u(i,j+1,k)-u(i,j,k))*dyi
          dudym = (u(i,j,k)-u(i,j-1,k))*dyi
          dudt(i,j,k) = dudt(i,j,k) + &
                        (dudxp-dudxm)*visc*dxi + &
                        (dudyp-dudym)*visc*dyi
        end do
      end do
    end do
  end subroutine momx_d_xy
  !
  subroutine momy_d_xy(nx,ny,nz,dxi,dyi,visc,v,dvdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,visc
    real(rp), dimension(0:,0:,0:), intent(in   ) :: v
    real(rp), dimension( :, :, :), intent(inout) :: dvdt
    real(rp) :: dvdxp,dvdxm,dvdyp,dvdym
    integer :: i,j,k
    !
    !$acc parallel     loop collapse(3) default(present) private(dvdxp,dvdxm,dvdyp,dvdym) async(1)
    !$omp target teams loop collapse(3)                  private(dvdxp,dvdxm,dvdyp,dvdym)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dvdxp = (v(i+1,j,k)-v(i,j,k))*dxi
          dvdxm = (v(i,j,k)-v(i-1,j,k))*dxi
          dvdyp = (v(i,j+1,k)-v(i,j,k))*dyi
          dvdym = (v(i,j,k)-v(i,j-1,k))*dyi
          dvdt(i,j,k) = dvdt(i,j,k) + &
                        (dvdxp-dvdxm)*visc*dxi + &
                        (dvdyp-dvdym)*visc*dyi
        end do
      end do
    end do
  end subroutine momy_d_xy
  !
  subroutine momz_d_xy(nx,ny,nz,dxi,dyi,visc,w,dwdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,visc
    real(rp), dimension(0:,0:,0:), intent(in   ) :: w
    real(rp), dimension( :, :, :), intent(inout) :: dwdt
    integer :: i,j,k
    real(rp) :: dwdxp,dwdxm,dwdyp,dwdym
    !
    !$acc parallel     loop collapse(3) default(present) private(dwdxp,dwdxm,dwdyp,dwdym) async(1)
    !$omp target teams loop collapse(3)                  private(dwdxp,dwdxm,dwdyp,dwdym)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dwdxp = (w(i+1,j,k)-w(i,j,k))*dxi
          dwdxm = (w(i,j,k)-w(i-1,j,k))*dxi
          dwdyp = (w(i,j+1,k)-w(i,j,k))*dyi
          dwdym = (w(i,j,k)-w(i,j-1,k))*dyi
          dwdt(i,j,k) = dwdt(i,j,k) + &
                        (dwdxp-dwdxm)*visc*dxi + &
                        (dwdyp-dwdym)*visc*dyi
        end do
      end do
    end do
  end subroutine momz_d_xy
  !
  subroutine cmpt_wallshear(n,is_cmpt,is_bound,l,dli,dzci,dzfi,visc,u,v,w,taux,tauy,tauz)
    use mod_param, only: cbcpre
    implicit none
    integer , intent(in ), dimension(3) :: n
    logical , intent(in ), dimension(    3) :: is_cmpt
    logical , intent(in ), dimension(0:1,3) :: is_bound
    real(rp), intent(in ), dimension(3)     :: l,dli
    real(rp), intent(in ), dimension(0:)    :: dzci,dzfi
    real(rp), intent(in )                   :: visc
    real(rp), intent(in ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out), dimension(0:1,3) :: taux,tauy,tauz
    real(rp) :: dudyp,dudym,dudzp,dudzm, &
                dvdxp,dvdxm,dvdzp,dvdzm, &
                dwdxp,dwdxm,dwdyp,dwdym
    !
    ! n.b.: replace scalars with reduction of tau(1:3,1:3) once the
    !       nvfortran bug for array reductions on Pascal architectures
    !       is solved; this subroutine is not used in production anyway
    !
    real(rp) :: tau21p,tau31p,tau12p,tau32p,tau13p,tau23p, &
                tau21m,tau31m,tau12m,tau32m,tau13m,tau23m
    integer :: i,j,k,nx,ny,nz
    real(rp) :: dxi,dyi,lx,ly,lz
    real(rp) :: tau(0:1,3,3)
    !
    nx = n(1); ny = n(2); nz = n(3)
    dxi = dli(1); dyi = dli(2)
    lx = l(1); ly = l(2); lz = l(3)
    tau21p = 0.; tau21m = 0.
    tau31p = 0.; tau31m = 0.
    if(is_cmpt(1)) then
      !$acc        data copy(      tau21p,tau21m,tau31p,tau31m) async(1)
      !$omp target data map(tofrom:tau21p,tau21m,tau31p,tau31m)
      if(is_bound(0,2).and.cbcpre(0,2)//cbcpre(1,2) /= 'PP') then
        !$acc parallel     loop collapse(2) default(present) private(dudyp) reduction(+:tau21p) async(1)
        !$omp target teams loop collapse(2)                  private(dudyp) reduction(+:tau21p)
        do k=1,nz
          do i=1,nx
            dudyp  = (u(i,1 ,k)-u(i,0   ,k))*dyi*visc
            tau21p = tau21p + dudyp/(dxi*dzfi(k)*lx*lz)
          end do
        end do
      end if
      if(is_bound(1,2).and.cbcpre(0,2)//cbcpre(1,2) /= 'PP') then
        !$acc parallel     loop collapse(2) default(present) private(dudym) reduction(+:tau21m) async(1)
        !$omp target teams loop collapse(2)                  private(dudym) reduction(+:tau21m)
        do k=1,nz
          do i=1,nx
            dudym  = (u(i,ny,k)-u(i,ny+1,k))*dyi*visc
            tau21m = tau21m + dudym/(dxi*dzfi(k)*lx*lz)
          end do
        end do
      end if
      if(is_bound(0,3).and.cbcpre(0,3)//cbcpre(1,3) /= 'PP') then
        !$acc parallel     loop collapse(2) default(present) private(dudzp) reduction(+:tau31p) async(1)
        !$omp target teams loop collapse(2)                  private(dudzp) reduction(+:tau31p)
        do j=1,ny
          do i=1,nx
            dudzp  = (u(i,j,1 )-u(i,j,0   ))*dzci(0)*visc
            tau31p = tau31p + dudzp/(dxi*dyi*lx*ly)
          end do
        end do
      end if
      if(is_bound(1,3).and.cbcpre(0,3)//cbcpre(1,3) /= 'PP') then
        !$acc parallel     loop collapse(2) default(present) private(dudzm) reduction(+:tau31m) async(1)
        !$omp target teams loop collapse(2)                  private(dudzm) reduction(+:tau31m)
        do j=1,ny
          do i=1,nx
            dudzm  = (u(i,j,nz)-u(i,j,nz+1))*dzci(nz)*visc
            tau31m = tau31m + dudzm/(dxi*dyi*lx*ly)
          end do
        end do
      end if
      !$omp end target data
      !$acc end        data
    end if
    !
    tau12p = 0.; tau12m = 0.
    tau32p = 0.; tau32m = 0.
    if(is_cmpt(2)) then
      !$acc        data copy(      tau12p,tau12m,tau32p,tau32m) async(1)
      !$omp target data map(tofrom:tau12p,tau12m,tau32p,tau32m)
      if(is_bound(0,1).and.cbcpre(0,1)//cbcpre(1,1) /= 'PP') then
        !$acc parallel     loop collapse(2) default(present) private(dvdxp) reduction(+:tau12p) async(1)
        !$omp target teams loop collapse(2)                  private(dvdxp) reduction(+:tau12p)
        do k=1,nz
          do j=1,ny
            dvdxp  = (v(1 ,j,k)-v(0  ,j,k))*dxi*visc
            tau12p = tau12p + dvdxp/(dyi*dzfi(k)*ly*lz)
          end do
        end do
      end if
      if(is_bound(1,1).and.cbcpre(0,1)//cbcpre(1,1) /= 'PP') then
        !$acc parallel     loop collapse(2) default(present) private(dvdxm) reduction(+:tau12m) async(1)
        !$omp target teams loop collapse(2)                  private(dvdxm) reduction(+:tau12m)
        do k=1,nz
          do j=1,ny
            dvdxm  = (v(nx,j,k)-v(nx+1,j,k))*dxi*visc
            tau12m = tau12m + dvdxm/(dyi*dzfi(k)*ly*lz)
          end do
        end do
      end if
      if(is_bound(0,3).and.cbcpre(0,3)//cbcpre(1,3) /= 'PP') then
        !$acc parallel     loop collapse(2) default(present) private(dvdzp) reduction(+:tau32p) async(1)
        !$omp target teams loop collapse(2)                  private(dvdzp) reduction(+:tau32p)
        do j=1,ny
          do i=1,nx
            dvdzp  = (v(i,j,1 )-v(i,j,0   ))*dzci(0)*visc
            tau32p = tau32p + dvdzp/(dxi*dyi*lx*ly)
          end do
        end do
      end if
      if(is_bound(1,3).and.cbcpre(0,3)//cbcpre(1,3) /= 'PP') then
        !$acc parallel     loop collapse(2) default(present) private(dvdzm) reduction(+:tau32m) async(1)
        !$omp target teams loop collapse(2)                  private(dvdzm) reduction(+:tau32m)
        do j=1,ny
          do i=1,nx
            dvdzm  = (v(i,j,nz)-v(i,j,nz+1))*dzci(nz)*visc
            tau32m = tau32m + dvdzm/(dxi*dyi*lx*ly)
          end do
        end do
      end if
      !$omp end target data
      !$acc end        data
    end if
    !
    tau13p = 0.; tau13m = 0.
    tau23p = 0.; tau23m = 0.
    if(is_cmpt(3)) then
      !$acc        data copy(      tau13p,tau13m,tau23p,tau23m) async(1)
      !$omp target data map(tofrom:tau13p,tau13m,tau23p,tau23m)
      if(is_bound(0,1).and.cbcpre(0,1)//cbcpre(1,1) /= 'PP') then
        !$acc parallel     loop collapse(2) default(present) private(dwdxp) reduction(+:tau13p) async(1)
        !$omp target teams loop collapse(2)                  private(dwdxp) reduction(+:tau13p)
        do k=1,nz
          do j=1,ny
            dwdxp  = (w(1 ,j,k)-w(0   ,j,k))*dxi*visc
            tau13p = tau13p + dwdxp/(dyi*dzfi(k)*ly*lz)
          end do
        end do
      end if
      if(is_bound(1,1).and.cbcpre(0,1)//cbcpre(1,1) /= 'PP') then
        !$acc parallel     loop collapse(2) default(present) private(dwdxm) reduction(+:tau13m) async(1)
        !$omp target teams loop collapse(2)                  private(dwdxm) reduction(+:tau13m)
        do k=1,nz
          do j=1,ny
            dwdxm  = (w(nx,j,k)-w(nx+1,j,k))*dxi*visc
            tau13m = tau13m + dwdxm/(dyi*dzfi(k)*ly*lz)
          end do
        end do
      end if
      if(is_bound(0,2).and.cbcpre(0,2)//cbcpre(1,2) /= 'PP') then
        !$acc parallel     loop collapse(2) default(present) private(dwdyp) reduction(+:tau23p) async(1)
        !$omp target teams loop collapse(2)                  private(dwdyp) reduction(+:tau23p)
        do k=1,nz
          do i=1,nx
            dwdyp  = (w(i,1 ,k)-w(i,0   ,k))*dyi*visc
            tau23p = tau23p + dwdyp/(dxi*dzfi(k)*lx*lz)
          end do
        end do
      end if
      if(is_bound(1,2).and.cbcpre(0,2)//cbcpre(1,2) /= 'PP') then
        !$acc parallel     loop collapse(2) default(present) private(dwdym) reduction(+:tau23m) async(1)
        !$omp target teams loop collapse(2)                  private(dwdym) reduction(+:tau23m)
        do k=1,nz
          do i=1,nx
            dwdym  = (w(i,ny,k)-w(i,ny+1,k))*dyi*visc
            tau23m = tau23m + dwdym/(dxi*dzfi(k)*lx*lz)
          end do
        end do
      end if
      !$omp end target data
      !$acc end        data
    end if
    !$acc wait(1)
    tau(:,:,:) = 0.
    tau(0:1,2,1) = [tau21m,tau21p]
    tau(0:1,3,1) = [tau31m,tau31p]
    tau(0:1,1,2) = [tau12m,tau12p]
    tau(0:1,3,2) = [tau32m,tau32p]
    tau(0:1,1,3) = [tau13m,tau13p]
    tau(0:1,2,3) = [tau23m,tau23p]
    call MPI_ALLREDUCE(MPI_IN_PLACE,tau(0,1,1),9*2,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    taux(:,:) = tau(:,:,1)
    tauy(:,:) = tau(:,:,2)
    tauz(:,:) = tau(:,:,3)
  end subroutine cmpt_wallshear
  !
  subroutine bulk_forcing(n,is_forced,f,u,v,w)
    implicit none
    integer , intent(in   ), dimension(3) :: n
    logical , intent(in   ), dimension(3) :: is_forced
    real(rp), intent(in   ), dimension(3) :: f
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    real(rp) :: ff
    integer  :: i,j,k
    if(is_forced(1)) then
      ff = f(1)
      !$acc parallel     loop collapse(3) default(present) async(1)
      !$omp target teams loop collapse(3)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            u(i,j,k) = u(i,j,k) + ff
          end do
        end do
      end do
    end if
    if(is_forced(2)) then
      ff = f(2)
      !$acc parallel     loop collapse(3) default(present) async(1)
      !$omp target teams loop collapse(3)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            v(i,j,k) = v(i,j,k) + ff
          end do
        end do
      end do
    end if
    if(is_forced(3)) then
      ff = f(3)
      !$acc parallel     loop collapse(3) default(present) async(1)
      !$omp target teams loop collapse(3)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            w(i,j,k) = w(i,j,k) + ff
          end do
        end do
      end do
    end if
  end subroutine bulk_forcing
  !
  subroutine mom_xyz_ad(nx,ny,nz,dxi,dyi,dzci,dzfi,visc,u,v,w,dudt,dvdt,dwdt,dudtd,dvdtd,dwdtd)
    use mod_param, only: is_impdiff,is_impdiff_1d
    !
    ! lump all r.h.s. of momentum terms (excluding pressure) into a single fast kernel
    !
    implicit none
    integer , intent(in   ) :: nx,ny,nz
    real(rp), intent(in   ) :: dxi,dyi
    real(rp), intent(in   ), dimension(0:) :: dzci,dzfi
    real(rp), intent(in   ) :: visc
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(inout), dimension( :, :, :) :: dudt,dvdt,dwdt
    real(rp), intent(inout), dimension( :, :, :), optional :: dudtd,dvdtd,dwdtd
    integer  :: i,j,k
    real(rp) :: u_ccm,u_pcm,u_cpm,u_cmc,u_pmc,u_mcc,u_ccc,u_pcc,u_mpc,u_cpc,u_cmp,u_mcp,u_ccp, &
                v_ccm,v_pcm,v_cpm,v_cmc,v_pmc,v_mcc,v_ccc,v_pcc,v_mpc,v_cpc,v_cmp,v_mcp,v_ccp, &
                w_ccm,w_pcm,w_cpm,w_cmc,w_pmc,w_mcc,w_ccc,w_pcc,w_mpc,w_cpc,w_cmp,w_mcp,w_ccp
    real(rp) :: uuip,uuim,vujp,vujm,wukp,wukm, &
                uvip,uvim,vvjp,vvjm,wvkp,wvkm, &
                uwip,uwim,vwjp,vwjm,wwkp,wwkm
    real(rp) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm, &
                dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm, &
                dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
    real(rp) :: dudt_s ,dvdt_s ,dwdt_s , &
                dudtd_s,dvdtd_s,dwdtd_s, &
                dudtd_xy_s,dudtd_z_s   , &
                dvdtd_xy_s,dvdtd_z_s   , &
                dwdtd_xy_s,dwdtd_z_s
    !
#if !defined(_LOOP_UNSWITCHING)
    !$acc parallel     loop collapse(3) default(present) async(1) &
    !$acc private(u_ccm,u_pcm,u_cpm,u_cmc,u_pmc,u_mcc,u_ccc,u_pcc,u_mpc,u_cpc,u_cmp,u_mcp,u_ccp) &
    !$acc private(v_ccm,v_pcm,v_cpm,v_cmc,v_pmc,v_mcc,v_ccc,v_pcc,v_mpc,v_cpc,v_cmp,v_mcp,v_ccp) &
    !$acc private(w_ccm,w_pcm,w_cpm,w_cmc,w_pmc,w_mcc,w_ccc,w_pcc,w_mpc,w_cpc,w_cmp,w_mcp,w_ccp) &
    !$acc private(uuip,uuim,vujp,vujm,wukp,wukm) &
    !$acc private(uvip,uvim,vvjp,vvjm,wvkp,wvkm) &
    !$acc private(uwip,uwim,vwjp,vwjm,wwkp,wwkm) &
    !$acc private(dudxp,dudxm,dudyp,dudym,dudzp,dudzm) &
    !$acc private(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm) &
    !$acc private(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm) &
    !$acc private(dudt_s ,dvdt_s ,dwdt_s ) &
    !$acc private(dudtd_s,dvdtd_s,dwdtd_s) &
    !$acc private(dudtd_xy_s,dudtd_z_s) &
    !$acc private(dvdtd_xy_s,dvdtd_z_s) &
    !$acc private(dwdtd_xy_s,dwdtd_z_s)
    !$omp target teams loop collapse(3) &
    !$omp private(u_ccm,u_pcm,u_cpm,u_cmc,u_pmc,u_mcc,u_ccc,u_pcc,u_mpc,u_cpc,u_cmp,u_mcp,u_ccp) &
    !$omp private(v_ccm,v_pcm,v_cpm,v_cmc,v_pmc,v_mcc,v_ccc,v_pcc,v_mpc,v_cpc,v_cmp,v_mcp,v_ccp) &
    !$omp private(w_ccm,w_pcm,w_cpm,w_cmc,w_pmc,w_mcc,w_ccc,w_pcc,w_mpc,w_cpc,w_cmp,w_mcp,w_ccp) &
    !$omp private(uuip,uuim,vujp,vujm,wukp,wukm) &
    !$omp private(uvip,uvim,vvjp,vvjm,wvkp,wvkm) &
    !$omp private(uwip,uwim,vwjp,vwjm,wwkp,wwkm) &
    !$omp private(dudxp,dudxm,dudyp,dudym,dudzp,dudzm) &
    !$omp private(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm) &
    !$omp private(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm) &
    !$omp private(dudt_s ,dvdt_s ,dwdt_s ) &
    !$omp private(dudtd_s,dvdtd_s,dwdtd_s) &
    !$omp private(dudtd_xy_s,dudtd_z_s) &
    !$omp private(dvdtd_xy_s,dvdtd_z_s) &
    !$omp private(dwdtd_xy_s,dwdtd_z_s)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          ! touch u,v,w sequentially
          !
          u_ccm = u(i  ,j  ,k-1)
          u_pcm = u(i+1,j  ,k-1)
          u_cpm = u(i  ,j+1,k-1)
          u_cmc = u(i  ,j-1,k  )
          u_pmc = u(i+1,j-1,k  )
          u_mcc = u(i-1,j  ,k  )
          u_ccc = u(i  ,j  ,k  )
          u_pcc = u(i+1,j  ,k  )
          u_mpc = u(i-1,j+1,k  )
          u_cpc = u(i  ,j+1,k  )
          u_cmp = u(i  ,j-1,k+1)
          u_mcp = u(i-1,j  ,k+1)
          u_ccp = u(i  ,j  ,k+1)
          !
          v_ccm = v(i  ,j  ,k-1)
          v_pcm = v(i+1,j  ,k-1)
          v_cpm = v(i  ,j+1,k-1)
          v_cmc = v(i  ,j-1,k  )
          v_pmc = v(i+1,j-1,k  )
          v_mcc = v(i-1,j  ,k  )
          v_ccc = v(i  ,j  ,k  )
          v_pcc = v(i+1,j  ,k  )
          v_mpc = v(i-1,j+1,k  )
          v_cpc = v(i  ,j+1,k  )
          v_cmp = v(i  ,j-1,k+1)
          v_mcp = v(i-1,j  ,k+1)
          v_ccp = v(i  ,j  ,k+1)
          !
          w_ccm = w(i  ,j  ,k-1)
          w_pcm = w(i+1,j  ,k-1)
          w_cpm = w(i  ,j+1,k-1)
          w_cmc = w(i  ,j-1,k  )
          w_pmc = w(i+1,j-1,k  )
          w_mcc = w(i-1,j  ,k  )
          w_ccc = w(i  ,j  ,k  )
          w_pcc = w(i+1,j  ,k  )
          w_mpc = w(i-1,j+1,k  )
          w_cpc = w(i  ,j+1,k  )
          w_cmp = w(i  ,j-1,k+1)
          w_mcp = w(i-1,j  ,k+1)
          w_ccp = w(i  ,j  ,k+1)
          !
          ! x diffusion
          !
          dudxp = (u_pcc-u_ccc)*dxi
          dudxm = (u_ccc-u_mcc)*dxi
          dudyp = (u_cpc-u_ccc)*dyi
          dudym = (u_ccc-u_cmc)*dyi
          dudzp = (u_ccp-u_ccc)*dzci(k  )
          dudzm = (u_ccc-u_ccm)*dzci(k-1)
          !
          ! x advection
          !
          uuip  = 0.25*(u_pcc+u_ccc)*(u_ccc+u_pcc)
          uuim  = 0.25*(u_mcc+u_ccc)*(u_ccc+u_mcc)
          vujp  = 0.25*(v_pcc+v_ccc)*(u_ccc+u_cpc)
          vujm  = 0.25*(v_pmc+v_cmc)*(u_ccc+u_cmc)
          wukp  = 0.25*(w_pcc+w_ccc)*(u_ccc+u_ccp)
          wukm  = 0.25*(w_pcm+w_ccm)*(u_ccc+u_ccm)
          !
          dudtd_xy_s = &
                         visc*(dudxp-dudxm)*dxi + &
                         visc*(dudyp-dudym)*dyi
          dudtd_z_s  = &
                         visc*(dudzp-dudzm)*dzfi(k)
          dudt_s     = &
                             -(uuip -uuim )*dxi - &
                              (vujp -vujm )*dyi - &
                              (wukp -wukm )*dzfi(k)
          !
          ! y diffusion
          !
          dvdxp = (v_pcc-v_ccc)*dxi
          dvdxm = (v_ccc-v_mcc)*dxi
          dvdyp = (v_cpc-v_ccc)*dyi
          dvdym = (v_ccc-v_cmc)*dyi
          dvdzp = (v_ccp-v_ccc)*dzci(k  )
          dvdzm = (v_ccc-v_ccm)*dzci(k-1)
          !
          ! y advection
          !
          uvip  = 0.25*(u_ccc+u_cpc)*(v_ccc+v_pcc)
          uvim  = 0.25*(u_mcc+u_mpc)*(v_ccc+v_mcc)
          vvjp  = 0.25*(v_ccc+v_cpc)*(v_ccc+v_cpc)
          vvjm  = 0.25*(v_ccc+v_cmc)*(v_ccc+v_cmc)
          wvkp  = 0.25*(w_ccc+w_cpc)*(v_ccc+v_ccp)
          wvkm  = 0.25*(w_ccm+w_cpm)*(v_ccc+v_ccm)
          !
          dvdtd_xy_s = &
                         visc*(dvdxp-dvdxm)*dxi + &
                         visc*(dvdyp-dvdym)*dyi
          dvdtd_z_s  = &
                         visc*(dvdzp-dvdzm)*dzfi(k)
          dvdt_s     = &
                             -(uvip -uvim )*dxi - &
                              (vvjp -vvjm )*dyi - &
                              (wvkp -wvkm )*dzfi(k)
          !
          ! z diffusion
          !
          dwdxp = (w_pcc-w_ccc)*dxi
          dwdxm = (w_ccc-w_mcc)*dxi
          dwdyp = (w_cpc-w_ccc)*dyi
          dwdym = (w_ccc-w_cmc)*dyi
          dwdzp = (w_ccp-w_ccc)*dzfi(k+1)
          dwdzm = (w_ccc-w_ccm)*dzfi(k  )
          !
          ! z advection
          !
          uwip  = 0.25*(u_ccc+u_ccp)*(w_ccc+w_pcc)
          uwim  = 0.25*(u_mcc+u_mcp)*(w_ccc+w_mcc)
          vwjp  = 0.25*(v_ccc+v_ccp)*(w_ccc+w_cpc)
          vwjm  = 0.25*(v_cmc+v_cmp)*(w_ccc+w_cmc)
          wwkp  = 0.25*(w_ccc+w_ccp)*(w_ccc+w_ccp)
          wwkm  = 0.25*(w_ccc+w_ccm)*(w_ccc+w_ccm)
          !
          dwdtd_xy_s =  &
                          visc*(dwdxp-dwdxm)*dxi + &
                          visc*(dwdyp-dwdym)*dyi
          dwdtd_z_s  =  &
                          visc*(dwdzp-dwdzm)*dzci(k)
          dwdt_s     =  &
                              -(uwip -uwim )*dxi - &
                               (vwjp -vwjm )*dyi - &
                               (wwkp -wwkm )*dzci(k)
          if(is_impdiff) then
            if(is_impdiff_1d) then
              dudt_s = dudt_s + dudtd_xy_s
              dvdt_s = dvdt_s + dvdtd_xy_s
              dwdt_s = dwdt_s + dwdtd_xy_s
              dudtd_s = dudtd_z_s
              dvdtd_s = dvdtd_z_s
              dwdtd_s = dwdtd_z_s
            else
              dudtd_s = dudtd_xy_s + dudtd_z_s
              dvdtd_s = dvdtd_xy_s + dvdtd_z_s
              dwdtd_s = dwdtd_xy_s + dwdtd_z_s
            end if
            dudt( i,j,k) = dudt_s
            dvdt( i,j,k) = dvdt_s
            dwdt( i,j,k) = dwdt_s
            dudtd(i,j,k) = dudtd_s
            dvdtd(i,j,k) = dvdtd_s
            dwdtd(i,j,k) = dwdtd_s
          else
            dudt_s = dudt_s + dudtd_xy_s + dudtd_z_s
            dvdt_s = dvdt_s + dvdtd_xy_s + dvdtd_z_s
            dwdt_s = dwdt_s + dwdtd_xy_s + dwdtd_z_s
            dudt(i,j,k) = dudt_s
            dvdt(i,j,k) = dvdt_s
            dwdt(i,j,k) = dwdt_s
          end if
        end do
      end do
    end do
#else
    if(.not.is_impdiff) then
      !$acc parallel     loop collapse(3) default(present) async(1) &
      !$acc private(u_ccm,u_pcm,u_cpm,u_cmc,u_pmc,u_mcc,u_ccc,u_pcc,u_mpc,u_cpc,u_cmp,u_mcp,u_ccp) &
      !$acc private(v_ccm,v_pcm,v_cpm,v_cmc,v_pmc,v_mcc,v_ccc,v_pcc,v_mpc,v_cpc,v_cmp,v_mcp,v_ccp) &
      !$acc private(w_ccm,w_pcm,w_cpm,w_cmc,w_pmc,w_mcc,w_ccc,w_pcc,w_mpc,w_cpc,w_cmp,w_mcp,w_ccp) &
      !$acc private(uuip,uuim,vujp,vujm,wukp,wukm) &
      !$acc private(uvip,uvim,vvjp,vvjm,wvkp,wvkm) &
      !$acc private(uwip,uwim,vwjp,vwjm,wwkp,wwkm) &
      !$acc private(dudxp,dudxm,dudyp,dudym,dudzp,dudzm) &
      !$acc private(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm) &
      !$acc private(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm) &
      !$acc private(dudt_s,dvdt_s,dwdt_s) &
      !$acc private(dudtd_s,dvdtd_s,dwdtd_s) &
      !$acc private(dudtd_xy_s,dudtd_z_s) &
      !$acc private(dvdtd_xy_s,dvdtd_z_s) &
      !$acc private(dwdtd_xy_s,dwdtd_z_s)
      !$omp target teams loop collapse(3) &
      !$omp private(u_ccm,u_pcm,u_cpm,u_cmc,u_pmc,u_mcc,u_ccc,u_pcc,u_mpc,u_cpc,u_cmp,u_mcp,u_ccp) &
      !$omp private(v_ccm,v_pcm,v_cpm,v_cmc,v_pmc,v_mcc,v_ccc,v_pcc,v_mpc,v_cpc,v_cmp,v_mcp,v_ccp) &
      !$omp private(w_ccm,w_pcm,w_cpm,w_cmc,w_pmc,w_mcc,w_ccc,w_pcc,w_mpc,w_cpc,w_cmp,w_mcp,w_ccp) &
      !$omp private(uuip,uuim,vujp,vujm,wukp,wukm) &
      !$omp private(uvip,uvim,vvjp,vvjm,wvkp,wvkm) &
      !$omp private(uwip,uwim,vwjp,vwjm,wwkp,wwkm) &
      !$omp private(dudxp,dudxm,dudyp,dudym,dudzp,dudzm) &
      !$omp private(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm) &
      !$omp private(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm) &
      !$omp private(dudt_s,dvdt_s,dwdt_s) &
      !$omp private(dudtd_s,dvdtd_s,dwdtd_s) &
      !$omp private(dudtd_xy_s,dudtd_z_s) &
      !$omp private(dvdtd_xy_s,dvdtd_z_s) &
      !$omp private(dwdtd_xy_s,dwdtd_z_s)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            u_ccm = u(i  ,j  ,k-1)
            u_pcm = u(i+1,j  ,k-1)
            u_cpm = u(i  ,j+1,k-1)
            u_cmc = u(i  ,j-1,k  )
            u_pmc = u(i+1,j-1,k  )
            u_mcc = u(i-1,j  ,k  )
            u_ccc = u(i  ,j  ,k  )
            u_pcc = u(i+1,j  ,k  )
            u_mpc = u(i-1,j+1,k  )
            u_cpc = u(i  ,j+1,k  )
            u_cmp = u(i  ,j-1,k+1)
            u_mcp = u(i-1,j  ,k+1)
            u_ccp = u(i  ,j  ,k+1)
            v_ccm = v(i  ,j  ,k-1)
            v_pcm = v(i+1,j  ,k-1)
            v_cpm = v(i  ,j+1,k-1)
            v_cmc = v(i  ,j-1,k  )
            v_pmc = v(i+1,j-1,k  )
            v_mcc = v(i-1,j  ,k  )
            v_ccc = v(i  ,j  ,k  )
            v_pcc = v(i+1,j  ,k  )
            v_mpc = v(i-1,j+1,k  )
            v_cpc = v(i  ,j+1,k  )
            v_cmp = v(i  ,j-1,k+1)
            v_mcp = v(i-1,j  ,k+1)
            v_ccp = v(i  ,j  ,k+1)
            w_ccm = w(i  ,j  ,k-1)
            w_pcm = w(i+1,j  ,k-1)
            w_cpm = w(i  ,j+1,k-1)
            w_cmc = w(i  ,j-1,k  )
            w_pmc = w(i+1,j-1,k  )
            w_mcc = w(i-1,j  ,k  )
            w_ccc = w(i  ,j  ,k  )
            w_pcc = w(i+1,j  ,k  )
            w_mpc = w(i-1,j+1,k  )
            w_cpc = w(i  ,j+1,k  )
            w_cmp = w(i  ,j-1,k+1)
            w_mcp = w(i-1,j  ,k+1)
            w_ccp = w(i  ,j  ,k+1)
            dudxp = (u_pcc-u_ccc)*dxi
            dudxm = (u_ccc-u_mcc)*dxi
            dudyp = (u_cpc-u_ccc)*dyi
            dudym = (u_ccc-u_cmc)*dyi
            dudzp = (u_ccp-u_ccc)*dzci(k)
            dudzm = (u_ccc-u_ccm)*dzci(k-1)
            uuip  = 0.25*(u_pcc+u_ccc)*(u_ccc+u_pcc)
            uuim  = 0.25*(u_mcc+u_ccc)*(u_ccc+u_mcc)
            vujp  = 0.25*(v_pcc+v_ccc)*(u_ccc+u_cpc)
            vujm  = 0.25*(v_pmc+v_cmc)*(u_ccc+u_cmc)
            wukp  = 0.25*(w_pcc+w_ccc)*(u_ccc+u_ccp)
            wukm  = 0.25*(w_pcm+w_ccm)*(u_ccc+u_ccm)
            dudtd_xy_s = visc*(dudxp-dudxm)*dxi + visc*(dudyp-dudym)*dyi
            dudtd_z_s  = visc*(dudzp-dudzm)*dzfi(k)
            dudt_s     = -(uuip-uuim)*dxi - (vujp-vujm)*dyi - (wukp-wukm)*dzfi(k)
            dvdxp = (v_pcc-v_ccc)*dxi
            dvdxm = (v_ccc-v_mcc)*dxi
            dvdyp = (v_cpc-v_ccc)*dyi
            dvdym = (v_ccc-v_cmc)*dyi
            dvdzp = (v_ccp-v_ccc)*dzci(k)
            dvdzm = (v_ccc-v_ccm)*dzci(k-1)
            uvip  = 0.25*(u_ccc+u_cpc)*(v_ccc+v_pcc)
            uvim  = 0.25*(u_mcc+u_mpc)*(v_ccc+v_mcc)
            vvjp  = 0.25*(v_ccc+v_cpc)*(v_ccc+v_cpc)
            vvjm  = 0.25*(v_ccc+v_cmc)*(v_ccc+v_cmc)
            wvkp  = 0.25*(w_ccc+w_cpc)*(v_ccc+v_ccp)
            wvkm  = 0.25*(w_ccm+w_cpm)*(v_ccc+v_ccm)
            dvdtd_xy_s = visc*(dvdxp-dvdxm)*dxi + visc*(dvdyp-dvdym)*dyi
            dvdtd_z_s  = visc*(dvdzp-dvdzm)*dzfi(k)
            dvdt_s     = -(uvip-uvim)*dxi - (vvjp-vvjm)*dyi - (wvkp-wvkm)*dzfi(k)
            dwdxp = (w_pcc-w_ccc)*dxi
            dwdxm = (w_ccc-w_mcc)*dxi
            dwdyp = (w_cpc-w_ccc)*dyi
            dwdym = (w_ccc-w_cmc)*dyi
            dwdzp = (w_ccp-w_ccc)*dzfi(k+1)
            dwdzm = (w_ccc-w_ccm)*dzfi(k)
            uwip  = 0.25*(u_ccc+u_ccp)*(w_ccc+w_pcc)
            uwim  = 0.25*(u_mcc+u_mcp)*(w_ccc+w_mcc)
            vwjp  = 0.25*(v_ccc+v_ccp)*(w_ccc+w_cpc)
            vwjm  = 0.25*(v_cmc+v_cmp)*(w_ccc+w_cmc)
            wwkp  = 0.25*(w_ccc+w_ccp)*(w_ccc+w_ccp)
            wwkm  = 0.25*(w_ccc+w_ccm)*(w_ccc+w_ccm)
            dwdtd_xy_s = visc*(dwdxp-dwdxm)*dxi + visc*(dwdyp-dwdym)*dyi
            dwdtd_z_s  = visc*(dwdzp-dwdzm)*dzci(k)
            dwdt_s     = -(uwip-uwim)*dxi - (vwjp-vwjm)*dyi - (wwkp-wwkm)*dzfi(k)
            dudt_s = dudt_s + dudtd_xy_s + dudtd_z_s
            dvdt_s = dvdt_s + dvdtd_xy_s + dvdtd_z_s
            dwdt_s = dwdt_s + dwdtd_xy_s + dwdtd_z_s
            dudt(i,j,k) = dudt_s
            dvdt(i,j,k) = dvdt_s
            dwdt(i,j,k) = dwdt_s
          end do
        end do
      end do
    else if(is_impdiff .and. is_impdiff_1d) then
      !$acc parallel     loop collapse(3) default(present) async(1) &
      !$acc private(u_ccm,u_pcm,u_cpm,u_cmc,u_pmc,u_mcc,u_ccc,u_pcc,u_mpc,u_cpc,u_cmp,u_mcp,u_ccp) &
      !$acc private(v_ccm,v_pcm,v_cpm,v_cmc,v_pmc,v_mcc,v_ccc,v_pcc,v_mpc,v_cpc,v_cmp,v_mcp,v_ccp) &
      !$acc private(w_ccm,w_pcm,w_cpm,w_cmc,w_pmc,w_mcc,w_ccc,w_pcc,w_mpc,w_cpc,w_cmp,w_mcp,w_ccp) &
      !$acc private(uuip,uuim,vujp,vujm,wukp,wukm) &
      !$acc private(uvip,uvim,vvjp,vvjm,wvkp,wvkm) &
      !$acc private(uwip,uwim,vwjp,vwjm,wwkp,wwkm) &
      !$acc private(dudxp,dudxm,dudyp,dudym,dudzp,dudzm) &
      !$acc private(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm) &
      !$acc private(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm) &
      !$acc private(dudt_s,dvdt_s,dwdt_s) &
      !$acc private(dudtd_s,dvdtd_s,dwdtd_s) &
      !$acc private(dudtd_xy_s,dudtd_z_s) &
      !$acc private(dvdtd_xy_s,dvdtd_z_s) &
      !$acc private(dwdtd_xy_s,dwdtd_z_s)
      !$omp target teams loop collapse(3) &
      !$omp private(u_ccm,u_pcm,u_cpm,u_cmc,u_pmc,u_mcc,u_ccc,u_pcc,u_mpc,u_cpc,u_cmp,u_mcp,u_ccp) &
      !$omp private(v_ccm,v_pcm,v_cpm,v_cmc,v_pmc,v_mcc,v_ccc,v_pcc,v_mpc,v_cpc,v_cmp,v_mcp,v_ccp) &
      !$omp private(w_ccm,w_pcm,w_cpm,w_cmc,w_pmc,w_mcc,w_ccc,w_pcc,w_mpc,w_cpc,w_cmp,w_mcp,w_ccp) &
      !$omp private(uuip,uuim,vujp,vujm,wukp,wukm) &
      !$omp private(uvip,uvim,vvjp,vvjm,wvkp,wvkm) &
      !$omp private(uwip,uwim,vwjp,vwjm,wwkp,wwkm) &
      !$omp private(dudxp,dudxm,dudyp,dudym,dudzp,dudzm) &
      !$omp private(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm) &
      !$omp private(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm) &
      !$omp private(dudt_s,dvdt_s,dwdt_s) &
      !$omp private(dudtd_s,dvdtd_s,dwdtd_s) &
      !$omp private(dudtd_xy_s,dudtd_z_s) &
      !$omp private(dvdtd_xy_s,dvdtd_z_s) &
      !$omp private(dwdtd_xy_s,dwdtd_z_s)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            u_ccm = u(i  ,j  ,k-1)
            u_pcm = u(i+1,j  ,k-1)
            u_cpm = u(i  ,j+1,k-1)
            u_cmc = u(i  ,j-1,k  )
            u_pmc = u(i+1,j-1,k  )
            u_mcc = u(i-1,j  ,k  )
            u_ccc = u(i  ,j  ,k  )
            u_pcc = u(i+1,j  ,k  )
            u_mpc = u(i-1,j+1,k  )
            u_cpc = u(i  ,j+1,k  )
            u_cmp = u(i  ,j-1,k+1)
            u_mcp = u(i-1,j  ,k+1)
            u_ccp = u(i  ,j  ,k+1)
            v_ccm = v(i  ,j  ,k-1)
            v_pcm = v(i+1,j  ,k-1)
            v_cpm = v(i  ,j+1,k-1)
            v_cmc = v(i  ,j-1,k  )
            v_pmc = v(i+1,j-1,k  )
            v_mcc = v(i-1,j  ,k  )
            v_ccc = v(i  ,j  ,k  )
            v_pcc = v(i+1,j  ,k  )
            v_mpc = v(i-1,j+1,k  )
            v_cpc = v(i  ,j+1,k  )
            v_cmp = v(i  ,j-1,k+1)
            v_mcp = v(i-1,j  ,k+1)
            v_ccp = v(i  ,j  ,k+1)
            w_ccm = w(i  ,j  ,k-1)
            w_pcm = w(i+1,j  ,k-1)
            w_cpm = w(i  ,j+1,k-1)
            w_cmc = w(i  ,j-1,k  )
            w_pmc = w(i+1,j-1,k  )
            w_mcc = w(i-1,j  ,k  )
            w_ccc = w(i  ,j  ,k  )
            w_pcc = w(i+1,j  ,k  )
            w_mpc = w(i-1,j+1,k  )
            w_cpc = w(i  ,j+1,k  )
            w_cmp = w(i  ,j-1,k+1)
            w_mcp = w(i-1,j  ,k+1)
            w_ccp = w(i  ,j  ,k+1)
            dudxp = (u_pcc-u_ccc)*dxi
            dudxm = (u_ccc-u_mcc)*dxi
            dudyp = (u_cpc-u_ccc)*dyi
            dudym = (u_ccc-u_cmc)*dyi
            dudzp = (u_ccp-u_ccc)*dzci(k)
            dudzm = (u_ccc-u_ccm)*dzci(k-1)
            uuip  = 0.25*(u_pcc+u_ccc)*(u_ccc+u_pcc)
            uuim  = 0.25*(u_mcc+u_ccc)*(u_ccc+u_mcc)
            vujp  = 0.25*(v_pcc+v_ccc)*(u_ccc+u_cpc)
            vujm  = 0.25*(v_pmc+v_cmc)*(u_ccc+u_cmc)
            wukp  = 0.25*(w_pcc+w_ccc)*(u_ccc+u_ccp)
            wukm  = 0.25*(w_pcm+w_ccm)*(u_ccc+u_ccm)
            dudtd_xy_s = visc*(dudxp-dudxm)*dxi + visc*(dudyp-dudym)*dyi
            dudtd_z_s  = visc*(dudzp-dudzm)*dzfi(k)
            dudt_s     = -(uuip-uuim)*dxi - (vujp-vujm)*dyi - (wukp-wukm)*dzfi(k)
            dvdxp = (v_pcc-v_ccc)*dxi
            dvdxm = (v_ccc-v_mcc)*dxi
            dvdyp = (v_cpc-v_ccc)*dyi
            dvdym = (v_ccc-v_cmc)*dyi
            dvdzp = (v_ccp-v_ccc)*dzci(k)
            dvdzm = (v_ccc-v_ccm)*dzci(k-1)
            uvip  = 0.25*(u_ccc+u_cpc)*(v_ccc+v_pcc)
            uvim  = 0.25*(u_mcc+u_mpc)*(v_ccc+v_mcc)
            vvjp  = 0.25*(v_ccc+v_cpc)*(v_ccc+v_cpc)
            vvjm  = 0.25*(v_ccc+v_cmc)*(v_ccc+v_cmc)
            wvkp  = 0.25*(w_ccc+w_cpc)*(v_ccc+v_ccp)
            wvkm  = 0.25*(w_ccm+w_cpm)*(v_ccc+v_ccm)
            dvdtd_xy_s = visc*(dvdxp-dvdxm)*dxi + visc*(dvdyp-dvdym)*dyi
            dvdtd_z_s  = visc*(dvdzp-dvdzm)*dzfi(k)
            dvdt_s     = -(uvip-uvim)*dxi - (vvjp-vvjm)*dyi - (wvkp-wvkm)*dzfi(k)
            dwdxp = (w_pcc-w_ccc)*dxi
            dwdxm = (w_ccc-w_mcc)*dxi
            dwdyp = (w_cpc-w_ccc)*dyi
            dwdym = (w_ccc-w_cmc)*dyi
            dwdzp = (w_ccp-w_ccc)*dzfi(k+1)
            dwdzm = (w_ccc-w_ccm)*dzfi(k)
            uwip  = 0.25*(u_ccc+u_ccp)*(w_ccc+w_pcc)
            uwim  = 0.25*(u_mcc+u_mcp)*(w_ccc+w_mcc)
            vwjp  = 0.25*(v_ccc+v_ccp)*(w_ccc+w_cpc)
            vwjm  = 0.25*(v_cmc+v_cmp)*(w_ccc+w_cmc)
            wwkp  = 0.25*(w_ccc+w_ccp)*(w_ccc+w_ccp)
            wwkm  = 0.25*(w_ccc+w_ccm)*(w_ccc+w_ccm)
            dwdtd_xy_s = visc*(dwdxp-dwdxm)*dxi + visc*(dwdyp-dwdym)*dyi
            dwdtd_z_s  = visc*(dwdzp-dwdzm)*dzci(k)
            dwdt_s     = -(uwip-uwim)*dxi - (vwjp-vwjm)*dyi - (wwkp-wwkm)*dzfi(k)
            dudt_s = dudt_s + dudtd_xy_s
            dvdt_s = dvdt_s + dvdtd_xy_s
            dwdt_s = dwdt_s + dwdtd_xy_s
            dudtd_s = dudtd_z_s
            dvdtd_s = dvdtd_z_s
            dwdtd_s = dwdtd_z_s
            dudt(i,j,k) = dudt_s
            dvdt(i,j,k) = dvdt_s
            dwdt(i,j,k) = dwdt_s
            dudtd(i,j,k) = dudtd_s
            dvdtd(i,j,k) = dvdtd_s
            dwdtd(i,j,k) = dwdtd_s
          end do
        end do
      end do
    else if(is_impdiff) then
      !$acc parallel     loop collapse(3) default(present) async(1) &
      !$acc private(u_ccm,u_pcm,u_cpm,u_cmc,u_pmc,u_mcc,u_ccc,u_pcc,u_mpc,u_cpc,u_cmp,u_mcp,u_ccp) &
      !$acc private(v_ccm,v_pcm,v_cpm,v_cmc,v_pmc,v_mcc,v_ccc,v_pcc,v_mpc,v_cpc,v_cmp,v_mcp,v_ccp) &
      !$acc private(w_ccm,w_pcm,w_cpm,w_cmc,w_pmc,w_mcc,w_ccc,w_pcc,w_mpc,w_cpc,w_cmp,w_mcp,w_ccp) &
      !$acc private(uuip,uuim,vujp,vujm,wukp,wukm) &
      !$acc private(uvip,uvim,vvjp,vvjm,wvkp,wvkm) &
      !$acc private(uwip,uwim,vwjp,vwjm,wwkp,wwkm) &
      !$acc private(dudxp,dudxm,dudyp,dudym,dudzp,dudzm) &
      !$acc private(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm) &
      !$acc private(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm) &
      !$acc private(dudt_s,dvdt_s,dwdt_s) &
      !$acc private(dudtd_s,dvdtd_s,dwdtd_s) &
      !$acc private(dudtd_xy_s,dudtd_z_s) &
      !$acc private(dvdtd_xy_s,dvdtd_z_s) &
      !$acc private(dwdtd_xy_s,dwdtd_z_s)
      !$omp target teams loop collapse(3) &
      !$omp private(u_ccm,u_pcm,u_cpm,u_cmc,u_pmc,u_mcc,u_ccc,u_pcc,u_mpc,u_cpc,u_cmp,u_mcp,u_ccp) &
      !$omp private(v_ccm,v_pcm,v_cpm,v_cmc,v_pmc,v_mcc,v_ccc,v_pcc,v_mpc,v_cpc,v_cmp,v_mcp,v_ccp) &
      !$omp private(w_ccm,w_pcm,w_cpm,w_cmc,w_pmc,w_mcc,w_ccc,w_pcc,w_mpc,w_cpc,w_cmp,w_mcp,w_ccp) &
      !$omp private(uuip,uuim,vujp,vujm,wukp,wukm) &
      !$omp private(uvip,uvim,vvjp,vvjm,wvkp,wvkm) &
      !$omp private(uwip,uwim,vwjp,vwjm,wwkp,wwkm) &
      !$omp private(dudxp,dudxm,dudyp,dudym,dudzp,dudzm) &
      !$omp private(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm) &
      !$omp private(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm) &
      !$omp private(dudt_s,dvdt_s,dwdt_s) &
      !$omp private(dudtd_s,dvdtd_s,dwdtd_s) &
      !$omp private(dudtd_xy_s,dudtd_z_s) &
      !$omp private(dvdtd_xy_s,dvdtd_z_s) &
      !$omp private(dwdtd_xy_s,dwdtd_z_s)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            u_ccm = u(i  ,j  ,k-1)
            u_pcm = u(i+1,j  ,k-1)
            u_cpm = u(i  ,j+1,k-1)
            u_cmc = u(i  ,j-1,k  )
            u_pmc = u(i+1,j-1,k  )
            u_mcc = u(i-1,j  ,k  )
            u_ccc = u(i  ,j  ,k  )
            u_pcc = u(i+1,j  ,k  )
            u_mpc = u(i-1,j+1,k  )
            u_cpc = u(i  ,j+1,k  )
            u_cmp = u(i  ,j-1,k+1)
            u_mcp = u(i-1,j  ,k+1)
            u_ccp = u(i  ,j  ,k+1)
            v_ccm = v(i  ,j  ,k-1)
            v_pcm = v(i+1,j  ,k-1)
            v_cpm = v(i  ,j+1,k-1)
            v_cmc = v(i  ,j-1,k  )
            v_pmc = v(i+1,j-1,k  )
            v_mcc = v(i-1,j  ,k  )
            v_ccc = v(i  ,j  ,k  )
            v_pcc = v(i+1,j  ,k  )
            v_mpc = v(i-1,j+1,k  )
            v_cpc = v(i  ,j+1,k  )
            v_cmp = v(i  ,j-1,k+1)
            v_mcp = v(i-1,j  ,k+1)
            v_ccp = v(i  ,j  ,k+1)
            w_ccm = w(i  ,j  ,k-1)
            w_pcm = w(i+1,j  ,k-1)
            w_cpm = w(i  ,j+1,k-1)
            w_cmc = w(i  ,j-1,k  )
            w_pmc = w(i+1,j-1,k  )
            w_mcc = w(i-1,j  ,k  )
            w_ccc = w(i  ,j  ,k  )
            w_pcc = w(i+1,j  ,k  )
            w_mpc = w(i-1,j+1,k  )
            w_cpc = w(i  ,j+1,k  )
            w_cmp = w(i  ,j-1,k+1)
            w_mcp = w(i-1,j  ,k+1)
            w_ccp = w(i  ,j  ,k+1)
            dudxp = (u_pcc-u_ccc)*dxi
            dudxm = (u_ccc-u_mcc)*dxi
            dudyp = (u_cpc-u_ccc)*dyi
            dudym = (u_ccc-u_cmc)*dyi
            dudzp = (u_ccp-u_ccc)*dzci(k)
            dudzm = (u_ccc-u_ccm)*dzci(k-1)
            uuip  = 0.25*(u_pcc+u_ccc)*(u_ccc+u_pcc)
            uuim  = 0.25*(u_mcc+u_ccc)*(u_ccc+u_mcc)
            vujp  = 0.25*(v_pcc+v_ccc)*(u_ccc+u_cpc)
            vujm  = 0.25*(v_pmc+v_cmc)*(u_ccc+u_cmc)
            wukp  = 0.25*(w_pcc+w_ccc)*(u_ccc+u_ccp)
            wukm  = 0.25*(w_pcm+w_ccm)*(u_ccc+u_ccm)
            dudtd_xy_s = visc*(dudxp-dudxm)*dxi + visc*(dudyp-dudym)*dyi
            dudtd_z_s  = visc*(dudzp-dudzm)*dzfi(k)
            dudt_s     = -(uuip-uuim)*dxi - (vujp-vujm)*dyi - (wukp-wukm)*dzfi(k)
            dvdxp = (v_pcc-v_ccc)*dxi
            dvdxm = (v_ccc-v_mcc)*dxi
            dvdyp = (v_cpc-v_ccc)*dyi
            dvdym = (v_ccc-v_cmc)*dyi
            dvdzp = (v_ccp-v_ccc)*dzci(k)
            dvdzm = (v_ccc-v_ccm)*dzci(k-1)
            uvip  = 0.25*(u_ccc+u_cpc)*(v_ccc+v_pcc)
            uvim  = 0.25*(u_mcc+u_mpc)*(v_ccc+v_mcc)
            vvjp  = 0.25*(v_ccc+v_cpc)*(v_ccc+v_cpc)
            vvjm  = 0.25*(v_ccc+v_cmc)*(v_ccc+v_cmc)
            wvkp  = 0.25*(w_ccc+w_cpc)*(v_ccc+v_ccp)
            wvkm  = 0.25*(w_ccm+w_cpm)*(v_ccc+v_ccm)
            dvdtd_xy_s = visc*(dvdxp-dvdxm)*dxi + visc*(dvdyp-dvdym)*dyi
            dvdtd_z_s  = visc*(dvdzp-dvdzm)*dzfi(k)
            dvdt_s     = -(uvip-uvim)*dxi - (vvjp-vvjm)*dyi - (wvkp-wvkm)*dzfi(k)
            dwdxp = (w_pcc-w_ccc)*dxi
            dwdxm = (w_ccc-w_mcc)*dxi
            dwdyp = (w_cpc-w_ccc)*dyi
            dwdym = (w_ccc-w_cmc)*dyi
            dwdzp = (w_ccp-w_ccc)*dzfi(k+1)
            dwdzm = (w_ccc-w_ccm)*dzfi(k)
            uwip  = 0.25*(u_ccc+u_ccp)*(w_ccc+w_pcc)
            uwim  = 0.25*(u_mcc+u_mcp)*(w_ccc+w_mcc)
            vwjp  = 0.25*(v_ccc+v_ccp)*(w_ccc+w_cpc)
            vwjm  = 0.25*(v_cmc+v_cmp)*(w_ccc+w_cmc)
            wwkp  = 0.25*(w_ccc+w_ccp)*(w_ccc+w_ccp)
            wwkm  = 0.25*(w_ccc+w_ccm)*(w_ccc+w_ccm)
            dwdtd_xy_s = visc*(dwdxp-dwdxm)*dxi + visc*(dwdyp-dwdym)*dyi
            dwdtd_z_s  = visc*(dwdzp-dwdzm)*dzci(k)
            dwdt_s     = -(uwip-uwim)*dxi - (vwjp-vwjm)*dyi - (wwkp-wwkm)*dzfi(k)
            dudtd_s = dudtd_xy_s + dudtd_z_s
            dvdtd_s = dvdtd_xy_s + dvdtd_z_s
            dwdtd_s = dwdtd_xy_s + dwdtd_z_s
            dudt(i,j,k) = dudt_s
            dvdt(i,j,k) = dvdt_s
            dwdt(i,j,k) = dwdt_s
            dudtd(i,j,k) = dudtd_s
            dvdtd(i,j,k) = dvdtd_s
            dwdtd(i,j,k) = dwdtd_s
          end do
        end do
      end do
    end if
#endif
  end subroutine mom_xyz_ad
end module mod_mom
