! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_scal
  use mod_types
  implicit none
  private
  public scal,cmpt_scalflux
  contains
  subroutine scal(nx,ny,nz,dxi,dyi,dzi,dzci,dzfi,visc,u,v,w,s,dsdt)
    !
    !
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,dzi,visc
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in) :: u,v,w,s
    real(rp), dimension(:,:,:), intent(out) :: dsdt
    integer :: i,j,k
    real(rp) :: usip,usim,vsjp,vsjm,wskp,wskm
    real(rp) :: dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm
    !
    !$acc parallel loop collapse(3) default(present) private(usip,usim,vsjp,vsjm,wskp,wskm,dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm) async(1)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(usip,usim,vsjp,vsjm,wskp,wskm) &
    !$OMP PRIVATE(dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzi,visc,u,v,w,s,dsdt,dzci,dzfi)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          usim  = 0.5*( s(i-1,j,k)+s(i,j,k) )*u(i-1,j,k)
          usip  = 0.5*( s(i+1,j,k)+s(i,j,k) )*u(i  ,j,k)
          vsjm  = 0.5*( s(i,j-1,k)+s(i,j,k) )*v(i,j-1,k)
          vsjp  = 0.5*( s(i,j+1,k)+s(i,j,k) )*v(i,j  ,k)
          wskm  = 0.5*( s(i,j,k-1)+s(i,j,k) )*w(i,j,k-1)
          wskp  = 0.5*( s(i,j,k+1)+s(i,j,k) )*w(i,j,k  )
          dsdxp = (s(i+1,j,k)-s(i  ,j,k))*dxi
          dsdxm = (s(i  ,j,k)-s(i-1,j,k))*dxi
          dsdyp = (s(i,j+1,k)-s(i,j  ,k))*dyi
          dsdym = (s(i,j  ,k)-s(i,j-1,k))*dyi
          dsdzp = (s(i,j,k+1)-s(i,j,k  ))*dzci(k  )
          dsdzm = (s(i,j,k  )-s(i,j,k-1))*dzci(k-1)
          !
          dsdt(i,j,k) = dxi*(     -usip + usim ) + (dsdxp-dsdxm)*visc*dxi + &
                        dyi*(     -vsjp + vsjm ) + (dsdyp-dsdym)*visc*dyi + &
                        dzfi(k)*( -wskp + wskm ) + (dsdzp-dsdzm)*visc*dzfi(k)
        end do
      end do
    end do
  end subroutine scal
  !
  subroutine cmpt_scalflux(n,is_bound,l,dli,dzci,dzfi,alpha,s,flux)
  use mpi
  use mod_param, only: cbcpre ! needs to be replaced with scal bc
  implicit none
  integer , intent(in ), dimension(3) :: n
  logical , intent(in ), dimension(0:1,3) :: is_bound
  real(rp), intent(in ), dimension(3)     :: l,dli
  real(rp), intent(in ), dimension(0:)    :: dzci,dzfi
  real(rp), intent(in )                   :: alpha
  real(rp), intent(in ), dimension(0:,0:,0:) :: s
  real(rp), intent(out), dimension(3) :: flux
  real(rp) :: dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm
  real(rp) :: flux_x,flux_y,flux_z
  integer :: ierr
  !
  integer :: i,j,k,nx,ny,nz
  real(rp) :: dxi,dyi,lx,ly,lz
  !
  nx = n(1); ny = n(2); nz = n(3)
  dxi = dli(1); dyi = dli(2)
  lx = l(1); ly = l(2); lz = l(3)
  flux_x = 0._rp
  !$acc data copy(flux_x) async(1)
  if(is_bound(0,1).and.cbcpre(0,1)//cbcpre(1,1) /= 'PP') then
    !$acc parallel loop collapse(2) default(present) private(dsdxp) reduction(+:flux_x) async(1)
    !$OMP PARALLEL DO DEFAULT(shared) private(dsdxp) reduction(+:flux_x)
    do k=1,nz
      do j=1,ny
        dsdxp = (s(1 ,j,k)-s(0   ,j,k))*dxi*alpha
        flux_x = flux_x + dsdxp/(dyi*dzfi(k)*ly*lz)
      end do
    end do
  end if
  if(is_bound(1,1).and.cbcpre(0,1)//cbcpre(1,1) /= 'PP') then
    !$acc parallel loop collapse(2) default(present) private(dsdxm) reduction(+:flux_x) async(1)
    !$OMP PARALLEL DO DEFAULT(shared) private(dsdxm) reduction(+:flux_x)
    do k=1,nz
      do j=1,ny
        dsdxm = (s(nx,j,k)-s(nx+1,j,k))*dxi*alpha
        flux_x = flux_x + dsdxm/(dyi*dzfi(k)*ly*lz)
      end do
    end do
  end if
  !$acc end data
  flux_y = 0._rp
  !$acc data copy(flux_y) async(1)
  if(is_bound(0,2).and.cbcpre(0,2)//cbcpre(1,2) /= 'PP') then
    !$acc parallel loop collapse(2) default(present) private(dsdyp) reduction(+:flux_y) async(1)
    !$OMP PARALLEL DO DEFAULT(shared) private(dsdyp) reduction(+:flux_y)
    do k=1,nz
      do i=1,nx
        dsdyp = (s(i,1 ,k)-s(i,0   ,k))*dyi*alpha
        flux_y = flux_y + dsdyp/(dxi*dzfi(k)*lx*lz)
      end do
    end do
  end if
  if(is_bound(1,2).and.cbcpre(0,2)//cbcpre(1,2) /= 'PP') then
    !$acc parallel loop collapse(2) default(present) private(dsdym) reduction(+:flux_y) async(1)
    !$OMP PARALLEL DO DEFAULT(shared) private(dsdym) reduction(+:flux_y)
    do k=1,nz
      do i=1,nx
        dsdym = (s(i,ny,k)-s(i,ny+1,k))*dyi*alpha
        flux_y = flux_y + dsdym/(dxi*dzfi(k)*lx*lz)
      end do
    end do
  end if
  !$acc end data
  flux_z = 0._rp
  !$acc data copy(flux_z) async(1)
  if(is_bound(0,3).and.cbcpre(0,3)//cbcpre(1,3) /= 'PP') then
    !$acc parallel loop collapse(2) default(present) private(dsdzp) reduction(+:flux_z) async(1)
    !$OMP PARALLEL DO DEFAULT(shared) private(dsdzp) reduction(+:flux_z)
    do j=1,ny
      do i=1,nx
        dsdzp = (s(i,j,1 )-s(i,j,0   ))*dzci(0)*alpha
        flux_z = flux_z + dsdzp/(dxi*dyi*lx*ly)
      end do
    end do
  end if
  if(is_bound(1,3).and.cbcpre(0,3)//cbcpre(1,3) /= 'PP') then
    !$acc parallel loop collapse(2) default(present) private(dsdzm) reduction(+:flux_z) async(1)
    !$OMP PARALLEL DO DEFAULT(shared) private(dsdzm) reduction(+:flux_z)
    do j=1,ny
      do i=1,nx
        dsdzm = (s(i,j,nz)-s(i,j,nz+1))*dzci(nz)*alpha
        flux_z = flux_z + dsdzm/(dxi*dyi*lx*ly)
      end do
    end do
  end if
  !$acc end data
  !$acc wait(1)
  flux(:) = [flux_x,flux_y,flux_z]
  call MPI_ALLREDUCE(MPI_IN_PLACE,flux(3),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
end subroutine cmpt_scalflux
!
end module mod_scal
