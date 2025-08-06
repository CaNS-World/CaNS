! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2025 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_scal
  use, intrinsic :: iso_c_binding, only: C_PTR
  use mod_types
  use mod_solve_helmholtz, only: rhs_bound
  implicit none
  private
  public scal,scalar,cmpt_scalflux,bulk_forcing_s,initialize_scalars
  !
  ! scalar derived type
  !
  type scalar
    real(rp), allocatable, dimension(:,:,:) :: val
    real(rp), allocatable, dimension(:,:,:) :: dsdtrko
    real(rp) :: alpha
    character(len=100) :: ini
    character(len=1), dimension(0:1,3) :: cbc
    real(rp),         dimension(0:1,3) :: bc
    real(rp) :: source
    logical  :: is_forced
    real(rp) :: scalf
    real(rp) :: f
    real(rp), dimension(0:1,3) :: fluxo
#if !defined(_OPENACC) || defined(_USE_HIP)
    type(C_PTR), dimension(2,2) :: arrplan
#else
    integer    , dimension(2,2) :: arrplan
#endif
    real(rp), allocatable, dimension(:,:) :: lambdaxy
    real(rp), allocatable, dimension(:) :: a,b,c
    real(rp) :: normfft
    type(rhs_bound) :: rhsb
  end type scalar
  !
  contains
  subroutine scal(nx,ny,nz,dxi,dyi,dzci,dzfi,visc,u,v,w,s,dsdt,dsdtd)
    use mod_param, only: is_impdiff,is_impdiff_1d
    !
    ! computes convective and diffusive fluxes
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,visc
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in) :: u,v,w,s
    real(rp), dimension(:,:,:), intent(out) :: dsdt
    real(rp), dimension(:,:,:), intent(out), optional :: dsdtd
    integer :: i,j,k
    real(rp) :: usip,usim,vsjp,vsjm,wskp,wskm
    real(rp) :: dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm
    real(rp) :: dsdtd_xy,dsdtd_z
    !
#if !defined(_LOOP_UNSWITCHING)
    !$acc parallel loop collapse(3) default(present) &
    !$acc private(usip,usim,vsjp,vsjm,wskp,wskm,dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm,dsdtd_xy,dsdtd_z) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared) &
    !$OMP PRIVATE(usip,usim,vsjp,vsjm,wskp,wskm,dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm,dsdtd_xy,dsdtd_z)
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
          dsdt(i,j,k) = dxi*(     -usip + usim ) + &
                        dyi*(     -vsjp + vsjm ) + &
                        dzfi(k)*( -wskp + wskm )
          dsdtd_xy = (dsdxp-dsdxm)*visc*dxi + &
                     (dsdyp-dsdym)*visc*dyi
          dsdtd_z  = (dsdzp-dsdzm)*visc*dzfi(k)
          if(is_impdiff) then
            if(is_impdiff_1d) then
              dsdt(i,j,k)  = dsdt(i,j,k) + dsdtd_xy
              dsdtd(i,j,k) = dsdtd_z
            else
              dsdtd(i,j,k) = dsdtd_xy + dsdtd_z
            end if
          else
            dsdt(i,j,k)  = dsdt(i,j,k) + dsdtd_xy + dsdtd_z
          end if
        end do
      end do
    end do
#else
    if(.not.is_impdiff) then
      !$acc parallel loop collapse(3) default(present) &
      !$acc private(usip,usim,vsjp,vsjm,wskp,wskm,dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm,dsdtd_xy,dsdtd_z) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared) &
      !$OMP PRIVATE(usip,usim,vsjp,vsjm,wskp,wskm,dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm,dsdtd_xy,dsdtd_z)
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
            dsdt(i,j,k) = dxi*(    -usip + usim ) + &
                          dyi*(    -vsjp + vsjm ) + &
                          dzfi(k)*( -wskp + wskm )
            dsdtd_xy = (dsdxp-dsdxm)*visc*dxi + &
                       (dsdyp-dsdym)*visc*dyi
            dsdtd_z  = (dsdzp-dsdzm)*visc*dzfi(k)
            dsdt(i,j,k) = dsdt(i,j,k) + dsdtd_xy + dsdtd_z
          end do
        end do
      end do
    else if(is_impdiff .and. is_impdiff_1d) then
      !$acc parallel loop collapse(3) default(present) &
      !$acc private(usip,usim,vsjp,vsjm,wskp,wskm,dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm,dsdtd_xy,dsdtd_z) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared) &
      !$OMP PRIVATE(usip,usim,vsjp,vsjm,wskp,wskm,dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm,dsdtd_xy,dsdtd_z)
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
            dsdt(i,j,k) = dxi*(    -usip + usim ) + &
                          dyi*(    -vsjp + vsjm ) + &
                          dzfi(k)*( -wskp + wskm )
            dsdtd_xy = (dsdxp-dsdxm)*visc*dxi + &
                       (dsdyp-dsdym)*visc*dyi
            dsdtd_z  = (dsdzp-dsdzm)*visc*dzfi(k)
            dsdt(i,j,k)  = dsdt(i,j,k) + dsdtd_xy
            dsdtd(i,j,k) = dsdtd_z
          end do
        end do
      end do
    else
      !$acc parallel loop collapse(3) default(present) &
      !$acc private(usip,usim,vsjp,vsjm,wskp,wskm,dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm,dsdtd_xy,dsdtd_z) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared) &
      !$OMP PRIVATE(usip,usim,vsjp,vsjm,wskp,wskm,dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm,dsdtd_xy,dsdtd_z)
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
            dsdt(i,j,k) = dxi*(    -usip + usim ) + &
                          dyi*(    -vsjp + vsjm ) + &
                          dzfi(k)*( -wskp + wskm )
            dsdtd_xy = (dsdxp-dsdxm)*visc*dxi + &
                       (dsdyp-dsdym)*visc*dyi
            dsdtd_z  = (dsdzp-dsdzm)*visc*dzfi(k)
            dsdtd(i,j,k) = dsdtd_xy + dsdtd_z
          end do
        end do
      end do
    end if
#endif
  end subroutine scal
  !
  subroutine cmpt_scalflux(n,is_bound,l,dli,dzci,dzfi,alpha,s,flux)
    use mpi
    use mod_param, only: cbcpre ! it is fine to use the pressure BC to check for periodicity
    implicit none
    integer , intent(in ), dimension(3) :: n
    logical , intent(in ), dimension(0:1,3) :: is_bound
    real(rp), intent(in ), dimension(3)     :: l,dli
    real(rp), intent(in ), dimension(0:)    :: dzci,dzfi
    real(rp), intent(in )                   :: alpha
    real(rp), intent(in ), dimension(0:,0:,0:) :: s
    real(rp), intent(out), dimension(0:1,1:3) :: flux
    real(rp) :: dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm
    real(rp) :: flux_xp,flux_xm,flux_yp,flux_ym,flux_zp,flux_zm
    integer :: ierr
    !
    integer :: i,j,k,nx,ny,nz
    real(rp) :: dxi,dyi,lx,ly,lz
    !
    nx = n(1); ny = n(2); nz = n(3)
    dxi = dli(1); dyi = dli(2)
    lx = l(1); ly = l(2); lz = l(3)
    flux_xp = 0.
    flux_xm = 0.
    !$acc data copy(flux_xp,flux_xm) async(1)
    if(is_bound(0,1).and.cbcpre(0,1)//cbcpre(1,1) /= 'PP') then
      !$acc parallel loop collapse(2) default(present) private(dsdxp) reduction(+:flux_xp) async(1)
      !$OMP PARALLEL DO   COLLAPSE(2) DEFAULT(shared ) private(dsdxp) reduction(+:flux_xp)
      do k=1,nz
        do j=1,ny
          dsdxp   = (s(1 ,j,k)-s(0   ,j,k))*dxi*alpha
          flux_xp = flux_xp + dsdxp/(dyi*dzfi(k)*ly*lz)
        end do
      end do
    end if
    if(is_bound(1,1).and.cbcpre(0,1)//cbcpre(1,1) /= 'PP') then
      !$acc parallel loop collapse(2) default(present) private(dsdxm) reduction(+:flux_xm) async(1)
      !$OMP PARALLEL DO   collapse(2) default(shared ) private(dsdxm) reduction(+:flux_xm)
      do k=1,nz
        do j=1,ny
          dsdxm   = (s(nx,j,k)-s(nx+1,j,k))*dxi*alpha
          flux_xm = flux_xm + dsdxm/(dyi*dzfi(k)*ly*lz)
        end do
      end do
    end if
    !$acc end data
    flux_yp = 0.
    flux_ym = 0.
    !$acc data copy(flux_yp,flux_ym) async(1)
    if(is_bound(0,2).and.cbcpre(0,2)//cbcpre(1,2) /= 'PP') then
      !$acc parallel loop collapse(2) default(present) private(dsdyp) reduction(+:flux_yp) async(1)
      !$OMP PARALLEL DO   collapse(2) default(shared ) private(dsdyp) reduction(+:flux_yp)
      do k=1,nz
        do i=1,nx
          dsdyp   = (s(i,1 ,k)-s(i,0   ,k))*dyi*alpha
          flux_yp = flux_yp + dsdyp/(dxi*dzfi(k)*lx*lz)
        end do
      end do
    end if
    if(is_bound(1,2).and.cbcpre(0,2)//cbcpre(1,2) /= 'PP') then
      !$acc parallel loop collapse(2) default(present) private(dsdym) reduction(+:flux_ym) async(1)
      !$OMP PARALLEL DO   collapse(2) default(shared ) private(dsdym) reduction(+:flux_ym)
      do k=1,nz
        do i=1,nx
          dsdym   = (s(i,ny,k)-s(i,ny+1,k))*dyi*alpha
          flux_ym = flux_ym + dsdym/(dxi*dzfi(k)*lx*lz)
        end do
      end do
    end if
    !$acc end data
    flux_zp = 0.
    flux_zm = 0.
    !$acc data copy(flux_zp,flux_zm) async(1)
    if(is_bound(0,3).and.cbcpre(0,3)//cbcpre(1,3) /= 'PP') then
      !$acc parallel loop collapse(2) default(present) private(dsdzp) reduction(+:flux_zp) async(1)
      !$OMP PARALLEL DO   collapse(2) default(shared ) private(dsdzp) reduction(+:flux_zp)
      do j=1,ny
        do i=1,nx
          dsdzp   = (s(i,j,1 )-s(i,j,0   ))*dzci(0)*alpha
          flux_zp = flux_zp + dsdzp/(dxi*dyi*lx*ly)
        end do
      end do
    end if
    if(is_bound(1,3).and.cbcpre(0,3)//cbcpre(1,3) /= 'PP') then
      !$acc parallel loop collapse(2) default(present) private(dsdzm) reduction(+:flux_zm) async(1)
      !$OMP PARALLEL DO   collapse(2) default(shared ) private(dsdzm) reduction(+:flux_zm)
      do j=1,ny
        do i=1,nx
          dsdzm   = (s(i,j,nz)-s(i,j,nz+1))*dzci(nz)*alpha
          flux_zm = flux_zm + dsdzm/(dxi*dyi*lx*ly)
        end do
      end do
    end if
    !$acc end data
    !$acc wait(1)
    !
    !flux(:) = [flux_x,flux_y,flux_z]
    !
    ! discern upper and lower boundary fluxes
    !
    flux(:,1) = [flux_xp,flux_xm]
    flux(:,2) = [flux_yp,flux_ym]
    flux(:,3) = [flux_zp,flux_zm]
    !
    call MPI_ALLREDUCE(MPI_IN_PLACE,flux,6,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  end subroutine cmpt_scalflux
  !
  subroutine bulk_forcing_s(n,is_forced,ff,p)
    implicit none
    integer , intent(in   ), dimension(3) :: n
    logical , intent(in   )               :: is_forced
    real(rp), intent(in   )               :: ff
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    integer :: i,j,k
    if(is_forced) then
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP parallel do   collapse(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            p(i,j,k) = p(i,j,k) + ff
          end do
        end do
      end do
    end if
  end subroutine bulk_forcing_s
  !
  subroutine initialize_scalars(scalars,nscal,n,n_z)
    use mod_param, only: alphai,iniscal,cbcscal,bcscal,ssource,is_sforced,scalf, &
                         is_impdiff
    !
    ! initializes/allocates members of an array of `nscal` scalar derived types
    !
    implicit none
    type(scalar), intent(inout), dimension(:) :: scalars
    integer     , intent(in   ) :: nscal,n(3),n_z(3)
    integer :: iscal
    do iscal=1,nscal
      allocate(scalars(iscal)%val(0:n(1)+1,0:n(2)+1,0:n(3)+1))
      allocate(scalars(iscal)%dsdtrko(n(1),n(2),n(3)))
      if(is_impdiff) then
        allocate(scalars(iscal)%lambdaxy(n_z(1),n_z(2)))
        allocate(scalars(iscal)%a(n_z(3)), &
                 scalars(iscal)%b(n_z(3)), &
                 scalars(iscal)%c(n_z(3)))
        allocate(scalars(iscal)%rhsb%x(n(2),n(3),0:1), &
                 scalars(iscal)%rhsb%y(n(1),n(3),0:1), &
                 scalars(iscal)%rhsb%z(n(1),n(2),0:1))
      end if
      scalars(iscal)%alpha     = alphai(iscal)**(-1)
      scalars(iscal)%ini       = iniscal(iscal)
      scalars(iscal)%cbc       = cbcscal(:,:,iscal)
      scalars(iscal)%bc        = bcscal(:,:,iscal)
      scalars(iscal)%source    = ssource(iscal)
      scalars(iscal)%is_forced = is_sforced(iscal)
      scalars(iscal)%scalf     = scalf(iscal)
    end do
  end subroutine initialize_scalars
end module mod_scal
