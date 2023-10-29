! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
#define _FAST_MOM_KERNELS
module mod_rk
  use mod_mom  , only: momx_a,momy_a,momz_a, &
                       momx_d,momy_d,momz_d, &
                       momx_p,momy_p,momz_p, cmpt_wallshear
#if defined(_IMPDIFF_1D)
  use mod_mom  , only: momx_d_xy,momy_d_xy,momz_d_xy, &
                       momx_d_z ,momy_d_z ,momz_d_z
#endif
#if defined(_FAST_MOM_KERNELS)
  use mod_mom  , only: mom_xyz_ad
#endif
  use mod_scal , only: scal,cmpt_scalflux
  use mod_utils, only: bulk_mean,swap
  use mod_types
  implicit none
  private
  public rk
  contains
  subroutine rk(rkpar,n,dli,dzci,dzfi,grid_vol_ratio_c,grid_vol_ratio_f,visc,dt,p, &
                is_forced,velf,bforce,u,v,w,f)
    !
    ! low-storage 3rd-order Runge-Kutta scheme
    ! for time integration of the momentum equations.
    !
    implicit none
    real(rp), intent(in), dimension(2) :: rkpar
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in) :: visc,dt
    real(rp), intent(in   ), dimension(3) :: dli
    real(rp), intent(in   ), dimension(0:) :: dzci,dzfi
    real(rp), intent(in   ), dimension(0:) :: grid_vol_ratio_c,grid_vol_ratio_f
    real(rp), intent(in   ), dimension(0:,0:,0:) :: p
    logical , intent(in   ), dimension(3)        :: is_forced
    real(rp), intent(in   ), dimension(3)        :: velf,bforce
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out), dimension(3) :: f
    real(rp), target     , allocatable, dimension(:,:,:), save :: dudtrk_t ,dvdtrk_t ,dwdtrk_t , &
                                                                  dudtrko_t,dvdtrko_t,dwdtrko_t
    real(rp), pointer    , contiguous , dimension(:,:,:), save :: dudtrk   ,dvdtrk   ,dwdtrk   , &
                                                                  dudtrko  ,dvdtrko  ,dwdtrko
    real(rp),              allocatable, dimension(:,:,:), save :: dudtrkd  ,dvdtrkd  ,dwdtrkd
    logical, save :: is_first = .true.
    real(rp) :: factor1,factor2,factor12
    integer :: i,j,k
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2
    !
    ! initialization
    !
    if(is_first) then ! leverage save attribute to allocate these arrays on the device only once
      is_first = .false.
      allocate(dudtrk_t( n(1),n(2),n(3)),dvdtrk_t( n(1),n(2),n(3)),dwdtrk_t( n(1),n(2),n(3)))
      allocate(dudtrko_t(n(1),n(2),n(3)),dvdtrko_t(n(1),n(2),n(3)),dwdtrko_t(n(1),n(2),n(3)))
      !$acc enter data create(dudtrk_t ,dvdtrk_t ,dwdtrk_t ) async(1)
      !$acc enter data create(dudtrko_t,dvdtrko_t,dwdtrko_t) async(1)
      !$acc kernels default(present) async(1) ! not really necessary
      dudtrko_t(:,:,:) = 0._rp
      dvdtrko_t(:,:,:) = 0._rp
      dwdtrko_t(:,:,:) = 0._rp
      !$acc end kernels
#if defined(_IMPDIFF)
      allocate(dudtrkd(n(1),n(2),n(3)),dvdtrkd(n(1),n(2),n(3)),dwdtrkd(n(1),n(2),n(3)))
      !$acc enter data create(dudtrkd,dvdtrkd,dwdtrkd) async(1)
#endif
      dudtrk  => dudtrk_t
      dvdtrk  => dvdtrk_t
      dwdtrk  => dwdtrk_t
      dudtrko => dudtrko_t
      dvdtrko => dvdtrko_t
      dwdtrko => dwdtrko_t
    end if
    !
#if defined(_FAST_MOM_KERNELS)
    call mom_xyz_ad(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,u,v,w,dudtrk,dvdtrk,dwdtrk,dudtrkd,dvdtrkd,dwdtrkd)
#else
    !$acc kernels default(present) async(1)
    !$OMP PARALLEL WORKSHARE
    dudtrk(:,:,:) = 0._rp
    dvdtrk(:,:,:) = 0._rp
    dwdtrk(:,:,:) = 0._rp
    !$OMP END PARALLEL WORKSHARE
    !$acc end kernels
#if !defined(_IMPDIFF)
    call momx_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,u,dudtrk)
    call momy_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,v,dvdtrk)
    call momz_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,w,dwdtrk)
#else
    !$acc kernels default(present) async(1)
    !$OMP PARALLEL WORKSHARE
    dudtrkd(:,:,:) = 0._rp
    dvdtrkd(:,:,:) = 0._rp
    dwdtrkd(:,:,:) = 0._rp
    !$OMP END PARALLEL WORKSHARE
    !$acc end kernels
#if !defined(_IMPDIFF_1D)
    call momx_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,u,dudtrkd)
    call momy_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,v,dvdtrkd)
    call momz_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,w,dwdtrkd)
#else
    call momx_d_xy(n(1),n(2),n(3),dli(1),dli(2),visc,u,dudtrk )
    call momy_d_xy(n(1),n(2),n(3),dli(1),dli(2),visc,v,dvdtrk )
    call momz_d_xy(n(1),n(2),n(3),dli(1),dli(2),visc,w,dwdtrk )
    call momx_d_z( n(1),n(2),n(3),dzci  ,dzfi  ,visc,u,dudtrkd)
    call momy_d_z( n(1),n(2),n(3),dzci  ,dzfi  ,visc,v,dvdtrkd)
    call momz_d_z( n(1),n(2),n(3),dzci  ,dzfi  ,visc,w,dwdtrkd)
#endif
    call momx_a(n(1),n(2),n(3),dli(1),dli(2),dzfi,u,v,w,dudtrk)
    call momy_a(n(1),n(2),n(3),dli(1),dli(2),dzfi,u,v,w,dvdtrk)
    call momz_a(n(1),n(2),n(3),dli(1),dli(2),dzci,u,v,w,dwdtrk)
#endif
#endif
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
#if !defined(_FAST_MOM_KERNELS)
          u(i,j,k) = u(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k)
          v(i,j,k) = v(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k)
          w(i,j,k) = w(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k)
#else
          u(i,j,k) = u(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k) + &
                                factor12*(bforce(1) - dli(1)*( p(i+1,j,k)-p(i,j,k)))
          !
          v(i,j,k) = v(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k) + &
                                factor12*(bforce(2) - dli(2)*( p(i,j+1,k)-p(i,j,k)))
          !
          w(i,j,k) = w(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k) + &
                                factor12*(bforce(3) - dzci(k)*(p(i,j,k+1)-p(i,j,k)))
#endif
#if defined(_IMPDIFF)
          u(i,j,k) = u(i,j,k) + factor12*dudtrkd(i,j,k)
          v(i,j,k) = v(i,j,k) + factor12*dvdtrkd(i,j,k)
          w(i,j,k) = w(i,j,k) + factor12*dwdtrkd(i,j,k)
#endif
        end do
      end do
    end do
    !
    ! swap d?dtrk <-> d?dtrko
    !
    call swap(dudtrk,dudtrko)
    call swap(dvdtrk,dvdtrko)
    call swap(dwdtrk,dwdtrko)
!#if 0 /*pressure gradient term treated explicitly later */
!    !$acc kernels
!    !$OMP PARALLEL WORKSHARE
!    dudtrk(:,:,:) = 0._rp
!    dvdtrk(:,:,:) = 0._rp
!    dwdtrk(:,:,:) = 0._rp
!    !$OMP END PARALLEL WORKSHARE
!    !$acc end kernels
!    call momx_p(n(1),n(2),n(3),dli(1),bforce(1),p,dudtrk)
!    call momy_p(n(1),n(2),n(3),dli(2),bforce(2),p,dvdtrk)
!    call momz_p(n(1),n(2),n(3),dzci  ,bforce(3),p,dwdtrk)
!    !$acc parallel loop collapse(3)
!    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
!    do k=1,n(3)
!      do j=1,n(2)
!        do i=1,n(1)
!          u(i,j,k) = u(i,j,k) + factor12*dudtrk(i,j,k)
!          v(i,j,k) = v(i,j,k) + factor12*dvdtrk(i,j,k)
!          w(i,j,k) = w(i,j,k) + factor12*dwdtrk(i,j,k)
!        end do
!      end do
!    end do
!#endif
#if !defined(_FAST_MOM_KERNELS)
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          u(i,j,k) = u(i,j,k) + factor12*(bforce(1) - dli(1)*( p(i+1,j,k)-p(i,j,k)))
          v(i,j,k) = v(i,j,k) + factor12*(bforce(2) - dli(2)*( p(i,j+1,k)-p(i,j,k)))
          w(i,j,k) = w(i,j,k) + factor12*(bforce(3) - dzci(k)*(p(i,j,k+1)-p(i,j,k)))
        end do
      end do
    end do
#endif
    !
    ! compute bulk velocity forcing
    !
    call cmpt_bulk_forcing(n,is_forced,velf,grid_vol_ratio_c,grid_vol_ratio_f,u,v,w,f)
    !
#if defined(_IMPDIFF)
    !
    ! compute rhs of Helmholtz equation
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          u(i,j,k) = u(i,j,k) - .5_rp*factor12*dudtrkd(i,j,k)
          v(i,j,k) = v(i,j,k) - .5_rp*factor12*dvdtrkd(i,j,k)
          w(i,j,k) = w(i,j,k) - .5_rp*factor12*dwdtrkd(i,j,k)
        end do
      end do
    end do
#endif
  end subroutine rk
  !
  subroutine rk_scal(rkpar,n,dli,l,dzci,dzfi,grid_vol_ratio_f,alpha,dt,is_bound,u,v,w, &
                     is_forced,scalf,ssource,fluxo,s,f)
    !
    ! low-storage 3rd-order Runge-Kutta scheme
    ! for time integration of the scalar field.
    !
    ! n.b.: since we leverage the `save` attribute for dsdtrk*, this subroutine only supports
    !       transport of a single scalar; extension to n arbtritrary scalars could be done,
    !       e.g., using a loop through an array of scalar derived types, with an ASSOCIATE
    !       statement to keep the same piece of code to for each scalar
    !
    implicit none
    logical , parameter :: is_cmpt_wallflux = .false.
    real(rp), intent(in   ), dimension(2) :: rkpar
    integer , intent(in   ), dimension(3) :: n
    real(rp), intent(in   ), dimension(3) :: dli,l
    real(rp), intent(in   ), dimension(0:) :: dzci,dzfi
    real(rp), intent(in   ), dimension(0:) :: grid_vol_ratio_f
    real(rp), intent(in   ) :: alpha,dt
    logical , intent(in   ), dimension(0:1,3)    :: is_bound
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u,v,w
    logical , intent(in   ) :: is_forced
    real(rp), intent(in   ) :: scalf,ssource
    real(rp), intent(inout), dimension(3) :: fluxo
    real(rp), intent(inout), dimension(0:,0:,0:) :: s
    real(rp), intent(out  ) :: f
    real(rp), target     , allocatable, dimension(:,:,:), save :: dsdtrk_t,dsdtrko_t
    real(rp), pointer    , contiguous , dimension(:,:,:), save :: dsdtrk  ,dsdtrko
    logical, save :: is_first = .true.
    real(rp) :: factor1,factor2,factor12
    real(rp), dimension(3) :: flux
    integer :: i,j,k
    real(rp) :: mean
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2
    if(is_first) then ! leverage save attribute to allocate these arrays on the device only once
      is_first = .false.
      allocate(dsdtrk_t(n(1),n(2),n(3)),dsdtrko_t(n(1),n(2),n(3)))
      !$acc enter data create(dsdtrk_t,dsdtrko_t) async(1)
      !$acc kernels default(present) async(1) ! not really necessary
      dsdtrko_t(:,:,:) = 0._rp
      !$acc end kernels
      dsdtrk  => dsdtrk_t
      dsdtrko => dsdtrko_t
    end if
    if(is_cmpt_wallflux) then
      call cmpt_scalflux(n,is_bound,l,dli,dzci,dzfi,alpha,s,flux)
      f = (factor1*sum(flux(:)/l(:)) + factor2*sum(fluxo(:)/l(:)))
      fluxo(:) = flux(:)
    end if
    call scal(n(1),n(2),n(3),dli(1),dli(2),dli(3),dzci,dzfi,alpha,u,v,w,s,dsdtrk)
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          s(i,j,k) = s(i,j,k) + factor1*dsdtrk(i,j,k) + factor2*dsdtrko(i,j,k) + factor12*ssource
        end do
      end do
    end do
    !
    ! swap d?dtrk <-> d?dtrko
    !
    call swap(dsdtrk,dsdtrko)
    !
    ! bulk scalar forcing
    !
    if(is_forced) then
      call bulk_mean(n,grid_vol_ratio_f,s,mean)
      f = scalf - mean
    end if
  end subroutine rk_scal
  !
  subroutine cmpt_bulk_forcing(n,is_forced,velf,grid_vol_ratio_c,grid_vol_ratio_f,u,v,w,f)
    implicit none
    integer , intent(in   ), dimension(3) :: n
    logical , intent(in   ), dimension(3) :: is_forced
    real(rp), intent(in   ), dimension(3) :: velf
    real(rp), intent(in   ), dimension(0:) :: grid_vol_ratio_c,grid_vol_ratio_f
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out  ), dimension(3) :: f
    real(rp) :: mean
    !
    ! bulk velocity forcing
    !
    f(:) = 0.
    if(is_forced(1)) then
      call bulk_mean(n,grid_vol_ratio_f,u,mean)
      f(1) = velf(1) - mean
    end if
    if(is_forced(2)) then
      call bulk_mean(n,grid_vol_ratio_f,v,mean)
      f(2) = velf(2) - mean
    end if
    if(is_forced(3)) then
      call bulk_mean(n,grid_vol_ratio_c,w,mean)
      f(3) = velf(3) - mean
    end if
  end subroutine cmpt_bulk_forcing
  !
  subroutine cmpt_bulk_forcing_alternative(rkpar,n,dli,l,dzci,dzfi,visc,dt,is_bound,is_forced,u,v,w,tauxo,tauyo,tauzo,f,is_first)
    !
    ! computes the pressure gradient to be added to the flow that perfectly balances the wall shear stresses
    ! this effectively prescribes zero net acceleration, which allows to sustain a constant mass flux
    !
    implicit none
    real(rp), intent(in), dimension(2) :: rkpar
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in) :: visc,dt
    real(rp), intent(in   ), dimension(3) :: dli,l
    real(rp), intent(in   ), dimension(0:) :: dzci,dzfi
    logical , intent(in   ), dimension(0:1,3)    :: is_bound
    logical , intent(in   ), dimension(3) :: is_forced
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(inout), dimension(3) :: tauxo,tauyo,tauzo
    real(rp), intent(inout), dimension(3) :: f
    real(rp), dimension(3) :: f_aux
    logical , intent(in   ) :: is_first
    real(rp), dimension(3) :: taux,tauy,tauz
    real(rp) :: factor1,factor2,factor12
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = (factor1 + factor2)/2.
    !
    taux(:) = 0._rp
    tauy(:) = 0._rp
    tauz(:) = 0._rp
    call cmpt_wallshear(n,is_forced,is_bound,l,dli,dzci,dzfi,visc,u,v,w,taux,tauy,tauz)
#if !defined(_IMPDIFF)
    if(is_first) then
      f(1) = (factor1*sum(taux(:)/l(:)) + factor2*sum(tauxo(:)/l(:)))
      f(2) = (factor1*sum(tauy(:)/l(:)) + factor2*sum(tauyo(:)/l(:)))
      f(3) = (factor1*sum(tauz(:)/l(:)) + factor2*sum(tauzo(:)/l(:)))
      tauxo(:) = taux(:)
      tauyo(:) = tauy(:)
      tauzo(:) = tauz(:)
    end if
#else
#if defined(_IMPDIFF_1D)
    f_aux(1) = factor12*taux(3)/l(3)
    f_aux(2) = factor12*tauy(3)/l(3)
    if(is_first) then
      f(1) = factor1*taux(2)/l(2) + factor2*tauxo(2)/l(2) + f_aux(1)
      f(2) = factor1*tauy(1)/l(1) + factor2*tauyo(1)/l(1) + f_aux(2)
      f(3) = factor1*sum(tauz(1:2)/l(1:2)) + factor2*sum(tauzo(1:2)/l(1:2))
      tauxo(1:2) = taux(1:2)
      tauyo(1:2) = tauy(1:2)
      tauzo(1:2) = tauz(1:2)
    else
      f(1) = f(1) + f_aux(1)
      f(2) = f(2) + f_aux(2)
    end if
#else
    f_aux(:) = factor12*[sum(taux(:)/l(:)), &
                         sum(tauy(:)/l(:)), &
                         sum(tauz(:)/l(:))]
    if(is_first) then
       f(:) = f_aux(:)
    else
       f(:) = f(:) + f_aux(:)
    endif
#endif
#endif
  end subroutine cmpt_bulk_forcing_alternative
end module mod_rk
