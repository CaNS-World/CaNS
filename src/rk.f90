! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2025 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_rk
  use mod_mom  , only: momx_a,momy_a,momz_a, &
                       momx_d,momy_d,momz_d, &
                       momx_p,momy_p,momz_p, &
                       cmpt_wallshear, &
                       momx_d_xy,momy_d_xy,momz_d_xy, &
                       momx_d_z ,momy_d_z ,momz_d_z, &
                       mom_xyz_ad
  use mod_param, only: is_impdiff,is_impdiff_1d,is_boussinesq_buoyancy,is_fast_mom_kernels
  use mod_scal , only: scal,cmpt_scalflux,scalar
  use mod_utils, only: bulk_mean
  use mod_types
  implicit none
  public rk,rk_scal
  contains
  subroutine rk(rkpar,n,dli,dzci,dzfi,grid_vol_ratio_c,grid_vol_ratio_f,visc,dt,p, &
                is_forced,velf,bforce,gacc,beta,scalars,dudtrko,dvdtrko,dwdtrko,u,v,w,f)
#if defined(_OPENACC)
    use mod_common_cudecomp, only: dudtrk_t => work, &
                                   dvdtrk_t => solver_buf_0, &
                                   dwdtrk_t => solver_buf_1
#endif
    !
    ! low-storage 3rd-order Runge-Kutta scheme
    ! for time integration of the momentum equations.
    !
    implicit none
    real(rp), intent(in   ), dimension(2)        :: rkpar
    integer , intent(in   ), dimension(3)        :: n
    real(rp), intent(in   )                      :: visc,dt
    real(rp), intent(in   ), dimension(3)        :: dli
    real(rp), intent(in   ), dimension(0:)       :: dzci,dzfi
    real(rp), intent(in   ), dimension(0:)       :: grid_vol_ratio_c,grid_vol_ratio_f
    real(rp), intent(in   ), dimension(0:,0:,0:) :: p
    logical , intent(in   ), dimension(3)        :: is_forced
    real(rp), intent(in   ), dimension(3)        :: velf,bforce
    real(rp), intent(in   ), dimension(3)        :: gacc
    real(rp), intent(in   )                      :: beta
    type(scalar), intent(in   ), target, dimension(:), optional :: scalars
    real(rp), intent(inout), dimension(1:,1:,1:) :: dudtrko,dvdtrko,dwdtrko
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out  ), dimension(3)        :: f
    real(rp), pointer, contiguous, dimension(:,:,:) :: s
#if !defined(_OPENACC)
    real(rp), target       , allocatable, dimension(:,:,:), save :: dudtrk_t,dvdtrk_t,dwdtrk_t
#endif
    real(rp), pointer      , contiguous , dimension(:,:,:), save :: dudtrk  ,dvdtrk  ,dwdtrk
    real(rp),                allocatable, dimension(:,:,:), save :: dudtrkd ,dvdtrkd ,dwdtrkd
    logical, save :: is_first = .true.
    real(rp) :: factor1,factor2,factor12
    logical :: is_buoyancy
    integer  :: i,j,k
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2
    !
    ! is_boussinesq_buoyancy = T will always imply present(scalars), as we do not allow for
    ! nscal = 0 with buoyancy on (see `param.f90`)
    !
    is_buoyancy = present(scalars).and.is_boussinesq_buoyancy
    if(is_buoyancy) then
      s => scalars(1)%val
    end if
    !
    ! initialization
    !
    if(is_first) then ! leverage save attribute to allocate these arrays on the device only once
      is_first = .false.
#if !defined(_OPENACC)
      allocate(dudtrk_t(n(1),n(2),n(3)),dvdtrk_t(n(1),n(2),n(3)),dwdtrk_t(n(1),n(2),n(3)))
#endif
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP parallel do   collapse(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            dudtrko(i,j,k) = 0._rp
            dvdtrko(i,j,k) = 0._rp
            dwdtrko(i,j,k) = 0._rp
          end do
        end do
      end do
      if(is_impdiff) then
        allocate(dudtrkd(n(1),n(2),n(3)),dvdtrkd(n(1),n(2),n(3)),dwdtrkd(n(1),n(2),n(3)))
        !$acc enter data create(dudtrkd,dvdtrkd,dwdtrkd) async(1)
      end if
#if defined(_OPENACC)
      dudtrk(1:n(1),1:n(2),1:n(3)) => dudtrk_t(1:product(n(:)))
      dvdtrk(1:n(1),1:n(2),1:n(3)) => dvdtrk_t(1:product(n(:)))
      dwdtrk(1:n(1),1:n(2),1:n(3)) => dwdtrk_t(1:product(n(:)))
#else
      dudtrk  => dudtrk_t
      dvdtrk  => dvdtrk_t
      dwdtrk  => dwdtrk_t
#endif
    end if
    !
    if(is_fast_mom_kernels) then
      call mom_xyz_ad(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,u,v,w,dudtrk,dvdtrk,dwdtrk,dudtrkd,dvdtrkd,dwdtrkd)
    else
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP parallel do   collapse(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            dudtrk(i,j,k) = 0._rp
            dvdtrk(i,j,k) = 0._rp
            dwdtrk(i,j,k) = 0._rp
          end do
        end do
      end do
      if(.not.is_impdiff) then
        call momx_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,u,dudtrk)
        call momy_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,v,dvdtrk)
        call momz_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,w,dwdtrk)
      else
        !$acc parallel loop collapse(3) default(present) async(1)
        !$OMP parallel do   collapse(3) DEFAULT(shared)
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              dudtrkd(i,j,k) = 0._rp
              dvdtrkd(i,j,k) = 0._rp
              dwdtrkd(i,j,k) = 0._rp
            end do
          end do
        end do
        if(.not.is_impdiff_1d) then
          call momx_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,u,dudtrkd)
          call momy_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,v,dvdtrkd)
          call momz_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,w,dwdtrkd)
        else
          call momx_d_xy(n(1),n(2),n(3),dli(1),dli(2),visc,u,dudtrk )
          call momy_d_xy(n(1),n(2),n(3),dli(1),dli(2),visc,v,dvdtrk )
          call momz_d_xy(n(1),n(2),n(3),dli(1),dli(2),visc,w,dwdtrk )
          call momx_d_z( n(1),n(2),n(3),dzci  ,dzfi  ,visc,u,dudtrkd)
          call momy_d_z( n(1),n(2),n(3),dzci  ,dzfi  ,visc,v,dvdtrkd)
          call momz_d_z( n(1),n(2),n(3),dzci  ,dzfi  ,visc,w,dwdtrkd)
        end if
        call momx_a(n(1),n(2),n(3),dli(1),dli(2),dzfi,u,v,w,dudtrk)
        call momy_a(n(1),n(2),n(3),dli(1),dli(2),dzfi,u,v,w,dvdtrk)
        call momz_a(n(1),n(2),n(3),dli(1),dli(2),dzci,u,v,w,dwdtrk)
      end if
    end if
    !
#if !defined(_LOOP_UNSWITCHING)
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          u(i,j,k) = u(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k) + &
                                factor12*(bforce(1) - dli(1)*( p(i+1,j,k)-p(i,j,k)))
          if(is_buoyancy) then
            u(i,j,k) = u(i,j,k) - factor12*gacc(1)*beta*0.5*(s(i+1,j,k)+s(i,j,k))
          end if
          !
          v(i,j,k) = v(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k) + &
                                factor12*(bforce(2) - dli(2)*( p(i,j+1,k)-p(i,j,k)))
          if(is_buoyancy) then
            v(i,j,k) = v(i,j,k) - factor12*gacc(2)*beta*0.5*(s(i,j+1,k)+s(i,j,k))
          end if
          !
          w(i,j,k) = w(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k) + &
                                factor12*(bforce(3) - dzci(k)*(p(i,j,k+1)-p(i,j,k)))
          if(is_buoyancy) then
            w(i,j,k) = w(i,j,k) - factor12*gacc(3)*beta*0.5*(s(i,j,k+1)+s(i,j,k))
          end if
          !
          if(is_impdiff) then
            u(i,j,k) = u(i,j,k) + factor12*dudtrkd(i,j,k)
            v(i,j,k) = v(i,j,k) + factor12*dvdtrkd(i,j,k)
            w(i,j,k) = w(i,j,k) + factor12*dwdtrkd(i,j,k)
          end if
        end do
      end do
    end do
#else
    if(.not.is_impdiff .and. .not.is_buoyancy) then
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            u(i,j,k) = u(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k) + &
                                  factor12*(bforce(1) - dli(1)*(p(i+1,j,k)-p(i,j,k)))
            v(i,j,k) = v(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k) + &
                                  factor12*(bforce(2) - dli(2)*(p(i,j+1,k)-p(i,j,k)))
            w(i,j,k) = w(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k) + &
                                  factor12*(bforce(3) - dzci(k)*(p(i,j,k+1)-p(i,j,k)))
          end do
        end do
      end do
    else if(is_impdiff .and. .not.is_buoyancy) then
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            u(i,j,k) = u(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k) + &
                                  factor12*(bforce(1) - dli(1)*(p(i+1,j,k)-p(i,j,k)) + &
                                            dudtrkd(i,j,k))
            v(i,j,k) = v(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k) + &
                                  factor12*(bforce(2) - dli(2)*(p(i,j+1,k)-p(i,j,k)) + &
                                            dvdtrkd(i,j,k))
            w(i,j,k) = w(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k) + &
                                  factor12*(bforce(3) - dzci(k)*(p(i,j,k+1)-p(i,j,k)) + &
                                            dwdtrkd(i,j,k))
          end do
        end do
      end do
    else if(.not.is_impdiff .and. is_buoyancy) then
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            u(i,j,k) = u(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k) + &
                                  factor12*(bforce(1) - dli(1)*(p(i+1,j,k)-p(i,j,k))) &
                                          - gacc(1)*beta*0.5*(s(i+1,j,k)+s(i,j,k))
            v(i,j,k) = v(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k) + &
                                  factor12*(bforce(2) - dli(2)*(p(i,j+1,k)-p(i,j,k))) &
                                          - gacc(2)*beta*0.5*(s(i,j+1,k)+s(i,j,k))
            w(i,j,k) = w(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k) + &
                                  factor12*(bforce(3) - dzci(k)*(p(i,j,k+1)-p(i,j,k))) &
                                          - gacc(3)*beta*0.5*(s(i,j,k+1)+s(i,j,k))
          end do
        end do
      end do
    else
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            u(i,j,k) = u(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k) + &
                                  factor12*(bforce(1) - dli(1)*(p(i+1,j,k)-p(i,j,k)) &
                                          - gacc(1)*beta*0.5*(s(i+1,j,k)+s(i,j,k)) + &
                                            dudtrkd(i,j,k))
            v(i,j,k) = v(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k) + &
                                  factor12*(bforce(2) - dli(2)*(p(i,j+1,k)-p(i,j,k)) &
                                          - gacc(2)*beta*0.5*(s(i,j+1,k)+s(i,j,k)) + &
                                            dvdtrkd(i,j,k))
            w(i,j,k) = w(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k) + &
                                  factor12*(bforce(3) - dzci(k)*(p(i,j,k+1)-p(i,j,k)) &
                                          - gacc(3)*beta*0.5*(s(i,j,k+1)+s(i,j,k)) + &
                                            dwdtrkd(i,j,k))
          end do
        end do
      end do
    end if
#endif
    !
    ! replaced previous pointer swap to save memory on GPUs by using already allocated
    ! buffers
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP parallel do   collapse(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          dudtrko(i,j,k) = dudtrk(i,j,k)
          dvdtrko(i,j,k) = dvdtrk(i,j,k)
          dwdtrko(i,j,k) = dwdtrk(i,j,k)
        end do
      end do
    end do
!#if 0 /*pressure gradient term treated explicitly above */
!    !$acc parallel loop collapse(3) default(present) async(1)
!    !$OMP parallel do   collapse(3) DEFAULT(shared)
!    do k=1,n(3)
!      do j=1,n(2)
!        do i=1,n(1)
!          dudtrk(i,j,k) = 0._rp
!          dvdtrk(i,j,k) = 0._rp
!          dwdtrk(i,j,k) = 0._rp
!        end do
!      end do
!    end do
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
    !
    ! compute bulk velocity forcing
    !
    call cmpt_bulk_forcing(n,is_forced,velf,grid_vol_ratio_c,grid_vol_ratio_f,u,v,w,f)
    !
    if(is_impdiff) then
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
    end if
  end subroutine rk
  !
  subroutine rk_scal(rkpar,n,dli,l,dzci,dzfi,grid_vol_ratio_f,alpha,dt,is_bound,u,v,w, &
                     is_forced,scalf,ssource,fluxo,dsdtrko,s,f)
#if defined(_OPENACC)
    use mod_common_cudecomp, only: dsdtrk_t => work
#endif
    !
    ! low-storage 3rd-order Runge-Kutta scheme
    ! for time integration of the scalar field.
    !
    implicit none
    logical , parameter :: is_cmpt_wallflux = .false.
    real(rp), intent(in   ), dimension(2) :: rkpar
    integer , intent(in   ), dimension(3) :: n
    real(rp), intent(in   ), dimension(3) :: dli,l
    real(rp), intent(in   ), dimension(0:) :: dzci,dzfi
    real(rp), intent(in   ), dimension(:) :: grid_vol_ratio_f
    real(rp), intent(in   ) :: alpha,dt
    logical , intent(in   ), dimension(0:1,3)    :: is_bound
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u,v,w
    logical , intent(in   ) :: is_forced
    real(rp), intent(in   ) :: scalf,ssource
    real(rp), intent(inout), dimension(0:1,3) :: fluxo
    real(rp), intent(inout), dimension(1:,1:,1:) :: dsdtrko
    real(rp), intent(inout), dimension(0:,0:,0:) :: s
    real(rp), intent(out  ) :: f
    !
#if !defined(_OPENACC)
    real(rp), target       , allocatable, dimension(:,:,:), save :: dsdtrk_t
#endif
    real(rp), pointer      , contiguous , dimension(:,:,:), save :: dsdtrk
    real(rp), target       , allocatable, dimension(:,:,:), save :: dsdtrkd
    logical, save :: is_first = .true.
    !
    real(rp) :: factor1,factor2,factor12
    real(rp), dimension(0:1,3) :: flux
    integer :: i,j,k
    real(rp) :: mean
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2
    if(is_first) then ! leverage save attribute to allocate these arrays on the device only once
      is_first = .false.
#if !defined(_OPENACC)
      allocate(dsdtrk_t(1:n(1),1:n(2),1:n(3)))
#endif
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP parallel do   collapse(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            dsdtrko(i,j,k) = 0._rp
          end do
        end do
      end do
      if(is_impdiff) then
        allocate(dsdtrkd(n(1),n(2),n(3)))
        !$acc enter data create(dsdtrkd) async(1)
        !$acc parallel loop collapse(3) default(present) async(1)
        !$OMP parallel do   collapse(3) DEFAULT(shared)
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              dsdtrkd(i,j,k) = 0._rp
            end do
          end do
        end do
      end if
    end if
#if defined(_OPENACC)
    dsdtrk(1:n(1),1:n(2),1:n(3)) => dsdtrk_t(1:product(n(:)))
#else
    dsdtrk => dsdtrk_t
#endif
    !
    call scal(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,alpha,u,v,w,s,dsdtrk,dsdtrkd)
#if !defined(_LOOP_UNSWITCHING)
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          s(i,j,k) = s(i,j,k) + factor1*dsdtrk(i,j,k) + factor2*dsdtrko(i,j,k) + factor12*ssource
          if(is_impdiff) then
            s(i,j,k) = s(i,j,k) + factor12*dsdtrkd(i,j,k)
          end if
        end do
      end do
    end do
#else
    if(.not.is_impdiff) then
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            s(i,j,k) = s(i,j,k) + factor1*dsdtrk(i,j,k) + factor2*dsdtrko(i,j,k) + factor12*ssource
          end do
        end do
      end do
    else
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            s(i,j,k) = s(i,j,k) + factor1*dsdtrk(i,j,k) + factor2*dsdtrko(i,j,k) + &
                                  factor12*(ssource + dsdtrkd(i,j,k))
          end do
        end do
      end do
    end if
#endif
    !
    ! compute wall scalar flux
    !
    if(is_cmpt_wallflux) then
      call cmpt_scalflux(n,is_bound,l,dli,dzci,dzfi,alpha,s(:,:,:),flux)
      f = (factor1*sum((flux( 0,:)+flux( 1,:))/l(:)) + &
           factor2*sum((fluxo(0,:)+fluxo(1,:))/l(:)))
      fluxo(:,:) = flux(:,:)
    end if
    !
    ! bulk scalar forcing
    !
    if(is_forced) then
      call bulk_mean(n,grid_vol_ratio_f,s(:,:,:),mean)
      f = scalf - mean
    end if
    if(is_impdiff) then
      !
      ! compute rhs of Helmholtz equation
      !
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            s(i,j,k) = s(i,j,k) - .5_rp*factor12*dsdtrkd(i,j,k)
          end do
        end do
      end do
    end if
    !
    ! replaced previous pointer swap to save memory on GPUs by using already allocated
    ! buffers
    !
    !$acc parallel loop collapse(3) default(present) async(1) ! not really necessary
    !$OMP parallel do   collapse(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          dsdtrko(i,j,k) = dsdtrk(i,j,k)
        end do
      end do
    end do
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
    real(rp), intent(inout), dimension(0:1,3) :: tauxo,tauyo,tauzo
    real(rp), intent(inout), dimension(3) :: f
    real(rp), dimension(3) :: f_aux
    logical , intent(in   ) :: is_first
    real(rp), dimension(0:1,3) :: taux,tauy,tauz
    real(rp), dimension(3) :: taux_tot,tauy_tot,tauz_tot,tauxo_tot,tauyo_tot,tauzo_tot
    real(rp) :: factor1,factor2,factor12
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = (factor1 + factor2)/2.
    !
    call cmpt_wallshear(n,is_forced,is_bound,l,dli,dzci,dzfi,visc,u,v,w,taux,tauy,tauz)
    taux_tot(:) = sum(taux(0:1,:),1); tauxo_tot(:) = sum(tauxo(0:1,:),1)
    tauy_tot(:) = sum(tauy(0:1,:),1); tauyo_tot(:) = sum(tauyo(0:1,:),1)
    tauz_tot(:) = sum(tauz(0:1,:),1); tauzo_tot(:) = sum(tauzo(0:1,:),1)
    if(.not.is_impdiff) then
      if(is_first) then
        f(1) = (factor1*sum(taux_tot(:)/l(:)) + factor2*sum(tauxo_tot(:)/l(:)))
        f(2) = (factor1*sum(tauy_tot(:)/l(:)) + factor2*sum(tauyo_tot(:)/l(:)))
        f(3) = (factor1*sum(tauz_tot(:)/l(:)) + factor2*sum(tauzo_tot(:)/l(:)))
        tauxo(:,:) = taux(:,:)
        tauyo(:,:) = tauy(:,:)
        tauzo(:,:) = tauz(:,:)
      end if
    else
      if(is_impdiff_1d) then
        f_aux(1) = factor12*taux_tot(3)/l(3)
        f_aux(2) = factor12*tauy_tot(3)/l(3)
        if(is_first) then
          f(1) = factor1*taux_tot(2)/l(2) + factor2*tauxo_tot(2)/l(2) + f_aux(1)
          f(2) = factor1*tauy_tot(1)/l(1) + factor2*tauyo_tot(1)/l(1) + f_aux(2)
          f(3) = factor1*sum(tauz_tot(1:2)/l(1:2)) + factor2*sum(tauzo_tot(1:2)/l(1:2))
          tauxo(:,1:2) = taux(:,1:2)
          tauyo(:,1:2) = tauy(:,1:2)
          tauzo(:,1:2) = tauz(:,1:2)
        else
          f(1) = f(1) + f_aux(1)
          f(2) = f(2) + f_aux(2)
        end if
      else
        f_aux(:) = factor12*[sum(taux_tot(:)/l(:)), &
                             sum(tauy_tot(:)/l(:)), &
                             sum(tauz_tot(:)/l(:))]
        if(is_first) then
           f(:) = f_aux(:)
        else
           f(:) = f(:) + f_aux(:)
        end if
      end if
    end if
  end subroutine cmpt_bulk_forcing_alternative
end module mod_rk
