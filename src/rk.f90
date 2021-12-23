module mod_rk
  use mod_debug, only: chk_mean
  use mod_mom  , only: momx_a,momy_a,momz_a, &
                       momx_d,momy_d,momz_d, &
                       momx_p,momy_p,momz_p, cmpt_wallshear
#if defined(_IMPDIFF) && defined(_IMPDIFF_1D)
  use mod_mom  , only: momx_d_xy,momy_d_xy,momz_d_xy, &
                 only: momx_d_z ,momy_d_z ,momz_d_z
#endif
  use mod_scal , only: scal
  use mod_types
  implicit none
  private
  public rk
  contains
  subroutine rk(rkpar,n,dli,l,dzci,dzfi,visc,dt,u,v,w,p,is_bound,is_forced,velf,bforce, &
                dudtrko,dvdtrko,dwdtrko,tauxo,tauyo,tauzo,up,vp,wp,f)
    !
    ! low-storage 3rd-order Runge-Kutta scheme
    ! for time integration of the momentum equations.
    !
    implicit none
    real(rp), intent(in), dimension(2) :: rkpar
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in) :: visc,dt
    real(rp), intent(in   ), dimension(3) :: dli,l
    real(rp), intent(in   ), dimension(0:) :: dzci,dzfi
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u,v,w,p
    logical , intent(in   ), dimension(0:1,3)    :: is_bound
    logical , intent(in   ), dimension(3)        :: is_forced
    real(rp), intent(in   ), dimension(3)        :: velf,bforce
    real(rp), intent(inout), dimension(:,:,:) :: dudtrko,dvdtrko,dwdtrko
    real(rp), intent(inout), dimension(3) :: tauxo,tauyo,tauzo
    real(rp), intent(out), dimension(0:,0:,0:) :: up,vp,wp
    real(rp), intent(out), dimension(3) :: f
    real(rp),              dimension(n(1),n(2),n(3)) :: dudtrk ,dvdtrk ,dwdtrk
#if defined(_IMPDIFF)
    real(rp),              dimension(n(1),n(2),n(3)) :: dudtrkd,dvdtrkd,dwdtrkd
#endif
    real(rp) :: factor1,factor2,factor12
    real(rp), dimension(3) :: taux,tauy,tauz
    integer :: i,j,k
    real(rp) :: mean,grid_vol_ratio(0:n(3)+1)
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2
    !
    !$OMP WORKSHARE
    dudtrk(:,:,:) = 0._rp
    dvdtrk(:,:,:) = 0._rp
    dwdtrk(:,:,:) = 0._rp
    !$OMP END WORKSHARE
#if !defined(_IMPDIFF)
    call momx_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,u,dudtrk)
    call momy_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,v,dvdtrk)
    call momz_d(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,visc,w,dwdtrk)
#else
    !$OMP WORKSHARE
    dudtrkd(:,:,:) = 0._rp
    dvdtrkd(:,:,:) = 0._rp
    dwdtrkd(:,:,:) = 0._rp
    !$OMP END WORKSHARE
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
#endif
    call momx_a(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,u,v,w,dudtrk)
    call momy_a(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,u,v,w,dvdtrk)
    call momz_a(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,u,v,w,dwdtrk)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
#if defined(_IMPDIFF)
    !$OMP SHARED(factor12,dudtrkd,dvdtrkd,dwdtrkd) &
#endif
    !$OMP SHARED(n,factor1,factor2,u,v,w,up,vp,wp,dudtrk,dvdtrk,dwdtrk,dudtrko,dvdtrko,dwdtrko)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          up(i,j,k) = u(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k)
          vp(i,j,k) = v(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k)
          wp(i,j,k) = w(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k)
#if defined(_IMPDIFF)
          up(i,j,k) = up(i,j,k) + factor12*dudtrkd(i,j,k)
          vp(i,j,k) = vp(i,j,k) + factor12*dvdtrkd(i,j,k)
          wp(i,j,k) = wp(i,j,k) + factor12*dwdtrkd(i,j,k)
#endif
          dudtrko(i,j,k) = dudtrk(i,j,k)
          dvdtrko(i,j,k) = dvdtrk(i,j,k)
          dwdtrko(i,j,k) = dwdtrk(i,j,k)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
!#if 0 /*pressure gradient term treated explicitly later */
!    !$OMP WORKSHARE
!    dudtrk(:,:,:) = 0._rp
!    dvdtrk(:,:,:) = 0._rp
!    dwdtrk(:,:,:) = 0._rp
!    !$OMP END WORKSHARE
!    call momx_p(n(1),n(2),n(3),dli(1),bforce(1),p,dudtrk)
!    call momy_p(n(1),n(2),n(3),dli(2),bforce(2),p,dvdtrk)
!    call momz_p(n(1),n(2),n(3),dzci  ,bforce(3),p,dwdtrk)
!    !$OMP PARALLEL DO DEFAULT(none) &
!    !$OMP PRIVATE(i,j,k) &
!    !$OMP SHARED(n,factor12,u,v,w,up,vp,wp,dudtrk,dvdtrk,dwdtrk)
!    do k=1,n(3)
!      do j=1,n(2)
!        do i=1,n(1)
!          up(i,j,k) = up(i,j,k) + factor12*dudtrk(i,j,k)
!          vp(i,j,k) = vp(i,j,k) + factor12*dvdtrk(i,j,k)
!          wp(i,j,k) = wp(i,j,k) + factor12*dwdtrk(i,j,k)
!        end do
!      end do
!    end do
!    !$OMP END PARALLEL DO
!#endif
    ! pressure gradient term in the loop below
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,factor12,dli,dzci,bforce,u,v,w,up,vp,wp,p)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          up(i,j,k) = up(i,j,k) + factor12*(bforce(1) - dli(1)*( p(i+1,j,k)-p(i,j,k)))
          vp(i,j,k) = vp(i,j,k) + factor12*(bforce(2) - dli(2)*( p(i,j+1,k)-p(i,j,k)))
          wp(i,j,k) = wp(i,j,k) + factor12*(bforce(3) - dzci(k)*(p(i,j,k+1)-p(i,j,k)))
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    !
    ! compute mean wall shear stresses
    !
    call cmpt_wallshear(n,is_bound,l,dli,dzci,dzfi,visc,u,v,w,taux,tauy,tauz)
#if !defined(_IMPDIFF)
    f(1) = (factor1*sum(taux(:)/l(:)) + factor2*sum(tauxo(:)/l(:)))
    f(2) = (factor1*sum(tauy(:)/l(:)) + factor2*sum(tauyo(:)/l(:)))
    f(3) = (factor1*sum(tauz(:)/l(:)) + factor2*sum(tauzo(:)/l(:)))
    tauxo(:) = taux(:)
    tauyo(:) = tauy(:)
    tauzo(:) = tauz(:)
#else
    f(:) = factor12*[sum(taux(:)/l(:)), &
                     sum(tauy(:)/l(:)), &
                     sum(tauz(:)/l(:))]
#endif
    !
    ! bulk velocity forcing
    !
    f(:) = 0.
    if(is_forced(1)) then
      grid_vol_ratio(:) = 1./(dzfi(:)*dli(1)*dli(2)*l(1)*l(2)*l(3))
      call chk_mean(n,grid_vol_ratio,up,mean)
      f(1) = velf(1) - mean
    end if
    if(is_forced(2)) then
      grid_vol_ratio(:) = 1./(dzfi(:)*dli(1)*dli(2)*l(1)*l(2)*l(3))
      call chk_mean(n,grid_vol_ratio,vp,mean)
      f(2) = velf(2) - mean
    end if
    if(is_forced(3)) then
      grid_vol_ratio(:) = 1./(dzci(:)*dli(1)*dli(2)*l(1)*l(2)*l(3))
      call chk_mean(n,grid_vol_ratio,wp,mean)
      f(3) = velf(3) - mean
    end if
#if defined(_IMPDIFF)
    !
    ! compute rhs of helmholtz equation
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,factor12,factor2,visc,up,vp,wp,dudtrkd,dvdtrkd,dwdtrkd)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          up(i,j,k) = up(i,j,k) - .5_rp*factor12*dudtrkd(i,j,k)
          vp(i,j,k) = vp(i,j,k) - .5_rp*factor12*dvdtrkd(i,j,k)
          wp(i,j,k) = wp(i,j,k) - .5_rp*factor12*dwdtrkd(i,j,k)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
#endif
  end subroutine rk
  subroutine rk_scal(rkpar,n,dli,dzci,dzfi,visc,dt,u,v,w,dsdtrko,s)
    !
    ! low-storage 3rd-order Runge-Kutta scheme
    ! for time integration of the scalar field.
    !
    implicit none
    real(rp), intent(in   ), dimension(2) :: rkpar
    integer , intent(in   ), dimension(3) :: n
    real(rp), intent(in   ), dimension(3) :: dli
    real(rp), intent(in   ), dimension(0:) :: dzci,dzfi
    real(rp), intent(in   ) :: visc,dt
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(inout), dimension(:,:,:) :: dsdtrko
    real(rp), intent(inout), dimension(0:,0:,0:) :: s
    real(rp),              dimension(n(1),n(2),n(3)) :: dsdtrk
    real(rp) :: factor1,factor2,factor12
    integer :: i,j,k
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2
    call scal(n(1),n(2),n(3),dli(1),dli(2),dli(3),dzci,dzfi,visc,u,v,w,s,dsdtrk)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,factor1,factor2,s,dsdtrk,dsdtrko)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          s(i,j,k) = s(i,j,k) + factor1*dsdtrk(i,j,k) + factor2*dsdtrko(i,j,k)
          dsdtrko(i,j,k) = dsdtrk(i,j,k)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine rk_scal
end module mod_rk
