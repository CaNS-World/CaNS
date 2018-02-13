module mod_rk
use mod_param, only: forceinx,forceiny,forceinz, &
                     velfx   ,velfy   ,velfz
use mod_debug, only: chkmean
#ifdef IMPDIFF
  use mod_momd , only: momxpd,momypd,momzpd,momxa,momya,momza
#else
  use mod_mom , only: momxad,momyad,momzad,momxp,momyp,momzp
#endif
  implicit none
  private
  public rk
  contains
  subroutine rk(rkpar,n,dli,dzci,dzfi,dzflzi,visc,dt,l,u,v,w,p,dudtrko,dvdtrko,dwdtrko,tauxo,tauyo,tauzo,up,vp,wp,f)
    !
    ! low-storage 3rd-order Runge-Kutta scheme 
    ! for time integration of the momentum equations.
    !
    implicit none
    real(8), intent(in), dimension(2) :: rkpar
    integer, intent(in), dimension(3) :: n
    real(8), intent(in) :: visc,dt
    real(8), intent(in   ), dimension(3) :: dli,l 
    real(8), intent(in   ), dimension(0:) :: dzci,dzfi,dzflzi
    real(8), intent(in   ), dimension(0:,0:,0:) :: u ,v ,w,p
    real(8), intent(inout), dimension(n(1),n(2),n(3)) :: dudtrko,dvdtrko,dwdtrko
    real(8), intent(inout), dimension(3) :: tauxo,tauyo,tauzo
    real(8), intent(out), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1) :: up,vp,wp
    real(8), intent(out), dimension(3) :: f
    real(8),              dimension(n(1),n(2),n(3)) ::          dudtrk, dvdtrk, dwdtrk
#ifdef IMPDIFF
    real(8),              dimension(n(1),n(2),n(3)) ::          dudtrkd, dvdtrkd, dwdtrkd
    real(8) :: alpha
#endif
    real(8) :: factor1,factor2,factor12
    real(8), dimension(3) :: taux,tauy,tauz
    integer :: i,j,k
    real(8) :: mean
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2
    !
#ifdef IMPDIFF
    call momxpd(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,dzflzi,visc,p,u,dudtrk,dudtrkd,taux)
    call momypd(n(1),n(2),n(3),dli(2),dli(2),dzci,dzfi,dzflzi,visc,p,v,dvdtrk,dvdtrkd,tauy)
    call momzpd(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,dzflzi,visc,p,w,dwdtrk,dwdtrkd,tauz)
    f(1) = factor12*sum(taux(:)/l(:))
    f(2) = factor12*sum(tauy(:)/l(:))
    f(3) = factor12*sum(tauz(:)/l(:))
    ! alternatively, calculate force from the mean velocity directly
#else
    call momxp(n(1),n(2),n(3),dli(1),p,dudtrk)
    call momyp(n(1),n(2),n(3),dli(2),p,dvdtrk)
    call momzp(n(1),n(2),n(3),dzci  ,p,dwdtrk)
#endif
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,factor12,u,v,w,up,vp,wp,dudtrk,dvdtrk,dwdtrk)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          up(i,j,k) = u(i,j,k) + factor12*dudtrk(i,j,k)
          vp(i,j,k) = v(i,j,k) + factor12*dvdtrk(i,j,k)
          wp(i,j,k) = w(i,j,k) + factor12*dwdtrk(i,j,k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
#ifdef IMPDIFF
    call momxa(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,u,v,w,dudtrk)
    call momya(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,u,v,w,dvdtrk)
    call momza(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,u,v,w,dwdtrk)
#else
    call momxad(n(1),n(2),n(3),dli(1),dli(2),dli(3),dzci,dzfi,dzflzi,visc,u,v,w,dudtrk,taux)
    call momyad(n(1),n(2),n(3),dli(1),dli(2),dli(3),dzci,dzfi,dzflzi,visc,u,v,w,dvdtrk,tauy)
    call momzad(n(1),n(2),n(3),dli(1),dli(2),dli(3),dzci,dzfi,dzflzi,visc,u,v,w,dwdtrk,tauz)
    f(1) = (factor1*sum(taux(:)/l(:)) + factor2*sum(tauxo(:)/l(:)))
    f(2) = (factor1*sum(tauy(:)/l(:)) + factor2*sum(tauyo(:)/l(:)))
    f(3) = (factor1*sum(tauz(:)/l(:)) + factor2*sum(tauzo(:)/l(:)))
    tauxo(:) = taux(:)
    tauyo(:) = tauy(:)
    tauzo(:) = tauz(:)
#endif
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,factor1,factor2,up,vp,wp,dudtrk,dvdtrk,dwdtrk,dudtrko,dvdtrko,dwdtrko)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          ! could be split in two loops, because factor2=0 for istep=1, but like this reads nicer
          up(i,j,k) = up(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k)
          vp(i,j,k) = vp(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k)
          wp(i,j,k) = wp(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k)
          dudtrko(i,j,k) = dudtrk(i,j,k)
          dvdtrko(i,j,k) = dvdtrk(i,j,k)
          dwdtrko(i,j,k) = dwdtrk(i,j,k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    !
    ! bulk velocity forcing
    !
    f(:) = 0.
    if(forceinx) then
      call chkmean(n,dzflzi,up,mean)
      f(1) = velfx - mean
    endif
    if(forceiny) then
      call chkmean(n,dzflzi,vp,mean)
      f(2) = velfy - mean
    endif
    if(forceinz) then
      call chkmean(n,dzflzi,wp,mean) ! here should be dzclzi
      f(3) = velfz - mean
    endif
#ifdef IMPDIFF
    ! compute rhs of helmholtz equation
    alpha = -1.d0/(.5d0*visc*factor12)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,factor12,factor2,visc,up,vp,wp,dudtrkd,dvdtrkd,dwdtrkd,alpha,f)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          up(i,j,k) = ( up(i,j,k) + f(1) - .5d0*factor12*dudtrkd(i,j,k) )*alpha
          vp(i,j,k) = ( vp(i,j,k) + f(2) - .5d0*factor12*dvdtrkd(i,j,k) )*alpha
          wp(i,j,k) = ( wp(i,j,k) + f(3) - .5d0*factor12*dwdtrkd(i,j,k) )*alpha
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
#endif
    return
  end subroutine rk
end module mod_rk
