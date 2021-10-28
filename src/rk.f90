module mod_rk
  use mod_param, only: is_forced,velf
  use mod_debug, only: chk_mean
  use mod_mom  , only: momxad,momyad,momzad,momxp,momyp,momzp
  use mod_momd , only: momxpd,momypd,momzpd,momxa,momya,momza
  use mod_scal , only: scal
  use mod_types
  implicit none
  private
  public rk,rk_id
  contains
  subroutine rk(rkpar,n,dli,dzci,dzfi,dzflzi,dzclzi,visc,dt,l,u,v,w,p,dudtrko,dvdtrko,dwdtrko,tauxo,tauyo,tauzo,up,vp,wp,f)
    !
    ! low-storage 3rd-order Runge-Kutta scheme
    ! for time integration of the momentum equations.
    !
    implicit none
    real(rp), intent(in), dimension(2) :: rkpar
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in) :: visc,dt
    real(rp), intent(in   ), dimension(3) :: dli,l
    real(rp), intent(in   ), dimension(0:) :: dzci,dzfi,dzflzi,dzclzi
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u,v,w,p
    real(rp), intent(inout), dimension(:,:,:) :: dudtrko,dvdtrko,dwdtrko
    real(rp), intent(inout), dimension(3) :: tauxo,tauyo,tauzo
    real(rp), intent(out), dimension(0:,0:,0:) :: up,vp,wp
    real(rp), intent(out), dimension(3) :: f
    real(rp),              dimension(n(1),n(2),n(3)) :: dudtrk,dvdtrk,dwdtrk
    real(rp) :: factor1,factor2,factor12
    real(rp), dimension(3) :: taux,tauy,tauz
    integer :: i,j,k
    real(rp) :: mean
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2
    !
    call momxp(n(1),n(2),n(3),dli(1),p,dudtrk)
    call momyp(n(1),n(2),n(3),dli(2),p,dvdtrk)
    call momzp(n(1),n(2),n(3),dzci  ,p,dwdtrk)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,factor12,u,v,w,up,vp,wp,dudtrk,dvdtrk,dwdtrk)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          up(i,j,k) = u(i,j,k) + factor12*dudtrk(i,j,k)
          vp(i,j,k) = v(i,j,k) + factor12*dvdtrk(i,j,k)
          wp(i,j,k) = w(i,j,k) + factor12*dwdtrk(i,j,k)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    call momxad(n(1),n(2),n(3),dli(1),dli(2),dli(3),dzci,dzfi,dzflzi,visc,u,v,w,dudtrk,taux)
    call momyad(n(1),n(2),n(3),dli(1),dli(2),dli(3),dzci,dzfi,dzflzi,visc,u,v,w,dvdtrk,tauy)
    call momzad(n(1),n(2),n(3),dli(1),dli(2),dli(3),dzci,dzfi,dzflzi,visc,u,v,w,dwdtrk,tauz)
    f(1) = (factor1*sum(taux(:)/l(:)) + factor2*sum(tauxo(:)/l(:)))
    f(2) = (factor1*sum(tauy(:)/l(:)) + factor2*sum(tauyo(:)/l(:)))
    f(3) = (factor1*sum(tauz(:)/l(:)) + factor2*sum(tauzo(:)/l(:)))
    tauxo(:) = taux(:)
    tauyo(:) = tauy(:)
    tauzo(:) = tauz(:)
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
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    !
    ! bulk velocity forcing
    !
    f(:) = 0.
    if(is_forced(1)) then
      call chk_mean(n,dzflzi*1./(dli(1)*dli(2)*l(1)*l(2)),up,mean)
      f(1) = velf(1) - mean
    end if
    if(is_forced(2)) then
      call chk_mean(n,dzflzi*1./(dli(1)*dli(2)*l(1)*l(2)),vp,mean)
      f(2) = velf(2) - mean
    end if
    if(is_forced(3)) then
      call chk_mean(n,dzclzi*1./(dli(1)*dli(2)*l(1)*l(2)),wp,mean)
      f(3) = velf(3) - mean
    end if
  end subroutine rk
  subroutine rk_id(rkpar,n,dli,dzci,dzfi,dzflzi,dzclzi,visc,dt,l,u,v,w,p,dudtrko,dvdtrko,dwdtrko,tauxo,tauyo,tauzo,up,vp,wp,f)
    !
    ! low-storage 3rd-order Runge-Kutta scheme
    ! for time integration of the momentum equations with implicit diffusion.
    !
    implicit none
    real(rp), intent(in), dimension(2) :: rkpar
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in) :: visc,dt
    real(rp), intent(in   ), dimension(3) :: dli,l
    real(rp), intent(in   ), dimension(0:) :: dzci,dzfi,dzflzi,dzclzi
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u,v,w,p
    real(rp), intent(inout), dimension(:,:,:) :: dudtrko,dvdtrko,dwdtrko
    real(rp), intent(inout), dimension(3) :: tauxo,tauyo,tauzo
    real(rp), intent(out), dimension(0:,0:,0:) :: up,vp,wp
    real(rp), intent(out), dimension(3) :: f
    real(rp),              dimension(n(1),n(2),n(3)) :: dudtrk , dvdtrk , dwdtrk
    real(rp),              dimension(n(1),n(2),n(3)) :: dudtrkd, dvdtrkd, dwdtrkd
    real(rp) :: factor1,factor2,factor12
    real(rp), dimension(3) :: taux,tauy,tauz
    integer :: i,j,k
    real(rp) :: mean
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2
    !
    call momxpd(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,dzflzi,visc,p,u,dudtrk,dudtrkd,taux)
    call momypd(n(1),n(2),n(3),dli(2),dli(2),dzci,dzfi,dzflzi,visc,p,v,dvdtrk,dvdtrkd,tauy)
    call momzpd(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,dzflzi,visc,p,w,dwdtrk,dwdtrkd,tauz)
    f(1) = factor12*sum(taux(:)/l(:))
    f(2) = factor12*sum(tauy(:)/l(:))
    f(3) = factor12*sum(tauz(:)/l(:))
    ! alternatively, calculate force from the mean velocity directly
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,factor12,u,v,w,up,vp,wp,dudtrk,dvdtrk,dwdtrk)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          up(i,j,k) = u(i,j,k) + factor12*dudtrk(i,j,k)
          vp(i,j,k) = v(i,j,k) + factor12*dvdtrk(i,j,k)
          wp(i,j,k) = w(i,j,k) + factor12*dwdtrk(i,j,k)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    call momxa(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,u,v,w,dudtrk)
    call momya(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,u,v,w,dvdtrk)
    call momza(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,u,v,w,dwdtrk)
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
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    !
    ! bulk velocity forcing
    !
    f(:) = 0.
    if(is_forced(1)) then
      call chk_mean(n,dzflzi*1./(dli(1)*dli(2)*l(1)*l(2)),up,mean)
      f(1) = velf(1) - mean
    end if
    if(is_forced(2)) then
      call chk_mean(n,dzflzi*1./(dli(1)*dli(2)*l(1)*l(2)),vp,mean)
      f(2) = velf(2) - mean
    end if
    if(is_forced(3)) then
      call chk_mean(n,dzclzi*1./(dli(1)*dli(2)*l(1)*l(2)),wp,mean)
      f(3) = velf(3) - mean
    end if
    !
    ! compute rhs of helmholtz equation
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,factor12,factor2,visc,up,vp,wp,dudtrkd,dvdtrkd,dwdtrkd)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          up(i,j,k) = up(i,j,k) - .5*factor12*dudtrkd(i,j,k)
          vp(i,j,k) = vp(i,j,k) - .5*factor12*dvdtrkd(i,j,k)
          wp(i,j,k) = wp(i,j,k) - .5*factor12*dwdtrkd(i,j,k)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine rk_id
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
