! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_initgrid
  use mod_param, only:pi
  use mod_types
  implicit none
  private
  public initgrid
  contains
  subroutine initgrid(gtype,n,gr,lz,dzc,dzf,zc,zf)
    !
    ! initializes the non-uniform grid along z
    !
    implicit none
    integer, parameter :: CLUSTER_TWO_END              = 1, &
                          CLUSTER_ONE_END              = 2, &
                          CLUSTER_ONE_END_R            = 3, &
                          CLUSTER_MIDDLE               = 4
    integer , intent(in ) :: gtype,n
    real(rp), intent(in ) :: gr,lz
    real(rp), intent(out), dimension(0:n+1) :: dzc,dzf,zc,zf
    real(rp) :: z0
    real(dp),              dimension(0:n+1) :: dzc_aux,dzf_aux,zc_aux,zf_aux
    integer :: k
    procedure (), pointer :: gridpoint => null()
    select case(gtype)
    case(CLUSTER_TWO_END)
      gridpoint => gridpoint_cluster_two_end
    case(CLUSTER_ONE_END)
      gridpoint => gridpoint_cluster_one_end
    case(CLUSTER_ONE_END_R)
      gridpoint => gridpoint_cluster_one_end_r
    case(CLUSTER_MIDDLE)
      gridpoint => gridpoint_cluster_middle
    case default
      gridpoint => gridpoint_cluster_two_end
    end select
    !
    ! step 1) determine coordinates of cell faces zf_aux
    !
    zf_aux(0) = 0.
    do k=1,n
      z0  = (k-0.)/(1.*n)
#if !defined(_GRIDPOINT_NATURAL_CHANNEL)
      call gridpoint(gr,z0,zf_aux(k))
#else
      call gridpoint_natural(k,n,zf_aux(k))
#endif
      zf_aux(k) = zf_aux(k)*lz
    end do
    !
    ! step 2) determine grid spacing between faces dzf_aux
    !
    do k=1,n
      dzf_aux(k) = zf_aux(k)-zf_aux(k-1)
    end do
    dzf_aux(0  ) = dzf_aux(1)
    dzf_aux(n+1) = dzf_aux(n)
    !
    ! step 3) determine grid spacing between centers dzc_aux
    !
    do k=0,n
      dzc_aux(k) = .5*(dzf_aux(k)+dzf_aux(k+1))
    end do
    dzc_aux(n+1) = dzc_aux(n)
    !
    ! step 4) compute coordinates of cell centers zc_aux and faces zf_aux
    !
    zc_aux(0)    = -dzc_aux(0)/2.
    zf_aux(0)    = 0.
    do k=1,n+1
      zc_aux(k) = zc_aux(k-1) + dzc_aux(k-1)
      zf_aux(k) = zf_aux(k-1) + dzf_aux(k)
    end do
    zf(:) = zf_aux(:)
    zc(:) = zc_aux(:)
    dzf(:) = dzf_aux(:)
    dzc(:) = dzc_aux(:)
  end subroutine initgrid
  !
  ! grid stretching functions
  ! see e.g., Fluid Flow Phenomena -- A Numerical Toolkit, by P. Orlandi
  !           Pirozzoli et al. JFM 788, 614–639 (commented)
  !
  subroutine gridpoint_cluster_two_end(alpha,z0,z)
    !
    ! clustered at the two sides
    !
    implicit none
    real(rp), intent(in) :: alpha,z0
    real(rp), intent(out) :: z
    if(alpha /= 0.) then
      z = 0.5*(1.+tanh((z0-0.5)*alpha)/tanh(alpha/2.))
      !z = 0.5*(1.+erf( (z0-0.5)*alpha)/erf( alpha/2.))
    else
      z = z0
    end if
  end subroutine gridpoint_cluster_two_end
  subroutine gridpoint_cluster_one_end(alpha,z0,z)
    !
    ! clustered at the lower side
    !
    implicit none
    real(rp), intent(in) :: alpha,z0
    real(rp), intent(out) :: z
    if(alpha /= 0.) then
      z = 1.0*(1.+tanh((z0-1.0)*alpha)/tanh(alpha/1.))
      !z = 1.0*(1.+erf( (z0-1.0)*alpha)/erf( alpha/1.))
    else
      z = z0
    end if
  end subroutine gridpoint_cluster_one_end
  subroutine gridpoint_cluster_one_end_r(alpha,r0,r)
    !
    ! clustered at the upper side
    !
    implicit none
    real(rp), intent(in ) :: alpha,r0
    real(rp), intent(out) :: r
    if(alpha /= 0._rp) then
      r = 1._rp-1.0_rp*(1._rp+tanh((1._rp-r0-1.0_rp)*alpha)/tanh(alpha/1._rp))
      !r = 1._rp-1.0_rp*(1._rp+erf( (1._rp-r0-1.0_rp)*alpha)/erf( alpha/1._rp))
    else
      r = r0
    end if
  end subroutine gridpoint_cluster_one_end_r
  subroutine gridpoint_cluster_middle(alpha,z0,z)
    !
    ! clustered in the middle
    !
    implicit none
    real(rp), intent(in) :: alpha,z0
    real(rp), intent(out) :: z
    if(alpha /= 0.) then
      if(     z0 <= 0.5) then
        z = 0.5*(1.-1.+tanh(2.*alpha*(z0-0.))/tanh(alpha))
        !z = 0.5*(1.-1.+erf( 2.*alpha*(z0-0.))/erf( alpha))
      else if(z0  > 0.5) then
        z = 0.5*(1.+1.+tanh(2.*alpha*(z0-1.))/tanh(alpha))
        !z = 0.5*(1.+1.+erf( 2.*alpha*(z0-1.))/erf( alpha))
      end if
    else
      z = z0
    end if
  end subroutine gridpoint_cluster_middle
  subroutine gridpoint_natural(kg,nzg,z,kb_a,alpha_a,c_eta_a,dyp_a)
    !
    ! a physics-based, 'natural' grid stretching function for wall-bounded turbulence
    ! see Pirozzoli & Orlandi, JCP 439 - 110408 (2021)
    !
    ! clustered at the two sides
    !
    implicit none
    real(dp), parameter :: kb_p     = 32._dp,    &
                           alpha_p  = pi/1.5_dp, &
                           c_eta_p  = 0.8_dp,    &
                           dyp_p    = 0.05_dp
    integer , intent(in ) :: kg,nzg
    real(dp), intent(out) :: z
    real(dp), intent(in ), optional :: kb_a,alpha_a,c_eta_a,dyp_a
    real(dp)                        :: kb  ,alpha  ,c_eta  ,dyp
    real(dp) :: retau,n,k
    !
    ! handle input parameters
    !
    kb    = kb_p   ; if(present(kb_a   )) kb    = kb_a
    alpha = alpha_p; if(present(alpha_a)) alpha = alpha_a
    c_eta = c_eta_p; if(present(c_eta_a)) c_eta = c_eta_a
    dyp   = dyp_p  ; if(present(dyp_a  )) dyp   = dyp_a
    !
    ! determine retau
    !
    n = nzg/2._dp
    retau = 1._dp/(1._dp+(n/kb)**2)*(dyp*n+(3._dp/4._dp*alpha*c_eta*n)**(4._dp/3._dp)*(n/kb)**2)
#if defined(_DEBUG)
    if(kg==1) print*,'Grid targeting Retau = ',retau
#endif
    k = 1._dp*min(kg,(nzg-kg))
    !
    ! dermine z/(2h)
    !
    z = 1._dp/(1._dp+(k/kb)**2)*(dyp*k+(3._dp/4._dp*alpha*c_eta*k)**(4._dp/3._dp)*(k/kb)**2)/(2._dp*retau)
    if( kg >= nzg-kg ) z = 1._dp-z
  end subroutine gridpoint_natural
end module mod_initgrid
