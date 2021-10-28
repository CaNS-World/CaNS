module mod_initgrid
  use mod_param, only:pi
  use mod_types
  implicit none
  private
  public initgrid
  contains
  subroutine initgrid(inivel,n,gr,lz,dzc,dzf,zc,zf)
    !
    ! initializes the non-uniform grid in z
    !
    implicit none
    character(len=3), intent(in) :: inivel
    integer , intent(in) :: n
    real(rp), intent(in) :: gr,lz
    real(rp), intent(out), dimension(0:n+1) :: dzc,dzf,zc,zf
    real(rp) :: z0
    integer :: k
    procedure (), pointer :: gridpoint => null()
    select case(inivel)
    case('zer','log','poi','cou')
      gridpoint => gridpoint_cluster_two_end
    case('hcl','hcp')
      gridpoint => gridpoint_cluster_one_end
    case default
      gridpoint => gridpoint_cluster_two_end
    end select
    !
    ! step 1) determine coordinates of cell faces zf
    !
    do k=1,n
      z0  = (k-0.)/(1.*n)
      call gridpoint(gr,z0,zf(k))
      zf(k) = zf(k)*lz
    end do
    zf(0) = 0.
    !
    ! step 2) determine grid spacing between faces dzf
    !
    do k=1,n
      dzf(k) = zf(k)-zf(k-1)
    end do
    dzf(0  ) = dzf(1)
    dzf(n+1) = dzf(n)
    !
    ! step 3) determine grid spacing between centers dzc
    !
    do k=0,n
      dzc(k) = .5*(dzf(k)+dzf(k+1))
    end do
    dzc(n+1) = dzc(n)
    !
    ! step 4) compute coordinates of cell centers zc and faces zf
    !
    zc(0)    = -dzc(0)/2.
    zf(0)    = 0.
    do k=1,n+1
      zc(k) = zc(k-1) + dzc(k-1)
      zf(k) = zf(k-1) + dzf(k)
    end do
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
    if(alpha.ne.0.) then
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
    if(alpha.ne.0.) then
      z = 1.0*(1.+tanh((z0-1.0)*alpha)/tanh(alpha/1.))
      !z = 1.0*(1.+erf( (z0-1.0)*alpha)/erf( alpha/1.))
    else
      z = z0
    end if
  end subroutine gridpoint_cluster_one_end
  subroutine gridpoint_cluster_middle(alpha,z0,z)
    !
    ! clustered in the middle
    !
    implicit none
    real(rp), intent(in) :: alpha,z0
    real(rp), intent(out) :: z
    if(alpha.ne.0.) then
      if(    z0 <= 0.5) then
        z = 0.5*(1.-1.+tanh(2.*alpha*(z0-0.))/tanh(alpha))
        !z = 0.5*(1.-1.+erf( 2.*alpha*(z0-0.))/erf( alpha))
      elseif(z0 > 0.5) then
        z = 0.5*(1.+1.+tanh(2.*alpha*(z0-1.))/tanh(alpha))
        !z = 0.5*(1.+1.+erf( 2.*alpha*(z0-1.))/erf( alpha))
      end if
    else
      z = z0
    end if
  end subroutine gridpoint_cluster_middle
end module mod_initgrid
