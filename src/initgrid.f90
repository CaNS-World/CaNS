module mod_initgrid
  use mod_param, only:pi
  implicit none
  private
  public initgrid
  contains
  subroutine initgrid(inivel,n,dzmin,gr,lz,dzc,dzf,zc,zf) ! ADD LZ
    implicit none
    character(len=3) :: inivel
    integer, intent(in) :: n
    real(8), intent(in) :: dzmin,gr,lz 
    real(8), intent(out), dimension(0:n+1) :: dzc,dzf,zc,zf
    real(8) :: z,z0,period
    integer :: k
    real(8) :: norm
    select case(inivel)
    case('zer','log','poi','cou') ! stretch till half of the channel
      z0 = .5d0
      period = 1.d0
    case('hcl','hcp')
    !  z0 = 1.d0
    !  period = 2.d0
    !  not working yet
      z0 = .5d0
      period = 1.d0
    case default
      z0 = .5d0
      period = 1.d0
    end select
    do k=1,n
      z  = (k-.5d0)/(1.d0*n)
      dzf(k) = gridpoint(gr,z0,period,dzmin,z)
    enddo
    dzf(0)   =  dzf(1)
    dzf(n+1) =  dzf(n)
    !
    do k=0,n
      dzc(k) = .5d0*(dzf(k+1)+dzf(k))
    enddo
    dzc(n+1) = dzc(n)
    zc(0) = -dzc(0)/2.d0
    zf(0) = 0.
    do k=0,n
      zc(k+1) = zc(k) + dzc(k)
      zf(k+1) = zf(k) + dzf(k)
    enddo
    norm = zf(n)
    dzc(:) = dzc(:)/norm*lz
    dzf(:) = dzf(:)/norm*lz
    zc(:)  =  zc(:)/norm*lz
    zf(:)  =  zf(:)/norm*lz
    return
  end subroutine initgrid
  !
  real(8) function gridpoint(alpha,z0,period,dz,z)
    implicit none
    real(8), intent(in) :: alpha,z0,period,dz,z
    gridpoint = dz*( 1.d0+.5d0*(alpha-1.d0)*( 1.d0+cos( 2.d0*pi*(z-z0)/period ) ) )
    return
  end function gridpoint
end module mod_initgrid
