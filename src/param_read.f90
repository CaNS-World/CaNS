module mod_param
implicit none
private read_params
public
real(8), parameter ::  :: pi = acos(-1.d0)
real(8), parameter :: small = 1.d-9
logical, parameter :: dimension(2,3) :: no_outflow = & 
    reshape((/.false.,.false.,   & ! no outflow in x lower,upper bound
              .false.,.false.,   & ! no outflow in y lower,upper bound
              .false.,.false./), & ! no outflow in z lower,upper bound
              shape(no_outflow))
character(len=100), parameter :: datadir = 'data/'
real(8), parameter, dimension(2,3) :: rkcoeff = reshape((/ 32.d0/60.d0,  0.d0       , &
                                                   25.d0/60.d0, -17.d0/60.d0, &
                                                   45.d0/60.d0, -25.d0/60.d0/), shape(rkcoeff))
real(8), parameter, dimension(3)   :: rkcoeff12 = rkcoeff(1,:)+rkcoeff(2,:)
!
integer :: itot,jtot,ktot,imax,jmax
real(8) :: lx,ly,lz,dx,dy,dz,dxi,dyi,dzi,gr
real(8) :: cfl
real(8) :: visc
!
character(len=*) :: inivel
logical :: is_wallturb
!
integer :: nstep
logical :: restart
integer :: icheck,iout0d,iout1d,iout2d,iout3d,isave
!
integer, dimension(2) :: dims = (/2,2/)
integer :: nthreadsmax = 4
!
character(len=1), dimension(0:1,3,3) ::  cbcvel
real(8)         , dimension(0:1,3,3) :: bcvel
character(len=1), dimension(0:1,3)   ::  cbcpre
real(8)         , dimension(0:1,3)   ::   bcpre
!
logical, dimension(3) :: is_forced
real(8), dimension(3) :: velf
logical, dimension(0:1,3) :: is_outflow
!
contains 
  subroutine read_input
  implicit none
    open(unit=99,file='dns.in')
      read(unit=99,*) itot,jtot,ktot
      read(unit=99,*) lx,ly,lz
      read(unit=99,*) gr
      read(unit=99,*) cfl
      read(unit=99,*) inivel
      read(unit=99,*) is_wallturb
      read(unit=99,*) nstep
      read(unit=99,*) restart
      read(unit=99,*) icheck,iout0d,iout1d,iout2d,iout3d,isave
      read(unit=99,*) dims(1),dims(2)
      read(unit=99,*) nthreadsmax
      read(unit=99,*) cbcvel(0,1,1),cbcvel(1,1,1),cbcvel(0,2,1),cbcvel(1,2,1),cbcvel(0,3,1),cbcvel(1,3,1)
      read(unit=99,*) cbcvel(0,1,2),cbcvel(1,1,2),cbcvel(0,2,2),cbcvel(1,2,2),cbcvel(0,3,2),cbcvel(1,3,2)
      read(unit=99,*) cbcvel(0,1,3),cbcvel(1,1,3),cbcvel(0,2,3),cbcvel(1,2,3),cbcvel(0,3,3),cbcvel(1,3,3)
      read(unit=99,*) cbcpre(0,1  ),cbcpre(1,1  ),cbcpre(0,2  ),cbcpre(1,2  ),cbcpre(0,3  ),cbcpre(1,3  )
      read(unit=99,*)  bcvel(0,1,1), bcvel(1,1,1), bcvel(0,2,1), bcvel(1,2,1), bcvel(0,3,1), bcvel(1,3,1)
      read(unit=99,*)  bcvel(0,1,2), bcvel(1,1,2), bcvel(0,2,2), bcvel(1,2,2), bcvel(0,3,2), bcvel(1,3,2)
      read(unit=99,*)  bcvel(0,1,3), bcvel(1,1,3), bcvel(0,2,3), bcvel(1,2,3), bcvel(0,3,3), bcvel(1,3,3)
      read(unit=99,*)  bcpre(0,1  ), bcpre(1,1  ), bcpre(0,2  ), bcpre(1,2  ), bcpre(0,3  ), bcpre(1,3  )
      read(unit=99,*)  is_forced(1),is_forced(2),is_forced(3)
      read(unit=99,*)  velf(1),velf(2),velf(3)
      read(unit=99,*)  is_outflow(0,1),is_outflow(1,1),is_outflow(0,2),is_outflow(1,2),is_outflow(0,3),is_outflow(1,3)
    close(99)
    dx = lx/(1.d0*itot), &
    dy = ly/(1.d0*jtot), &
    dz = lz/(1.d0*ktot), &
    dxi    = dx**(-1),   &
    dyi    = dy**(-1),   &
    dzi    = dz**(-1),   &
    imax = itot/dims(1)
    jmax = jtot/dims(2)
  return
  end subroutine read_params
end module mod_param
