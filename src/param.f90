module mod_param
implicit none
public
!
! parameters
!
real(8), parameter :: pi = acos(-1.d0)
real(8), parameter :: small = 1.d-9
logical, parameter, dimension(2,3) :: no_outflow = & 
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
! variables to be determined from the input file 'dns.in'
!
integer :: itot,jtot,ktot,imax,jmax
real(8) :: lx,ly,lz,dx,dy,dz,dxi,dyi,dzi,gr
real(8) :: cfl
real(8) :: visc
!
character(len=100) :: inivel ! DON'T FORGET TO ADD A TRIM IN THE SWITCH COMMAND!
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
integer, dimension(3) :: ng
integer, dimension(3) :: n
real(8), dimension(3) :: l
real(8), dimension(3) :: dl
real(8), dimension(3) :: dli
!
contains 
  subroutine read_input(myid)
  use mpi
  implicit none
  integer, intent(in) :: myid
  integer :: iunit,ierr
    open(newunit=iunit,file='dns.in',status='old',action='read',iostat=ierr)
      if( ierr.eq.0 ) then
        read(iunit,*) itot,jtot,ktot
        read(iunit,*) lx,ly,lz
        read(iunit,*) gr
        read(iunit,*) cfl
        read(iunit,*) visc
        read(iunit,*) inivel
        read(iunit,*) is_wallturb
        read(iunit,*) nstep
        read(iunit,*) restart
        read(iunit,*) icheck,iout0d,iout1d,iout2d,iout3d,isave
        read(iunit,*) dims(1),dims(2)
        read(iunit,*) nthreadsmax
        read(iunit,*) cbcvel(0,1,1),cbcvel(1,1,1),cbcvel(0,2,1),cbcvel(1,2,1),cbcvel(0,3,1),cbcvel(1,3,1)
        read(iunit,*) cbcvel(0,1,2),cbcvel(1,1,2),cbcvel(0,2,2),cbcvel(1,2,2),cbcvel(0,3,2),cbcvel(1,3,2)
        read(iunit,*) cbcvel(0,1,3),cbcvel(1,1,3),cbcvel(0,2,3),cbcvel(1,2,3),cbcvel(0,3,3),cbcvel(1,3,3)
        read(iunit,*) cbcpre(0,1  ),cbcpre(1,1  ),cbcpre(0,2  ),cbcpre(1,2  ),cbcpre(0,3  ),cbcpre(1,3  )
        read(iunit,*)  bcvel(0,1,1), bcvel(1,1,1), bcvel(0,2,1), bcvel(1,2,1), bcvel(0,3,1), bcvel(1,3,1)
        read(iunit,*)  bcvel(0,1,2), bcvel(1,1,2), bcvel(0,2,2), bcvel(1,2,2), bcvel(0,3,2), bcvel(1,3,2)
        read(iunit,*)  bcvel(0,1,3), bcvel(1,1,3), bcvel(0,2,3), bcvel(1,2,3), bcvel(0,3,3), bcvel(1,3,3)
        read(iunit,*)  bcpre(0,1  ), bcpre(1,1  ), bcpre(0,2  ), bcpre(1,2  ), bcpre(0,3  ), bcpre(1,3  )
        read(iunit,*)  is_forced(1),is_forced(2),is_forced(3)
        read(iunit,*)  velf(1),velf(2),velf(3)
        read(iunit,*)  is_outflow(0,1),is_outflow(1,1),is_outflow(0,2),is_outflow(1,2),is_outflow(0,3),is_outflow(1,3)
      else
        if(myid.eq.0) print*, 'Error reading the input file' 
        if(myid.eq.0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        call exit
    endif
    close(iunit)
    dx = lx/(1.d0*itot)
    dy = ly/(1.d0*jtot)
    dz = lz/(1.d0*ktot)
    dxi = dx**(-1)
    dyi = dy**(-1)
    dzi = dz**(-1)
    imax = itot/dims(1)
    jmax = jtot/dims(2)
    !
    ng  = (/itot,jtot,ktot/)
    n   = (/imax,jmax,ktot/)
    l   = (/lx,ly,lz/)
    dl  = (/dx,dy,dz/)
    dli = (/dxi,dyi,dzi/)
  return
  end subroutine read_input
end module mod_param
