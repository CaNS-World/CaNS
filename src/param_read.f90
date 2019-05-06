module mod_param
implicit none
public
   real(8) :: pi = acos(-1.d0)
   real(8) :: small = 1.d-9
   logical, dimension(2,3) :: no_outflow = & 
       reshape((/.false.,.false.,   & ! no outflow in x lower,upper bound
                 .false.,.false.,   & ! no outflow in y lower,upper bound
                 .false.,.false./), & ! no outflow in z lower,upper bound
                 shape(no_outflow))
   character(len=100) :: datadir = 'data/'
   real(8), dimension(2,3) :: rkcoeff = reshape((/ 32.d0/60.d0,  0.d0       , &
                                                   25.d0/60.d0, -17.d0/60.d0, &
                                                   45.d0/60.d0, -25.d0/60.d0/), shape(rkcoeff))
   real(8), dimension(3)   :: rkcoeff12 = rkcoeff(1,:)+rkcoeff(2,:)
   integer :: itot,jtot,ktot
   real(8) :: lx, ly, lz
   real(8) :: gr     = 0.0d0 ! grid stretching parameter (0. -> no stretching;
                             !                            see initgrid.f90)
   real(8) :: cfl
   real(8) :: visc
   !
   character(len=*) :: inivel = 'log'
   !
   ! is_wallturb: if .true., generates an initial condition that triggers 
   ! transition to turbulence in wall-bounded flows (assumes streamwise direction x)
   !
   logical :: is_wallturb = .false.
   !
   integer :: nstep = 100000 ! number of time steps
   logical :: restart = .false. ! restart or not from a checkpoint
   !
   ! -> every *icheck* time steps compute the new time step size dt
   ! according to the new stability criterion and cfl (above)
   ! -> every *iout0d* time steps update the history files with global scalar variables;
   ! currently the forcing pressure gradient and time step history are reported
   ! -> every *iout1d* time steps write 1d profiles (velocity and its moments)
   ! to a file
   ! -> every *iout2d* time steps write a 2d slice of a 3d scalar field to a file
   ! -> every *iout3d* time steps write a 3d scalar field into a file
   ! -> every *isave*  time steps write a checkpoint file
   ! note: in order to edit the outputs based on the specific flow case, edit 
   ! main.f90 or even output.f90 accordingly. Currently these assume the z to be
   ! an inhomogeneous direction.
   !
   integer :: icheck, iout0d, iout1d,iout2d,iout3d,isave
   ! grid of computational subdomains
   integer, dimension(2) :: dims = (/2,2/)
   ! x and y sizes of local arrays in the basic 2D z-pencil decomposition
   integer :: imax = itot/dims(1), jmax = jtot/dims(2)
   !
   ! number of OpenMP threads
   !
   integer :: nthreadsmax = 4
   !--------------------------------------------------------------
   ! boundary conditions
   ! P -> periodic, D -> Dirichlet, N -> Neumann
   !--------------------------------------------------------------
   character(len=1), dimension(0:1,3,3) ::  cbcvel
   real(8)         , dimension(0:1,3,3) :: bcvel
   character(len=1), dimension(0:1,3)   ::  cbcpre
   real(8)         , dimension(0:1,3)   ::   bcpre
   !-------------------------------------------------------------------
   ! forcing the flow with a pressure gradient
   ! that balances the total wall shear 
   ! (e.g. for a pressure-driven channel) 
   !
   logical, dimension(3) :: is_forced
   !
   ! desired values of bulk velocities 
   ! (only relevant if the corresponding boolean
   !  above is .true.)
   !
   real(8), dimension(3) :: velf
   !
   ! outflow boundary condition
   ! if is_outflow(i,j) is true, an outflow condition is prescribed for the
   ! face-centered velocity at that boundary
   ! the outflow BC is determined from the condition of zero divergence
   !
   logical, dimension(0:1,3) :: is_outflow
end module mod_param
