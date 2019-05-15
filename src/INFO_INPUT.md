## input file `dns.in`
Hereby the input file dns.in is explained. Consider the following file as example (corresponds to a turbulent plane channel flow with $Re_\tau = 180$):

~~~
64 64 64                 ! itot, jtot, ktot
3. 1.5 1.                ! lx, ly, lz
0.                       ! gr
.95                      ! cfl
0.00001                  ! visc
log                      ! inivel
F                        ! is_wallturb
10000                    ! nstep
F                        ! restart
10 10 20 5000 10000 2000 ! iout0d, iout1d, iout2d, iout3d, isave
2 2                      ! dims(1:2)
4                        ! numthreadsmax
P P  P P  D D            ! cbcvel(0:1,1:3,1)
P P  P P  D D            ! cbcvel(0:1,1:3,2)
P P  P P  D D            ! cbcvel(0:1,1:3,3)
P P  P P  N N            ! cbcpre(0:1,1:3,1)
0. 0.  0. 0.  0. 0.      ! cbcvel(0:1,1:3,1)
0. 0.  0. 0.  0. 0.      ! cbcvel(0:1,1:3,2)
0. 0.  0. 0.  0. 0.      ! cbcvel(0:1,1:3,3)
0. 0.  0. 0.  0. 0.      ! cbcpre(0:1,1:3,1)
T F F                    ! is_forced(1:3)
1. 0. 0.                 ! velf(1:3)
F F  F F  F F            ! is_outflow(0:1,1:3)
~~~

The following lines set up the computational grid:

~~~
64 64 64                 ! itot, jtot, ktot
3. 1.5 1.                ! lx, ly, lz
0.                       ! gr
~~~

where `itot, jtot, ktot ` and `lx, ly, lz` are the number of points (integer) and domain size (real) in each direction. `gr` is a grid stretching parameter (real) that tweaks the non-uniform grid in the third direction; `gr` equal to zero implies no stretching. For more details see `initgrid.f90`.

~~~
.95                      ! cfl
~~~

this line controls the simulation time step, set to be equal to `cfl` times the maximum allowable timestep for the low-storage RK3 temporal integration.

~~~
0.00001                  ! visc
~~~

this line sets the fluid kinematic viscosity.

~~~
log                      ! inivel
~~~

this line stets up the initial velocity field. The following options are available:

(add table here).