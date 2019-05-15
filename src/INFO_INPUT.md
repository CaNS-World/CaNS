# about `dns.in`

Consider the following input file as example (corresponds to a turbulent plane channel flow with friction Reynolds number = 180):

~~~
64 64 64                 ! itot, jtot, ktot
3. 1.5 1.                ! lx, ly, lz
0.                       ! gr
.95                      ! cfl
1. 1. 5600.              ! uref, lref, visc
poi                      ! inivel
T                        ! is_wallturb
10000                    ! nstep
F                        ! restart
10 10 20 5000 10000 2000 ! icheck, iout0d, iout1d, iout2d, iout3d, isave
P P  P P  D D            ! cbcvel(0:1,1:3,1); u BC type
P P  P P  D D            ! cbcvel(0:1,1:3,2); v BC type
P P  P P  D D            ! cbcvel(0:1,1:3,3); w BC type
P P  P P  N N            ! cbcpre(0:1,1:3,1); p BC type
0. 0.  0. 0.  0. 0.      ! cbcvel(0:1,1:3,1); u BC type
0. 0.  0. 0.  0. 0.      ! cbcvel(0:1,1:3,2); v BC type
0. 0.  0. 0.  0. 0.      ! cbcvel(0:1,1:3,3); w BC type
0. 0.  0. 0.  0. 0.      ! cbcpre(0:1,1:3,1); p BC type
T F F                    ! is_forced(1:3)
1. 0. 0.                 ! velf(1:3)
F F  F F  F F            ! is_outflow(0:1,1:3)
2 2                      ! dims(1:2)
4                        ! numthreadsmax
~~~

---
---

~~~
64 64 64                 ! itot, jtot, ktot
3. 1.5 1.                ! lx, ly, lz
0.                       ! gr
~~~

These lines set the computational grid.

`itot, jtot, ktot ` (integer) and `lx, ly, lz` (real) are the number of points  and domain size (real) in each direction.

`gr` (real) is a grid stretching parameter that tweaks the non-uniform grid in the third direction; zero `gr` implies no stretching. fFor more details see `initgrid.f90`.

-

~~~
.95                      ! cfl
~~~

this line controls the simulation time step, set to be equal to `cfl` (real) times the maximum allowable timestep for the low-storage RK3 temporal integration.

-

~~~
1. 1. 5640.              ! uref, lref, rey
~~~

this line sets the.

-

~~~
log                      ! inivel
T                        ! is_wallturb
~~~

these lines set up the initial velocity field. The following options for variable `initvel` (character) are available:

* `zer`: zero velocity field
* `cou`: plane Couette flow profile with symmetric wall velocities equal to `uref/2`; streamwise direction in `x`
* `poi`: plane Poiseuille flow profile with mean velocity `uref`                    ; streamwise direction in `x`
* `log`: logarithmic profile with mean velocity `uref`                              ; streamwise direction in `x`
* `hcp`: half channel with plane poiseuille profile and mean velocity `uref`        ; streamwise direction in `x`
* `hcp`: half channel with logarithmic profile and mean velocity `uref`             ; streamwise direction in `x`
* `tgv`: Taylor-Green vortex

`is_wallturb` (logical) is a boolean variable that, if `.true.`, superimposes to the velocity field with an high amplitude disturbance that effectively triggers transition to turbulence in a wall-bounded shear flow; see `initflow.f90` for more details.

-

~~~
P P  P P  D D          ! cbcvel(0:1,1:3,1); u BC type
P P  P P  D D          ! cbcvel(0:1,1:3,2); v BC type
P P  P P  D D          ! cbcvel(0:1,1:3,3); w BC type
P P  P P  N N          ! cbcpre(0:1,1:3,1); p BC type
0. 0.  0. 0.  0. 0.    ! cbcvel(0:1,1:3,1); u BC value
0. 0.  0. 0.  0. 0.    ! cbcvel(0:1,1:3,2); v BC value
0. 0.  0. 0.  0. 0.    ! cbcvel(0:1,1:3,3); w BC value
0. 0.  0. 0.  0. 0.    ! cbcpre(0:1,1:3,1); p BC value
~~~

The *type* of boundary conditions (BC) for each field variable are set by a row of six characters: `X0 X1  Y0 Y1  Z0 Z1`

where `X0 X1` set the type of BC the field variable for the lower and upper boundaries in `x`; ditto for `Y0 Y1` in `y` and `Z0 Z1` in `z`. 

The four rows correspond to the `x`, `y` and `z` velocity components, and pressure, i.e. `u`, `v`, `w`, and `p`.

The following options are available:

* `P` periodic
* `D` Dirichlet
* `N` Neumann

The last for rows follow the same logic, but now for the BC *values* (reals; dummy for a periodic direction).

-
~~~
T F F                    ! is_forced(1:3)
1. 0. 0.                 ! velf(1:3)
F F  F F  F F            ! is_outflow(0:1,1:3)
~~~

-

~~~
10000                    ! nstep
F                        ! restart
~~~

-

~~~
10 10 20 5000 10000 2000 ! icheck, iout0d, iout1d, iout2d, iout3d, isave
~~~


* every `icheck` time steps compute the new time step size dt according to the new stability criterion and cfl (above)
* every `iout0d` time steps update the history files with global scalar variables; currently the forcing pressure gradient and time step history are reported
* every `iout1d` time steps write 1d profiles (velocity and its moments) to a file
* every `iout2d` time steps write a 2d slice of a 3d scalar field to a file
* every `iout3d` time steps write a 3d scalar field into a file
* every `isave`  time steps write a checkpoint file

 note: in order to edit the outputs based on the specific flow case, edit
 main.f90 or even output.f90 accordingly. Currently these assume the z to be
 an inhomogeneous direction.

-

~~~
2 2                      ! dims(1:2)
4                        ! numthreadsmax
~~~

grid of computational subdomains
