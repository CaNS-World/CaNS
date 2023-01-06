# about `dns.in`

Consider the following input file as example (corresponds to a turbulent plane channel flow):

~~~
512 256 144              ! itot, jtot, ktot
6. 3. 1.                 ! lx, ly, lz
0.                       ! gr
.95 1.0e5                ! cfl, dtmin
5640.                    ! visci
poi                      ! inivel
T                        ! is_wallturb
100000 100. 0.1          ! nstep, time_max, tw_max
T F F                    ! stop_type(1:3)
F T 0                    ! restart,is_overwrite_save,nsaves_max
10 10 100 500 10000 5000 ! icheck, iout0d, iout1d, iout2d, iout3d, isave
P P  P P  D D            ! cbcvel(0:1,1:3,1) [u BC type]
P P  P P  D D            ! cbcvel(0:1,1:3,2) [v BC type]
P P  P P  D D            ! cbcvel(0:1,1:3,3) [w BC type]
P P  P P  N N            ! cbcpre(0:1,1:3  ) [p BC type]
0. 0.  0. 0.  0. 0.      !  bcvel(0:1,1:3,1) [u BC value]
0. 0.  0. 0.  0. 0.      !  bcvel(0:1,1:3,2) [v BC value]
0. 0.  0. 0.  0. 0.      !  bcvel(0:1,1:3,3) [w BC value]
0. 0.  0. 0.  0. 0.      !  bcpre(0:1,1:3  ) [p BC value]
0. 0. 0.                 ! bforce(1:3)
T F F                    ! is_forced(1:3)
1. 0. 0.                 ! velf(1:3)
2 2                      ! dims(1:2)
~~~

---
---

~~~
512 256 144              ! itot, jtot, ktot
6. 3. 1.                 ! lx, ly, lz
0.                       ! gr
~~~

These lines set the computational grid.

`itot, jtot, ktot ` and `lx, ly, lz` are the **number of points**  and **domain length** in each direction.

`gr` is a **grid stretching parameter** that tweaks the non-uniform grid in the third direction; zero `gr` implies no stretching. See `initgrid.f90` for more details.

---

~~~
.95 1.0e5                ! cfl, dtmin
~~~

This line controls the simulation time step.

The time step is set to be equal to `min(cfl*dtmax,dtmin)`, i.e. the minimum value between `dtmin` and `cfl` times the maximum allowable time step `dtmax` (computed every `ickeck` time steps; see below).
`dtmin` is therefore used when a constant time step, smaller than `cfl*dtmax`, is required. If not, it should be set to a high value so that the time step is dynamically adjusted to `cfl*dtmax`.

---

~~~
5640.                    ! visci
~~~

This line defines the inverse of the fluid viscosity, `visci`, meaning that the viscosity is `visc = visci**(-1)`. Note that, for a setup defined with unit reference length and velocity scales, `visci` has the same value as the flow Reynolds number.

---

~~~
poi                      ! inivel
T                        ! is_wallturb
~~~

These lines set the initial velocity field.

`initvel` **chooses the initial velocity field**. The following options are available:

* `zer`: zero velocity field
* `uni`: uniform velocity field equal to `uref`                                     ; streamwise direction in `x`
* `cou`: plane Couette flow profile with symmetric wall velocities equal to `uref/2`; streamwise direction in `x`
* `poi`: plane Poiseuille flow profile with mean velocity `uref`                    ; streamwise direction in `x`
* `tbl`: temporal boundary layer profile with wall velocity `uref`                  ; streamwise direction in `x`
* `log`: logarithmic profile with mean velocity `uref`                              ; streamwise direction in `x`
* `hcp`: half channel with plane Poiseuille profile and mean velocity `uref`        ; streamwise direction in `x`
* `hcl`: half channel with logarithmic profile and mean velocity `uref`             ; streamwise direction in `x`
* `tgv`: three-dimensional Taylor-Green vortex
* `tgw`: two-dimensional   Taylor-Green vortex

`is_wallturb`, if true, **superimposes a high amplitude disturbance on the initial velocity field** that effectively triggers transition to turbulence in a wall-bounded shear flow.

See `initflow.f90` for more details.

---

~~~
100000 100. 0.1          ! nstep, time_max, tw_max
T F F                    ! stop_type(1:3)
F T 0                    ! restart,is_overwrite_save,nsaves_max
~~~

These lines set the simulation termination criteria and whether the simulation should be restarted from a checkpoint file.

`nstep` is the **total number of time steps**.

`time_max` is the **maximum physical time**.

`tw_max` is the **maximum total simulation wall-clock time**.

`stop_type` sets which criteria for terminating the simulation are to be used (more than one can be selected, and at least one of them must be `T`)

* `stop_type(1)`, if true (`T`), the simulation will terminate after `nstep` time steps have been simulated;
* `stop_type(2)`, if true (`T`), the simulation will terminate after `time_max` physical time units have been reached;
* `stop_type(3)`, if true (`T`), the simulation will terminate after `tw_max` simulation wall-clock time (in hours) has been reached;

a checkpoint file `fld.bin` will be saved before the simulation is terminated.

`restart`, if true, **restarts the simulation** from a previously saved checkpoint file, named `fld.bin`.

`is_overwrite_save`, if true, overwrites the checkpoint file `fld.bin` at every save; if false, a symbolic link is created which makes `fld.bin` point to the last checkpoint file with name `fld_???????.bin` (with `???????` denoting the corresponding time step number). In the latter case, to restart a run from a different checkpoint one just has to point the file `fld.bin` to the right file, e.g.: ` ln -sf fld_0000100.bin fld.bin`.

`nsaves_max` limits the number of saved checkpoint files, if `is_over_write_save` is false; a value of `0` or any negative integer corresponds to no limit, and the code uses the file format described above; otherwise, only `nsaves_max` checkpoint files are saved, with the oldest save being overwritten when the number of saved checkpoints exceeds this threshold; in this case, files with a format `fld_????.bin` are saved (with `????` denoting the saved file number), with `fld.bin` pointing to the last checkpoint file as described above; moreover, a file `log_checkpoints.out` records information about the time step number and physical time corresponding to each saved file number.

---

~~~
10 10 20 5000 10000 2000 ! icheck, iout0d, iout1d, iout2d, iout3d, isave
~~~

These lines set the frequency of time step checking and output:

* every `icheck` time steps **the new time step size** is computed according to the new stability criterion and cfl (above);
* every `iout0d` time steps **history files with global scalar variables** are appended; currently the forcing pressure gradient and time step history are reported;
* every `iout1d` time steps **1d profiles** are written (e.g. velocity and its moments) to a file;
* every `iout2d` time steps **2d slices of a 3d scalar field** are written to a file;
* every `iout3d` time steps **3d scalar fields** are written to a file;
* every `isave`  time steps a **checkpoint file** is written (`fld_???????.bin`), and a symbolic link for the restart file, `fld.bin`, will point to this last save so that, by default, the last saved checkpoint file is used to restart the simulation.

1d, 2d and 3d outputs can be tweaked modifying files `out?d.h90`, and re-compiling the source. See also `output.f90` for more details.

---

~~~
P P  P P  D D          ! cbcvel(0:1,1:3,1) [u BC type]
P P  P P  D D          ! cbcvel(0:1,1:3,2) [v BC type]
P P  P P  D D          ! cbcvel(0:1,1:3,3) [w BC type]
P P  P P  N N          ! cbcpre(0:1,1:3  ) [p BC type]
0. 0.  0. 0.  0. 0.    !  bcvel(0:1,1:3,1) [u BC value]
0. 0.  0. 0.  0. 0.    !  bcvel(0:1,1:3,2) [v BC value]
0. 0.  0. 0.  0. 0.    !  bcvel(0:1,1:3,3) [w BC value]
0. 0.  0. 0.  0. 0.    !  bcpre(0:1,1:3  ) [p BC value]
~~~

These lines set the boundary conditions (BC).

The **type** (BC) for each field variable are set by a row of six characters, `X0 X1  Y0 Y1  Z0 Z1` where,

* `X0` `X1` set the type of BC the field variable for the **lower** and **upper** boundaries in `x`;
* `Y0` `Y1` set the type of BC the field variable for the **lower** and **upper** boundaries in `y`;
* `Z0` `Z1` set the type of BC the field variable for the **lower** and **upper** boundaries in `z`.

The four rows correspond to the three velocity components, and pressure, i.e. `u`, `v`, `w`, and `p`.

The following options are available:

* `P` periodic;
* `D` Dirichlet;
* `N` Neumann.

The **last four rows** follow the same logic, but now for the BC **values** (dummy for a periodic direction).

---

~~~
0. 0. 0.                 ! bforce(1:3)
T F F                    ! is_forced(1:3)
1. 0. 0.                 ! velf(1:3)
~~~

These lines set the flow forcing.

`bforce`, is a constant **body force density term** in the direction in question (e.g., the negative of a constant pressure gradient) that can be added to the right-hand-side of the momentum equation. The three values correspond to three domain directions.

`is_forced`, if true in the direction in question, **forces the flow** with a pressure gradient that balances the total wall shear (e.g., for a pressure-driven channel). The three boolean values correspond to three domain directions.

`velf`, is the **target bulk velocity** in the direction in question (where `is_forced` is true). The three values correspond to three domain directions.

---

~~~
2 2                      ! dims(1:2)
~~~

This line set the grid of computational subdomains.

`dims` is the **number of computational** subdomains in `x` and `y`.

# about `cudecomp.in`

In addition to the `dns.in` file, there is an **optional** input file to set some runtime configurations for the *cuDecomp* library. Consider the following file, which corresponds to the default options in case the file is not provided:

~~~
0 T T ! cudecomp_t_comm_backend,cudecomp_is_t_enable_nccl,cudecomp_is_t_enable_nvshmem
0 T T ! cudecomp_t_comm_backend,cudecomp_is_t_enable_nccl,cudecomp_is_t_enable_nvshmem
~~~

The first line sets the configuration for the transpose communication backend autotuning. Here `cudecomp_t_comm_backend` can be one of:

* `1` -> `CUDECOMP_TRANSPOSE_COMM_MPI_P2P`
* `2` -> `CUDECOMP_TRANSPOSE_COMM_MPI_P2P_PL`
* `3` -> `CUDECOMP_TRANSPOSE_COMM_MPI_A2A`
* `4` -> `CUDECOMP_TRANSPOSE_COMM_NCCL`
* `5` -> `CUDECOMP_TRANSPOSE_COMM_NCCL_PL`
* `6` -> `CUDECOMP_TRANSPOSE_COMM_NVSHMEM`
* `7` -> `CUDECOMP_TRANSPOSE_COMM_NVSHMEM_PL`
* any other value -> enable runtime transpose backend autotuning

The other two boolean values, enable/disable the NCCL (`cudecomp_is_t_enable_nccl`) and NVSHMEM (`cudecomp_is_t_enable_nvshmem`) options for *transpose* communication backend autotuning.

The second line is analogous to the first one, but for halo communication backend autotuning. Here `cudecomp_h_comm_backend` can be one of:

* `1` -> `CUDECOMP_HALO_COMM_MPI`
* `2` -> `CUDECOMP_HALO_COMM_MPI_BLOCKING`
* `3` -> `CUDECOMP_HALO_COMM_NCCL`
* `4` -> `CUDECOMP_HALO_COMM_NVSHMEM`
* `5` -> `CUDECOMP_HALO_COMM_NVSHMEM_BLOCKING`
* any other value -> enable runtime halo backend autotuning

The other two boolean values, enable/disable the NCCL (`cudecomp_is_h_enable_nccl`) and NVSHMEM (`cudecomp_is_h_enable_nvshmem`) options for *halo* communication backend autotuning.

Finally, it is worth recalling that passing `dims(1:2) = [0,0]` under `dns.in` will trigger the *processor grid* autotuning, so there is no need to provide that option in the `cudecomp.in` input file.
