# about the input file `input.nml`

Consider the following input file as example (corresponds to a turbulent plane channel flow). `&dns` defines a so-called Fortran namelist containing all the necessary physical and computational parameters to set a case.


```fortran
&dns
ng(1:3) = 512, 256, 144
l(1:3) = 6., 3., 1.
gtype = 1, gr = 0.
cfl = 0.95, dtmax = 1.e5, dt_f = -1.
visci = 5640.
inivel = 'poi'
is_wallturb = T
nstep = 100000, time_max = 100., tw_max = 0.1
stop_type(1:3) = T, F, F
restart = F, is_overwrite_save = T, nsaves_max = 0
icheck = 10, iout0d = 10, iout1d = 100, iout2d = 500, iout3d = 10000, isave = 5000
cbcvel(0:1,1:3,1) = 'P','P',  'P','P',  'D','D'
cbcvel(0:1,1:3,2) = 'P','P',  'P','P',  'D','D'
cbcvel(0:1,1:3,3) = 'P','P',  'P','P',  'D','D'
cbcpre(0:1,1:3)   = 'P','P',  'P','P',  'N','N'
bcvel(0:1,1:3,1) =  0.,0.,   0.,0.,   0.,0.
bcvel(0:1,1:3,2) =  0.,0.,   0.,0.,   0.,0.
bcvel(0:1,1:3,3) =  0.,0.,   0.,0.,   0.,0.
bcpre(0:1,1:3  ) =  0.,0.,   0.,0.,   0.,0.
bforce(1:3) = 0., 0., 0.
is_forced(1:3) = T, F, F
velf(1:3) = 1., 0., 0.
dims(1:2) = 2, 2, ipecil_axis = 1
/
```
<details>

<summary>Tip for vim/nvim users</summary>
Consider adding the following lines in your `.vimrc` file for syntax highlighting of the namelist file:

```vim
if has("autocmd")
  au BufNewFile,BufRead *.nml set filetype=fortran
  au BufNewFile,BufRead *.namelist set filetype=fortran
endif
```

</details>

---
---

```fortran
ng(1:3) = 512, 256, 144
l(1:3) = 6., 3., 1.
gtype = 1, gr = 0.
```

These lines set the computational grid.

`ng(1:3)` and `l(1:3)` are the **number of points**  and **domain length** in each direction.

`gtype` and `gr` are the **grid stretching type** and **grid stretching parameter** that tweak the non-uniform grid in the third direction; zero `gr` implies no stretching. See `initgrid.f90` for more details. The following options are available for `gtype`:

* `1`: grid clustered towards both ends (default)
* `2`: grid clustered towards the lower end
* `3`: grid clustered towards the upper end
* `4`: grid clustered towards the middle

---

```fortran
cfl = 0.95, dtmax = 1.e5, dt_f = -1.
```

This line controls the simulation time step size.

The time step size is set to be equal to `min(cfl*dt_cfl,dtmax)` if `dt_f < 0`, and to `dt_f` otherwise. In the former case, the code prescribes the minimum value between `dtmax` and `cfl` times the maximum allowable time step `dt_cfl` (computed every `ickeck` time steps; see below). `dtmax` is therefore used when a constant time step, smaller than `cfl*dt_cfl`, is required. If not, it should be set to a high value so that the time step is dynamically adjusted to `cfl*dt_cfl`. Finally, a constant time step size time step may be forced, irrespective of the temporal stability evaluation through `dt_f`.

---

```fortran
visci = 5640.
```

This line defines the inverse of the fluid viscosity, `visci`, meaning that the viscosity is `visc = visci**(-1)`. Note that, for a setup defined with unit reference length and velocity scales, `visci` has the same value as the flow Reynolds number.

---

```fortran
inivel = 'poi'
is_wallturb = T
```

These lines set the initial velocity field.

`initvel` **chooses the initial velocity field**. The following options are available:

* `zer`: zero velocity field
* `uni`: uniform velocity field equal to `uref`                                     ; streamwise direction in `x`
* `cou`: plane Couette flow profile with symmetric wall velocities equal to `uref/2`; streamwise direction in `x`
* `poi`: plane Poiseuille flow profile with mean velocity `uref`                    ; streamwise direction in `x`
* `tbl`: temporal boundary layer profile with wall velocity `uref`                  ; streamwise direction in `x`
* `pdc`: plane Poiseuille flow profile with constant pressure gradient              ; streamwise direction in `x`
* `log`: logarithmic profile with mean velocity `uref`                              ; streamwise direction in `x`
* `hcp`: half channel with plane Poiseuille profile and mean velocity `uref`        ; streamwise direction in `x`
* `hcl`: half channel with logarithmic profile and mean velocity `uref`             ; streamwise direction in `x`
* `hdc`: half plane Poiseuille flow profile with constant pressure gradient         ; streamwise direction in `x`
* `tgv`: three-dimensional Taylor-Green vortex
* `tgw`: two-dimensional   Taylor-Green vortex
* `ant`: three-dimensional Antuono vortex

`is_wallturb`, if true, **superimposes a high amplitude disturbance on the initial velocity field** that effectively triggers transition to turbulence in a wall-bounded shear flow.

See `initflow.f90` for more details.

---

```fortran
nstep = 100000, time_max = 100., tw_max = 0.1
stop_type(1:3) = T, F, F
restart = F, is_overwrite_save = T, nsaves_max = 0
```

These lines set the simulation termination criteria and whether the simulation should be restarted from a checkpoint file.

`nstep` is the **total number of time steps**.

`time_max` is the **maximum physical time**.

`tw_max` is the **maximum total simulation wall-clock time**.

`stop_type` sets which criteria for terminating the simulation are to be used (more than one can be selected, and at least one of them must be `T`)

* `stop_type(1)`, if true (`T`), the simulation will terminate after `nstep` time steps have been simulated
* `stop_type(2)`, if true (`T`), the simulation will terminate after `time_max` physical time units have been reached
* `stop_type(3)`, if true (`T`), the simulation will terminate after `tw_max` simulation wall-clock time (in hours) has been reached

a checkpoint file `fld.bin` will be saved before the simulation is terminated.

`restart`, if true, **restarts the simulation** from a previously saved checkpoint file, named `fld.bin`.

`is_overwrite_save`, if true, overwrites the checkpoint file `fld.bin` at every save; if false, a symbolic link is created which makes `fld.bin` point to the last checkpoint file with name `fld_???????.bin` (with `???????` denoting the corresponding time step number). In the latter case, to restart a run from a different checkpoint one just has to point the file `fld.bin` to the right file, e.g.: ` ln -sf fld_0000100.bin fld.bin`.

`nsaves_max` limits the number of saved checkpoint files, if `is_over_write_save` is false; a value of `0` or any negative integer corresponds to no limit, and the code uses the file format described above; otherwise, only `nsaves_max` checkpoint files are saved, with the oldest save being overwritten when the number of saved checkpoints exceeds this threshold; in this case, files with a format `fld_????.bin` are saved (with `????` denoting the saved file number), with `fld.bin` pointing to the last checkpoint file as described above; moreover, a file `log_checkpoints.out` records information about the time step number and physical time corresponding to each saved file number.

---

```fortran
icheck = 10, iout0d = 10, iout1d = 100, iout2d = 500, iout3d = 10000, isave = 5000
```

These lines set the frequency of time step checking and output:

* every `icheck` time steps **the new time step size** is computed according to the new stability criterion and cfl (above)
* every `iout0d` time steps **history files with global scalar variables** are appended; currently the forcing pressure gradient and time step history are reported
* every `iout1d` time steps **1d profiles** are written (e.g. velocity and its moments) to a file
* every `iout2d` time steps **2d slices of a 3d scalar field** are written to a file
* every `iout3d` time steps **3d scalar fields** are written to a file
* every `isave`  time steps a **checkpoint file** is written (`fld_???????.bin`), and a symbolic link for the restart file, `fld.bin`, will point to this last save so that, by default, the last saved checkpoint file is used to restart the simulation

1d, 2d and 3d outputs can be tweaked modifying files `out?d.h90`, and re-compiling the source. See also `output.f90` for more details. _Set any of these variables to `0` to skip the corresponding operation._

---

```fortran
cbcvel(0:1,1:3,1) = 'P','P',  'P','P',  'D','D'
cbcvel(0:1,1:3,2) = 'P','P',  'P','P',  'D','D'
cbcvel(0:1,1:3,3) = 'P','P',  'P','P',  'D','D'
cbcpre(0:1,1:3)   = 'P','P',  'P','P',  'N','N'
bcvel(0:1,1:3,1) =  0.,0.,   0.,0.,   0.,0.
bcvel(0:1,1:3,2) =  0.,0.,   0.,0.,   0.,0.
bcvel(0:1,1:3,3) =  0.,0.,   0.,0.,   0.,0.
bcpre(0:1,1:3  ) =  0.,0.,   0.,0.,   0.,0.
```

These lines set the boundary conditions (BC).

The **type** (BC) for each field variable are set by a row of six characters, `X0 X1  Y0 Y1  Z0 Z1` where,

* `X0` `X1` set the type of BC the field variable for the **lower** and **upper** boundaries in `x`
* `Y0` `Y1` set the type of BC the field variable for the **lower** and **upper** boundaries in `y`
* `Z0` `Z1` set the type of BC the field variable for the **lower** and **upper** boundaries in `z`

The four rows correspond to the three velocity components, and pressure, i.e. `u`, `v`, `w`, and `p`.

The following options are available:

* `P` periodic
* `D` Dirichlet
* `N` Neumann

The **last four rows** follow the same logic, but now for the BC **values** (dummy for a periodic direction).

---

```fortran
bforce(1:3) = 0., 0., 0.
is_forced(1:3) = T, F, F
velf(1:3) = 1., 0., 0.
```

These lines set the flow forcing.

`bforce`, is a constant **body force density term** in the direction in question (e.g., the negative of a constant pressure gradient) that can be added to the right-hand-side of the momentum equation. The three values correspond to three domain directions.

`is_forced`, if true in the direction in question, **forces the flow** with a pressure gradient that balances the total wall shear (e.g., for a pressure-driven channel with zero net acceleration/constant bulk velocity). The three boolean values correspond to three domain directions.

`velf`, is the **target bulk velocity** in the direction in question (where `is_forced` is true). The three values correspond to three domain directions.

---

```fortran
dims(1:2) = 2, 2, ipencil_axis = 1
```

This line set the domain decomposition and orientation of the computational subdomains.

`dims` is the **processor grid**, the number of domain partitions along the first and second decomposed directions (which depend on the selected default pencil orientation). `dims(1)*dims(2)` corresponds therefore to the total number of computational subdomains. Setting `dims(:) = [0,0]` will trigger a runtime autotuning step to find the processor grid that minimizes transpose times. Note, however, that other components of the algorithm (e.g., collective I/O) may also be affected by the choice of processor grid.

`ipencil_axis` sets the **orientation of the computational subdomains** (or pencils), being one of [1,2,3] for [X,Y,Z]-aligned pencils. X-aligned is the default if this option is not set, and should be optimal for all cases except for Z-implicit diffusion, where using Z-pencils are recommended if `dims(2) > 1` in the input file; see the description of the `&numerics` namelist below.

# about the `&cudecomp` namelist under `input.nml`

In addition to the `&dns` namelist in the input file, there is an **optional** namelist to set some runtime configurations for the *cuDecomp* library. Consider the following `&cudecomp` namelist, which corresponds to the default options in case the file is not provided:

```fortran
&cudecomp
cudecomp_t_comm_backend = 0, cudecomp_is_t_enable_nccl = T, cudecomp_is_t_enable_nvshmem = T
cudecomp_h_comm_backend = 0, cudecomp_is_h_enable_nccl = T, cudecomp_is_h_enable_nvshmem = T
```

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

Finally, it is worth recalling that passing `dims(1:2) = [0,0]` under `&dns` will trigger the *processor grid* autotuning, so there is no need to provide that option in the `&cudecomp` namelist.

# about the `&scalar` namelist under `input.nml`

This namelist is **optional** and defines the parameters needed to solve the transport equations associated with an arbitrary number of scalars.

The **number of scalars** `nscal` has to be set in the `&dns` namelist, e.g., for a single scalar:
```fortran
&dns
! (...)
nscal = 1
\
```
The default value of `nscal` is zero.

---

The following example namelist defines the parameters needed for the differentially heated cavity case that may be found under `examples/`.
```fortran
&scalar
alphai(1)              = 842.61498
beta                   = 1
iniscal(1)             = 'dhc'
cbcscal(0:1,1:3,1)     = 'D'  ,'D' ,  'P','P',  'N','N'
bcscal(0:1,1:3,1)      =  -0.5,0.5 ,   0.,0. ,   0.,0.
ssource(1)             = 0.
is_sforced(1)          = F
scalf(1)               = 0.
is_boussinesq_buoyancy = T
/
```

---

```fortran
alphai(1)        = 842.61498
beta             = 1
```
These lines define the **inverse of the diffusivity** `alphai` of the scalar(s) and, in this case, the **thermal expansion coefficient** `beta`.

`alphai` is an array with size `nscal`. To define several scalars, the inverse diffusivities can be defined as:
```fortran
alphai(:)        = alphai_1, alphai_2, alphai_3, ..., alphai_nscal
```

---

```fortran
iniscal(1)       = 'dhc'
```
This line sets the **initial scalar field(s)**. `iniscal` is an array with size `nscal`. The following options are available:
* `zer`: uniform scalar field equal to zero
* `uni`: uniform scalar field equal to `sref`
* `cou`: linearly varying scalar field from bottom (`z = 0`) to top   (`z = l(3)`)
* `dhc`: linearly varying scalar field from left   (`x = 0`) to right (`x = l(1)`)
* `tbl`: temporal boundary layer profile with scalar bc `sref` at the bottom wall

---

```fortran
cbcscal(0:1,1:3,1) = 'D'  ,'D' ,  'P','P',  'N','N'
bcscal(0:1,1:3,1)  =  -0.5,0.5 ,   0.,0. ,   0.,0.
```
These lines set the **boundary conditions** for each scalar, just like the velocity components. The last dimension corresponds to the scalar index.

---

```fortran
ssource(1)       = 0.
is_sforced(1)    = F
scalf(1)         = 0.
```
These lines set the **scalar forcing**.

`ssource` is an array of size `nscal` defining a **uniform volumetric source term** added to the righ-hand-side of each scalar equation (analogous to `source` in momentum).

`is_sforced` is a logical array of size `nscal` that **triggers the bulk forcing** of the corresponding scalar field (analogous to `is_forced` in momentum).

`scalf` is an array of size `nscal` defining **the target bulk mean** for each scalar (where `is_sforced` is true; analogous to `velf` in momentum).

---

**Buoyancy effects**

The input value `beta` above is the thermal expansion coefficient that relates density to temperature in the so-called _Boussinesq approximation_, thereby triggering buoyancy-induced momentum transport. Its default value is zero. When buoyancy is activated, **the first scalar is always associated with the temperature difference**. Moreover, in that case, the value for the **gravitational acceleration vector** `gacc` should be provided in the `&dns` namelist:
```fortran
&dns
! (...)
gacc(1:3) = 0., 0., -1.
nscal = 1
/
```
which, in this example, defines a negative gravitational acceleration along `z` with unit magnitude. Its default value is `[0., 0., 0.]`.

Finally, **activate the buoyancy term** with `is_boussinesq_buoyancy = T` in the `&scalar` namelist.

# about the `&numerics` namelist under `input.nml`

This namelist defines parameters related to the numerical discretization and computational method. The values below are the default ones, in case the namelist is not specified in the input file.

```fortran
&numerics
is_impdiff = F, is_impdiff_1d = F
is_poisson_pcr_tdma = F
/
```

In these lines, `is_impdiff` and `is_impdiff_1d` enable the (semi-) **implicit temporal integration of diffusion terms**:

* `is_impdiff`, if `.true.`, the diffusion term of the Navier-Stokes and scalar equations is integrated in time implicitly, which may improve the stability of the numerical algorithm for viscous-dominated flows.
* `is_impdiff_1d`, is similar to `is_impdiff`, but with implicit diffusion *only* along Z, which may be advantageous when the grid along Z is much finer than along the other directions; *for optimal parallel performance, the domain should not be decomposed along Z* (`ipencil_axis=3`, or `ipencil_axis = 1/2` with `dims(2) = 1`)

Finally, `is_poisson_pcr_tdma`, if `.true.`, allows for solving the Poisson/Helmhotlz equations along Z with a parallel cyclic reduction--tridiagonal matrix algorithm (PCR-TDMA) method. This approach may result in major gains in scalability for pencil-distributed simulations at scale, on many GPUs.

# about the `&other_options` namelist under `input.nml`

This namelist defines other parameters related to the monitoring, debugging, or benchmarking of a computation. The values below are the default ones, in case the namelist is not specified in the input file.

```fortran
&other_options
is_debug = T, is_timing = T
/
```

In these lines, `is_debug` performs some **sanity checks for debugging purposes** that do not introduce computational overhead, and `is_timing` reports the wall-clock time per step.

Note: other parameters for `&numerics` and `&other_options` are not exposed here, as they are meant for very specific, advanced use or developers. They can be found in `param.f90`.
