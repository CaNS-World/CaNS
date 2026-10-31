# Coding style

This document describes the coding style of *CaNS* and *CaNS-EIGEN*, for contributors and coding agents. The rules apply to new and edited code; existing code that deviates from them does not need to be reformatted, and unrelated code should not be touched. Where this document is silent, follow the style of the surrounding routine.

## Quick checklist

- six-line SPDX banner; one `mod_<file>` module per file, with `implicit none` and `private` plus a `public` list;
- declarations as `real(rp), intent(in), dimension(...) :: x`, with types padded (`integer ,`) and `::` aligned; intents padded only when they differ;
- 2-space indentation, no tabs, no blank lines (use a bare `!`), at most 132 characters per line;
- related lines aligned vertically: declarations, groups of assignments and calls, stencils, and paired directives;
- no spaces in `if(`, `do k=1,n(3)`, `call f(a,b)`, or around `.and.`/`.or.`/`.not.`; spaces around ` = `, ` == `, ` /= `;
- real literals bare only if exactly representable (`0.`, `.5`, `2.`), otherwise with `_rp`;
- field kernels as explicit `k`-`j`-`i` loops with an `!$acc parallel loop ... async(1)` directive and, right below it, the aligned `!$OMP PARALLEL DO` counterpart;
- messages prefixed with `ERROR: `, `WARNING: `, or wrapped in `*** ... ***`, printed by rank 0;
- only OpenACC constructs that have an OpenMP target-offload counterpart (see [OpenMP port](#keeping-main-portable-to-the-openmp-port)).

## Files and modules

- Every source file (`.f90` and `.h90`) starts with this banner:
  ```fortran
  ! -
  !
  ! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
  ! SPDX-License-Identifier: MIT
  !
  ! -
  ```
  An optional description of the file follows as a `!` block.
- One module per file, named `mod_<file>` (`bound.f90` holds `mod_bound`). The module skeleton is:
  ```fortran
  module mod_bound
    use mpi
    use mod_common_mpi, only: ierr,halo,myid
    use mod_param     , only: ipencil_axis
    use mod_types
    implicit none
    private
    integer, parameter :: HALO_LOWER = -1, HALO_NONE = 0, HALO_UPPER = 1
    public boundp,bounduvw,updt_rhs_b,set_bc,bc_dist
    contains
    subroutine bounduvw(...)
      ...
    end subroutine bounduvw
    !
    subroutine boundp(...)
  ```
- `use` statements come in this order: intrinsic modules (`use, intrinsic :: ...`), `mpi`, `decomp_2d`, the `mod_*` modules, and `mod_types` last. All modules are used with `only:`, except `mpi`, `decomp_2d`, and `mod_types`. In a group of `use` statements, the module names are padded so that the commas line up; long `only:` lists continue aligned under the first name.
- `implicit none` appears in the module and again in every procedure.
- Modules are `private` with an explicit `public` list; a bare `public` is only for modules that hold global state (`mod_param`, `mod_common_mpi`).
- `contains` is indented like the module body (2 spaces); procedures are indented 2 spaces and their bodies 4.
- Procedures end with their name: `end subroutine bounduvw`, `end function is_output_due`, `end module mod_bound`.
- Consecutive procedures are separated by a line with a bare `!`.
- Internal procedures go under a `contains` at the end of the host procedure, indented 2 more spaces.
- Procedures without arguments are declared and called without parentheses: `subroutine abortit`, `call abortit`.

## Declarations

- Kinds come from `mod_types`: `rp` for reals (`sp` or `dp`, depending on `SINGLE_PRECISION`), `dp` only where double precision is required regardless of `rp` (e.g., timers), `i8` for large integer sizes, and `MPI_REAL_RP` for MPI.
- Attributes go in this order: type, `intent`, `target`/`pointer`/`contiguous`/`allocatable`, `dimension(...)`, `optional`, `save`. Shapes are given with `dimension(...)`, not after the name.
- Types are padded to the same width, so that `integer ,` and `logical ,` line up with `real(rp),`, and `::` is aligned within a group:
  ```fortran
  integer , intent(in), dimension(3) :: n
  integer , intent(in)               :: nh,nh_g
  ```
- When all dummy arguments of a procedure are `intent(in)`, write `intent(in)`. When intents differ, pad all of them to the longest one present (`intent(in   )`, `intent(inout)`, `intent(out  )`, or `intent(in )`, `intent(out)`), as in `updatep.f90` (abridged):
  ```fortran
  integer , intent(in   ), dimension(3) :: n
  integer , intent(in   )               :: nh,nh_g
  real(rp), intent(in   ), dimension(1-nh_g:) :: dzci,dzfi
  real(rp), intent(in   ), dimension(1-nh:,1-nh:,1-nh:) :: pp
  real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
  ```
- Array dummy arguments are assumed-shape with explicit lower bounds that encode the halo: `dimension(1-nh:,1-nh:,1-nh:)` for fields, `dimension(1-nh_g:)` for grid arrays, `dimension(lo(1)-nh:,lo(2)-nh:,lo(3)-nh:)` for globally indexed arrays. Arrays of boundary-condition data have the shape `(0:1,3[,ivar])`: side, direction, and variable.
- Dummy arguments are declared in the order of the argument list, then the locals: arrays and `save` variables first, then scalars, then loop indices.
- Arrays that are computed once are `save` variables, allocated on the first call:
  ```fortran
  real(rp), save :: dlmin
  logical , save :: is_first = .true.
  ...
  if(is_first) then ! calculate dlmin only once
    is_first = .false.
    ...
  end if
  ```
  or guarded by `if(.not.allocated(x))`, followed by `!$acc enter data create(x)` for device arrays.
- Input parameters live in `mod_param`, are `protected`, and get their default values before the input file is read.

## Literals and arithmetic

- A real literal may be written without kind only if its value is exactly representable: `0.`, `1.`, `2.`, `.5`, `.25`. Any other value carries `_rp`, so that it is not evaluated in single precision in double-precision builds: `1._rp/3.`, `1.65_rp/12.`, `acos(-1._rp)`. Arguments of inquiry intrinsics carry `_rp` too: `epsilon(1._rp)`, `huge(1._rp)`.
- Literals are written `0.`, `.5`, `1.e-12`, never `0.0`, `0.5d0`, or `1.d-12`.
- Integers are converted to reals with `1.*n` (`dl(:) = l(:)/(1.*ng(:))`); reciprocals are written `1./x` (`x**(-1)` is also accepted).

## Naming

- Everything is lowercase, with words separated by underscores where needed. MPI routines and constants are uppercase (`call MPI_ALLREDUCE(MPI_IN_PLACE,...)`), external APIs keep their own case (`cudecompTransposeXtoY`), and named constants that act as enumerators or flags are uppercase (`HALO_LOWER`, `MPI_REAL_RP`).
- Names are terse, in the style of the existing ones. Routines start with a short verb: `chk*`/`chk_` (checks: `chkdt`, `chk_bc`), `cmpt_` (compute: `cmpt_wallshear`), `updt_`/`upd` (update: `updt_rhs_b`, `updatep`), `init*` (initialization), `out*` (output). Variants get a suffix: `_gpu` for the GPU version, `_s` for the version for scalars, and letters for parts of a term (`momx_a`, `momx_d`, `momx_p` for advection, diffusion, pressure).
- Variables follow this vocabulary:

  | name | meaning |
  |---|---|
  | `n`, `ng`, `lo`, `hi` | local sizes, global sizes, global index range of the local block (`dimension(3)`) |
  | `n_z`, `lo_z`, `hi_z` | the same, for the z-pencil layout of the Poisson solver |
  | `nh`, `nh_v`, `nh_p`, `nh_s`, `nh_g` | halo width of the arrays of a routine; of the velocity, pressure, scalars, and grid arrays |
  | `l`, `dl`, `dli` | domain lengths, uniform grid spacings, and their inverses |
  | `zc`, `zf`, `dzc`, `dzf`, `dzci`, `dzfi` | cell-center and face coordinates, spacings, and inverse spacings |
  | `<name>_g` | global (non-distributed) copy of an array, e.g. `dzci_g` |
  | `cbc<var>`, `bc<var>` | boundary condition types (`'P'`, `'D'`, `'N'`) and values, e.g. `cbcvel`, `bcpre` |
  | `nb`, `is_bound` | neighboring tasks, and whether a side of the block is a domain boundary |
  | `c_or_f` | location of a variable (`'c'`: cell center, `'f'`: face) |
  | `u`, `v`, `w`, `p`, `pp`, `s` | velocity components, pressure, pressure correction, scalar |
  | `dudt`, `dudtrk`, `dudtrko` | right-hand sides; `o` marks a value of the previous substage |
  | `...p`, `...m` | stencil terms at the plus/minus neighbor, e.g. `uuip`, `dudxm` |
  | `is_<what>` | logical variables, e.g. `is_forced`, `is_first` |
  | `i<what>`, `n<what>` | integer identifiers and counts, e.g. `iunit`, `istep`, `nscal`, `nstep` |
  | `<name>_loc` | local part of a result, e.g. `passed_loc` |

- Loop indices: `i`, `j`, `k` along x, y, z; `idir` for the direction, `ib` for the side (0 or 1), `ivel`, `iscal`, `irk` for velocity components, scalars, and Runge-Kutta substages.
- New logical variables start with `is_`.
- An optional argument with a default value keeps its public name, and its local copy with the default applied gets a trailing underscore:
  ```fortran
  logical :: is_dtdma_update_
  ...
  is_dtdma_update_ = .true.
  if(present(is_dtdma_update)) is_dtdma_update_ = is_dtdma_update
  ```

## Argument order

Arguments come in this order: mode or boundary-condition selectors, sizes (`n`, `ng`, `lo`, `hi`), halo widths (field widths `nh_v`, `nh_p`, `nh_s`, or `nh`, then `nh_g`), boundary-condition values and topology (`bc`, `nb`, `is_bound`), grid (`l`, `dl`, `dli`, `dzc`, `dzf`, ...), scalar parameters (`dt`, `visc`, `alpha`), input arrays, input/output fields, outputs, and optional arguments:

```fortran
subroutine bounduvw(cbc,n,nh,nh_g,bc,nb,is_bound,dl,dzc,dzf,u,v,w,keep_norm_values,halos_norm,halos_tang)
subroutine correc(n,nh_v,nh_p,nh_g,dli,dzci,dt,p,u,v,w)
```

## Layout

- Indentation is 2 spaces per level; no tabs, no trailing spaces.
- Lines are at most 132 characters long.
- There are no blank lines: a line with a bare `!` separates blocks of code, declarations from statements, and procedures.
- Continuation lines end with ` &` (with a space) and align with the opening parenthesis, the first argument, or the right-hand side of the assignment. A leading `&` is used only to continue a character string:
  ```fortran
  print*, 'ERROR: `is_bulk_forcing_wall_shear = T` requires no-slip walls for the forced velocity components, &
          &no fully implicit diffusion with walls along x or y, and no gravity along a forced direction with &
          &Boussinesq buoyancy.'
  ```
- Spacing:
  - no space between a keyword and its parenthesis: `if(`, `select case(`, `do while(`, `allocate(`;
  - no spaces in loop headers and argument lists: `do k=1,n(3)`, `call set_bc(cbc(0,1,1),0,1,nh,.false.,bc(0,1,1),drf(:,0,1),u)`;
  - spaces around the assignment ` = `, but not around `=` in keyword arguments: `iostat=ierr`, `mold=dzc`;
  - spaces around top-level `+` and `-`, none around `*`, `/`, `**`, or in index arithmetic: `u(i,j,k) = u(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k)`;
  - spaces around relational operators, which are always symbolic (`==`, `/=`, `<`, `>=`, never `.eq.`): `if(myid == 0)`;
  - no spaces around the logical operators `.and.`, `.or.`, `.not.`: `if(myid == 0.and.(.not.passed_loc))`, `is_bound(0,1).and.update_halo(0)`;
  - `print*, 'text', x` with a space after the comma.
- Related lines are aligned vertically; see [Vertical alignment](#vertical-alignment).
- Short, related statements may share a line with semicolons: `nx = n(1); ny = n(2); nz = n(3)`.
- `else if`, `end if`, `end do`, `end select` are written as two words; `case` is indented like its `select case`.
- Single-statement conditionals are written on one line: `if(present(keep_norm_values)) keep_norm = keep_norm_values`.
- Errors terminate the program with `error stop`, never `stop`.

## Vertical alignment

Lines that belong together are padded with spaces so that corresponding parts start in the same column. This shows the structure of the code and makes differences between similar lines stand out.

- Declarations: types are padded (`integer ,`, `logical ,`), intents are padded when they differ (see [Declarations](#declarations)), and `dimension(...)` lists are padded so that `::` lines up within a group:
  ```fortran
  real(rp), intent(in), dimension(0:1,3,3) :: bc
  integer , intent(in), dimension(0:1,3  ) :: nb
  logical , intent(in), dimension(0:1,3  ) :: is_bound
  ```
- `use` statements: module names are padded so that the commas line up (see [Files and modules](#files-and-modules)).
- Groups of assignments: the `=` signs line up:
  ```fortran
  beta          = 0.
  ssource(:)    = 0.
  is_sforced(:) = .false.
  scalf(:)      = 0.
  ```
- Stencils: indices are padded so that the offsets line up across terms and lines:
  ```fortran
  uuip = 0.25*( u(i,j  ,k  )+u(i+1,j  ,k  ) )*( u(i,j,k)+u(i+1,j,k) )
  vujm = 0.25*( v(i,j-1,k  )+v(i+1,j-1,k  ) )*( u(i,j,k)+u(i,j-1,k) )
  wukm = 0.25*( w(i,j  ,k-1)+w(i+1,j  ,k-1) )*( u(i,j,k)+u(i,j,k-1) )
  ```
- Groups of calls: arguments line up; a shorter argument is padded (`.true. ,` under `.false.,`), and a missing leading part (e.g., a condition) is replaced by spaces:
  ```fortran
  if(apply_norm_bc(0,1)) call set_bc(cbc(0,1,1),0,1,nh,.false.,bc(0,1,1),drf(:,0,1),u)
                         call set_bc(cbc(0,1,2),0,1,nh,.true. ,bc(0,1,2),drc(:,0,1),v)
  ```
- Lists of names on consecutive lines, e.g., in directive clauses, are padded so that related names line up:
  ```fortran
  !$acc private(dudt_s ,dvdt_s ,dwdt_s ) &
  !$acc private(dudtd_s,dvdtd_s,dwdtd_s) &
  ```
- Continuation lines align with the opening parenthesis, the first argument, or the right-hand side (see [Layout](#layout)); the arrows of argument descriptions line up (see [Comments](#comments)).
- Paired directives: when two directive lines describe the same construct (an OpenACC directive and its CPU OpenMP counterpart in main, or its OpenMP target-offload counterpart in the `openmp-port` branch), each token of one line starts in the same column as the corresponding token of the other:
  - the shorter keyword or clause is padded with spaces, on whichever line it is (`!$acc parallel loop ` above `!$OMP PARALLEL DO   `);
  - when the clause arguments differ only in a prefix, the padding goes inside the parentheses, so that the variable lists line up (`create(   p1d)` above `map(alloc:p1d)`, `device(gg)` above `to(    gg)`);
  - a clause that has no counterpart is replaced by spaces on the other line (e.g., `default(present) ` or `gang vector `), unless it is the last clause (e.g., `async(1)`), which is simply omitted.

## Comments

- Comments are lowercase and terse, usually without a final period; identifiers and file names go in backticks.
- Block comments are sandwiched between bare `!` lines, at the indentation of the code:
  ```fortran
  !
  ! mirror distances of the centered (drc) and non-centered (drf) components, computed once
  !
  ```
- A procedure starts with a short description of what it does, right after its first line, e.g., `! imposes velocity boundary conditions`. Arguments, when documented, are listed as `! name -> meaning`, aligned:
  ```fortran
  ! fname -> name of the file
  ! n     -> number of entries
  ! var   -> input array of real values
  ```
- Caveats are introduced with `! note:`; references are cited inline, e.g., `! see M. Antuono, JFM 890, A23 (2020)`.
- Commented-out code is allowed sparingly, when it documents an option (e.g., optional outputs in `out1d.h90`).

## Loops and kernels

- Kernels on fields are explicit triple loops, with `k` outermost and `i` innermost.
- Array syntax, with an explicit `(:)` on the left-hand side, is for setup and one-dimensional work on the host: `dl(:) = l(:)/(1.*ng(:))`.
- In code that runs on the GPU, use explicit loops: no array syntax and no `!$acc kernels`.

## OpenACC

- Compute loops use one form, with the clauses in this order:
  ```fortran
  !$acc parallel loop collapse(3) default(present) private(div) reduction(+:divtot) reduction(max:divmax) async(1)
  ```
  Long clause lists continue with ` &` and a new `!$acc` line per group of clauses (`!$acc private(...) &`).
- Arrays must be present on the device (`default(present)`); scalars are not listed (they are firstprivate by default). Small arrays used in kernels (e.g., `lo`, `hi`, `n`, `dl`, `dli`) are copied to the device once, in the main program.
- Device arrays are created right after they are allocated, with `!$acc enter data create(...)` or `copyin(...)`; for derived types, the parent first, then its components: `!$acc enter data copyin(rhsbp,rhsbp%x,rhsbp%y,rhsbp%z)`. Short-lived device arrays are paired with `!$acc exit data delete(...)`.
- All device work in the time loop (kernels, data regions, `enter`/`exit data`) runs on queue 1, `async(1)`, whose CUDA stream is also used by cuFFT and cuDecomp. Wait with `!$acc wait(1)` before the host uses the data or before MPI. Setup, I/O, and diagnostics may be synchronous. Do not use a bare `async` in new code.
- A reduction to a host scalar follows this pattern (from `chkdiv.f90`, abridged):
  ```fortran
  divtot = 0.
  !$acc data copy(divtot) async(1)
  !$acc parallel loop collapse(3) default(present) private(div) reduction(+:divtot) async(1)
  !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  PRIVATE(div) REDUCTION(+:divtot)
  do k=lo(3),hi(3)
    ...
  end do
  !$acc end data
  !$acc wait(1)
  call MPI_ALLREDUCE(MPI_IN_PLACE,divtot,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  ```
- Device pointers are passed to libraries (cuFFT, cuDecomp) through `!$acc host_data use_device(...)`, guarded for diezDecomp, which takes host arrays:
  ```fortran
  #if !defined(_USE_DIEZDECOMP)
      !$acc host_data use_device(px,py,work)
  #endif
      istat = cudecompTransposeXtoY(ch,gd,px,py,work,dtype_rp,stream=istream)
  #if !defined(_USE_DIEZDECOMP)
      !$acc end host_data
  #endif
  ```
- Line-by-line solvers (e.g., tridiagonal systems along z) use `!$acc parallel loop gang vector collapse(2) ... async(1)` over the other two directions, with `!$acc loop seq` along the solve.

## CPU OpenMP

- Each OpenACC compute loop that also runs on the CPU has its OpenMP counterpart right below, in uppercase, padded so that the clauses line up with those of the OpenACC line (see [paired directives](#vertical-alignment)):
  ```fortran
  !$acc parallel loop collapse(3) default(present) private(div) reduction(+:divtot) async(1)
  !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  PRIVATE(div) REDUCTION(+:divtot)
  ```
  Long lists continue with ` &` and a new `!$OMP PRIVATE(...)` line.
- Loops that run only on the CPU use single spaces: `!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(shared) PRIVATE(den,pivot_tol)`.
- Statements needed only with OpenMP use the conditional sentinel: `!$ use omp_lib`, `!$ call omp_set_num_threads(...)`.

## Preprocessor

- Conditional compilation uses `#if defined(_X)`, `#if !defined(_X)`, `#elif defined(_X)`, and `#endif`, at column 1 regardless of the indentation of the code; no `#ifdef`.
- Macros are uppercase with a leading underscore (`_USE_HIP`, `_SINGLE_PRECISION`). Build options become macros through `CUSTOM_DEFINES` in `configs/defaults/flags-default.mk` (a `build.conf` variable `X=1` defines `_X`); `_OPENACC` comes from the compiler (`GPU=1`).
- CPU and GPU code share one source by default: the same loop, with both directives. A separate GPU version (a `_gpu` procedure, or a GPU module wrapped entirely in `#if defined(_OPENACC)` and renamed at the `use`, e.g. `use mod_solver_gpu, only: solver => solver_gpu`) is only for algorithms or library calls that differ. Small differences go in `#if !defined(_OPENACC)` ... `#else` ... `#endif` blocks.

## MPI, errors, and messages

- Use `use mpi`, with `myid` and `ierr` from `mod_common_mpi`, and `MPI_REAL_RP` for reals of kind `rp`. Global reductions are done in place: `call MPI_ALLREDUCE(MPI_IN_PLACE,x,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)`.
- Only rank 0 prints, one guarded statement per line: `if(myid == 0) print*, '*** Beginning of simulation ***'`.
- Messages start with `ERROR: ` or `WARNING: `; status lines are wrapped in `*** ... ***`.
- A fatal error is reported by rank 0, followed by a clean exit:
  ```fortran
  if(myid == 0) print*, 'ERROR: reading the input file: ', trim(c_iomsg)
  if(myid == 0) print*, 'Aborting...'
  call MPI_FINALIZE(ierr)
  error stop
  ```
  Checks of the input in `sanity.f90` call `abortit` instead.
- Checks of the input are `chk_<what>(...,passed)` routines in `sanity.f90`, each made of tests of the form:
  ```fortran
  passed = .true.
  passed_loc = all(ng(:) >= 1)
  if(myid == 0.and.(.not.passed_loc)) &
    print*, 'ERROR: ng(:) must be positive in every direction.'
  passed = passed.and.passed_loc
  ```
  and the caller does `call chk_<what>(...,passed); if(.not.passed) call abortit`.
- Input namelists are read in `param.f90`: defaults first, then `rewind(iunit)` and `read(iunit,nml=<name>,iostat=ierr,iomsg=c_iomsg)`. Mandatory namelists abort on any error; optional ones accept `ierr == iostat_end`.

## Keeping main portable to the OpenMP port

The OpenMP target-offload version of the code lives in the `openmp-port` branch. There, each OpenACC directive is followed by its OpenMP counterpart, the CPU `!$OMP` directives are removed, and the guards `#if defined(_OPENACC)` become `#if defined(_OPENACC) || defined(_OPENMP)`. To keep this port mechanical, code in main uses only OpenACC constructs that have a counterpart in this table:

| OpenACC (main) | OpenMP (`openmp-port`) |
|---|---|
| `parallel loop [gang vector] collapse(N) default(present) private(...) reduction(...) async(1)` | `target teams loop collapse(N) private(...) reduction(...)` |
| `loop collapse(N) reduction(...)` (inside a `parallel loop gang`) | `loop collapse(N) reduction(...)` |
| `parallel default(present) async(1)` ... `end parallel` (serial code) | `target` ... `end target` |
| `data copy(x)` / `copyin(x)` / `copyout(x)` | `target data map(tofrom:x)` / `map(to:x)` / `map(from:x)` |
| `enter data create(x)` / `copyin(x)` | `target enter data map(alloc:x)` / `map(to:x)` |
| `exit data delete(x)` / `copyout(x)` | `target exit data map(delete:x)` / `map(from:x)` |
| `update self(x)` / `update device(x)` | `target update from(x)` / `target update to(x)` |
| `host_data use_device(x)` ... `end host_data` | `target data use_device_addr(x)` ... `end target data` |
| `wait`, `wait(1)` | no counterpart (OpenMP target constructs without `nowait` are synchronous) |

In the port, the OpenMP directive goes right below the OpenACC one when a construct opens, and right above it when the construct closes, so that the two nest properly. The two lines follow the rule for [paired directives](#vertical-alignment): corresponding keywords, clauses, and variable lists start in the same column, the shorter one padded with spaces, and a clause without counterpart (e.g., `default(present)`, `gang vector`) leaves a gap:

```fortran
!$acc parallel     loop collapse(3) default(present) private(i) async(1)
!$omp target teams loop collapse(3)                  private(i)
!$acc parallel     loop gang vector collapse(2) default(present) private(den,pivot_tol) async(1)
!$omp target teams loop             collapse(2)                  private(den,pivot_tol)
!$acc        enter data create(   p1d)
!$omp target enter data map(alloc:p1d)
!$acc        update device(gg(:,irk)) async(1)
!$omp target update to(    gg(:,irk))
!$acc        data copy(      divtot,divmax) async(1)
!$omp target data map(tofrom:divtot,divmax)
...
!$omp end target data
!$acc end        data
!$acc   host_data use_device(     p,work)
!$omp target data use_device_addr(p,work)
...
!$omp end target data
!$acc end   host_data
```

Hence, in main:

- use only the constructs above; in particular no `!$acc kernels` and no array syntax in device code (the port rewrites them as explicit loops);
- give every device array a lifetime that maps to `enter`/`exit data` or a structured `data` region;
- keep calls to the OpenACC runtime API (`acc_*`) in guarded setup code (e.g., `initmpi.f90`, `workspaces.f90`, `common_cudecomp.f90`).

## Syncing `openmp-port` with main

1. Create a branch off `openmp-port` (e.g., `updt-omp-port`) and merge `main` into it.
2. Resolve each conflict by keeping the code of main and applying the port's transformation to the new or changed lines:
   - add the OpenMP counterpart of each OpenACC directive, following the table and padding above;
   - widen the `_OPENACC` guards;
   - drop the CPU `!$OMP` directives;
   - keep code that exists only in the port.
3. Check that the merged sources contain:
   - no CPU `!$OMP PARALLEL` directives;
   - no OpenACC directive without its OpenMP counterpart that was not there before;
   - no new guard on `_OPENACC` alone.
4. Build with `GPU_BACKEND=OMP` and with the default OpenACC backend, and compare:
   - results with the OpenMP build of the port before the merge;
   - results with the OpenACC build of the same sources;
   - the CI tests run with the OpenMP build.

   Cases without global sums (e.g., cavities) agree bitwise. With nvfortran 26.9, OpenMP GPU sum reductions are not reproducible from run to run, so cases that depend on them (e.g., bulk forcing) agree only to round-off. When several tasks share one GPU, disable NCCL and NVSHMEM in `&cudecomp`.
5. Open a pull request into `openmp-port` and run the CI workflow by hand (`run-all`, `workflow_dispatch`), since it does not run for pull requests into that branch. Merge it with a merge commit, and delete the branch.
6. Do the same for *CaNS-EIGEN*.
