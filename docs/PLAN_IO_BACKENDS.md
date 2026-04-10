# I/O Backend Implementation Plan

## Scope

This document lays out a clean implementation plan for extending the CaNS I/O stack with:

- ADIOS2 checkpoint I/O
- HDF5-style checkpoint I/O
- HDF5-compatible NetCDF-style time-series I/O for visualization
- small updates to the existing XDMF helpers under `utils/`
- placement of a serial NumPy writer utility under `utils/`

The main functional requirement across all of these paths is:

- write arrays while excluding halos cleanly and consistently

## Goals

### ADIOS2

1. Checkpoint I/O for `u,v,w,p` with optional grid and metadata.
2. Checkpoint I/O for a single field with optional grid and metadata.
3. Configurable `compression_level`.
4. A visualization-oriented ADIOS2 mode with a suitable compression choice.
5. Configurable output precision independent of the run precision.

### HDF5 / NetCDF-style visualization

1. Standard checkpoint I/O for `u,v,w,p` with optional grid and metadata.
2. Standard checkpoint I/O for a single field with optional grid and metadata.
3. Compression support.
4. An extendible visualization format for `u,v,w,p` plus grid and time metadata, following the current HDF5-compatible NetCDF-style path that already worked with `ncview`.
5. Minimal updates to the XDMF helpers under `utils/` so the same visualization outputs remain easy to inspect in ParaView/VisIt.
6. Configurable output precision independent of the run precision.

### Other

1. Place the Fortran NumPy writer under `utils/`.
2. Keep it serial-only.

## Design Principles

### 1. Follow the current `load.f90` philosophy

As much as possible, the new I/O paths should follow the existing structure already familiar in `load.f90`:

- `load_one`
- `load_all`

That means the new backend work should aim for parallel concepts such as:

- one-field checkpoint operations
- multi-field checkpoint operations

instead of inventing a completely separate public calling style for every backend.

Internally, these routines can still use saved state and backend-specific helpers.

### 2. Split setup from data movement

For both ADIOS2 and HDF5, the cleaner structure is:

- `*_init(...)`
- `*_write(...)`
- `*_read(...)`
- `*_close(...)`

This avoids repeating expensive setup work every write and makes restart support much easier.

### 3. One initializer per schema

Initialization can be reused as long as the schema does not change. In practice the schema is defined by:

- file/stream identity
- variable names
- number of variables
- types
- global shapes
- compression/operator setup

This suggests separate logical modes, for example:

- full checkpoint: `u,v,w,p`
- single-field checkpoint
- time-series visualization output

These can share internal helper routines, but should keep separate saved state.

### 4. Halo exclusion must stay explicit

Halo exclusion should be implemented as part of the backend API, not as a caller-side convention.

- For ADIOS2: use `adios2_set_memory_selection` for contiguous interior data.
- For HDF5: use file/memory hyperslabs.

This keeps the interface safe and avoids accidental halo writes.

### 5. Preserve decomposition independence

All files should be defined in terms of global array shapes and global offsets so they can be read back with a different MPI decomposition.

### 6. Separate native-field output from visualization-friendly output

For staggered velocity fields, the plan should explicitly distinguish:

- native/raw solver storage locations
- visualization-friendly colocated fields

Policy:

- raw staggered velocity output remains supported
- XDMF/visualization can assume the existing nearest-neighbor-style interpretation when using raw components
- optionally, a cell-centered interpolated velocity field can be generated for visualization

This keeps checkpoint/restart faithful while still supporting cleaner vector visualization when requested.

### 7. Precision should be a user-controlled I/O property

By default, all backends should write using the run precision.

Optionally, all backends should allow writing in a chosen output precision, e.g.:

- `sp`
- `dp`

This is useful for:

- smaller visualization files
- reduced storage cost
- experimentation with lower-precision checkpoints

Checkpointing in reduced precision should be supported as an explicit option, not as the default.

## Format Decision: keep the working `ncview`-compatible path

Based on the current tests, the HDF5-compatible NetCDF-style writer that was added already worked with `ncview`, so this should be treated as the visualization path to preserve and extend.

This means we should not plan around arbitrary plain HDF5 for `ncview`. Instead, we should:

- keep the current NetCDF-style structure that `ncview` already accepts
- extend that format for time series
- make the XDMF generator understand that layout with only small changes

More generally, the structured visualization-oriented HDF5 family should aim to stay netCDF-compatible where practical, so that one backend layout can serve multiple consumers.

Plan implication:

- checkpoint HDF5 can remain a separate backend
- the visualization time-series path should follow the current `ncview`-compatible layout
- regular structured visualization outputs, including single-field outputs where practical, should also prefer netCDF-compatible organization
- XDMF should be updated to target that layout rather than introducing a second visualization-only file structure

Note:

The `ncview` visualization path should be handled as netCDF-style HDF5, not generic HDF5. Checkpoint HDF5 may use a freer internal layout, but visualization files intended for `ncview` should stay within netCDF-4-compatible conventions for dimensions, coordinate metadata, and extendible time-series structure.

## Proposed API Structure

The public shape should stay close to the `load_one` / `load_all` philosophy.

## Code Organization

### `load.f90`

`load.f90` should own the backend workflow and low-level I/O engine for:

- checkpoint/restart I/O
- generic structured subset-output I/O
- optional unstructured/point-cloud I/O helpers

This is where the implementation should keep:

- schema setup
- backend state
- halo exclusion logic
- precision conversion
- compression logic
- generic direct-selection logic
- specialized packing kernels when needed
- actual backend calls

### `output.f90`

`output.f90` should own the wrappers and output policies for visualization-style outputs.

In particular, `output.f90` is the natural place for wrappers that request:

- line outputs
- plane outputs
- volume outputs
- coarsened structured outputs

These wrappers should build a descriptor and call the reusable engine in `load.f90`.

### Structured vs unstructured

The generic subset-output API should cover only structured cases:

- line
- plane
- volume
- coarsened structured field
- coarsened structured subvolume

The structured subset descriptor should be centered on:

- global lower bounds `nmin`
- global upper bounds `nmax`
- skip/coarsening vector `nskip`

The non-coarsened case is simply:

- `nskip = [1,1,1]`

So the same user-facing API should handle both:

- full-resolution subsets
- coarsened subsets

Unstructured point-cloud output should stay in a separate API because:

- the packing model is different
- the metadata/XDMF model is different
- it is better treated as point-based rather than extent-based I/O

## Structured subset-output backends

The first structured subset-output backends should be:

- raw binary
- HDF5

These two cover the main needs well:

- raw binary for lightweight MPI-I/O and XDMF workflows
- HDF5 for self-describing subset files and Python-friendly access

ADIOS2 subset-output support can be added later, but it should not block the first implementation.

However, the API and file-layout choices should remain open enough that ADIOS2 can later support:

- line outputs
- plane outputs
- volume/subvolume outputs
- coarsened structured outputs

This is worth keeping open because ADIOS2 also offers attractive compression and ParaView-oriented metadata options.

The user-facing API should be unified, but the backend execution path does not have to be identical for all cases.

Recommended rule:

- if `nskip = [1,1,1]`, prefer the direct compact structured write path
- if `nskip /= [1,1,1]`, allow backend-specific choice between:
  - direct strided selection
  - packed write path

The packed path does not need to be "generic" in the same sense as the public API.

It is fine if the internal implementation uses:

- one generic direct-selection path
- a small number of specialized packed kernels for common cases

for example:

- line packing
- plane packing
- volume/subvolume packing
- coarsened structured packing

This preserves a clean interface while leaving room for performance-driven implementation choices.

## Unstructured point-cloud output

For point-cloud output, the first target should be:

- packed raw binary written in parallel
- Python XDMF generation in `utils/`

As a baseline, it is reasonable to include or document a small Fortran subroutine implementing this style of unstructured output for a particle-tracking code.

That baseline should be treated as:

- an example of the metadata/layout conventions
- guidance for how the unstructured XDMF path should look in practice

Suggested write pattern:

1. each task counts locally selected points
2. compute global displacements with `MPI_Exscan`
3. pack local point coordinates and attributes
4. write the packed global arrays collectively

HDF5 support for point clouds can be added later if needed.

## ADIOS2 API

Implementation note:

- when adding or updating ADIOS2 bindings, check the official ADIOS2 documentation directly for the Fortran API signatures and expected usage
- do not infer the CaNS implementation only from third-party wrappers or helper libraries
- local dependency code can still be useful as a secondary reference for patterns, but the ADIOS2 docs should be the primary source of truth

### Full checkpoint

- `adios2_load_all_init(...)`
- `adios2_load_all(io, filename, comm, ng, nh, lo, hi, nscal, u, v, w, p, s, time, istep, ...)`
- `adios2_load_all_close()`

### Single-field checkpoint

- `adios2_load_one_init(...)`
- `adios2_load_one(io, filename, comm, ng, nh, lo, hi, var, time, istep, varname, ...)`
- `adios2_load_one_close()`

### Compression mode

Use a small internal convention:

- `compression_level = 0` means no compression
- `compression_level = 1..9` means lossless Blosc compression level

For a visualization-oriented mode:

- default to modest compression, e.g. `3` to `5`
- keep checkpoints conservative, defaulting to `0` unless a user explicitly enables compression
- allow combining compression with chosen output precision

### Precision mode

Use a small internal convention such as:

- default: same as solver precision
- optional: explicit single-precision output
- optional: explicit double-precision output

For ADIOS2 in particular, lossy compression plus single precision is a sensible visualization-oriented combination, but should not be the default checkpoint path.

## HDF5 API

### Full checkpoint

- `hdf5_load_all_init(...)`
- `hdf5_load_all(io, filename, comm, ng, nh, lo, hi, nscal, u, v, w, p, s, time, istep, ...)`
- `hdf5_load_all_close()`

### Single-field checkpoint

- `hdf5_load_one_init(...)`
- `hdf5_load_one(io, filename, comm, ng, nh, lo, hi, var, time, istep, varname, ...)`
- `hdf5_load_one_close()`

For visualization-oriented single-field outputs, it is reasonable to use the same HDF5 backend while preferring a netCDF-compatible layout when the field is a regular structured array.

### Time-series visualization output

- `visu_series_init_uvwp(filename, ng, nh, lo, hi, x_g, y_g, z_g, compression_level)`
- `visu_series_append_uvwp(u, v, w, p, time)`
- `visu_series_close_uvwp()`

The time-series file should use the current `ncview`-compatible layout style and:

- extendible `time` dimension
- extendible datasets for `u,v,w,p`
- optional `istep(time)` integer metadata aligned with the time dimension
- one static grid definition
- chunked layout, which is required for both extension and compression
- optional raw staggered component storage
- optional cell-centered interpolated velocity output for visualization

## Internal State Strategy

Use backend-specific saved module state rather than repeated all-in-one routines.

For each logical mode, store:

- initialized flag
- file/engine handle
- IO/group/dataset handles
- variable handles
- schema description
- whether grid/meta are enabled
- current compression settings

This is preferable to repeatedly opening/closing/redefining everything on every write.

## Implementation Breakdown

## Phase 1: Refactor ADIOS2 into stateful `load_one` / `load_all` style routines

### Tasks

1. Extract current ADIOS2 checkpoint logic from the monolithic routine into stateful routines.
2. Support both:
  - `load_all`-style `u,v,w,p`
  - `load_one`-style single-field mode
3. Keep halo exclusion via `adios2_set_memory_selection`.
4. Add read support for restart.
5. Make decomposition-independent reads explicit by setting selections from `lo/hi`.
6. Add output-precision selection.

### Notes

- ADIOS2 is a strong candidate for flexible checkpoints and appendable outputs.
- Full-field checkpoint and single-field checkpoint should be separate schemas.
- Compression should remain optional and off by default for checkpoints until it is tested robustly on the target build.
- ADIOS2 should support visualization-oriented combinations such as reduced precision plus compression without making them the default restart path.

## Phase 2: Refactor HDF5 into matching stateful `load_one` / `load_all` style routines

### Tasks

1. Refactor the current HDF5 routines into the same pattern:
  - init
  - write/read
  - close
2. Add a single-field mode that mirrors the ADIOS2 API.
3. Keep halo exclusion via memory/file hyperslabs.
4. Add HDF5 compression support for chunked datasets.
5. Add output-precision selection.

### Compression note

For HDF5 compression:

- use chunked datasets
- expose a simple user-facing compression level
- hide filter details internally

We should expect:

- compression only for chunked datasets
- some performance tradeoff
- a cleaner path for full checkpoints than for very frequent output

### Precision note

HDF5 paths should support:

- native run precision by default
- optional selected output precision for smaller files or visualization-oriented workflows

## Phase 3: Extend the current `ncview`-compatible visualization path

### Tasks

1. Create a dedicated extendible writer for `u,v,w,p` following the same format family that already worked with `ncview`.
2. Store:
  - `x,y,z`
  - `time`
  - `istep`
  - `u,v,w,p`
3. Append steps along an extendible time dimension.
4. Ensure the layout is easy to target from XDMF.
5. Support:
  - raw staggered velocity components
  - optional cell-centered interpolated velocity for cleaner vector visualization
6. Support chosen output precision.
7. Where practical, keep the same netCDF-compatible organization for regular structured single-field visualization outputs as well.

### XDMF strategy

The XDMF generator can be useful in two cases:

1. stitching together separate files
2. describing the whole time series in one place

Plan recommendation:

- support both
- keep the stitched-separate-files case because it matches the current workflow already used in CaNS
- also support generating XDMF directly from a single extendible time-series file if the file layout stays regular
- optionally generate a lightweight full-domain bounding-box helper XDMF for spatial reference

### Optional `viewbox.xmf`

For subset, slice, and unstructured visualization workflows, it can be useful to also generate a small XDMF file describing only the full computational-domain box.

Recommendation:

- generate this from the Python XDMF helper path
- fetch the domain extents/geometry from `geometry.out` when that file is available, so the helper can mirror the existing binary-output metadata path
- use `viewbox.xmf` as the default output name
- allow the user to override that file name if desired
- treat it as a companion visualization aid rather than the primary data container
- an existing older domain-box helper can be used as a baseline/reference when implementing this

### Optional ADIOS2 `vtk.xml`

For ADIOS2 visualization-oriented outputs, a ParaView-friendly `vtk.xml` description should also remain an open option.

Recommendation:

- do not require the solver to emit this metadata in-core
- allow `vtk.xml` to be generated offline by helper scripts/utilities, in the same spirit as the current XDMF workflow
- treat this as a companion metadata layer for ADIOS2 BP outputs
- keep the underlying ADIOS2 array layout regular enough that offline metadata generation remains straightforward

### Opinion on whole-series XDMF

Generating one XDMF file for the whole extendible time series should be feasible and is probably desirable once the layout is fixed.

Why:

- one grid definition can be shared
- one temporal collection can point to all saved times
- each step can reference a hyperslab or step-specific view of the same file

So the intended evolution should be:

1. preserve the current stitched-separate-files workflow
2. add whole-series XDMF support for the single extendible file once the final dataset layout is agreed

This gives both:

- backwards-compatible visualization for existing separate-file workflows
- a cleaner one-file workflow for the new extendible series output

## Phase 4: Update XDMF helpers under `utils/`

Current helpers already include HDF5-oriented code, notably:

- `utils/visualize_fields/gen_xdmf_easy/write_xdmf_hdf5.py`

The desired change should be minimal:

1. teach the script to recognize the new visualization layout
2. support:
  - stitched separate files
  - extendible time-series layout from a single file
3. preserve current usage style as much as possible
4. handle both:
  - raw staggered velocity components
  - optional interpolated cell-centered velocity fields
5. optionally generate a domain-bounding `viewbox.xmf` helper, using `viewbox.xmf` as the default output name and allowing user override

This keeps XDMF as an optional companion description layer on top of the same netCDF-compatible HDF5 visualization data.

In the same spirit, helper utilities may later be extended to generate ADIOS2 `vtk.xml` companion metadata offline for ParaView-oriented workflows.

It is also worth keeping open a small ADIOS2-side inspection utility under `utils/`, in the same spirit as the existing `read_single_field_binary.py` and `read_restart_file.py` helpers. Even if the first version is simple, a lightweight reader/dumper for BP checkpoints could make validation, regression tests, and offline metadata generation much easier.

### Minimal-change strategy

- keep dataset naming stable (`fields/u`, `fields/v`, etc. or another fixed convention)
- keep grid datasets named consistently
- add only the logic needed to detect:
  - separate files that must be stitched
  - a single extendible time-series file

This minimizes churn in `utils/visualize_fields`.

## Metadata logging and `out0d`

The current `out0d` signature is:

- `out0d(fname, n, var)`

where `n` is the number of entries to write.

For new metadata logging, it would be cleaner if `out0d` also supported the natural usage:

- `out0d(fname, var)`

with `n = size(var)` inferred internally.

### Recommendation

- keep the current routine for backwards compatibility
- add a small overload or wrapper so the new code can just pass an array literal or allocatable/vector and let the routine infer its length

That makes metadata logging easier for:

- grid sizes
- skip factors
- time-series bookkeeping
- file-family descriptors

without having to explicitly thread `n` everywhere.

## Phase 5: NumPy writer evaluation

### Placement

If a Fortran NumPy writer is kept, it should live under `utils/`, not under `src/`, unless it becomes part of the production I/O stack.

Suggested locations:

- `utils/write_numpy_data/`
- or `utils/read_binary_data/` if it is tightly coupled to existing binary helpers

### MPI-I/O extension

This should not be pursued in this plan.

Recommendation:

- keep a NumPy writer as a serial convenience/debugging utility
- place it under `utils/`
- do not implement MPI-parallel `.npy` support

## Halo Exclusion Strategy by Backend

## ADIOS2

Use:

- global variable shape = interior/global field shape
- per-rank offsets = global offsets of the interior block
- memory selection = skip halo layers in local arrays

This is good for contiguous interior writes and avoids writing halos.

## HDF5

Use:

- global dataset shape = full global field shape
- file hyperslab = local interior block in global coordinates
- memory hyperslab = interior subset of halo-padded local array

This matches the current HDF5 logic and should remain the standard pattern.

## Sampling / skipping grid points

Strided outputs (every 2nd or every 4th point) should use the same structured subset API, with:

- `nmin`, `nmax` defining the rectangular region of interest
- `nskip` defining the sampling/coarsening

That means:

- arbitrary rectangular subvolumes are first-class outputs
- the non-coarsened case is `nskip = [1,1,1]`
- coarsened output is the same interface with `nskip /= [1,1,1]`

Implementation can still branch internally for performance:

- compact direct write path for `nskip = [1,1,1]`
- direct strided selection or pack-then-write for `nskip /= [1,1,1]`

Packing should therefore be treated as a performance fallback for non-unit skips, not as the defining API model.

## Other I/O modes worth keeping in mind

These are not all first-priority, but they are sensible extensions to leave in the plan.

### 1. Plane/slice outputs

- 2D fields extracted from 3D data
- useful for fast visualization and reduced storage
- natural fit for raw binary and HDF5 structured subset outputs
- conceptually a subvolume with one extent collapsed to size 1

### 2. Sampled/coarsened 3D fields

- write every `iskip/jskip/kskip` point
- especially useful for visualization output
- should be treated as a natural specialization of the generic structured subset API

### 2a. Subvolumes

- structured 3D blocks selected from the full domain
- natural fit for the generic structured subset API
- includes both full-resolution and coarsened rectangular selections

### 2b. Lines

- 1D fields extracted from 3D data
- conceptually a subvolume with two extents collapsed to size 1
- should reuse the same descriptor model, with only output-shape/metadata handling differing

### 3. Restart family management

- latest checkpoint alias
- rolling checkpoint windows
- append-vs-overwrite control

### 4. Scalar/metadata series

- a lightweight log for time, step, CFL, grid extents, compression level, skip factors, etc.
- this can be kept very simple and may naturally reuse the `out0d` path

### 4a. Unstructured point clouds

- selected points/particles with coordinates and optional fields
- should use a separate unstructured-output API
- first implementation should target packed raw binary + Python XDMF

### 5. Precision-aware output modes

- native-precision checkpoint
- reduced-precision checkpoint
- reduced-precision visualization output

These should be explicit user choices, not implicit backend behavior.

## Testing Plan

## Core correctness

1. Write/read round-trip for non-trivial values.
2. Confirm no halo contamination.
3. Test full checkpoint and single-field checkpoint.
4. Test optional metadata on/off.
5. Test optional grid on/off.

## Decomposition independence

1. Write with one decomposition.
2. Read with another decomposition.
3. Compare max errors against reference.

This is especially important for ADIOS2 restart support.

## Compression

1. Test `compression_level = 0`.
2. Test at least one non-zero level.
3. Compare file size and runtime.
4. Validate readback correctness after compression.

## Precision modes

1. Test native run precision output.
2. Test forced single-precision output from double-precision runs.
3. Test forced double-precision output from single-precision runs if supported.
4. Confirm metadata/grid paths stay consistent with chosen output precision.

## Staggered velocity visualization

1. Test raw staggered component output.
2. Test XDMF/visualization path using raw components.
3. Test optional cell-centered interpolated velocity output.
4. Compare the interpolated visualization field against the raw staggered field interpretation for sanity.

## Time-series HDF5

1. Append at least two snapshots.
2. Confirm time dimension grows correctly.
3. Confirm XDMF helper can visualize all saved steps.

## Visualization helpers

1. Run updated XDMF generator on new HDF5 outputs.
2. Open resulting `.xmf` in ParaView.
3. Confirm grid orientation and variable naming are correct.

## Integration Strategy in CaNS

There are two reasonable ways to expose the new I/O paths.

### Build-time vs runtime control

#### CPP macros

CPP macros should be used only for build-time backend availability, for example:

- `_USE_HDF5`
- `_USE_ADIOS2`

They should not be the main user-facing mechanism for selecting runtime output modes.

#### Namelist

The namelist should be the main user-facing control surface for standard workflows, including:

- backend choice
- precision choice
- compression level
- grid/meta options
- raw staggered vs interpolated visualization velocity
- subset bounds
- coarsening/skip values

### Option A: keep them out of the main CaNS feature set initially

This means:

- backend routines live in `load.f90`
- small call sites are tested manually or through isolated blocks
- no new public user-facing mode is added immediately

This is useful for stabilizing the backend before deciding how it should surface in the code.

### Option B: add a new namelist-controlled I/O mode

This would mean:

- a new I/O choice in the namelist
- the code selects backend/format from input
- easier user-facing workflows once stable

This is probably the cleaner long-term interface.

### Recommendation

Do this in two stages:

1. first implement the backends behind `load_one` / `load_all` style routines without making them a full CaNS feature yet
2. once stable, add a new namelist I/O selector and wire the supported options there

### Role of `out2d.h90` / `out3d.h90`

`out2d.h90` and `out3d.h90` should remain useful for advanced/custom extraction logic, but should reuse the common backend routines rather than re-implement low-level I/O.

Recommendation:

- standard workflows are driven by the namelist
- `out2d.h90` / `out3d.h90` remain the place for custom structured-output logic
- those wrappers should call the reusable structured subset engine in `load.f90`

### About `useful_fortran_blocks`

`utils/useful_fortran_blocks/` is good for reusable snippets, but it is probably not the right long-term home for full backend logic.

Recommendation:

- do not place the core new I/O implementation there
- use it only if a tiny reusable call example or integration snippet is helpful
- if useful, keep a small baseline example there showing a particle-tracking-style Fortran subroutine for unstructured XDMF-oriented output, as an extendable reference rather than production backend logic

## Recommended Order of Work

1. ADIOS2 refactor into `load_all`-style init/write/read/close
2. ADIOS2 `load_one`-style single-field mode
3. HDF5 refactor into `load_all`-style init/write/read/close
4. HDF5 `load_one`-style single-field mode
5. Extend the current `ncview`-compatible visualization writer into a time-series path
6. Add the generic structured subset engine in `load.f90`
7. Add line/plane/volume/coarsened wrappers in `output.f90`
8. XDMF helper updates, starting with stitched separate files
9. Add whole-series XDMF generation for the single extendible file
10. Add backend-wide output-precision controls
11. Add optional cell-centered interpolated velocity output for visualization
12. Add namelist-driven runtime selection for standard workflows
13. Add `out0d` convenience overload for metadata logging
14. Add the separate unstructured point-cloud output path
15. Place the NumPy utility under `utils/`

## Backlog / TODO List

- Implement ADIOS2 `read` path for restart.
- Refactor ADIOS2 full checkpoint into persistent `load_all`-style routines.
- Add ADIOS2 single-field `load_one`-style mode.
- Validate ADIOS2 compression levels on the local build and target machines.
- Add backend-wide output-precision selection.
- Allow visualization-oriented ADIOS2 modes such as lossy compression plus reduced precision.
- Refactor HDF5 checkpoint routines into persistent `load_all`-style routines.
- Add HDF5 single-field `load_one`-style mode.
- Extend the current `ncview`-compatible visualization writer to time series.
- Add a generic structured subset-writing engine in `load.f90`.
- Add wrappers in `output.f90` for line/plane/volume/coarsened structured outputs.
- Update `utils/visualize_fields/gen_xdmf_easy/write_xdmf_hdf5.py` for the new visualization layouts.
- Add a small ADIOS2 checkpoint inspection helper under `utils/`, analogous to the current binary restart/single-field readers.
- Document expected dataset naming conventions for XDMF generation.
- Support XDMF generation for both:
  - stitched separate files
  - a single extendible time-series file
- Keep raw staggered velocity output for fidelity.
- Add optional cell-centered interpolated velocity output for visualization.
- Keep structured subset backends focused on raw binary and HDF5 first.
- Add a separate unstructured point-cloud output path using packed raw binary + Python XDMF.
- Add an `out0d` convenience path that infers `n = size(var)` for metadata logging.
- Use CPP macros only for build-time backend availability.
- Use namelist options for runtime I/O mode selection and options.
- Place the Fortran NumPy writer under `utils/`.
- Do not implement MPI-enabled `.npy` support.

## Recommendation Summary

The cleanest path is:

- stateful backend APIs with `init/write/read/close`
- a public style that mirrors `load_one` / `load_all`
- backend workflow/engine in `load.f90`
- structured-output wrappers and policies in `output.f90`
- one schema per I/O mode
- halo exclusion embedded in backend logic
- the current `ncview`-compatible visualization format kept and extended
- XDMF updated for both stitched separate files and the whole extendible time series
- raw staggered velocity kept available, with optional cell-centered interpolated visualization output
- output precision made a first-class I/O option across backends
- a generic structured subset API for line/plane/volume/coarsened outputs
- a separate unstructured point-cloud API
- raw binary and HDF5 as the first structured subset-output backends
- `out0d` simplified for metadata logging by inferring the vector length when convenient
- namelist-driven standard workflows, with `out2d.h90` / `out3d.h90` reusing common backend routines
- NumPy writer treated as a serial utility under `utils/`, not as a core parallel I/O backend