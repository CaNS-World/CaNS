# Updated I/O Plan

## Purpose

This document replaces the earlier plan that mixed two different directions:

- checkpoint/backend work that is still relevant
- a NetCDF-style HDF5 visualization path built around `ncview`

The `ncview` / NetCDF-style direction is no longer part of the plan.

The goal now is simpler:

- keep checkpoint I/O robust and backend-aware
- keep HDF5 snapshots self-describing in a practical grouped layout
- keep room for richer structured I/O modes such as planes, lines, subvolumes, and coarsened outputs
- treat XDMF as a later visualization layer, not as a driver of the file layout

## Current Decisions

### Backends

- MPI-I/O remains the baseline checkpoint path
- HDF5 is supported for self-describing checkpoint files
- ADIOS2 is supported for checkpoint files, including optional compression

### HDF5 checkpoint layout

HDF5 snapshots should use the grouped layout:

- `fields/<var>`
- `grid/x`, `grid/y`, `grid/z`
- `meta/time`, `meta/istep`

This layout is now preferred over any NetCDF-style flat layout.

### ADIOS2 checkpoint layout

ADIOS2 currently uses flat variable naming:

- fields: `u`, `v`, `w`, `p`, `s_###`
- metadata: `time`, `istep`
- grid: `x`, `y`, `z`

This is acceptable for now, but it should be treated as an explicit schema choice and documented.

### Visualization policy

- do not shape checkpoint HDF5 around `ncview`
- do not pursue a dedicated NetCDF-style HDF5 writer at this stage
- do not pursue an extendible time-series HDF5 writer at this stage
- leave visualization metadata generation to XDMF helpers later

## Scope That Still Matters

The useful remaining roadmap is not about `ncview`. It is about broader I/O capability and keeping the backend layer clean.

The main families of work are:

1. checkpoint/backend stabilization
2. schema cleanup and documentation
3. structured subset output
4. later XDMF support on top of stable layouts

## Immediate Cleanup

### 1. Freeze the checkpoint schemas

Document the intended on-disk conventions clearly:

- MPI-I/O file structure and metadata behavior
- HDF5 grouped snapshot structure
- ADIOS2 flat variable naming

This avoids future churn in helper scripts and validation tools.

### 2. Confirm the public API

Review and explicitly confirm the intended interface shape in `src/load.f90`:

- `load_all`: non-optional `time`, `istep`; optional `x_g`, `y_g`, `z_g`
- `load_one`: optional `time`, `istep`; optional `x_g`, `y_g`, `z_g`

The goal is to separate intended API decisions from experimental drift.

### 3. Decide whether ADIOS2 naming should remain flat

Two valid options exist:

- keep ADIOS2 flat and just document it
- make ADIOS2 naming mirror HDF5 more closely, for example with names like `fields/u`, `grid/x`, `meta/time`

This is not urgent, but it is the main remaining schema-level design question.

## Main Next Workstream: Structured Subset Output

The next substantial feature area, apart from XDMF, is structured subset I/O.

This includes:

- planes
- lines
- rectangular subvolumes
- coarsened structured outputs

These should be treated as one conceptual family rather than separate ad hoc features.

### Descriptor model

Use a compact structured descriptor built around:

- `nmin`
- `nmax`
- `nskip`

Interpretation:

- `nmin`, `nmax` define the global rectangular region
- `nskip = [1,1,1]` means full-resolution output
- `nskip /= [1,1,1]` means coarsened/sampled output

Planes and lines are just special cases:

- a plane is a subvolume with one extent collapsed
- a line is a subvolume with two extents collapsed

### Backend priority

The first structured subset-output backends should be:

- raw binary
- HDF5

ADIOS2 subset output can remain a later extension.

### Backend rules

- halo exclusion must stay inside the backend implementation
- file layouts must remain decomposition-independent
- packed writes are acceptable as an internal optimization, not as the public API

## Suggested Architecture

### `src/load.f90`

Keep low-level backend mechanics here:

- schema setup
- selections/hyperslabs
- packing when needed
- metadata writing
- backend-specific read/write calls

### `src/output.f90`

Use this for policy-style wrappers that request:

- plane outputs
- line outputs
- subvolume outputs
- coarsened outputs

Those wrappers should call reusable lower-level routines from `src/load.f90`.

## Optional Future Work

These remain reasonable but are not required immediately.

### Precision-aware I/O

Potential future option:

- write in run precision by default
- optionally write in single or double precision independent of the run precision

Useful for:

- smaller visualization-oriented outputs
- storage tradeoff experiments
- explicit reduced-precision checkpoint studies

### Unstructured point output

Keep this as a separate future path:

- packed raw binary
- later Python-side XDMF generation

This should not be mixed into the structured subset API.

### Metadata helpers

A small convenience overload for `out0d` could still be useful for metadata logging, but it is not on the critical path.

## Deferred On Purpose

These are explicitly not part of the current plan:

- NetCDF-style HDF5 layout work
- `ncview` compatibility work
- HDF5 dimension scales for checkpoint files
- extendible time-series HDF5 files
- whole-series XDMF for a single appendable HDF5 file

## Later Visualization Phase

XDMF should be treated as a later companion layer built on top of stable data layouts.

When that work starts, the first target should be simple and practical:

- generate XDMF for existing snapshot/checkpoint outputs
- support the grouped HDF5 layout directly
- keep the file layout driven by checkpoint/data needs, not by a visualization tool

## Recommended Order

1. Freeze and document the current checkpoint schemas.
2. Confirm the intended `load_one` / `load_all` public API.
3. Decide whether ADIOS2 naming stays flat.
4. Design the structured subset descriptor using `nmin`, `nmax`, `nskip`.
5. Implement subvolume output first.
6. Extend that to planes and lines as special cases.
7. Add coarsened structured output support.
8. Add backend wrappers in `src/output.f90`.
9. Only after that, update or extend the XDMF helpers.

## Summary

The plan is no longer centered on making HDF5 look like NetCDF.

The plan is now:

- robust checkpoint backends
- stable and documented schemas
- a reusable structured subset-output capability
- later XDMF generation on top of those stable outputs
