# I/O Coding Style Notes

## Purpose

This note captures a first pass at style rules for the upcoming I/O work so that new code stays consistent with the current CaNS Fortran codebase, especially `src/load.f90` and `src/output.f90`.

It is not meant to introduce a new coding style. The goal is to make the new routines look like natural extensions of the existing code.

## General principle

Prefer consistency with the current code over stylistic reinvention.

In practice:

- follow the existing module structure and routine naming patterns
- preserve the current formatting style where possible
- keep new abstractions lightweight
- avoid introducing a second error-handling or logging style just for the new I/O code

## Module-level style

Follow the current module layout:

1. SPDX header comment
2. `module ...`
3. `use ...` statements
4. `implicit none`
5. `private`
6. `public ...`
7. `contains`

Guidelines:

- keep `public` exports explicit
- gate backend-specific public routines with the same preprocessor style already used in `load.f90`
- avoid exposing internal helpers unless they are part of the real module API

## Routine organization

New backend routines should follow the same overall style already present in the code:

- short comment header before the routine when helpful
- `implicit none` at the top of each routine
- declarations grouped near the top
- executable logic after declarations

For the new backend work, prefer routine families such as:

- `*_init(...)`
- `*_write(...)`
- `*_read(...)`
- `*_close(...)`

but keep the naming aligned with the existing CaNS vocabulary, especially the `load_one` / `load_all` style already used in `src/load.f90`.

## Naming

Prefer names that match current conventions:

- module names such as `mod_load`, `mod_output`
- short, descriptive routine names such as `load_all`, `load_one`, `io_field`
- lowercase identifiers
- concise local names such as `ng`, `nh`, `lo`, `hi`, `disp`, `fh`, `dset`

Recommendations:

- keep argument names consistent across backends when the concepts are the same
- use the same names for global sizes, halo sizes, local extents, time, and step counters across ADIOS2, HDF5, and raw MPI-I/O paths
- prefer descriptive backend-specific names only where needed, such as `engine`, `io`, `filespace`, `memspace`, `plist_id`

## Arguments and interfaces

Match existing interface style as closely as practical:

- put the main operation selector or main identifying arguments first
- keep array shape arguments explicit
- use `intent(in)`, `intent(out)`, and `intent(inout)` consistently
- keep assumed-shape and lower-bound-aware array declarations in the same style as existing routines

Recommendations:

- preserve the existing argument order where extending an existing conceptual API
- avoid unnecessary optional arguments in low-level kernels
- use wrappers when a friendlier public call shape is needed

## Formatting

Keep formatting close to the existing files:

- two-space indentation inside modules and routines
- aligned declaration style where it improves readability
- short continuation blocks using `&`
- compact `select case` and `if` structures
- in Fortran comment sections, prefer the existing block style with standalone `!` lines rather than visually blank separator lines

Avoid reformatting surrounding legacy code just to impose a new visual style.

## Comments

The current code uses brief, practical comments. Keep that tone.

Comment layout should also match the current codebase closely.

Preferred style:

- use comment blocks such as:
- `!`
- `! comment`
- `!`
- do not leave empty separator lines inside comment blocks unless they also begin with `!`

Use comments for:

- short routine purpose summaries
- non-obvious MPI/HDF5/ADIOS2 behavior
- layout assumptions
- decomposition or halo-handling details

Avoid comments that merely restate the next line of code.

## Error handling

Stay consistent with the current codebase rather than introducing a separate framework.

Guidelines:

- continue using `ierr` in the same style as the rest of the code
- keep collective failure behavior explicit in parallel code
- use the project's current fatal-error pattern where a condition is truly unrecoverable
- keep rank-0-only reporting consistent with existing practice

For now, consistency is more important than designing a perfect new abstraction for backend errors.

## Backend-specific code

Keep backend code localized.

Recommendations:

- keep low-level backend mechanics in `src/load.f90`
- keep visualization and user-facing wrappers in `src/output.f90`
- keep backend-specific helper state private unless it must be shared
- use preprocessor guards in the same style already present for optional backends

The new code should not spread HDF5- or ADIOS2-specific details through unrelated parts of the code.

## Structured subset I/O style

For the planned structured subset work:

- keep one common descriptor vocabulary, centered on `nmin`, `nmax`, and `nskip`
- treat lines, planes, and subvolumes as the same conceptual family
- allow internal specialization without changing the public API

That means it is fine to have:

- one direct generic path for MPI-I/O or HDF5 selections
- specialized packed kernels for lines, planes, subvolumes, or coarsened outputs

The style goal is a clean interface outside and pragmatic specialization inside.

## Saved state

The codebase already leans toward direct procedural routines. For new backend state:

- keep saved state small and purposeful
- use one state object per logical schema or mode when needed
- avoid hidden cross-coupling between unrelated output modes
- make setup and teardown explicit through `init` and `close` routines

Saved state is fine, but it should stay easy to reason about.

## What to avoid

- introducing a brand-new naming system only for the new I/O stack
- mixing user-facing workflow logic deeply into `main.f90`
- duplicating the same backend logic across `main.f90`, `output.f90`, and `load.f90`
- exposing too many internal helpers as public entry points
- over-generalizing too early when a small specialized helper is clearer

## Initial review checklist

When adding a new I/O routine, check:

1. Does it look like it belongs next to the existing routines in `src/load.f90` or `src/output.f90`?
2. Are the argument names and ordering consistent with nearby code?
3. Is the backend-specific logic localized?
4. Are halo handling, decomposition assumptions, and metadata layout explicit?
5. Did we avoid creating a second coding style by accident?

## Status

This is a starting note only.

Before substantial implementation begins, this document should be refined against the exact patterns already present in the target files we will edit most heavily.
