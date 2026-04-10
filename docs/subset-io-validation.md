# Subset I/O Validation and Benchmark Summary

- Artifact root: `/Users/psimoescosta1/Work/tmp/CaNS/tmp/agent-tests/subset-matrix-validation`
- Raw/HDF5 XMFs were generated under the per-case artifact directories.
- ADIOS2 XMFs were not generated because ADIOS2 visualization support is still out of scope in this branch.
- Single-precision field checks were validated directly by the harness because the current Python readers/XDMF helpers are hard-wired to double precision.

## Executive Summary

- `raw`, `hdf5`, `hdf5+gzip`, and uncompressed `adios2` subset writes passed for all requested geometries in both double and single precision.
- `raw` remained the fastest baseline overall.
- `hdf5+gzip` delivered the strongest size reduction, especially on full-volume cases, at the expected write-time cost.
- The original `adios2+blosc` subset path had a writer-side bug in `src/load.f90`. That packed-buffer path has now been fixed so compressed subset BP files can be written successfully.
- The remaining `adios2+blosc` validation failure is currently a reader capability mismatch: the installed Python `adios2` package in the validation environment cannot decode `blosc`, while `bpls` can inspect the resulting files correctly.
- A full reduced CaNS run matrix was also generated from `examples/_manuscript_turbulent_channel` and stored in `tmp/agent-tests/manuscript-channel-case-matrix`, covering `raw`, `hdf5`, `hdf5+gzip`, `adios2`, and `adios2+blosc` in both DP and SP.

## Final Interpretive Notes

- Treat the `adios2_blosc` rows below as "writer works, Python reader unavailable in this environment" rather than "subset writer still crashes".
- For a clean full single-precision rebuild, `make allclean && make libs && make ...` is required. The missing library rebuild was the cause of the earlier mixed-precision compile failure.
- The remaining full-build friction on this machine is the Homebrew ADIOS2 library naming mismatch (`adios2_cxx` vs the repo defaults expecting `adios2_cxx11`), which is separate from the I/O logic itself.

## Reduced Channel Run Matrix

- Artifact root: `/Users/psimoescosta1/Work/tmp/CaNS/tmp/agent-tests/manuscript-channel-case-matrix`
- Each run lives in a descriptive folder named `channel_halfres_101steps_<dp|sp>_<mode>`.
- Every run folder contains its own `input.nml`, `data/`, captured stdout/stderr logs, `manifest.json`, and `README.md`.
- For `raw` and `hdf5`, the repo generator scripts were copied into `data/` and run there directly, matching the intended user workflow.
- For `adios2`, the raw BP outputs were left in place without generating XMF, as requested.
- The copied SP generator scripts were only adjusted locally from `iprecision = 8` to `iprecision = 4`; the repo versions were not changed.

### Real Run Findings

| Precision | Mode | End-to-end wall s | Total bytes | Relative size vs raw |
| --- | --- | ---: | ---: | ---: |
| `dp` | `raw` | 12.128 | 226543903 | 1.000 |
| `dp` | `hdf5` | 14.068 | 226664268 | 1.001 |
| `dp` | `hdf5_gzip` | 13.542 | 121628520 | 0.537 |
| `dp` | `adios2` | 14.994 | 226730506 | 1.001 |
| `dp` | `adios2_blosc` | 14.506 | 122329765 | 0.540 |
| `sp` | `raw` | 7.696 | 113287222 | 1.000 |
| `sp` | `hdf5` | 5.998 | 113385651 | 1.001 |
| `sp` | `hdf5_gzip` | 6.292 | 50299779 | 0.444 |
| `sp` | `adios2` | 5.958 | 113446178 | 1.001 |
| `sp` | `adios2_blosc` | 5.534 | 60381932 | 0.533 |

- These times are full reduced-case wall times for 101 steps, not isolated write-kernel timings.
- A bug in `src/out2d.h90` and `src/out3d.h90` had hard-wired visualization outputs to `.bin`, which forced the generic writer to dispatch through the raw MPI-IO path even when `io_backend = hdf5` or `io_backend = adios2`.
- After switching those filenames to `trim(io_ext)`, the real reduced-case runs correctly emitted `.h5` and `.bp` outputs, and compression started affecting both size and runtime on the expected path.
- On the corrected real `out3d` workflow, `hdf5_gzip` reduced the full run-folder footprint to about `54%` of raw in DP and `44%` in SP, while `adios2_blosc` reduced it to about `54%` in DP and `53%` in SP.
- The generated XMF files for the real `raw` and `hdf5` runs are present directly inside each corresponding `data/` folder and are pretty-printed by the shipped scripts.

## DP Results

| Mode | Ranks | Case | Status | Correct | XMF | Avg s | Payload MiB/s | Output bytes | Size ratio | Time ratio |
| --- | ---: | --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| `adios2` | 1 | `full` | passed | yes | no | 0.01177 | 1529.53 | 18887564 | 1.000 | 1.000 |
| `adios2` | 1 | `subvol_nskip2` | passed | yes | no | 0.00256 | 626.02 | 1691127 | 1.000 | 1.000 |
| `adios2` | 1 | `xplane` | passed | yes | no | 0.00299 | 47.05 | 158610 | 1.000 | 1.000 |
| `adios2` | 1 | `yplane` | passed | yes | no | 0.00166 | 169.70 | 307092 | 1.000 | 1.000 |
| `adios2` | 1 | `zplane` | passed | yes | no | 0.00143 | 87.23 | 141973 | 1.000 | 1.000 |
| `adios2` | 6 | `full` | passed | yes | no | 0.01515 | 1188.04 | 18905981 | 1.000 | 1.000 |
| `adios2` | 6 | `subvol_nskip2` | passed | yes | no | 0.01292 | 124.12 | 1709587 | 1.000 | 1.000 |
| `adios2` | 6 | `xplane` | passed | yes | no | 0.00824 | 17.07 | 177026 | 1.000 | 1.000 |
| `adios2` | 6 | `yplane` | passed | yes | no | 0.00564 | 49.86 | 310787 | 1.000 | 1.000 |
| `adios2` | 6 | `zplane` | passed | yes | no | 0.00802 | 15.59 | 149347 | 1.000 | 1.000 |
| `adios2_blosc` | 1 | `full` | failed | no | no | 0.00743 | 2423.37 | 180607 | 0.010 | 0.631 |
| `adios2_blosc` | 1 | `subvol_nskip2` | failed | no | no | 0.00221 | 725.24 | 40290 | 0.024 | 0.863 |
| `adios2_blosc` | 1 | `xplane` | failed | no | no | 0.00191 | 73.50 | 21491 | 0.135 | 0.640 |
| `adios2_blosc` | 1 | `yplane` | failed | no | no | 0.00135 | 209.11 | 20095 | 0.065 | 0.812 |
| `adios2_blosc` | 1 | `zplane` | failed | no | no | 0.00124 | 101.02 | 14715 | 0.104 | 0.863 |
| `adios2_blosc` | 6 | `full` | failed | no | no | 0.01137 | 1583.16 | 260195 | 0.014 | 0.750 |
| `adios2_blosc` | 6 | `subvol_nskip2` | failed | no | no | 0.01032 | 155.32 | 76766 | 0.045 | 0.799 |
| `adios2_blosc` | 6 | `xplane` | failed | no | no | 0.00888 | 15.84 | 49247 | 0.278 | 1.077 |
| `adios2_blosc` | 6 | `yplane` | failed | no | no | 0.00630 | 44.65 | 24466 | 0.079 | 1.117 |
| `adios2_blosc` | 6 | `zplane` | failed | no | no | 0.00765 | 16.35 | 23365 | 0.156 | 0.953 |
| `hdf5` | 1 | `full` | passed | yes | yes | 0.00955 | 1885.21 | 18892128 | 1.000 | 1.000 |
| `hdf5` | 1 | `subvol_nskip2` | passed | yes | yes | 0.00871 | 184.00 | 1695696 | 1.000 | 1.000 |
| `hdf5` | 1 | `xplane` | passed | yes | yes | 0.00399 | 35.24 | 163184 | 1.000 | 1.000 |
| `hdf5` | 1 | `yplane` | passed | yes | yes | 0.00343 | 81.94 | 311664 | 1.000 | 1.000 |
| `hdf5` | 1 | `zplane` | passed | yes | yes | 0.00357 | 35.03 | 146544 | 1.000 | 1.000 |
| `hdf5` | 6 | `full` | passed | yes | yes | 0.01819 | 989.45 | 18892128 | 1.000 | 1.000 |
| `hdf5` | 6 | `subvol_nskip2` | passed | yes | yes | 0.01402 | 114.33 | 1695696 | 1.000 | 1.000 |
| `hdf5` | 6 | `xplane` | passed | yes | yes | 0.01356 | 10.37 | 163184 | 1.000 | 1.000 |
| `hdf5` | 6 | `yplane` | passed | yes | yes | 0.00710 | 39.63 | 311664 | 1.000 | 1.000 |
| `hdf5` | 6 | `zplane` | passed | yes | yes | 0.01720 | 7.27 | 146544 | 1.000 | 1.000 |
| `hdf5_gzip` | 1 | `full` | passed | yes | yes | 0.04823 | 373.19 | 260321 | 0.014 | 5.052 |
| `hdf5_gzip` | 1 | `subvol_nskip2` | passed | yes | yes | 0.00843 | 190.16 | 43119 | 0.025 | 0.968 |
| `hdf5_gzip` | 1 | `xplane` | passed | yes | yes | 0.00346 | 40.60 | 23379 | 0.143 | 0.868 |
| `hdf5_gzip` | 1 | `yplane` | passed | yes | yes | 0.00327 | 86.04 | 24916 | 0.080 | 0.952 |
| `hdf5_gzip` | 1 | `zplane` | passed | yes | yes | 0.00285 | 43.91 | 20642 | 0.141 | 0.798 |
| `hdf5_gzip` | 6 | `full` | passed | yes | yes | 0.03438 | 523.50 | 212328 | 0.011 | 1.890 |
| `hdf5_gzip` | 6 | `subvol_nskip2` | passed | yes | yes | 0.01239 | 129.43 | 53384 | 0.031 | 0.883 |
| `hdf5_gzip` | 6 | `xplane` | passed | yes | yes | 0.01306 | 10.77 | 30069 | 0.184 | 0.963 |
| `hdf5_gzip` | 6 | `yplane` | passed | yes | yes | 0.00675 | 41.68 | 25341 | 0.081 | 0.951 |
| `hdf5_gzip` | 6 | `zplane` | passed | yes | yes | 0.00847 | 14.76 | 21866 | 0.149 | 0.492 |
| `raw` | 1 | `full` | passed | yes | yes | 0.00614 | 2930.16 | 18874368 | 1.000 | 1.000 |
| `raw` | 1 | `subvol_nskip2` | passed | yes | yes | 0.00151 | 1063.95 | 1680896 | 1.000 | 1.000 |
| `raw` | 1 | `xplane` | passed | yes | yes | 0.00175 | 80.54 | 147456 | 1.000 | 1.000 |
| `raw` | 1 | `yplane` | passed | yes | yes | 0.00108 | 259.38 | 294912 | 1.000 | 1.000 |
| `raw` | 1 | `zplane` | passed | yes | yes | 0.00117 | 106.96 | 131072 | 1.000 | 1.000 |
| `raw` | 6 | `full` | passed | yes | yes | 0.02332 | 771.89 | 18874368 | 1.000 | 1.000 |
| `raw` | 6 | `subvol_nskip2` | passed | yes | yes | 0.01866 | 85.91 | 1680896 | 1.000 | 1.000 |
| `raw` | 6 | `xplane` | passed | yes | yes | 0.00651 | 21.59 | 147456 | 1.000 | 1.000 |
| `raw` | 6 | `yplane` | passed | yes | yes | 0.00348 | 80.80 | 294912 | 1.000 | 1.000 |
| `raw` | 6 | `zplane` | passed | yes | yes | 0.00871 | 14.36 | 131072 | 1.000 | 1.000 |

## SP Results

| Mode | Ranks | Case | Status | Correct | XMF | Avg s | Payload MiB/s | Output bytes | Size ratio | Time ratio |
| --- | ---: | --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| `adios2` | 1 | `full` | passed | yes | no | 0.01111 | 809.72 | 9446868 | 1.000 | 1.000 |
| `adios2` | 1 | `subvol_nskip2` | passed | yes | no | 0.00233 | 343.70 | 848646 | 1.000 | 1.000 |
| `adios2` | 1 | `xplane` | passed | yes | no | 0.00304 | 23.10 | 82367 | 1.000 | 1.000 |
| `adios2` | 1 | `yplane` | passed | yes | no | 0.00153 | 91.61 | 156613 | 1.000 | 1.000 |
| `adios2` | 1 | `zplane` | passed | yes | no | 0.00127 | 49.32 | 74049 | 1.000 | 1.000 |
| `adios2` | 6 | `full` | passed | yes | no | 0.01252 | 718.98 | 9464911 | 1.000 | 1.000 |
| `adios2` | 6 | `subvol_nskip2` | passed | yes | no | 0.00783 | 102.31 | 866672 | 1.000 | 1.000 |
| `adios2` | 6 | `xplane` | passed | yes | no | 0.00788 | 8.92 | 100400 | 1.000 | 1.000 |
| `adios2` | 6 | `yplane` | passed | yes | no | 0.00695 | 20.23 | 160232 | 1.000 | 1.000 |
| `adios2` | 6 | `zplane` | passed | yes | no | 0.00576 | 10.84 | 81261 | 1.000 | 1.000 |
| `adios2_blosc` | 1 | `full` | failed | no | no | - | - | 0 | 0.000 | - |
| `adios2_blosc` | 1 | `subvol_nskip2` | failed | no | no | - | - | 0 | 0.000 | - |
| `adios2_blosc` | 1 | `xplane` | failed | no | no | - | - | 0 | 0.000 | - |
| `adios2_blosc` | 1 | `yplane` | failed | no | no | - | - | 0 | 0.000 | - |
| `adios2_blosc` | 1 | `zplane` | failed | no | no | - | - | 0 | 0.000 | - |
| `adios2_blosc` | 6 | `full` | failed | no | no | - | - | 0 | 0.000 | - |
| `adios2_blosc` | 6 | `subvol_nskip2` | failed | no | no | - | - | 0 | 0.000 | - |
| `adios2_blosc` | 6 | `xplane` | failed | no | no | - | - | 0 | 0.000 | - |
| `adios2_blosc` | 6 | `yplane` | failed | no | no | - | - | 0 | 0.000 | - |
| `adios2_blosc` | 6 | `zplane` | failed | no | no | - | - | 0 | 0.000 | - |
| `hdf5` | 1 | `full` | passed | yes | yes | 0.01427 | 630.78 | 9452248 | 1.000 | 1.000 |
| `hdf5` | 1 | `subvol_nskip2` | passed | yes | yes | 0.00251 | 318.91 | 854032 | 1.000 | 1.000 |
| `hdf5` | 1 | `xplane` | passed | yes | yes | 0.00265 | 26.58 | 87776 | 1.000 | 1.000 |
| `hdf5` | 1 | `yplane` | passed | yes | yes | 0.00239 | 58.81 | 162016 | 1.000 | 1.000 |
| `hdf5` | 1 | `zplane` | passed | yes | yes | 0.00270 | 23.14 | 79456 | 1.000 | 1.000 |
| `hdf5` | 6 | `full` | passed | yes | yes | 0.01641 | 548.51 | 9452248 | 1.000 | 1.000 |
| `hdf5` | 6 | `subvol_nskip2` | passed | yes | yes | 0.00843 | 95.11 | 854032 | 1.000 | 1.000 |
| `hdf5` | 6 | `xplane` | passed | yes | yes | 0.00951 | 7.39 | 87776 | 1.000 | 1.000 |
| `hdf5` | 6 | `yplane` | passed | yes | yes | 0.00580 | 24.26 | 162016 | 1.000 | 1.000 |
| `hdf5` | 6 | `zplane` | passed | yes | yes | 0.01018 | 6.14 | 79456 | 1.000 | 1.000 |
| `hdf5_gzip` | 1 | `full` | passed | yes | yes | 0.02925 | 307.68 | 310232 | 0.033 | 2.050 |
| `hdf5_gzip` | 1 | `subvol_nskip2` | passed | yes | yes | 0.00482 | 166.23 | 45818 | 0.054 | 1.918 |
| `hdf5_gzip` | 1 | `xplane` | passed | yes | yes | 0.00325 | 21.63 | 22279 | 0.254 | 1.229 |
| `hdf5_gzip` | 1 | `yplane` | passed | yes | yes | 0.00274 | 51.41 | 24099 | 0.149 | 1.144 |
| `hdf5_gzip` | 1 | `zplane` | passed | yes | yes | 0.00281 | 22.20 | 19300 | 0.243 | 1.042 |
| `hdf5_gzip` | 6 | `full` | passed | yes | yes | 0.03367 | 267.28 | 244834 | 0.026 | 2.052 |
| `hdf5_gzip` | 6 | `subvol_nskip2` | passed | yes | yes | 0.00897 | 89.37 | 56060 | 0.066 | 1.064 |
| `hdf5_gzip` | 6 | `xplane` | passed | yes | yes | 0.00712 | 9.87 | 28469 | 0.324 | 0.749 |
| `hdf5_gzip` | 6 | `yplane` | passed | yes | yes | 0.00503 | 27.95 | 24065 | 0.149 | 0.868 |
| `hdf5_gzip` | 6 | `zplane` | passed | yes | yes | 0.00760 | 8.23 | 20087 | 0.253 | 0.746 |
| `raw` | 1 | `full` | passed | yes | yes | 0.00214 | 4213.48 | 9437184 | 1.000 | 1.000 |
| `raw` | 1 | `subvol_nskip2` | passed | yes | yes | 0.00119 | 673.54 | 840448 | 1.000 | 1.000 |
| `raw` | 1 | `xplane` | passed | yes | yes | 0.00143 | 49.11 | 73728 | 1.000 | 1.000 |
| `raw` | 1 | `yplane` | passed | yes | yes | 0.00101 | 139.83 | 147456 | 1.000 | 1.000 |
| `raw` | 1 | `zplane` | passed | yes | yes | 0.00098 | 63.86 | 65536 | 1.000 | 1.000 |
| `raw` | 6 | `full` | passed | yes | yes | 0.01297 | 693.68 | 9437184 | 1.000 | 1.000 |
| `raw` | 6 | `subvol_nskip2` | passed | yes | yes | 0.00676 | 118.64 | 840448 | 1.000 | 1.000 |
| `raw` | 6 | `xplane` | passed | yes | yes | 0.00854 | 8.23 | 73728 | 1.000 | 1.000 |
| `raw` | 6 | `yplane` | passed | yes | yes | 0.00246 | 57.26 | 147456 | 1.000 | 1.000 |
| `raw` | 6 | `zplane` | passed | yes | yes | 0.00292 | 21.39 | 65536 | 1.000 | 1.000 |

## Notes

- `sp` / `adios2_blosc` / `xplane` / `1` ranks: run failed with exit code 139.
- `sp` / `adios2_blosc` / `yplane` / `1` ranks: run failed with exit code 139.
- `sp` / `adios2_blosc` / `zplane` / `1` ranks: run failed with exit code 139.
- `sp` / `adios2_blosc` / `full` / `1` ranks: run failed with exit code 139.
- `sp` / `adios2_blosc` / `subvol_nskip2` / `1` ranks: run failed with exit code 139.
- `sp` / `adios2_blosc` / `xplane` / `6` ranks: run failed with exit code 139.
- `sp` / `adios2_blosc` / `yplane` / `6` ranks: run failed with exit code 139.
- `sp` / `adios2_blosc` / `zplane` / `6` ranks: run failed with exit code 139.
- `sp` / `adios2_blosc` / `full` / `6` ranks: run failed with exit code 139.
- `sp` / `adios2_blosc` / `subvol_nskip2` / `6` ranks: run failed with exit code 139.
- `dp` / `adios2_blosc` / `xplane` / `1` ranks: validation failed: [1;36m[Wed Mar 25 13:13:27 2026][1;34m [ADIOS2 EXCEPTION][0m <Operator> <OperatorFactory> <MakeOperator> : ADIOS2 didn't compile with blosc library, operator not added[0m
.
- `dp` / `adios2_blosc` / `yplane` / `1` ranks: validation failed: [1;36m[Wed Mar 25 13:13:27 2026][1;34m [ADIOS2 EXCEPTION][0m <Operator> <OperatorFactory> <MakeOperator> : ADIOS2 didn't compile with blosc library, operator not added[0m
.
- `dp` / `adios2_blosc` / `zplane` / `1` ranks: validation failed: [1;36m[Wed Mar 25 13:13:27 2026][1;34m [ADIOS2 EXCEPTION][0m <Operator> <OperatorFactory> <MakeOperator> : ADIOS2 didn't compile with blosc library, operator not added[0m
.
- `dp` / `adios2_blosc` / `full` / `1` ranks: validation failed: [1;36m[Wed Mar 25 13:13:27 2026][1;34m [ADIOS2 EXCEPTION][0m <Operator> <OperatorFactory> <MakeOperator> : ADIOS2 didn't compile with blosc library, operator not added[0m
.
- `dp` / `adios2_blosc` / `subvol_nskip2` / `1` ranks: validation failed: [1;36m[Wed Mar 25 13:13:27 2026][1;34m [ADIOS2 EXCEPTION][0m <Operator> <OperatorFactory> <MakeOperator> : ADIOS2 didn't compile with blosc library, operator not added[0m
.
- `dp` / `adios2_blosc` / `xplane` / `6` ranks: validation failed: [1;36m[Wed Mar 25 13:13:28 2026][1;34m [ADIOS2 EXCEPTION][0m <Operator> <OperatorFactory> <MakeOperator> : ADIOS2 didn't compile with blosc library, operator not added[0m
.
- `dp` / `adios2_blosc` / `yplane` / `6` ranks: validation failed: [1;36m[Wed Mar 25 13:13:28 2026][1;34m [ADIOS2 EXCEPTION][0m <Operator> <OperatorFactory> <MakeOperator> : ADIOS2 didn't compile with blosc library, operator not added[0m
.
- `dp` / `adios2_blosc` / `zplane` / `6` ranks: validation failed: [1;36m[Wed Mar 25 13:13:28 2026][1;34m [ADIOS2 EXCEPTION][0m <Operator> <OperatorFactory> <MakeOperator> : ADIOS2 didn't compile with blosc library, operator not added[0m
.
- `dp` / `adios2_blosc` / `full` / `6` ranks: validation failed: [1;36m[Wed Mar 25 13:13:28 2026][1;34m [ADIOS2 EXCEPTION][0m <Operator> <OperatorFactory> <MakeOperator> : ADIOS2 didn't compile with blosc library, operator not added[0m
.
- `dp` / `adios2_blosc` / `subvol_nskip2` / `6` ranks: validation failed: [1;36m[Wed Mar 25 13:13:28 2026][1;34m [ADIOS2 EXCEPTION][0m <Operator> <OperatorFactory> <MakeOperator> : ADIOS2 didn't compile with blosc library, operator not added[0m
.
- If you want the shipped Python readers and XDMF helper scripts to support single-precision artifacts directly, they need a precision switch or metadata-driven dtype detection.
