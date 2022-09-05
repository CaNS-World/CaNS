# `CaNS 2.0` is finally released! :tada:

_This is the most significant revision of our toolkit so far._

Co-authored by Pedro Costa, Massimiliano Fatica, and Josh Romero.

## Summary

This release marks the ending of a fresh porting effort for massively parallel simulations on modern architectures, **_from one to thousands of GPUs_ with a focus on performance while ensuring a flexible and sustainable implementation that is easy to extend for more complex physical problems**. We used OpenACC directives to accelerate loops and for host/device data transfer, interoperated with NVIDIA's *cuFFT* and the new *cuDecomp* domain decomposition library.

*cuDecomp* is the heart of the multi-GPU implementation, ensuring the solver's performance by bringing a novel, hardware-adaptive parallelization of the transposes in the Poisson/Helmholtz solver, and of the halo-exchange operations.

Although quite performant, the implementation is also flexible, allowing for an easy change of solver profiles, such as X-aligned default pencils, which are optimal for a fully explicit time integration, or Z-aligned default pencils, which are optimal for a Z-implicit time integration for wall flows.

Finally, another noteworthy (optional) feature is `CaNS`' new mixed-precision mode, where the pressure Poisson equation is solved in lower precision. **This mode makes a huge difference in performance for many-GPU calculations across multiple nodes**.

In addition to these big-picture changes, there have been many impactful changes that make the solver more versatile and robust. All relevant changes are summarized below.

## Changes:
* **GPU acceleration** using OpenACC directives for loops and data movement, which is interfaced with CUDA whenever needed
* **Hardware-adaptive multi-GPU implementation** using the [**_cuDecomp_**](https://github.com/NVIDIA/cuDecomp) library for transposes (_seven_ possible communication backends) and halo exchanges (_five_ possible communication backends), with different flavors of MPI, NCCL and NVSHMEM implementations
* **Lean memory footprint on GPUs**, which can be made even leaner by exploiting *cuDecomp*'s in-place transposes
* **Mixed-precision mode** implemented on both CPUs and GPUs
* Hybrid MPI-OpenMP parallelization is still supported on CPUs
* Any default pencil orientation is supported, on both CPUs and GPUs
* A _fast-kernel_ mode is used by default to speed up the calculation of the prediction velocity, on both CPUs and GPUs
* The [*2DECOMP*](https://github.com/xcompact3d/2decomp-fft) library is still used for the many-CPU parallelization of the Poisson solver, and some of the parallel data I/O
* Build process made much simpler and more robust, with the dependencies determined automatically
* Refactoring of the FFT-based Fourier, cosine, and sine transforms on GPUs, together with the Gauss elimination kernels, with improvements both in terms of speed and maintainability
* Support for uneven decompositions and odd numbers along any direction; perhaps surprisingly, at times setups with odd numbers near the desired resolution may result in a more efficient FFT computation
* External domain decomposition libraries, *cuDecomp* and *2DECOMP*, loaded as Submodules
* **_Many_ changes for improved performance and robustness**, with a focus on minimizing the memory footprint and computation intensity while keeping the tool versatile

## Acknowledgements

`CaNS 2.0` has been tested in several GPU-accelerated systems such as Marconi 100, Meluxina, Perlmutter, Selene, Summit and Vega. We acknowledge the support from [CoE RAISE](https://www.coe-raise.eu), [NERSC](https://www.nersc.gov) and [EuroHPC](https://eurohpc-ju.europa.eu), which enabled thorough testing of `CaNS 2.0` in these state-of-the-art supercomputers.
