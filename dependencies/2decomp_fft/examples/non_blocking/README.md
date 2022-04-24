non_blocking
------------

This test contains two sample applications to compute multiple independent FFTs. The first application `blocking.f90` uses the standard blocking version of MPI communication code to transpose the data among different stages of the computation. The second application `non_blocking.f90` performs the same computation using the non-blocking communication routines supplied by 2DECOMP&FFT. 

These two applications are using FFTW APIs directly. Please compile them separately using the `Makefile` in this directory after building 2DECOMP&FFT with the FFTW engine.

Non-blocking collective communication is part of MPI 3 standard and it is now widely supported. Earlier users of 2DECOMP&FFT may remember the use of libNBC (http://www.unixer.de/research/nbcoll/libnbc/), a library implementing non-blocking MPI collectives (such as IALLTOALL) with existing MPI 1 functions. libNBC is now obsolete and reference to it has been removed from the source codes. 

To demonstrate the idea of overlap communication and computation, the 3D FFT is implemented using loops of 1D FFTs (without using the advanced interface of FFTW) so that MPI_TEST calls can be easily inserted in the computational part of the code. This is required because the communication has to be explicitly progressed when running on the same thread as the computation.

The two applications should produce identical results.

End users are responsible for identifying opportunities in their applications to overlap communication with computation.
