p3dfft
------

This test program performs the following tasks:

* It crosschecks 2DECOMP&FFT results against P3DFFT results.
* It compares the performance of the two codes. 

**How to set up this test?**

Due to external dependency, this test has to be built separately after building 2DECOMP&FFT.

First, download P3DFFT version 2.7.9, the most recent 2.x release. Install P3DFFT at a directory denoted as `$P3DFFT_HOME` and set this path properly in the `Makefile`. To build P3DFFT :
```
FC=mpif90 CC=mpicc LDFLAGS="-lm" ./configure --prefix=$HOME/software/build/p3dfft-2.7.9-dimsc \
--enable-gnu --enable-fftw --with-fftw=$HOME/software/build/fftw-3.3.9 --enable-openmpi --enable-dimsc
make
make install
```
This instruction is on a workstation with gcc 8.4.1 and OpenMPI 4.0.5. Adapt this accordlingly. Note that P3DFFT is built with `-DDIMS_C` flag (enabling same decomposition as 2DECOMP&FFT for a fair comparison). P3DFFT needs to link to *libm* but somehow its build system doesn't handle this correctly, thus the added `LDFLAGS` setting being required.

Both P3DFFT and 2DECOMP&FFT are built in double precision. P3DFFT is built against FFTW (provide its path in `Makefile`); 2DECOMP&FFT can use any FFT engine.

**What to expect:**

* Results from P3DFFT and 2DECOMP&FFT should be almost identical, even when different FFT engines are used.
* Each library should recover its input after one forward and one backward transform with errors reported close to machine accuracy.
