## API for Parallel Three-dimensional FFTs

#### Initialisation

To use the FFT programming interface, first of all, one additional Fortran module has to be used:
```
      use decomp_2d_fft
```

The FFT interface is built on top of the 2D decomposition library which, naturally, needs to be initialised first:
```
      call decomp_2d_init(nx, ny, nz, P_row, P_col)
```
where *nx\*ny\*nz* is the 3D domain size and *P_row \* P_col* is the 2D processor grid. 

Next one needs to initialise the FFT interface by:
```
      call decomp_2d_fft_init
```

The initialisation routine handles planing for the underlying FFT engine (if supported) and defines global data structures (such as temporary work spaces) for the computations. By default, it assumes that physical-space data is distributed in X-pencil format. The corresponding spectral-space data is stored in transposed Z-pencil format after the FFT. To give applications more flexibility, the library also supports the opposite direction, if an optional parameter is passed to the initialisation routine:
```
      call decomp_2d_fft_init(PHYSICAL_IN_Z)
```

Physical-space data in Y-pencil is not an option as it would require additional expensive transpositions which does not make economical sense. There is a third and the most flexible form of the initialisation routine:
```
      call decomp_2d_fft_init(pencil, n1, n2, n3)
```
It allows applications to initialise FFT computations using an arbitrary problem size *n1\*n2\*n3*, which can be different from the main domain size *nx\*ny\*nz*.

#### Complex-to-complex Transforms

The library supports three-dimensional FFTs whose data is distributed as 2D pencils and stored in ordinary ijk-ordered 3D arrays across processors. For complex-to-complex (c2c) FFTs, the user interface is:
```
      call decomp_2d_fft_3d(in, out, direction)
```
where direction can be either `DECOMP_2D_FFT_FORWARD` (-1) for forward transforms, or `DECOMP_2D_FFT_BACKWARD` (1) for backward transforms. The input array `in` and output array `out` are both complex and have to be either a X-pencil/Z-pencil combination or vice versa, depending on the direction of FFT and how the FFT interface is initialised earlier (`PHYSICAL_IN_X`, the optional default, or `PHYSICAL_IN_Z`).

#### Real-to-complex & Complex-to-Real Transforms

While the c2c interface is already in the simplest possible form, for r2c and c2r transforms, the 3D FFT interface can be used in a more compact form:
```
      call decomp_2d_fft_3d(in, out)
```
Here if `in` is a real array and `out` a complex array, then a forward FFT is implied. Similarly a backward FFT is computed if `in` is a complex array and `out` a real array.

When real input is involved, the corresponding complex output satisfies so-called ***Hermitian redundancy*** - i.e. some output values are complex conjugates of others. Taking advantage of this, FFT algorithms can normally compute r2c and c2r transforms twice as fast as c2c transforms while only using about half of the memory. Unfortunately, the price to pay is that application's data structures have to become slightly more complex. For a 3D real input data set of size nx*ny*nz, the complex output can be held in an array of size *(nx/2+1)\*ny\*nz*, with the first dimension being cut roughly in half<a href="#note1" id="note1ref"><sup>1</sup></a>. Applications can either rely on the advanced interface described in the [decomposition API](api_decomposition.md), or use the following utility routine to distribute the complex output as 2D pencils:
```
      call decomp_2d_fft_get_size(start,end,size)
```

Here all three arguments are 1D array of three elements, returning to the caller the starting index, ending index and size of the sub-domain held by the current processor - information very similar to the *start/end/size* variables defined in the main decomposition library.

Note that the complex output arrays obtained from X-pencil and Z-pencil input do not contain identical information (see the output of the fft_test_r2c [sample application](samples.md)). However, if 'Hermitian redundancy' is taken into account, no physical information is lost and the real input can be fully recovered through the corresponding inverse FFT from either complex array.

Also note that 2DECOMP&FFT does not scale the transforms. So a forward transform followed by a backward transform will not recover the input unless applications normalise the results by the sizes of the transforms.

#### Finalisation

Finally, to release the memory used by the FFT interface:
```
      call decomp_2d_fft_finalize
```

It is possible to re-initialise the FFT interface in the same application at the later stage after it has been finalised, if this becomes necessary.

To obtain first-hand experience on the FFT interface, users are advised to examine the [sample applications](samples.md)) distributed with the library.

<hr size="1">

<a id="note1" href="#note1ref"><sup>1</sup></a>The storage is for Fortran. In C/C++, the last dimension has to be cut in half due to different memory pattern. For Z-pencil input, the complex output is of size *nx\*ny\*(nz/2+1)* instead. Also note that the integer division is rounded down.