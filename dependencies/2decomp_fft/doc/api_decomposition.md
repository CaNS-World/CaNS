## 2D Pencil Decomposition API

This page explains the key public interfaces of the 2D decomposition library. After reading this section, users should be able to easily build applications using this domain decomposition strategy. The library interface is designed to be very simple. One can refer to the [sample applications](samples.md) for a quick start.

The 2D Pencil Decomposition API is defined in one Fortran module which should be used by applications:

```
      use decomp_2d
```

#### Global Variables

Following is the list of global variables defined by the library that can be used in applications. Obviously these names should not be redefined in applications to avoid conflict. Also note that some variables contain duplicate or redundant information just to simplify the programming.

- *mytype* - use this variable to define the KIND of floating-point data, e.g. `real(mytype) :: var` or `complex(mytype) :: cvar`. This makes it easy to switch between single precision and double precision (more details).
- *real_type, complex_type* - these are the proper MPI datatypes to use (for real and complex numbers, respectively) if applications need to invoke MPI library routines directly.
- *nx_global, ny_global, nz_global* - size of the global domain.
- *nproc* - the total number of MPI processes.
- *nrank* - the rank of the current MPI process.
- *xsize(i), ysize(i), zsize(i), i=1,2,3* - sizes of the sub-domains held by the current process. The first letter refers to the pencil orientation and the three 1D array elements contain the sub-domain sizes in X, Y and Z directions, respectively. In a 2D pencil decomposition, there is always one dimension which completely resides in local memory. So by definition *xsize(1)==nx_global*, *ysize(2)==ny_global* and *zsize(3)==nz_global*.
- *xstart(i), ystart(i), zstart(i), xend(i), yend(i), zend(i), i=1,2,3* - the starting and ending indices for each sub-domain, as in the global coordinate system. Obviously, it can be seen that *xsize(i)=xend(i)-xstart(i)+1*. It may be convenient for certain applications to use global coordinate. For example when extracting a 2D plane from a 3D dataset, it is easier to know which process owns the plane if global index is used.

It is recommended that an application define its major data structures using these global variables in its main program immediately after initialising the decomposition library. Using allocatable arrays is preferred, as shown in the following examples:

```
      ! allocate a X-pencil array
      allocate(ux(xsize(1), xsize(2), xsize(3))) 
      ! allocate a Y-pencil array using global indices
      allocate(uy(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
```

There are also a number of [utility routines](#utils) to help allocate 3D arrays.

For debugging or otherwise writing information to the standard output, the following is recommended:

```
      if (nrank==0) then
         write(*,*) ! list of variables ......
      end if
```

#### Basic 2D Decomposition API

All the global variables defined above are initialised by:
```
     call decomp_2d_init(nx, ny, nz, P_row, P_col)
```
where *nx*, *ny* and *nz* are the size of 3D global data to be distributed over a 2D processor grid *P_row \* P_col*. Note that these dimensions do not need to be divisible by *P_row* or *P_col*, i.e. the library can handle non-evenly distributed data. However, choose the numbers smartly to avoid significant load-imbalance, which could lead to poor parallel performance. Also note that constraints imposed by the algorithm are: *P_row <= min(nx, ny)* and *P_col <= min(ny, nz)*.

An optional parameter may be passed to this initialisation routine:

```
      call decomp_2d_init(nx, ny, nz, P_row, P_col, periodic_bc)
```

Here periodic_bc is a 1D array containing 3 logical values that specify whether periodic boundary condition should apply in certain dimensions. Note this is only applicable if [halo-cell communication](api_halo.md) is to be used.

To obtain optimal performance, it is very important to choose the best possible processor grid *P_row \* P_col*. Please follow [this note](pgrid.md) for more instructions.

A key part of this library is a set of communication routines that actually perform the data transpositions. As mentioned, one needs to perform 4 global transpositions to go through all 3 pencil orientations. Correspondingly, the library provides 4 communication routines:

```
      call transpose_x_to_y(in, out)
      call transpose_y_to_z(in, out)
      call transpose_z_to_y(in, out)
      call transpose_y_to_x(in, out)
```
The input array *in* and output array *out* should have been defined to contain distributed data for the correct pencil orientations.

Note that the library is written using Fortran's generic interface so different data types are supported without user intervention. That means in and out above can be either real arrays or complex arrays, the latter being useful for FFT-type of applications.

As seen, the communication details are packed within black boxes. From a user's perspective, it is not necessary to understand the internal logic of these transposition routines. From the developer's perspective, implementations can be changed without breaking user codes.

It is, however, noted that the communication routines are expensive, especially when running with large distributed datasets over a large number of processors. So applications should try to minimize the number of calls to them by adjusting the algorithms in use, even sometimes at the cost of duplicating computations.

Finally, before exit, applications should clean up the memory by:
```
      call decomp_2d_finalize
```

#### Advanced 2D Decomposition API

While the basic decomposition API is very user-friendly, there may be situations in which applications need to handle more complex data structures. There are quite a few examples:

- While performing real-to-complex FFTs, applications need to store both the real input (say, of global size *nx\*ny\*nz*) and the corresponding complex output (of smaller global size - such as *(nx/2+1)\*ny\*nz* - where roughly half the output is dropped due to conjugate symmetry).
- Many CFD applications use a staggered mesh system which requires different storage for global quantities (e.g. cell-centred vs. cell-interface variables).
- In applications using spectral method, for anti-aliasing purpose, it is a common practice to enlarge the spatial domain before applying the Fourier transforms.

In all these examples, there are ***multiple global sizes*** and applications need to distributed different global datasets as 2D pencils. 2DECOMP&FFT provides a powerful and flexible programming interface to handle this scenario:
```
      TYPE(DECOMP_INFO) :: decomp
      call decomp_info_init(n1, n2, n3, decomp)
```
Here `decomp` is an instance of Fortran derived data type `DECOMP_INFO` encapsulating the 2D decomposition information associated with one particular global size *n1\*n2\*n3*. The decomposition object can be initialised using the `decomp_info_init` routine. This object then can be passed to the communication routines defined in the basic interface as a third parameter. For example:
```
     call transpose_x_to_y(in, out, decomp)
```
The `decomp` parameter is optional and can be safely ignored by users using the basic interface only. When it is in use, the global transposition will be applied to data set of the associated global size instead of the default global size *nx\*ny\*nz*. The decomp object can also be used in many other routines where arbitrary global data size is required.

The following code snippet demonstrates the use of the advanced API in 2DECOMP&FFT's FFT implementations:

```
      ! Two decomposition objects used by physical- and spectral-space variables, respectively.
      TYPE(DECOMP_INFO), save :: ph  ! physical space
      TYPE(DECOMP_INFO), save :: sp  ! spectral space
      ! ......
      ! Initialise the decomposition objects. 
      call decomp_info_init(nx, ny, nz, ph)
      call decomp_info_init(nx/2+1, ny, nz, sp)  ! spectral space of roughly half size
```
The `DECOMP_INFO` type defines the following public components that can be used by applications:

- *xsz(i), ysz(i), zsz(i), i=1,2,3*
- *xst(i), yst(i), zst(i), i=1,2,3*
- *xen(i), yen(i), zen(i), i=1,2,3*

These describe the sizes, starting indices and ending indices of the sub-domains held by the current process - exactly the same meaning as the *size/start/end* variables defined earlier, except that they apply to the global size associated with one particular decomposition object, rather than the default global size. Here is more sample code:

```
      real(mytype), allocatable, dimension(:,:,:) : wk1
      complex(mytype), allocatable, dimension(:,:,:) : wk2

      ! allocate memory space
      allocate(wk1(ph%xsz(1),ph%xsz(2),ph%xsz(3)))  ! X-pencil, for physical-space data
      allocate(wk2(sp%zsz(1),sp%zsz(2),sp%zsz(3)))  ! Z-pencil, for spectral-space data

      ! to loop through the allocated array in physical space
      do k = 1, sp%xsz(3)
        do j = 1, sp%xsz(2)
          do i = 1, sp%xsz(1)
            wk1(i,j,k) = ......
          end do
        end do
      end do
```

The `DECOMP_INFO` type also defines other components that are mainly used internally in the library, such as parameters describing the memory addresses, sizes and displacements of the MPI buffers to be used by the communication code.

If a decomposition object is no longer needed in an application, one can release its associated memory by:
```
      call decomp_info_finalize(decomp)
```

#### Utility Routines to Help Allocate 3D Arrays<a name="utils"></a>

If a large number of distributed arrays are to be used in an application, using the global variables defined above for memory allocation can become a bit awkward. For this reason, a number of utility routines have been introduced in version 1.5 of the library to simplify such operation. In the simplest form, these utility routines look like:

```
      call alloc_x(var)
      call alloc_y(var)
      call alloc_z(var)
```
The three routines are used to create X-pencil, Y-pencil and Z-pencil distributed arrays, respectively. Here var is a three-dimensional real or complex allocatable array declared in the calling program and to be allocated in the subroutines. The following two lines of code are equivalent:
```
      call alloc_x(var)
      allocate(var(xsize(1), xsize(2), xsize(3)))
```
These routines also accept additional parameters for more general cases:

```
      call alloc_x(var, decomp, global)
```
Here the optional parameter `decomp` of type `DECOMP_INFO`, if present, specifies the global size of the array to be allocated; the optional logical parameter `global`, if set to `.true.`, indicates that the array should be allocated using the global coordinates. So the above statement is equivalent to:
```
      allocate(var(decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)))   ! if global==.false.
      allocate(var(decomp%xst(1):decomp%xen(1), decomp%xst(2):decomp%xen(2), &
                   decomp%xst(3):decomp%xen(3)))                   ! if global==.true.
```

As can be seen, these utility routines are in a much more compact form. They will also report an error if the memory allocation is somehow not successful.