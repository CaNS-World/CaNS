## API for Halo-cell Support

While most of the communications using the 2D decomposition are via the global transposition calls, it may become necessary for neighbouring blocks to exchange data explicitly. One such scenario is in CFD applications performing large-eddy simulations (LES). While most spatial derivatives are computed using the implicit formulation to achieve a high order of accuracy, some derivatives may be evaluated quickly using local stencils and explicit formulae, such as those used by sub-grid scale models (a model by definition does not require higher-order of accuracy).

The halo-cell support API provides data structures and nearest-neighbour communication routines that support explicit message passing between neighbouring pencils. As with the rest of the 2DECOMP&FFT library, the API is designed to be very user-friendly:

```
      call update_halo(var, var_halo, level)
```
Here the first parameter `var`, a 3D input array, contains the normal pencil-distributed data as defined by the decomposition library. After invoking the routine, the second parameter `var_halo`, an output, returns all original data plus halo data from the neighbouring processes. One can imagine that pencils are now fatter and overlap with the neighbouring pencils. The third parameter `level` defines how many layers of overlapping is required. `var_halo` should be defined from the calling routine as either a 3D allocatable array or pointer. Its memory space will be calculated and allocated by the library. When the routine returns, `var_halo` can be referenced by the calling program using the normal *i,j,k* indices.

As with the rest of the 2DECOMP&FFT library, a more general form of the routine is available (implemented using Fortran optional arguments):
```
      call update_halo(var, var_halo, level, opt_decomp, opt_global)
```
This supports halo-cell communications among pencils with arbitrary global sizes, as described by `opt_decomp`, the decomposition object. The last optional parameter `opt_global` is required (to be set to `.true.`) if global coordinate is used to define the pencils, i.e. the input array `var` is defined using the *start/end* variables rather than the *size* variables. This ensures the coordinate systems used by `var` and `var_halo` are consistent.

To demonstrate the use of this API, here is an example that computes spatial derivatives:

```
      ! to calculate dv/dy, assume that variables are stored in X-pencil
      
      real, allocatable, dimension(:,:,:) :: v, v_halo, dvdy
      
      allocate(v(xsize(1), xsize(2), xsize(3)))
      allocate(dvdy(xsize(1), xsize(2), xsize(3)))
      
      call update_halo(v,v_halo,level=1)
      
      ! compute derivatives
      do k=1,xsize(3)
         do j=1,xsize(2)
            do i=1,xsize(1)
               dvdy(i,j,k) = (v_halo(i,j+1,k)-v_halo(i,j-1,k)) / dy
            end do
         end do
      end do
```

As seen, the variables are stored in X-pencil and derivatives are to be evaluated over distributed data along Y direction using a central finite difference scheme. This is the perfect situation to use the halo-cell support API. Using global transpositions would be unnecessarily too expensive for this type of local/explicit calculations. After the call to `update_halo`, it is safe to refer to the *j+1* and *j-1* indices on array `v_halo` in order to compute the derivatives.

Note that for the pencils bordering the computational domain, it is up to the application to handle the physical boundary conditions. The library does support periodic condition, i.e. for processes near the boundary of the computational domain, a call to the update_halo routine will fill the halo cells of one side with values from the other side of the domain, when periodic condition is required. To specify periodic condition, one need to initialise the decomposition library with additional information:
```
      call decomp_2d_init(nx, ny, nz, P_row, P_col, periodic_bc)
```
The extra parameter `periodic_bc` is a 1D array containing 3 logical values that specify which dimension should be periodic. This parameter is optional and is only used with the halo-cell API. The domain decomposition should otherwise behaves exactly as normal.

Like the rest of 2DECOMP&FFT, the halo-cell support API is implemented in a black-box fashion. The library internally handles the communications between neighbouring blocks using the standard MPI non-blocking point-to-point communications.