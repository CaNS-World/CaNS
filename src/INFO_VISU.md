# how to visualize the output binary files form `CaNS`

1. after the simulation has run, copy the contents of `utils/visualize_fields/gen_xdmf_non_uniform_grid/` to the simulation `data` folder;
2. edit the `param.h90` file according to the flow physical and computational parameters, and output frequency;
3. run the `./genview.sh`, which generates a file `viewfld.xmf`. The script uses `gfortran` to compile `gen_xdmf.f90`;
4. the file `viewfld.xmf` can be use to visualize the contents of the raw binary files with e.g. paraview.

For more details see the first lines of `gen_xdmf.f90`.
