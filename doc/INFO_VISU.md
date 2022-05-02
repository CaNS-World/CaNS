# how to visualize the output binary files form *CaNS*

## \[*NEW*\] the easy way

in addition to the binary files for visualization, *CaNS* now generates a log file that contains information about the saved data (see `out2d.h90` and `out3d.h90` for more details); this new approach uses that log file to generate the `Xdmf` visualization file.

the steps are as follows:

1. after the simulation has run, copy the contents of `utils/visualize_fields/gen_xdmf_easy/write_xdmf.py` to the simulation `data` folder;
2. run the file with `python write_xdmf.py`.
3. load the generated Xdmf (`*.xmf`) file using paraview/visit or other visualization software.

## example: how to visualize the default binary output

### 3D fields

when running the script `write_xdmf.py` we get the following prompts:

~~~
 $ python write_xdmf.py
 Name of the log file written by CaNS [log_visu_3d.out]:
 Name to be appended to the grid files to prevent overwriting []:
 Name of the output file [viewfld_DNS.xmf]:
~~~

* the first value is the name of the file that logged the saved data;
* the second is a name to append to the grid files that are generated, which should change for different log files to prevent conflicts;
* the third is the name of the visualization file.

by pressing <kbd>enter</kbd> three times, the default values in the square brackets are assumed by the script; these correspond to the default steps required for visualizing 3D field data.

### 2D fields

the procedure for visualizing 2D field data that is saved by *CaNS* in `out2d.h90` is exactly the same; it is just that the correct log file should be selected. *CaNS* saves by default field data in a plane of constant `y=ly/2`, and logs the saves to a file named `log_visu_2d_slice_1.out`. If more planes are saved, the user should make sure that one log file per plane is saved by *CaNS* (e.g. if another plane is saved, the log file written in `out2d.h90` could be named `log_visu_2d_slice_2.out`); see `out2d.h90` for more details. The corresponding steps to generate the Xdmf file would be, for instance:

~~~
 $ python write_xdmf.py
 Name of the log file written by CaNS [log_visu_3d.out]: log_visu_2d_slice_1.out
 Name to be appended to the grid files to prevent overwriting []: 2d
 Name of the output file [viewfld_DNS.xmf]: viewfld_DNS_2d.xmf
~~~

### checkpoint files

A similar script also located in `utils/visualize_fields/gen_xdmf_easy/`, named `write_xdmf_restart.py`, can be used to generate metadata that allows to visualize the field data contained in all saved checkpoint files:

~~~
 $ python write_xdmf_restart.py
 Name of the pattern of the restart files to be visualized [fld?*.bin]:
 Names of stored variables [VEX VEY VEZ PRE]:
 Name to be appended to the grid files to prevent overwriting [_fld]:
 Name of the output file [viewfld_DNS_fld.xmf]:
~~~

## the older way

1. after the simulation has run, copy the contents of `utils/visualize_fields/gen_xdmf_non_uniform_grid/` to the simulation `data` folder;
2. edit the `param.h90` file according to the flow physical and computational parameters, and output frequency;
3. run the `./genview.sh`, which generates a file `viewfld_DNS.xmf`. The script uses `gfortran` to compile `gen_xdmf.f90`;
4. the file `viewfld_DNS.xmf` can be used to visualize the contents of the raw binary files with e.g. paraview.

For more details see the first lines of `gen_xdmf.f90`.
