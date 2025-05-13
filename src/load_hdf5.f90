! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_load_hdf5
  use iso_fortran_env
#if defined(_DECOMP_X)
#undef _DECOMP_X_IO
#endif
  use mpi
  use mod_common_mpi, only: myid,ierr
  use mod_types
  use mod_utils, only: f_sizeof
  use mod_scal, only: scalar
  use mod_param, only: ipencil_axis, l
  implicit none
  private
  public load_one
  interface load_one
    procedure :: io_field_hdf5
  end interface load_one
  contains
#if defined(_USE_HDF5)
  subroutine io_field_hdf5(io,filename,c_io_vars,comm,ng,nh,lo,hi,io_vars,time,istep, nvar)
    use hdf5
    use h5ds 
    !
    ! collective single field data I/O using HDF5
    !
    ! written with the help of Josh Romero,
    ! with the field data I/O inspired from the AFiD code
    !
    implicit none
    character(len=1), intent(in) :: io
    character(len=*), intent(in) :: filename
    integer         , intent(in) :: comm, nvar
    integer         , intent(in), dimension(3)   :: ng,nh,lo,hi
    ! Must statically define kind real64 for hdf5 routine,
    ! either duplicate for other precisions (and add to interface) or do some Macro wizardry
    !real(real64), intent(inout),contiguous, dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: var
    type(arr_ptr)    , allocatable, intent(in) ::   io_vars(:)
    character(len=10), allocatable, intent(in) :: c_io_vars(:)
    real(real64), intent(inout) :: time
    integer(int32) , intent(inout) :: istep

    real(real64), pointer, dimension(:,:,:,:) :: var
    real(real64), allocatable :: axis_x(:), axis_y(:), axis_z(:)
    
    !
    ! HDF5 variables
    !
    integer :: ndims, ierr, i
    integer(HID_T) :: file_id
    integer(HID_T) :: dspace_id, dset_id, attr_id
    integer(HID_T) :: time_dspace_id, time_id
    integer(HID_T) :: filespace, slabspace, memspace, timespace
    !
    integer(HSIZE_T) :: dims(4), chunk(4), ext_dims(4)
    !
    integer(HID_T) :: xfer_pid, file_pid, dset_pid, time_pid
    integer(HSIZE_T) , dimension(4) :: data_count
    integer(HSIZE_T) , dimension(4) :: data_halo_count
    integer(HSSIZE_T), dimension(4) :: data_offset
    integer(HSSIZE_T), dimension(4) :: halo_offset
    integer(HSSIZE_T), dimension(4) :: cur_dims
    integer(HSSIZE_T), dimension(1) :: ntime, time_ext
    logical :: file_exists, dset_exists
    character(len=20):: varname

    !
    ndims = 4 ! fourth dimension is time
    dims(1:3) = ng(:)
    dims(4) = 1
    ext_dims = dims
    ext_dims(4) = -1 !same as H5S_UNLIMITED_F but that is FUBAR: https://stackoverflow.com/questions/17703396/
    data_count(1:3) = hi(:) - lo(:) + 1
    data_count(4) = 1
    data_halo_count(1:3) = data_count(1:3) + 2*nh(:)
    data_halo_count(4) = 1
    data_offset(1:3) = lo(:) - 1 ! starts from 0
    halo_offset(1:3) = nh(:)
    halo_offset(4) = 0
    ntime = 1
    !

    allocate(axis_x(dims(1)))
    allocate(axis_y(dims(2)))
    allocate(axis_z(dims(3)))
    do i = 1, dims(1)
      axis_x(i) = (i-1)*l(1)/(1.*dims(1)-1)
    end do
    do i = 1, dims(2)
      axis_y(i) = (i-1)*l(2)/(1.*dims(2)-1)
    end do
    do i = 1, dims(3)
      axis_z(i) = (i-1)*l(3)/(1.*dims(3)-1)
    end do

    ! Common operations
    call h5open_f(ierr)

    call h5pcreate_f(h5p_dataset_xfer_f, xfer_pid, ierr)
    call h5pset_dxpl_mpio_f(xfer_pid, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5pcreate_f(h5p_file_access_f, file_pid, ierr)
    call h5pset_fapl_mpio_f(file_pid, comm, MPI_INFO_NULL, ierr)
    call h5pcreate_f(h5p_dataset_create_f, dset_pid, ierr)

    ! Sets the other dimensions of the chunk, 1 make a chunk of every independent pencil
    chunk = 1
    !Change chunking axis by editing the following line
    chunk(ipencil_axis) = ng(ipencil_axis)
    !Turn chunks on/off by commenting out the following line
    call h5pset_chunk_f(dset_pid, ndims, chunk , ierr)
    !Turn compression on/off by toggeling the following line, or change the level (1 least/fast, 9 most/slow)
    call h5pset_deflate_f(dset_pid, 1, ierr) 

    call h5pcreate_f(h5p_dataset_create_f, time_pid, ierr)
    call h5pset_chunk_f(time_pid, 1, (/100_HSIZE_T/) , ierr)

    call h5screate_simple_f(ndims, data_halo_count, memspace, ierr) 
    call h5sselect_hyperslab_f (memspace, H5S_SELECT_SET_F, halo_offset, data_count, ierr)

    select case(io)
    case('r')
      error stop "this extendible dataset implementation is not meant to be read from (although possible if you really want)" 

    case('w')
      inquire(file=filename,exist=file_exists)
      if (file_exists) then 
        call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, ierr, access_prp=file_pid)

        call h5dopen_f(file_id, "Time", time_id, ierr)
        call h5dget_space_f(time_id, time_dspace_id, ierr)
        call h5sget_simple_extent_dims_f(time_dspace_id, ntime, time_ext, ierr)
        ntime = ntime + 1
        call h5dset_extent_f(time_id, ntime, ierr)
      else
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr, access_prp=file_pid)

        ! Xaxis
        call h5screate_simple_f(1, dims(1:1), slabspace, ierr)
        call h5dcreate_f(file_id, "Xaxis", H5T_IEEE_F64LE, slabspace, dset_id, ierr)
        call h5dwrite_f(dset_id,H5T_IEEE_F64LE, axis_x, dims(1:1), ierr)
        call h5dclose_f(dset_id, ierr)
        call h5sclose_f(slabspace, ierr)

        ! Yaxis
        call h5screate_simple_f(1, dims(2:2), slabspace, ierr)
        call h5dcreate_f(file_id, "Yaxis", H5T_IEEE_F64LE, slabspace, dset_id, ierr)
        call h5dwrite_f(dset_id,H5T_IEEE_F64LE, axis_y, dims(2:2), ierr)
        call h5dclose_f(dset_id, ierr)
        call h5sclose_f(slabspace, ierr)

        ! Zaxis
        call h5screate_simple_f(1, dims(3:3), slabspace, ierr)
        call h5dcreate_f(file_id, "Zaxis", H5T_IEEE_F64LE, slabspace, dset_id, ierr)
        call h5dwrite_f(dset_id,H5T_IEEE_F64LE, axis_z, dims(3:3), ierr)
        call h5dclose_f(dset_id, ierr)
        call h5sclose_f(slabspace, ierr)

        call h5screate_simple_f(1, (/1_HSIZE_T/), time_dspace_id, ierr, (/-1_HSIZE_T/))
        call h5dcreate_f(file_id, "Time", H5T_IEEE_F64LE, time_dspace_id, time_id, ierr, &
                           dcpl_id=time_pid)
        call h5dsset_scale_f(time_id,ierr, "time")
      endif

      call h5screate_simple_f(1, ntime, slabspace, ierr)
      call h5screate_simple_f(1, (/1_HSIZE_T/), timespace, ierr)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, ntime-1, (/1_HSIZE_T/), ierr)
      call h5dwrite_f(time_id, H5T_IEEE_F64LE, time, (/1_HSIZE_T/), ierr, file_space_id=slabspace, mem_space_id=timespace)
      call h5sclose_f(slabspace, ierr)
      call h5sclose_f(timespace, ierr)
      call h5sclose_f(time_dspace_id, ierr)

      do i = 1, nvar
        varname = c_io_vars(i)
        call h5lexists_f(file_id, varname, dset_exists, ierr)
        if (dset_exists) then
          call h5dopen_f(file_id, varname, dset_id, ierr)
          call h5dget_space_f(dset_id, dspace_id, ierr)
          call h5sget_simple_extent_dims_f(dspace_id, dims, ext_dims, ierr)
          dims(4) = dims(4) + 1
          call h5dset_extent_f(dset_id, dims, ierr)
        else
          call h5screate_simple_f(ndims, dims, dspace_id, ierr, ext_dims)
          call h5dcreate_f(file_id, varname, H5T_IEEE_F64LE, dspace_id, dset_id, ierr, &
                             dcpl_id=dset_pid)

          call h5dopen_f(file_id, "Xaxis", attr_id, ierr) 
          call h5dsattach_scale_f(dset_id, attr_id, 4, ierr)
          call h5dclose_f(attr_id, ierr)

          call h5dopen_f(file_id, "Yaxis", attr_id, ierr) 
          call h5dsattach_scale_f(dset_id, attr_id, 3, ierr)
          call h5dclose_f(attr_id, ierr)

          call h5dopen_f(file_id, "Zaxis", attr_id, ierr) 
          call h5dsattach_scale_f(dset_id, attr_id, 2, ierr)
          call h5dclose_f(attr_id, ierr)

          call h5dsattach_scale_f(dset_id, time_id, 1, ierr)
        end if
        call h5screate_simple_f(ndims, dims, slabspace, ierr)
        data_offset(4) = dims(4)-1
        call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, data_offset, data_count, ierr)

        var(1:data_halo_count(1), 1:data_halo_count(2), 1:data_halo_count(3), 1:1) => io_vars(i)%arr

        call h5dwrite_f(dset_id, H5T_IEEE_F64LE, var, dims, ierr, &
                        file_space_id=slabspace, mem_space_id=memspace, xfer_prp=xfer_pid)
        call h5dclose_f(dset_id, ierr)
        call h5sclose_f(slabspace, ierr)
        call h5sclose_f(dspace_id, ierr)
      end do

      call h5dclose_f(time_id, ierr)
      call h5fclose_f(file_id, ierr)
    end select
10 continue
    call h5sclose_f(memspace, ierr)
    call h5pclose_f(file_pid, ierr)
    call h5pclose_f(time_pid, ierr)
    call h5pclose_f(xfer_pid, ierr)
    call h5close_f(ierr)
  end subroutine io_field_hdf5
#endif
end module mod_load_hdf5
