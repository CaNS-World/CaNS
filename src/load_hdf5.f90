! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
#if defined(_USE_HDF5)
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
  use mod_param, only: ipencil_axis, compression_level, chunk_checkpoint
  implicit none
  private
  public load_one
  interface load_one
    procedure :: io_field_hdf5
  end interface load_one
  contains
  subroutine io_field_hdf5(io,filename,c_io_vars,comm,ng,nh,lo,hi,io_vars,time,istep, nvar)
    use hdf5
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

    real(real64), pointer, dimension(:,:,:) :: var
    integer , dimension(3) :: n
    integer , dimension(3) :: subsizes, starts
    !
    ! HDF5 variables
    !
    integer :: ndims, ierr, i
    integer(HID_T) :: file_id, group_id
    integer(HID_T) :: dspace_id, dset_id, attr_id
    integer(HID_T) :: filespace
    integer(HID_T) :: slabspace
    integer(HID_T) :: memspace
    !
    integer(HSIZE_T) :: dims(3), chunk(3)
    !
    integer(HID_T) :: xfer_pid, file_pid, dset_pid
    integer(HSIZE_T) , dimension(3) :: data_count
    integer(HSSIZE_T), dimension(3) :: data_offset
    integer(HSSIZE_T), dimension(3) :: halo_offset
    logical :: file_exists, group_exists, dset_exists
    character(len=20):: name, varname
    integer(hsize_t) :: data_size

    !
    n(:)        = hi(:)-lo(:)+1
    subsizes(:) = n(:)
    starts(:)   = lo(:) - 1 ! starts from 0
    !
    ndims = 3
    dims(:) = ng(:)
    data_count(:) = subsizes(:)
    data_offset(:) = starts(:)
    halo_offset(:) = nh(:)+1
    !

    ! Common operations
    call h5open_f(ierr)

    call h5pcreate_f(h5p_dataset_xfer_f, xfer_pid, ierr)
    call h5pset_dxpl_mpio_f(xfer_pid, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5pcreate_f(h5p_file_access_f, file_pid, ierr)
    call h5pset_fapl_mpio_f(file_pid, comm, MPI_INFO_NULL, ierr)
    call h5pcreate_f(h5p_dataset_create_f, dset_pid, ierr)

    if(chunk_checkpoint) then
      ! Sets the other dimensions of the chunk, 1 make a chunk of every independent pencil
      chunk = 1
      !Change chunking axis by editing the following line
      chunk(ipencil_axis) = ng(ipencil_axis)
      !Turn chunks on/off by commenting out the following line
      call h5pset_chunk_f(dset_pid, ndims, chunk , ierr)

      !Turn compression on/off by toggeling the following line, or change the level (1 least/fast, 9 most/slow)
      call h5pset_deflate_f(dset_pid, compression_level, ierr) 
    endif

    call h5screate_simple_f(ndims, data_count+2*nh(:), memspace, ierr) 
    call h5screate_simple_f(ndims, dims, slabspace, ierr)

    call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, data_offset, data_count, ierr)
    call h5sselect_hyperslab_f (memspace, H5S_SELECT_SET_F, halo_offset, data_count, ierr)

    select case(io)
    case('r')

      inquire(file=filename,exist=file_exists)
      if (file_exists) then 
        call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, ierr, access_prp=file_pid)
      else
        error stop "Checkpoint file "//filename//"  does not exist"
      endif

      ! find latest group
      call h5lget_name_by_idx_f(file_id, ".", H5_INDEX_NAME_F, H5_ITER_DEC_F, 0_HSIZE_T, name, ierr, data_size)

      call h5aopen_by_name_f(file_id, trim(name), 'time', attr_id, ierr)
      call h5aread_f(attr_id, H5T_IEEE_F64LE, time, (/1_HSIZE_T/), ierr)
      call h5aclose_f(attr_id, ierr)

      call h5aopen_by_name_f(file_id, trim(name), 'istep', attr_id, ierr)
      call h5aread_f(attr_id, H5T_STD_I32LE, istep, (/1_HSIZE_T/), ierr)
      call h5aclose_f(attr_id, ierr)

      do i = 1, nvar
        varname = c_io_vars(i)
        var(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) => io_vars(i)%arr
        call h5dopen_f(file_id, trim(name)//'/'//varname, dset_id, ierr)
        call h5dread_f(dset_id,H5T_IEEE_F64LE,var,dims,ierr, &
                       file_space_id=slabspace,mem_space_id=memspace,xfer_prp=xfer_pid)
        call h5dclose_f(dset_id, ierr)
      end do

      call h5fclose_f(file_id, ierr)

    case('w')
      inquire(file=filename,exist=file_exists)
      if (file_exists) then 
        call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, ierr, access_prp=file_pid)
      else
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr, access_prp=file_pid)
      endif

      write(name,'(I8.8)') istep
      name = 'istep_'//trim(adjustl(name))
      call h5lexists_f(file_id, name, group_exists, ierr)
      if (group_exists) then
        call h5gopen_f(file_id, name, group_id, ierr)
      else
        call h5gcreate_f(file_id, name, group_id, ierr)

        call h5screate_simple_f(1, (/1_HSIZE_T/), dspace_id, ierr)
        call h5acreate_f(group_id, "time", H5T_IEEE_F64LE, dspace_id, attr_id, ierr)
        call h5awrite_f(attr_id, H5T_IEEE_F64LE, time, (/1_HSIZE_T/), ierr)
        call h5aclose_f(attr_id, ierr)

        call h5acreate_f(group_id, "istep", H5T_STD_I32LE, dspace_id, attr_id, ierr)
        call h5awrite_f(attr_id, H5T_STD_I32LE, istep, (/1_HSIZE_T/), ierr)
        call h5aclose_f(attr_id, ierr)
        call h5sclose_f(dspace_id, ierr)
      end if

      call h5screate_simple_f(ndims, dims, dspace_id, ierr)

      do i = 1, nvar
        varname = c_io_vars(i)
        call h5lexists_f(group_id, varname, dset_exists, ierr)
        if (dset_exists) then
          call h5dopen_f(group_id, varname, dset_id, ierr)
        else
          call h5dcreate_f(group_id, varname, H5T_IEEE_F64LE, dspace_id, dset_id, ierr, &
                           dcpl_id=dset_pid)
        end if


        var(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) => io_vars(i)%arr
        call h5dwrite_f(dset_id, H5T_IEEE_F64LE, var, dims, ierr, &
                        file_space_id=slabspace, mem_space_id=memspace, xfer_prp=xfer_pid)
        call h5dclose_f(dset_id, ierr)
      end do
      call h5sclose_f(dspace_id, ierr)
      call h5gclose_f(group_id, ierr)
      call h5fclose_f(file_id, ierr)
    end select
10 continue
    call h5sclose_f(memspace, ierr)
    call h5sclose_f(slabspace, ierr)
    call h5pclose_f(file_pid, ierr)
    call h5pclose_f(xfer_pid, ierr)
    call h5close_f(ierr)
  end subroutine io_field_hdf5
end module mod_load_hdf5
#endif
