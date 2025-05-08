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
  implicit none
  private
  public load_one
  interface load_one
    procedure :: io_field_hdf5
  end interface load_one
  contains
#if defined(_USE_HDF5)
  subroutine io_field_hdf5(io,filename,varname,comm,ng,nh,lo,hi,var,time,istep)
    use hdf5
    !
    ! collective single field data I/O using HDF5
    !
    ! written with the help of Josh Romero,
    ! with the field data I/O inspired from the AFiD code
    !
    implicit none
    character(len=1), intent(in) :: io
    character(len=*), intent(in) :: filename,varname
    integer         , intent(in) :: comm
    integer         , intent(in), dimension(3)   :: ng,nh,lo,hi
    ! Must statically define kind for hdf5 routine, either duplicate for other precisions (and add to interface) or do some Macro wizardry
    real(real64), intent(inout),contiguous, dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: var
    real(real64), intent(inout) :: time
    integer , intent(inout) :: istep

    integer , dimension(3) :: n
    integer , dimension(3) :: subsizes, starts
    !
    ! HDF5 variables
    !
    integer :: ndims, ierr
    integer(HID_T) :: file_id, group_id
    integer(HID_T) :: dspace_id, dset_id
    integer(HID_T) :: filespace
    integer(HID_T) :: slabspace
    integer(HID_T) :: memspace
    !
    integer(HSIZE_T) :: dims(3)
    !
    integer(HID_T) :: xfer_plist_id, file_plist_id
    integer(HSIZE_T) , dimension(3) :: data_count
    integer(HSSIZE_T), dimension(3) :: data_offset
    integer(HSSIZE_T), dimension(3) :: halo_offset
    logical :: file_exists, group_exists
    character(len=20):: name
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

    call h5pcreate_f(h5p_dataset_xfer_f, xfer_plist_id, ierr)
    call h5pset_dxpl_mpio_f(xfer_plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5pcreate_f(h5p_file_access_f, file_plist_id, ierr)
    call h5pset_fapl_mpio_f(file_plist_id, comm, MPI_INFO_NULL, ierr)

    call h5screate_simple_f(ndims, data_count+2*nh(:), memspace, ierr) 
    call h5screate_simple_f(ndims, dims, slabspace, ierr)

    call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, data_offset, data_count, ierr)
    call h5sselect_hyperslab_f (memspace, H5S_SELECT_SET_F, halo_offset, data_count, ierr)

    select case(io)
    case('r')

      inquire(file=filename,exist=file_exists)
      if (file_exists) then 
        call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, ierr, access_prp=file_plist_id)
      else
        error stop "Checkpoint file "//filename//"  does not exist"
      endif

      ! find latest group
      call h5lget_name_by_idx_f(file_id, ".", H5_INDEX_NAME_F, H5_ITER_DEC_F, 0_HSIZE_T, name, ierr, data_size)
      call h5dopen_f(file_id, trim(name)//'/'//varname, dset_id, ierr)


      call h5dread_f(dset_id,H5T_IEEE_F64LE,var,dims,ierr, &
                     file_space_id=slabspace,mem_space_id=memspace,xfer_prp=xfer_plist_id)

      call h5dclose_f(dset_id, ierr)
      call h5fclose_f(file_id, ierr)

      ! set time from group name (simpler then using h5 attributes, for now)
      read(name, '(6x,I8)') istep

    case('w')
      inquire(file=filename,exist=file_exists)
      if (file_exists) then 
        call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, ierr, access_prp=file_plist_id)
      else
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr, access_prp=file_plist_id)
      endif

      write(name,'(I8.8)') istep
      name = 'istep_'//trim(adjustl(name))
      call h5lexists_f(file_id, name, group_exists, ierr)
      if (group_exists) then
        call h5gopen_f(file_id, name, group_id, ierr)
      else
        call h5gcreate_f(file_id, name, group_id, ierr)
      end if

      call h5screate_simple_f(ndims, dims, dspace_id, ierr)

      call h5dcreate_f(group_id, varname, H5T_IEEE_F64LE, dspace_id, dset_id, ierr)
      call h5dwrite_f(dset_id, h5t_ieee_f64le, var, dims, ierr, &
                      file_space_id=slabspace, mem_space_id=memspace, xfer_prp=xfer_plist_id)
      call h5dclose_f(dset_id, ierr)
      call h5sclose_f(dspace_id, ierr)
      call h5gclose_f(group_id, ierr)
      call h5fclose_f(file_id, ierr)
    end select
10 continue
    call h5sclose_f(memspace, ierr)
    call h5sclose_f(slabspace, ierr)
    call h5pclose_f(file_plist_id, ierr)
    call h5pclose_f(xfer_plist_id, ierr)
    call h5close_f(ierr)
  end subroutine io_field_hdf5
#endif
end module mod_load_hdf5
