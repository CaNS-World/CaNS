! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_load
#if defined(_DECOMP_X)
#undef _DECOMP_X_IO
#endif
  use mpi
#if defined(_USE_ADIOS2)
  use adios2
#endif
  use mod_common_mpi, only: myid,ierr
  use mod_types
  use mod_utils, only: f_sizeof
  use mod_scal, only: scalar
  implicit none
  private
  integer, parameter :: FILETYPE_MPIIO = 0, FILETYPE_HDF5 = 1, FILETYPE_ADIOS2 = 2
#if defined(_USE_ADIOS2)
  interface io_field_adios2_1d
    module procedure io_field_adios2_1d_real
    module procedure io_field_adios2_1d_int
  end interface io_field_adios2_1d
  type :: adios2_compression_params
    !
    ! name  : compression algorithm name
    ! keys  : parameter names passed to the compressor
    ! values: parameter values paired with keys
    !
    character(len=64) :: name = ''
    character(len=64), dimension(10) :: keys = ''
    character(len=64), dimension(10) :: values = ''
  end type adios2_compression_params
#endif
  public :: load_all,load_one,io_field,io_write_subset
contains
  subroutine load_all(io,filename,comm,ng,nh,lo,hi,nscal,u,v,w,p,s,time,istep,x_g,y_g,z_g)
    !
    ! dispatches restart I/O to the selected backend
    !
    implicit none
    character(len=1), intent(in) :: io
    character(len=*), intent(in) :: filename
    integer         , intent(in) :: comm
    integer , intent(in), dimension(3) :: ng,nh,lo,hi
    integer , intent(in)               :: nscal
    real(rp), intent(inout), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: u,v,w,p
    type(scalar), intent(inout), dimension(:) :: s
    real(rp), intent(inout) :: time
    integer , intent(inout) :: istep
    real(rp), intent(inout), dimension(1-nh(1):), optional :: x_g
    real(rp), intent(inout), dimension(1-nh(2):), optional :: y_g
    real(rp), intent(inout), dimension(1-nh(3):), optional :: z_g
    integer :: file_type
    !
    file_type = fetch_file_type(filename)
    select case(file_type)
    case(FILETYPE_MPIIO)
      call load_all_mpiio(io,filename,comm,ng,nh,lo,hi,nscal,u,v,w,p,s,time,istep,x_g,y_g,z_g)
#if defined(_USE_HDF5)
    case(FILETYPE_HDF5)
      call load_all_hdf5(io,filename,comm,ng,nh,lo,hi,nscal,u,v,w,p,s,time,istep,x_g,y_g,z_g)
#endif
#if defined(_USE_ADIOS2)
    case(FILETYPE_ADIOS2)
      call load_all_adios2(io,filename,comm,ng,nh,lo,hi,nscal,u,v,w,p,s,time,istep,x_g,y_g,z_g)
#endif
    case default
      call load_all_mpiio(io,filename,comm,ng,nh,lo,hi,nscal,u,v,w,p,s,time,istep,x_g,y_g,z_g)
    end select
  end subroutine load_all
  !
  subroutine load_all_mpiio(io,filename,comm,ng,nh,lo,hi,nscal,u,v,w,p,s,time,istep,x_g,y_g,z_g)
    !
    ! reads/writes a restart file
    !
    implicit none
    character(len=1), intent(in) :: io
    character(len=*), intent(in) :: filename
    integer         , intent(in) :: comm
    integer , intent(in), dimension(3) :: ng,nh,lo,hi
    integer , intent(in)               :: nscal
    real(rp), intent(inout), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: u,v,w,p
    type(scalar), intent(inout), dimension(:) :: s
    real(rp), intent(inout) :: time
    integer , intent(inout) :: istep
    real(rp), intent(inout), dimension(1-nh(1):), optional :: x_g
    real(rp), intent(inout), dimension(1-nh(2):), optional :: y_g
    real(rp), intent(inout), dimension(1-nh(3):), optional :: z_g
    real(rp), dimension(2) :: fldinfo
    integer :: iscal
    integer :: fh
    integer :: nreals_myid
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp,good
    !
    select case(io)
    case('r')
      call MPI_FILE_OPEN(comm,filename, &
                         MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
      !
      ! check file size first
      !
      call MPI_FILE_GET_SIZE(fh,filesize,ierr)
      good = (product(int(ng(:),MPI_OFFSET_KIND))*(4+nscal)+2)*f_sizeof(1._rp)
      if(filesize /= good) then
        if(myid == 0) print*, ''
        if(myid == 0) print*, '*** Simulation aborted due to a checkpoint file with incorrect size ***'
        if(myid == 0) print*, '    file: ', filename, ' | expected size: ', good, '| actual size: ', filesize
        call MPI_FINALIZE(ierr)
        error stop
      end if
      !
      ! read
      !
      disp = 0_MPI_OFFSET_KIND
#if !defined(_DECOMP_X_IO)
      call io_field(io,fh,ng,nh,lo,hi,disp,u)
      call io_field(io,fh,ng,nh,lo,hi,disp,v)
      call io_field(io,fh,ng,nh,lo,hi,disp,w)
      call io_field(io,fh,ng,nh,lo,hi,disp,p)
      do iscal=1,nscal
        call io_field(io,fh,ng,nh,lo,hi,disp,s(iscal)%val)
      end do
#else
      block
        !
        ! I/O over x-aligned pencils
        !
        use decomp_2d
        use mod_param, only: ipencil => ipencil_axis
        real(rp), allocatable, dimension(:,:,:) :: tmp_x,tmp_y,tmp_z
        select case(ipencil)
        case(1)
          allocate(tmp_x(0,0,0),tmp_y(0,0,0),tmp_z(0,0,0))
        case(2)
          allocate(tmp_x(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), &
                   tmp_y(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)), &
                   tmp_z(0,0,0))
        case(3)
          allocate(tmp_x(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), &
                   tmp_y(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)), &
                   tmp_z(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
        end select
        call io_field(io,fh,ng,[0,0,0],lo,hi,disp,tmp_x)
        call transpose_to_or_from_x(io,ipencil,nh,u,tmp_x,tmp_y,tmp_z)
        call io_field(io,fh,ng,[0,0,0],lo,hi,disp,tmp_x)
        call transpose_to_or_from_x(io,ipencil,nh,v,tmp_x,tmp_y,tmp_z)
        call io_field(io,fh,ng,[0,0,0],lo,hi,disp,tmp_x)
        call transpose_to_or_from_x(io,ipencil,nh,w,tmp_x,tmp_y,tmp_z)
        call io_field(io,fh,ng,[0,0,0],lo,hi,disp,tmp_x)
        call transpose_to_or_from_x(io,ipencil,nh,p,tmp_x,tmp_y,tmp_z)
        do iscal=1,nscal
          call io_field(io,fh,ng,[0,0,0],lo,hi,disp,tmp_x)
          call transpose_to_or_from_x(io,ipencil,nh,s(iscal)%val,tmp_x,tmp_y,tmp_z)
        end do
        deallocate(tmp_x,tmp_y,tmp_z)
      end block
#endif
      call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,MPI_REAL_RP,'native',MPI_INFO_NULL,ierr)
      nreals_myid = 0
      if(myid == 0) nreals_myid = 2
      call MPI_FILE_READ(fh,fldinfo,nreals_myid,MPI_REAL_RP,MPI_STATUS_IGNORE,ierr)
      call MPI_FILE_CLOSE(fh,ierr)
      call MPI_BCAST(fldinfo,2,MPI_REAL_RP,0,comm,ierr)
      time  =      fldinfo(1)
      istep = nint(fldinfo(2))
    case('w')
      !
      ! write
      !
      call MPI_FILE_OPEN(comm,filename, &
                         MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
      filesize = 0_MPI_OFFSET_KIND
      call MPI_FILE_SET_SIZE(fh,filesize,ierr)
      disp = 0_MPI_OFFSET_KIND
#if !defined(_DECOMP_X_IO)
      call io_field(io,fh,ng,nh,lo,hi,disp,u)
      call io_field(io,fh,ng,nh,lo,hi,disp,v)
      call io_field(io,fh,ng,nh,lo,hi,disp,w)
      call io_field(io,fh,ng,nh,lo,hi,disp,p)
      do iscal=1,nscal
        call io_field(io,fh,ng,nh,lo,hi,disp,s(iscal)%val)
      end do
#else
      block
        !
        ! I/O over x-aligned pencils
        !
        use decomp_2d
        use mod_param, only: ipencil => ipencil_axis
        real(rp), allocatable, dimension(:,:,:) :: tmp_x,tmp_y,tmp_z
        select case(ipencil)
        case(1)
          allocate(tmp_x(0,0,0),tmp_y(0,0,0),tmp_z(0,0,0))
        case(2)
          allocate(tmp_x(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), &
                   tmp_y(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)), &
                   tmp_z(0,0,0))
        case(3)
          allocate(tmp_x(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), &
                   tmp_y(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)), &
                   tmp_z(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
        end select
        call transpose_to_or_from_x(io,ipencil,nh,u,tmp_x,tmp_y,tmp_z)
        call io_field(io,fh,ng,[0,0,0],lo,hi,disp,tmp_x)
        call transpose_to_or_from_x(io,ipencil,nh,v,tmp_x,tmp_y,tmp_z)
        call io_field(io,fh,ng,[0,0,0],lo,hi,disp,tmp_x)
        call transpose_to_or_from_x(io,ipencil,nh,w,tmp_x,tmp_y,tmp_z)
        call io_field(io,fh,ng,[0,0,0],lo,hi,disp,tmp_x)
        call transpose_to_or_from_x(io,ipencil,nh,p,tmp_x,tmp_y,tmp_z)
        call io_field(io,fh,ng,[0,0,0],lo,hi,disp,tmp_x)
        do iscal=1,nscal
          call transpose_to_or_from_x(io,ipencil,nh,s(iscal)%val,tmp_x,tmp_y,tmp_z)
          call io_field(io,fh,ng,[0,0,0],lo,hi,disp,tmp_x)
        end do
        deallocate(tmp_x,tmp_y,tmp_z)
      end block
#endif
      call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,MPI_REAL_RP,'native',MPI_INFO_NULL,ierr)
      fldinfo = [time,1._rp*istep]
      nreals_myid = 0
      if(myid == 0) nreals_myid = 2
      call MPI_FILE_WRITE(fh,fldinfo,nreals_myid,MPI_REAL_RP,MPI_STATUS_IGNORE,ierr)
      call MPI_FILE_CLOSE(fh,ierr)
    end select
  end subroutine load_all_mpiio
  !
  subroutine io_field(io,fh,ng,nh,lo,hi,disp,var)
    implicit none
    character(len=1), intent(in)                 :: io
    integer , intent(in)                         :: fh
    integer , intent(in), dimension(3)           :: ng,nh,lo,hi
    integer(kind=MPI_OFFSET_KIND), intent(inout) :: disp
    real(rp), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: var ! best skip intent() attribute here
    integer , dimension(3) :: n
    integer , dimension(3) :: sizes,subsizes,starts
    integer :: type_glob,type_loc
    n(:)        = hi(:)-lo(:)+1
    sizes(:)    = ng(:)
    subsizes(:) = n(:)
    starts(:)   = lo(:) - 1 ! starts from 0
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,MPI_REAL_RP,type_glob,ierr)
    call MPI_TYPE_COMMIT(type_glob,ierr)
    sizes(:)    = n(:) + 2*nh(:)
    subsizes(:) = n(:)
    starts(:)   = 0 + nh(:)
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,MPI_REAL_RP,type_loc ,ierr)
    call MPI_TYPE_COMMIT(type_loc,ierr)
    select case(io)
    case('r')
      call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,type_glob,'native',MPI_INFO_NULL,ierr)
      call MPI_FILE_READ_ALL(fh,var,1,type_loc,MPI_STATUS_IGNORE,ierr)
    case('w')
      call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,type_glob,'native',MPI_INFO_NULL,ierr)
      call MPI_FILE_WRITE_ALL(fh,var,1,type_loc,MPI_STATUS_IGNORE,ierr)
    end select
    disp = disp+product(int(ng(:),MPI_OFFSET_KIND))*f_sizeof(1._rp)
    call MPI_TYPE_FREE(type_glob,ierr)
    call MPI_TYPE_FREE(type_loc ,ierr)
  end subroutine io_field
  !
#if defined(_DECOMP_X_IO)
  subroutine transpose_to_or_from_x(io,ipencil_axis,nh,var,var_x,var_y,var_z)
    !
    ! transpose arrays for I/O over x-aligned pencils
    !
    use decomp_2d
    implicit none
    character(len=1), intent(in) :: io
    integer , intent(in) :: ipencil_axis,nh(3)
    real(rp), dimension(1-nh(1):,1-nh(2):,1-nh(3):) :: var
    real(rp), dimension(:,:,:) :: var_x,var_y,var_z
    integer, dimension(3) :: n
    n(:) = shape(var) - 2*nh(:)
    select case(ipencil_axis)
    case(1)
    case(2)
      select case(io)
      case('r')
        call transpose_x_to_y(var_x,var_y)
        !$OMP PARALLEL WORKSHARE
        var(1:n(1),1:n(2),1:n(3)) = var_y(:,:,:)
        !$OMP END PARALLEL WORKSHARE
      case('w')
        !$OMP PARALLEL WORKSHARE
        var_y(:,:,:) = var(1:n(1),1:n(2),1:n(3))
        !$OMP END PARALLEL WORKSHARE
        call transpose_y_to_x(var_y,var_x)
      end select
    case(3)
      select case(io)
      case('r')
        call transpose_x_to_y(var_x,var_y)
        call transpose_y_to_z(var_y,var_z)
        !$OMP PARALLEL WORKSHARE
        var(1:n(1),1:n(2),1:n(3)) = var_z(:,:,:)
        !$OMP END PARALLEL WORKSHARE
      case('w')
        !$OMP PARALLEL WORKSHARE
        var_z(:,:,:) = var(1:n(1),1:n(2),1:n(3))
        !$OMP END PARALLEL WORKSHARE
        call transpose_z_to_y(var_z,var_y)
        call transpose_y_to_x(var_y,var_x)
      end select
    end select
  end subroutine transpose_to_or_from_x
  !
#if defined(_OPENACC)
  subroutine transpose_to_or_from_x_gpu(io,ipencil_axis,nh,var_io,var)
    !
    ! transpose arrays for I/O over x-aligned pencils on GPUs
    !
    ! note: the Poisson solver buffers are being recycled here, meaning that
    ! I/O should use the same precision as these buffers if this routine
    ! is used in the future;
    ! alternatively one could *temporarily* (i.e., during I/O) offload
    ! device memory and allocate larger buffers, and use a dedicated
    ! cuDecomp grid descriptor without axis-contiguous layout
    !
    use cudecomp
    use mod_common_cudecomp, only: buf => solver_buf_0, work, &
                                   dtype_rp => cudecomp_real_rp, &
                                   ap_x   => ap_x_poi, &
                                   ap_y   => ap_y_poi, &
                                   ap_z   => ap_z_poi, &
                                   ap_x_0 => ap_x    , &
                                   ch => handle,gd => gd_poi
    implicit none
    character(len=1), intent(in) :: io
    integer , intent(in) :: ipencil_axis,nh(3)
    real(rp), dimension(1      :,1      :,       :) :: var_io
    real(rp), dimension(1-nh(1):,1-nh(2):,1-nh(3):) :: var
    real(rp), pointer, contiguous, dimension(:,:,:) :: var_x,var_y,var_z
    integer, dimension(3) :: n,n_x,n_y,n_z,n_x_0
    integer :: i,j,k
    integer :: istat
    !
    !$acc wait
    !
    n(:) = shape(var) - 2*nh(:)
    !
    n_x(:) = ap_x%shape(:)
    n_y(:) = ap_y%shape(:)
    n_z(:) = ap_z%shape(:)
    n_x_0(:) = ap_x_0%shape(:)
    !
    var_x(1:n_x(1),1:n_x(2),1:n_x(3)) => buf(1:product(n_x(:)))
    var_y(1:n_y(1),1:n_y(2),1:n_y(3)) => buf(1:product(n_y(:)))
    var_z(1:n_z(1),1:n_z(2),1:n_z(3)) => buf(1:product(n_z(:)))
    !
    select case(ipencil_axis)
    case(1)
    case(2)
      select case(io)
      case('r')
        !$acc data copyin(var_io) copyout(var)
        !$acc kernels default(present)
        var_x(1:n_x_0(1),1:n_x_0(2),1:n_x_0(3)) = var_io(:,:,:)
        !$acc end kernels
        !$acc host_data use_device(var_x,var_y,work)
        istat = cudecompTransposeXtoY(ch,gd,var_x,var_y,work,dtype_rp)
        !$acc end host_data
        !$acc kernels loop collapse(3) default(present)
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              var(i,j,k) = var_y(j,k,i) ! axis-contiguous layout along y
            end do
          end do
        end do
        !$acc end data
      case('w')
        !$acc data copyin(var) copyout(var_io) ! var already present, copyin will be ignored
        !$acc kernels loop collapse(3) default(present)
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              var_y(j,k,i) = var(i,j,k) ! axis-contiguous layout along y
            end do
          end do
        end do
        !$acc host_data use_device(var_y,var_x,work)
        istat = cudecompTransposeYtoX(ch,gd,var_y,var_x,work,dtype_rp)
        !$acc end host_data
        !$acc kernels default(present)
        var_io(:,:,:) = var_x(1:n_x_0(1),1:n_x_0(2),1:n_x_0(3))
        !$acc end kernels
        !$acc end data
      end select
    case(3)
      select case(io)
      case('r')
        !$acc data copyin(var_io) copyout(var)
        !$acc kernels default(present)
        var_x(1:n_x_0(1),1:n_x_0(2),1:n_x_0(3)) = var_io(:,:,:)
        !$acc end kernels
        !$acc host_data use_device(var_x,var_y,var_z,work)
        istat = cudecompTransposeXtoY(ch,gd,var_x,var_y,work,dtype_rp)
        istat = cudecompTransposeYtoZ(ch,gd,var_y,var_z,work,dtype_rp)
        !$acc end host_data
        !$acc kernels default(present)
        var(1:n(1),1:n(2),1:n(3)) = var_z(1:n(1),1:n(2),1:n(3))
        !$acc end kernels
        !$acc end data
      case('w')
        !$acc data copyin(var) copyout(var_io) ! var already present, copyin will be ignored
        !$acc kernels default(present)
        var_z(1:n(1),1:n(2),1:n(3)) = var(1:n(1),1:n(2),1:n(3))
        !$acc end kernels
        !$acc host_data use_device(var_z,var_y,var_x,work)
        istat = cudecompTransposeZtoY(ch,gd,var_z,var_y,work,dtype_rp)
        istat = cudecompTransposeYtoX(ch,gd,var_y,var_x,work,dtype_rp)
        !$acc end host_data
        !$acc kernels default(present)
        var_io(:,:,:) = var_x(1:n_x_0(1),1:n_x_0(2),1:n_x_0(3))
        !$acc end kernels
        !$acc end data
      end select
    end select
  end subroutine transpose_to_or_from_x_gpu
#endif
#endif
  subroutine load_one(io,filename,comm,ng,nh,lo,hi,p,varname,time,istep,x_g,y_g,z_g)
    !
    ! dispatches single-field restart I/O to the selected backend
    !
    implicit none
    character(len=1), intent(in) :: io
    character(len=*), intent(in) :: filename
    integer         , intent(in) :: comm
    integer , intent(in), dimension(3) :: ng,nh,lo,hi
    real(rp), intent(inout), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: p
    character(len=*), intent(in) :: varname
    real(rp), intent(inout), optional :: time
    integer , intent(inout), optional :: istep
    real(rp), intent(inout), dimension(1-nh(1):), optional :: x_g
    real(rp), intent(inout), dimension(1-nh(2):), optional :: y_g
    real(rp), intent(inout), dimension(1-nh(3):), optional :: z_g
    integer :: file_type
    !
    file_type = fetch_file_type(filename)
    select case(file_type)
    case(FILETYPE_MPIIO)
      call load_one_mpiio(io,filename,comm,ng,nh,lo,hi,p,varname,time,istep,x_g,y_g,z_g)
#if defined(_USE_HDF5)
    case(FILETYPE_HDF5)
      call load_one_hdf5(io,filename,comm,ng,nh,lo,hi,p,varname,time,istep,x_g,y_g,z_g)
#endif
#if defined(_USE_ADIOS2)
    case(FILETYPE_ADIOS2)
      call load_one_adios2(io,filename,comm,ng,nh,lo,hi,p,varname,time,istep,x_g,y_g,z_g)
#endif
    case default
      call load_one_mpiio(io,filename,comm,ng,nh,lo,hi,p,varname,time,istep,x_g,y_g,z_g)
    end select
  end subroutine load_one
  !
  subroutine load_one_mpiio(io,filename,comm,ng,nh,lo,hi,p,varname,time,istep,x_g,y_g,z_g)
    !
    ! reads/writes a restart file for a single field
    !
    implicit none
    character(len=1), intent(in) :: io
    character(len=*), intent(in) :: filename
    integer         , intent(in) :: comm
    integer , intent(in), dimension(3) :: ng,nh,lo,hi
    real(rp), intent(inout), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: p
    character(len=*), intent(in) :: varname
    real(rp), intent(inout), optional :: time
    integer , intent(inout), optional :: istep
    real(rp), intent(inout), dimension(1-nh(1):), optional :: x_g
    real(rp), intent(inout), dimension(1-nh(2):), optional :: y_g
    real(rp), intent(inout), dimension(1-nh(3):), optional :: z_g
    real(rp), dimension(2) :: fldinfo
    integer :: fh
    integer :: nreals_myid
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp,good
    !
    ! silence compiler warnings for backend-agnostic arguments unused by MPI-IO
    !
    associate(dummy_varname => varname)
    end associate
    !
    select case(io)
    case('r')
      call MPI_FILE_OPEN(comm,filename, &
                         MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
      !
      ! check file size first
      !
      call MPI_FILE_GET_SIZE(fh,filesize,ierr)
      good = (product(int(ng(:),MPI_OFFSET_KIND))*1+2)*f_sizeof(1._rp)
      if(filesize /= good) then
        if(myid == 0) print*, ''
        if(myid == 0) print*, '*** Simulation aborted due to a checkpoint file with incorrect size ***'
        if(myid == 0) print*, '    file: ', filename, ' | expected size: ', good, '| actual size: ', filesize
        call MPI_FINALIZE(ierr)
        error stop
      end if
      !
      ! read
      !
      disp = 0_MPI_OFFSET_KIND
#if !defined(_DECOMP_X_IO)
      call io_field(io,fh,ng,nh,lo,hi,disp,p)
#else
      block
        !
        ! I/O over x-aligned pencils
        !
        use decomp_2d
        use mod_param, only: ipencil => ipencil_axis
        real(rp), allocatable, dimension(:,:,:) :: tmp_x,tmp_y,tmp_z
        select case(ipencil)
        case(1)
          allocate(tmp_x(0,0,0),tmp_y(0,0,0),tmp_z(0,0,0))
        case(2)
          allocate(tmp_x(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), &
                   tmp_y(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)), &
                   tmp_z(0,0,0))
        case(3)
          allocate(tmp_x(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), &
                   tmp_y(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)), &
                   tmp_z(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
        end select
        call io_field(io,fh,ng,[0,0,0],lo,hi,disp,tmp_x)
        call transpose_to_or_from_x(io,ipencil,nh,p,tmp_x,tmp_y,tmp_z)
        deallocate(tmp_x,tmp_y,tmp_z)
      end block
#endif
      call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,MPI_REAL_RP,'native',MPI_INFO_NULL,ierr)
      nreals_myid = 0
      if(myid == 0) nreals_myid = 2
      call MPI_FILE_READ(fh,fldinfo,nreals_myid,MPI_REAL_RP,MPI_STATUS_IGNORE,ierr)
      call MPI_FILE_CLOSE(fh,ierr)
      call MPI_BCAST(fldinfo,2,MPI_REAL_RP,0,comm,ierr)
      if(present(time) .and. present(istep)) then
        time  =      fldinfo(1)
        istep = nint(fldinfo(2))
      end if
    case('w')
      !
      ! write
      !
      call MPI_FILE_OPEN(comm,filename, &
                         MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
      filesize = 0_MPI_OFFSET_KIND
      call MPI_FILE_SET_SIZE(fh,filesize,ierr)
      disp = 0_MPI_OFFSET_KIND
#if !defined(_DECOMP_X_IO)
      call io_field(io,fh,ng,nh,lo,hi,disp,p)
#else
      block
        !
        ! I/O over x-aligned pencils
        !
        use decomp_2d
        use mod_param, only: ipencil => ipencil_axis
        real(rp), allocatable, dimension(:,:,:) :: tmp_x,tmp_y,tmp_z
        select case(ipencil)
        case(1)
          allocate(tmp_x(0,0,0),tmp_y(0,0,0),tmp_z(0,0,0))
        case(2)
          allocate(tmp_x(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), &
                   tmp_y(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)), &
                   tmp_z(0,0,0))
        case(3)
          allocate(tmp_x(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), &
                   tmp_y(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)), &
                   tmp_z(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
        end select
        call transpose_to_or_from_x(io,ipencil,nh,p,tmp_x,tmp_y,tmp_z)
        call io_field(io,fh,ng,[0,0,0],lo,hi,disp,tmp_x)
        deallocate(tmp_x,tmp_y,tmp_z)
      end block
#endif
      if(present(time) .and. present(istep)) then
        call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,MPI_REAL_RP,'native',MPI_INFO_NULL,ierr)
        fldinfo = [time,1._rp*istep]
        nreals_myid = 0
        if(myid == 0) nreals_myid = 2
        call MPI_FILE_WRITE(fh,fldinfo,nreals_myid,MPI_REAL_RP,MPI_STATUS_IGNORE,ierr)
      end if
      call MPI_FILE_CLOSE(fh,ierr)
    end select
  end subroutine load_one_mpiio
  !
  subroutine io_write_subset(filename,varname,comm,ng,nh,lo,hi,lo_out,hi_out,nskip,var,time,istep,x_g,y_g,z_g,is_pack)
    !
    ! writes a structured subset of a 3D field
    !
    implicit none
    character(len=*), intent(in) :: filename,varname
    integer         , intent(in) :: comm
    integer         , intent(in), dimension(3) :: ng,nh,lo,hi,lo_out,hi_out,nskip
    real(rp), intent(in), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: var
    real(rp), intent(in), optional :: time
    integer , intent(in), optional :: istep
    real(rp), intent(in), dimension(1-nh(1):), optional :: x_g
    real(rp), intent(in), dimension(1-nh(2):), optional :: y_g
    real(rp), intent(in), dimension(1-nh(3):), optional :: z_g
    logical , intent(in), optional :: is_pack
    integer :: file_type
    !
    file_type = fetch_file_type(filename)
    select case(file_type)
    case(FILETYPE_MPIIO)
      call write_subset_mpiio(filename,varname,comm,ng,nh,lo,hi,lo_out,hi_out,nskip,var,time,istep,x_g,y_g,z_g,is_pack)
#if defined(_USE_HDF5)
    case(FILETYPE_HDF5)
      call write_subset_hdf5(filename,varname,comm,ng,nh,lo,hi,lo_out,hi_out,nskip,var,time,istep,x_g,y_g,z_g,is_pack)
#endif
#if defined(_USE_ADIOS2)
    case(FILETYPE_ADIOS2)
      call write_subset_adios2(filename,varname,comm,ng,nh,lo,hi,lo_out,hi_out,nskip,var,time,istep,x_g,y_g,z_g,is_pack)
#endif
    case default
      call write_subset_mpiio(filename,varname,comm,ng,nh,lo,hi,lo_out,hi_out,nskip,var,time,istep,x_g,y_g,z_g,is_pack)
    end select
    call MPI_BARRIER(comm,ierr)
  end subroutine io_write_subset
  !
  subroutine subset_get_local_extents(lo,hi,lo_out,hi_out,nskip,lo_sel,hi_sel,n_sel,lo_pack)
    !
    ! computes the local portion of the requested structured subset
    !
    implicit none
    integer, intent(in), dimension(3) :: lo,hi,lo_out,hi_out,nskip
    integer, intent(out), dimension(3) :: lo_sel,hi_sel,n_sel,lo_pack
    integer :: idir,gmin,gmax,first,last,rel
    !
    lo_sel(:) = 0
    hi_sel(:) = -1
    n_sel(:) = 0
    lo_pack(:) = 0
    do idir=1,3
      gmin = max(lo(idir),lo_out(idir))
      gmax = min(hi(idir),hi_out(idir))
      if(gmin > gmax) then
        exit
      end if
      rel = gmin - lo_out(idir)
      first = lo_out(idir) + ((rel + nskip(idir) - 1)/nskip(idir))*nskip(idir)
      if(first > gmax) then
        exit
      end if
      last = lo_out(idir) + ((gmax - lo_out(idir))/nskip(idir))*nskip(idir)
      lo_sel(idir) = first
      hi_sel(idir) = last
      n_sel(idir) = (last-first)/nskip(idir) + 1
      lo_pack(idir) = (first-lo_out(idir))/nskip(idir) + 1
    end do
  end subroutine subset_get_local_extents
  !
  subroutine subset_stage_field(lo,hi,nh,lo_sel,hi_sel,lo_pack,hi_pack,nskip,var,buf)
    !
    ! stages a local compact subset block, using direct copy for unit skip
    ! and explicit packing for strided selections
    !
    implicit none
    integer, intent(in), dimension(3) :: lo,hi,nh,lo_sel,hi_sel,lo_pack,hi_pack,nskip
    real(rp), intent(in), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: var
    real(rp), allocatable, intent(out), dimension(:,:,:) :: buf
    integer :: i,j,k,ii,jj,kk
    !
    allocate(buf(lo_pack(1):hi_pack(1),lo_pack(2):hi_pack(2),lo_pack(3):hi_pack(3)))
    if(all(nskip(:) == 1)) then
      buf(:,:,:) = var(lo_sel(1):hi_sel(1),lo_sel(2):hi_sel(2),lo_sel(3):hi_sel(3))
    else
      kk = lo_pack(3)-1
      do k=lo_sel(3),hi_sel(3),nskip(3)
        kk = kk + 1
        jj = lo_pack(2)-1
        do j=lo_sel(2),hi_sel(2),nskip(2)
          jj = jj + 1
          ii = lo_pack(1)-1
          do i=lo_sel(1),hi_sel(1),nskip(1)
            ii = ii + 1
            buf(ii,jj,kk) = var(i,j,k)
          end do
        end do
      end do
    end if
  end subroutine subset_stage_field
  !
  subroutine write_subset_mpiio(filename,varname,comm,ng,nh,lo,hi,lo_out,hi_out,nskip,var,time,istep,x_g,y_g,z_g,is_pack)
    !
    ! writes a structured subset to a raw binary file in parallel
    !
    implicit none
    character(len=*), intent(in) :: filename,varname
    integer         , intent(in) :: comm
    integer         , intent(in), dimension(3) :: ng,nh,lo,hi,lo_out,hi_out,nskip
    real(rp), intent(in), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: var
    real(rp), intent(in), optional :: time
    integer , intent(in), optional :: istep
    real(rp), intent(in), dimension(1-nh(1):), optional :: x_g
    real(rp), intent(in), dimension(1-nh(2):), optional :: y_g
    real(rp), intent(in), dimension(1-nh(3):), optional :: z_g
    logical , intent(in), optional :: is_pack
    real(rp), allocatable, dimension(:,:,:) :: buf
    real(rp), allocatable, dimension(:) :: x_sel,y_sel,z_sel
    integer, dimension(3) :: lo_sel,hi_sel,n_sel,lo_pack,hi_pack,ng_out
    integer :: fh,color,write_comm
    integer(MPI_OFFSET_KIND) :: filesize
    logical :: has_data,is_pack_,is_full_field
    !
    ng_out(:) = (hi_out(:)-lo_out(:))/nskip(:) + 1
    call subset_get_local_extents(lo,hi,lo_out,hi_out,nskip,lo_sel,hi_sel,n_sel,lo_pack)
    has_data = all(n_sel(:) > 0)
    is_pack_ = .not. all(nskip(:) == 1)
    if(present(is_pack)) is_pack_ = is_pack
    is_full_field = all(nskip(:) == 1) .and. all(lo_out(:) == 1) .and. all(hi_out(:) == ng(:))
    color = merge(1,MPI_UNDEFINED,has_data)
    call MPI_COMM_SPLIT(comm,color,0,write_comm,ierr)
    if(has_data) then
      hi_pack(:) = lo_pack(:) + n_sel(:) - 1
      call MPI_FILE_OPEN(write_comm,filename,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
      filesize = 0_MPI_OFFSET_KIND
      call MPI_FILE_SET_SIZE(fh,filesize,ierr)
      filesize = 0_MPI_OFFSET_KIND
      if(is_pack_) then
        call subset_stage_field(lo,hi,nh,lo_sel,hi_sel,lo_pack,hi_pack,nskip,var,buf)
        call io_field('w',fh,ng_out,[0,0,0],lo_pack,hi_pack,filesize,buf)
      else if(is_full_field) then
        call io_field('w',fh,ng_out,nh,lo,hi,filesize,var)
      else
        call io_field('w',fh,ng_out,[0,0,0],lo_pack,hi_pack,filesize, &
                      var(lo_sel(1):hi_sel(1),lo_sel(2):hi_sel(2),lo_sel(3):hi_sel(3)))
      end if
      call MPI_FILE_CLOSE(fh,ierr)
      call MPI_COMM_FREE(write_comm,ierr)
    end if
  end subroutine write_subset_mpiio
  !
#if defined(_USE_HDF5)
  subroutine write_subset_hdf5(filename,varname,comm,ng,nh,lo,hi,lo_out,hi_out,nskip,var,time,istep,x_g,y_g,z_g,is_pack)
    use hdf5
    use mod_param, only: is_use_compression
    !
    ! writes a structured subset to an HDF5 file in parallel
    !
    implicit none
    character(len=*), intent(in) :: filename,varname
    integer         , intent(in) :: comm
    integer         , intent(in), dimension(3) :: ng,nh,lo,hi,lo_out,hi_out,nskip
    real(rp), intent(in), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: var
    real(rp), intent(in), optional :: time
    integer , intent(in), optional :: istep
    real(rp), intent(in), dimension(1-nh(1):), optional :: x_g
    real(rp), intent(in), dimension(1-nh(2):), optional :: y_g
    real(rp), intent(in), dimension(1-nh(3):), optional :: z_g
    logical , intent(in), optional :: is_pack
    real(rp), allocatable, dimension(:,:,:) :: buf
    real(rp), allocatable, dimension(:) :: x_sel,y_sel,z_sel
    integer, dimension(3) :: ng_out,lo_sel,hi_sel,n_sel,lo_pack,hi_pack
    real(rp) :: time_
    integer(HID_T) :: file_id,group_id,dset,space
    integer(HSIZE_T), dimension(1) :: dims
    integer :: istep_,color,write_comm,myid_loc
    integer :: ierr_local
    logical :: has_data,is_pack_,is_full_field
    !
    ng_out(:) = (hi_out(:)-lo_out(:))/nskip(:) + 1
    call subset_get_local_extents(lo,hi,lo_out,hi_out,nskip,lo_sel,hi_sel,n_sel,lo_pack)
    has_data = all(n_sel(:) > 0)
    is_pack_ = .not. all(nskip(:) == 1)
    if(present(is_pack)) is_pack_ = is_pack
    is_full_field = all(nskip(:) == 1) .and. all(lo_out(:) == 1) .and. all(hi_out(:) == ng(:))
    color = merge(1,MPI_UNDEFINED,has_data)
    call MPI_COMM_SPLIT(comm,color,0,write_comm,ierr)
    if(has_data) then
      call MPI_COMM_RANK(write_comm,myid_loc,ierr)
      hi_pack(:) = lo_pack(:) + n_sel(:) - 1
      time_ = 0._rp
      istep_ = 0
      if(present(time)) time_ = time
      if(present(istep)) istep_ = istep
      if(present(x_g)) then
        allocate(x_sel(1:ng_out(1)))
        x_sel(1:ng_out(1)) = x_g(lo_out(1):hi_out(1):nskip(1))
      end if
      if(present(y_g)) then
        allocate(y_sel(1:ng_out(2)))
        y_sel(1:ng_out(2)) = y_g(lo_out(2):hi_out(2):nskip(2))
      end if
      if(present(z_g)) then
        allocate(z_sel(1:ng_out(3)))
        z_sel(1:ng_out(3)) = z_g(lo_out(3):hi_out(3):nskip(3))
      end if
      if(is_pack_) call subset_stage_field(lo,hi,nh,lo_sel,hi_sel,lo_pack,hi_pack,nskip,var,buf)
      if(present(x_g) .and. present(y_g) .and. present(z_g)) then
        if(is_pack_) then
          call io_field_hdf5('w',filename,trim(varname),write_comm,ng_out,[0,0,0],lo_pack,hi_pack,buf,time_,istep_, &
                             x_sel,y_sel,z_sel,first_write=.true.,is_compress=is_use_compression)
        else if(is_full_field) then
          call io_field_hdf5('w',filename,trim(varname),write_comm,ng_out,nh,lo,hi,var,time_,istep_, &
                             x_sel,y_sel,z_sel,first_write=.true.,is_compress=is_use_compression)
        else
          call io_field_hdf5('w',filename,trim(varname),write_comm,ng_out,[0,0,0],lo_pack,hi_pack, &
                             var(lo_sel(1):hi_sel(1),lo_sel(2):hi_sel(2),lo_sel(3):hi_sel(3)),time_,istep_, &
                             x_sel,y_sel,z_sel,first_write=.true.,is_compress=is_use_compression)
        end if
      else
        if(is_pack_) then
          call io_field_hdf5('w',filename,trim(varname),write_comm,ng_out,[0,0,0],lo_pack,hi_pack,buf,time_,istep_, &
                             first_write=.true.,is_compress=is_use_compression)
        else if(is_full_field) then
          call io_field_hdf5('w',filename,trim(varname),write_comm,ng_out,nh,lo,hi,var,time_,istep_, &
                             first_write=.true.,is_compress=is_use_compression)
        else
          call io_field_hdf5('w',filename,trim(varname),write_comm,ng_out,[0,0,0],lo_pack,hi_pack, &
                             var(lo_sel(1):hi_sel(1),lo_sel(2):hi_sel(2),lo_sel(3):hi_sel(3)),time_,istep_, &
                             first_write=.true.,is_compress=is_use_compression)
        end if
      end if
      if(myid_loc == 0) then
        dims(1) = 3_HSIZE_T
        call h5open_f(ierr_local)
        call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,ierr_local)
        call h5gopen_f(file_id,'meta',group_id,ierr_local)
        call h5screate_simple_f(1,dims,space,ierr_local)
        call h5dcreate_f(group_id,'lo',H5T_NATIVE_INTEGER,space,dset,ierr_local)
        call h5dwrite_f(dset,H5T_NATIVE_INTEGER,lo_out,dims,ierr_local)
        call h5dclose_f(dset,ierr_local)
        call h5dcreate_f(group_id,'hi',H5T_NATIVE_INTEGER,space,dset,ierr_local)
        call h5dwrite_f(dset,H5T_NATIVE_INTEGER,hi_out,dims,ierr_local)
        call h5dclose_f(dset,ierr_local)
        call h5dcreate_f(group_id,'nskip',H5T_NATIVE_INTEGER,space,dset,ierr_local)
        call h5dwrite_f(dset,H5T_NATIVE_INTEGER,nskip,dims,ierr_local)
        call h5dclose_f(dset,ierr_local)
        call h5sclose_f(space,ierr_local)
        call h5gclose_f(group_id,ierr_local)
        call h5fclose_f(file_id,ierr_local)
        call h5close_f(ierr_local)
      end if
      call MPI_COMM_FREE(write_comm,ierr)
    end if
  end subroutine write_subset_hdf5
  !
#endif
#if defined(_USE_ADIOS2)
  subroutine write_subset_adios2(filename,varname,comm,ng,nh,lo,hi,lo_out,hi_out,nskip,var,time,istep,x_g,y_g,z_g,is_pack)
    use mod_param, only: is_use_compression
    !
    ! writes a structured subset to an ADIOS2 file in parallel
    !
    implicit none
    character(len=*), intent(in) :: filename,varname
    integer         , intent(in) :: comm
    integer         , intent(in), dimension(3) :: ng,nh,lo,hi,lo_out,hi_out,nskip
    real(rp), intent(in), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: var
    real(rp), intent(in), optional :: time
    integer , intent(in), optional :: istep
    real(rp), intent(in), optional, dimension(1-nh(1):) :: x_g
    real(rp), intent(in), optional, dimension(1-nh(2):) :: y_g
    real(rp), intent(in), optional, dimension(1-nh(3):) :: z_g
    logical , intent(in), optional :: is_pack
    real(rp), allocatable, dimension(:,:,:) :: buf
    real(rp), allocatable, dimension(:) :: x_sel,y_sel,z_sel
    type(adios2_adios) :: adios
    type(adios2_io) :: io_handle
    type(adios2_engine) :: engine
    type(adios2_operator) :: compress
    type(adios2_compression_params) :: compress_params
    integer, dimension(3) :: ng_out,lo_sel,hi_sel,n_sel,lo_pack,hi_pack
    real(rp) :: meta_rp(1)
    integer :: meta_i4(1)
    integer :: meta_lo(3),meta_hi(3),meta_skip(3)
    real(rp) :: time_
    integer :: istep_,color,write_comm,myid_loc
    logical :: is_compress,has_data,is_pack_,is_full_field
    !
    ng_out(:) = (hi_out(:)-lo_out(:))/nskip(:) + 1
    call subset_get_local_extents(lo,hi,lo_out,hi_out,nskip,lo_sel,hi_sel,n_sel,lo_pack)
    has_data = all(n_sel(:) > 0)
    is_pack_ = .not. all(nskip(:) == 1)
    if(present(is_pack)) is_pack_ = is_pack
    is_full_field = all(nskip(:) == 1) .and. all(lo_out(:) == 1) .and. all(hi_out(:) == ng(:))
    color = merge(1,MPI_UNDEFINED,has_data)
    call MPI_COMM_SPLIT(comm,color,0,write_comm,ierr)
    if(.not. has_data) return
    hi_pack(:) = lo_pack(:) + n_sel(:) - 1
    if(is_pack_) call subset_stage_field(lo,hi,nh,lo_sel,hi_sel,lo_pack,hi_pack,nskip,var,buf)
    time_ = 0._rp
    istep_ = 0
    if(present(time)) time_ = time
    if(present(istep)) istep_ = istep
    if(present(x_g)) then
      allocate(x_sel(1:ng_out(1)))
      x_sel(1:ng_out(1)) = x_g(lo_out(1):hi_out(1):nskip(1))
    end if
    if(present(y_g)) then
      allocate(y_sel(1:ng_out(2)))
      y_sel(1:ng_out(2)) = y_g(lo_out(2):hi_out(2):nskip(2))
    end if
    if(present(z_g)) then
      allocate(z_sel(1:ng_out(3)))
      z_sel(1:ng_out(3)) = z_g(lo_out(3):hi_out(3):nskip(3))
    end if
     call adios2_open_engine('w',adios,io_handle,engine,filename,write_comm)
    is_compress = is_use_compression
    if(is_compress) then
      compress_params = adios_get_default_compression_params('blosc')
      is_pack_ = .true.
      if(.not. allocated(buf)) call subset_stage_field(lo,hi,nh,lo_sel,hi_sel,lo_pack,hi_pack,nskip,var,buf)
      call adios2_define_compress(adios,compress,compress_params)
      !
      ! Feed ADIOS2 compression with a canonical packed array starting at 1
      ! instead of relying on memory selections over an offset buffer due to
      ! a previously-described ADIOS2 2.11.0 limitation that will be fixed soon
      !
      call io_field_adios2('w',engine,io_handle,trim(varname),ng_out,[0,0,0],lo_pack,hi_pack,buf,.true.,compress,compress_params)
    else
      if(is_pack_) then
        call io_field_adios2('w',engine,io_handle,trim(varname),ng_out,[0,0,0],lo_pack,hi_pack,buf)
      else if(is_full_field) then
        call io_field_adios2('w',engine,io_handle,trim(varname),ng_out,nh,lo,hi,var,.false.)
      else
        call io_field_adios2('w',engine,io_handle,trim(varname),ng_out,[0,0,0],lo_pack,hi_pack, &
                             var(lo_sel(1):hi_sel(1),lo_sel(2):hi_sel(2),lo_sel(3):hi_sel(3)))
      end if
    end if
    meta_rp(1) = time_
    meta_i4(1) = istep_
    meta_lo(:) = lo_out(:)
    meta_hi(:) = hi_out(:)
    meta_skip(:) = nskip(:)
    call io_field_adios2_1d('w',engine,io_handle,'time',meta_rp,write_comm)
    call io_field_adios2_1d('w',engine,io_handle,'istep',meta_i4,write_comm)
    call io_field_adios2_1d('w',engine,io_handle,'lo',meta_lo,write_comm)
    call io_field_adios2_1d('w',engine,io_handle,'hi',meta_hi,write_comm)
    call io_field_adios2_1d('w',engine,io_handle,'nskip',meta_skip,write_comm)
    if(present(x_g)) call io_field_adios2_1d('w',engine,io_handle,'x',x_sel,write_comm)
    if(present(y_g)) call io_field_adios2_1d('w',engine,io_handle,'y',y_sel,write_comm)
    if(present(z_g)) call io_field_adios2_1d('w',engine,io_handle,'z',z_sel,write_comm)
    call adios2_close_engine(adios,engine)
    call MPI_COMM_FREE(write_comm,ierr)
  end subroutine write_subset_adios2
  !
#endif
  subroutine load_all_local(io,filename,n,nh,u,v,w,p,time,istep)
    !
    ! reads/writes a restart file
    !
    implicit none
    character(len=1), intent(in) :: io
    character(len=*), intent(in) :: filename
    integer , intent(in), dimension(3) :: n,nh
    real(rp), intent(inout), dimension(1-nh(1):,1-nh(2):,1-nh(3):) :: u,v,w,p
    real(rp), intent(inout) :: time
    integer , intent(inout) :: istep
    real(rp), dimension(2) :: fldinfo
    integer :: iunit
    !
    select case(io)
    case('r')
      open(newunit=iunit,file=filename,action='read' ,access='stream',form='unformatted',status='old'    )
      read(iunit) u(1:n(1),1:n(2),1:n(3)), &
                  v(1:n(1),1:n(2),1:n(3)), &
                  w(1:n(1),1:n(2),1:n(3)), &
                  p(1:n(1),1:n(2),1:n(3))
      read(iunit) fldinfo(1:2)
      close(iunit)
      time = fldinfo(1)
      istep = nint(fldinfo(2))
    case('w')
      !
      ! write
      !
      fldinfo = [time,1._rp*istep]
      open(newunit=iunit,file=filename,action='write',access='stream',form='unformatted',status='replace')
      write(iunit) u(1:n(1),1:n(2),1:n(3)), &
                   v(1:n(1),1:n(2),1:n(3)), &
                   w(1:n(1),1:n(2),1:n(3)), &
                   p(1:n(1),1:n(2),1:n(3))
      write(iunit) fldinfo(1:2)
      close(iunit)
    end select
  end subroutine load_all_local
  !
  subroutine load_one_local(io,filename,n,nh,p,time,istep)
    !
    ! reads/writes a restart file
    !
    implicit none
    character(len=1), intent(in) :: io
    character(len=*), intent(in) :: filename
    integer , intent(in), dimension(3) :: n,nh
    real(rp), intent(inout), dimension(1-nh(1):,1-nh(2):,1-nh(3):) :: p
    real(rp), intent(inout) :: time
    integer , intent(inout) :: istep
    real(rp), dimension(2) :: fldinfo
    integer :: iunit
    !
    select case(io)
    case('r')
      open(newunit=iunit,file=filename,action='read' ,access='stream',form='unformatted',status='old'    )
      read(iunit) p(1:n(1),1:n(2),1:n(3))
      read(iunit) fldinfo(1:2)
      close(iunit)
      time = fldinfo(1)
      istep = nint(fldinfo(2))
    case('w')
      !
      ! write
      !
      fldinfo = [time,1._rp*istep]
      open(newunit=iunit,file=filename,action='write',access='stream',form='unformatted',status='replace')
      write(iunit) p(1:n(1),1:n(2),1:n(3))
      write(iunit) fldinfo(1:2)
      close(iunit)
    end select
  end subroutine load_one_local
  !
  integer function fetch_file_type(filename) result(file_type)
    !
    ! returns the backend implied by the file suffix, defaulting to raw binary
    !
    implicit none
    character(len=*), intent(in) :: filename
    integer :: n,idot
    character(len=6) :: ext
    n = len_trim(filename)
    file_type = FILETYPE_MPIIO
    if(n > 0) then
      idot = scan(filename(1:n),'.',back=.true.)
      if(idot > 0 .and. idot < n) then
        ext = filename(idot:n)
        select case(trim(ext))
        case('.bin')
          file_type = FILETYPE_MPIIO
        case('.h5','.hdf','.hdf5','.nc')
          file_type = FILETYPE_HDF5
        case('.bp')
          file_type = FILETYPE_ADIOS2
        end select
      end if
    end if
  end function fetch_file_type
  !
#if defined(_USE_HDF5)
  subroutine load_all_hdf5(io,filename,comm,ng,nh,lo,hi,nscal,u,v,w,p,s,time,istep,x_g,y_g,z_g)
    use mod_param, only: is_use_compression
    !
    ! reads/writes a restart file for all fields using HDF5
    !
    implicit none
    character(len=1), intent(in) :: io
    character(len=*), intent(in) :: filename
    integer         , intent(in) :: comm
    integer , intent(in), dimension(3) :: ng,nh,lo,hi
    integer , intent(in)               :: nscal
    real(rp), intent(inout), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: u,v,w,p
    type(scalar), intent(inout), dimension(:) :: s
    real(rp), intent(inout) :: time
    integer , intent(inout) :: istep
    real(rp), intent(inout), dimension(1-nh(1):), optional :: x_g
    real(rp), intent(inout), dimension(1-nh(2):), optional :: y_g
    real(rp), intent(inout), dimension(1-nh(3):), optional :: z_g
    character(len=5) :: scalnum
    integer :: iscal
    !
    select case(io)
    case('r')
      call io_field_hdf5(io,filename,'u',comm,ng,nh,lo,hi,u,time,istep,x_g,y_g,z_g)
      call io_field_hdf5(io,filename,'v',comm,ng,nh,lo,hi,v)
      call io_field_hdf5(io,filename,'w',comm,ng,nh,lo,hi,w)
      call io_field_hdf5(io,filename,'p',comm,ng,nh,lo,hi,p)
      do iscal=1,nscal
        write(scalnum,'(i3.3)') iscal
        call io_field_hdf5(io,filename,'s_'//scalnum,comm,ng,nh,lo,hi,s(iscal)%val)
      end do
    case('w')
      call io_field_hdf5(io,filename,'u',comm,ng,nh,lo,hi,u,time,istep,x_g,y_g,z_g, &
                         first_write=.true.,is_compress=is_use_compression)
      call io_field_hdf5(io,filename,'v',comm,ng,nh,lo,hi,v,first_write=.false.,is_compress=is_use_compression)
      call io_field_hdf5(io,filename,'w',comm,ng,nh,lo,hi,w,first_write=.false.,is_compress=is_use_compression)
      call io_field_hdf5(io,filename,'p',comm,ng,nh,lo,hi,p,first_write=.false.,is_compress=is_use_compression)
      do iscal=1,nscal
        write(scalnum,'(i3.3)') iscal
        call io_field_hdf5(io,filename,'s_'//scalnum,comm,ng,nh,lo,hi,s(iscal)%val, &
                           first_write=.false.,is_compress=is_use_compression)
      end do
    end select
  end subroutine load_all_hdf5
  !
  subroutine load_one_hdf5(io,filename,comm,ng,nh,lo,hi,p,varname,time,istep,x_g,y_g,z_g)
    use mod_param, only: is_use_compression
    !
    ! reads/writes a restart file for a single field using HDF5
    !
    implicit none
    integer, parameter :: NAME_LEN_MAX = 256
    character(len=1), intent(in) :: io
    character(len=*), intent(in) :: filename
    integer         , intent(in) :: comm
    integer , intent(in), dimension(3) :: ng,nh,lo,hi
    real(rp), intent(inout), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: p
    character(len=*), intent(in) :: varname
    real(rp), intent(inout), optional :: time
    integer , intent(inout), optional :: istep
    real(rp), intent(inout), dimension(1-nh(1):), optional :: x_g
    real(rp), intent(inout), dimension(1-nh(2):), optional :: y_g
    real(rp), intent(inout), dimension(1-nh(3):), optional :: z_g
    character(len=NAME_LEN_MAX) :: field_name
    !
    field_name = trim(varname)
    select case(io)
    case('r')
      call io_field_hdf5(io,filename,field_name,comm,ng,nh,lo,hi,p,time,istep,x_g,y_g,z_g)
    case('w')
      ! single-field writes always start from a fresh file
      call io_field_hdf5(io,filename,field_name,comm,ng,nh,lo,hi,p,time,istep,x_g,y_g,z_g, &
                         first_write=.true.,is_compress=is_use_compression)
    end select
  end subroutine load_one_hdf5
  !
  subroutine hdf5_checkpoint_chunk_dims(comm,n,ng,chunk_dims)
    use hdf5, only: HSIZE_T
    !
    ! choose chunk sizes aligned with the current pencil layout
    !
    implicit none
    integer(i8), parameter :: MEGABYTE = 1024_i8*1024_i8
    integer(i8), parameter :: HDF5_CHUNK_MAX = 2._i8*MEGABYTE ! 2 MiB
    integer         , intent(in) :: comm
    integer         , intent(in),  dimension(3) :: n,ng
    integer(HSIZE_T), intent(out), dimension(3) :: chunk_dims
    integer, dimension(3) :: nn
    integer(i8) :: chunk_size
    integer :: idir
    logical :: is_reduced
    !
    call MPI_ALLREDUCE(n,nn,3,MPI_INTEGER,MPI_MIN,comm,ierr)
    !
    chunk_size = product(int(nn,i8))*int(f_sizeof(1._rp),i8)
    do while(chunk_size > HDF5_CHUNK_MAX)
      is_reduced = .false.
      do idir=1,3
        if(nn(idir) > 1) then
          nn(idir) = max(1,nn(idir)/2)
          is_reduced = .true.
          chunk_size = product(int(nn,i8))*int(f_sizeof(1._rp),i8)
          if(chunk_size <= HDF5_CHUNK_MAX) exit
        end if
      end do
      if(.not. is_reduced) exit
    end do
    chunk_dims(:) = nn(:)
  end subroutine hdf5_checkpoint_chunk_dims
  !
  subroutine io_field_hdf5(io,filename,varname,comm,ng,nh,lo,hi,var,time,istep,x_g,y_g,z_g,first_write,is_compress)
    use hdf5
    !
    ! collective single field data I/O using HDF5
    !
    ! written with the help of Josh Romero,
    ! with the field data I/O inspired from the AFiD code
    !
    implicit none
    character(len=*), parameter :: HDF5_COMPRESSION = 'gzip'
    integer         , parameter :: HDF5_COMPRESSION_LEVEL = 4
    character(len=1), intent(in) :: io
    character(len=*), intent(in) :: filename,varname
    integer         , intent(in) :: comm
    integer         , intent(in), dimension(3)   :: ng,nh,lo,hi
    real(rp), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: var
    real(rp), intent(inout), optional :: time
    integer , intent(inout), optional :: istep
    real(rp), intent(inout), dimension(:), optional :: x_g
    real(rp), intent(inout), dimension(:), optional :: y_g
    real(rp), intent(inout), dimension(:), optional :: z_g
    logical , intent(in), optional :: first_write
    logical , intent(in), optional :: is_compress
    integer , dimension(3) :: n
    integer , dimension(3) :: sizes,subsizes,starts
    logical :: first_write_
    real(rp) :: meta_time(1)
    integer  :: meta_istep(1)
    !
    ! HDF5 variables
    !
    integer :: ndims,ierr,myid_loc
    integer(HID_T) :: file_id,group_id
    integer(HID_T) :: filespace
    integer(HID_T) :: slabspace
    integer(HID_T) :: memspace
    !
    integer(HID_T) :: dset
    integer(HID_T) :: dtype_rp,dtype_int
    integer(HID_T) :: dcpl_id
    !
    integer(HSIZE_T) :: dims(3)
    integer(HSIZE_T) :: chunk_dims(3)
    !
    integer(HID_T) :: plist_id
    integer(HSIZE_T) , dimension(3) :: data_count
    integer(HSIZE_T) , dimension(3) :: mem_dims
    integer(HSSIZE_T), dimension(3) :: data_offset
    integer(HSSIZE_T), dimension(3) :: halo_offset
    logical :: use_dcpl
    logical :: is_compress_
    !
    n(:)        = hi(:)-lo(:)+1
    sizes(:)    = ng(:)
    subsizes(:) = n(:)
    starts(:)   = lo(:) - 1 ! starts from 0
    !
    ndims = 3
    dims(:) = ng(:)
    data_count(:) = subsizes(:)
    mem_dims(:) = subsizes(:)+2*nh(:)
    data_offset(:) = starts(:)
    halo_offset(:) = nh(:)
    first_write_ = .true.
    if(present(first_write)) first_write_ = first_write
    is_compress_ = .false.
    if(present(is_compress)) is_compress_ = is_compress
    call MPI_COMM_RANK(comm,myid_loc,ierr)
    call h5open_f(ierr)
    dtype_rp = HDF5_REAL_RP()
    dtype_int = H5T_NATIVE_INTEGER
    !
    select case(io)
    case('r')
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,ierr)
      call h5pset_fapl_mpio_f(plist_id,comm,MPI_INFO_NULL,ierr)
      call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,ierr,access_prp=plist_id)
      call h5pclose_f(plist_id,ierr)
      !
      call h5dopen_f(file_id,'fields/'//trim(varname),dset,ierr)
      call h5screate_simple_f(3,mem_dims,memspace,ierr)
      !
      call h5dget_space_f(dset,slabspace,ierr)
      call h5sselect_hyperslab_f(slabspace,H5S_SELECT_SET_F,data_offset,data_count,ierr)
      call h5sselect_hyperslab_f(memspace,H5S_SELECT_SET_F,halo_offset,data_count,ierr)
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,ierr)
      call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,ierr)
      !
      call h5dread_f(dset,dtype_rp,var,mem_dims,ierr,file_space_id=slabspace,mem_space_id=memspace,xfer_prp=plist_id)
      !
      call h5pclose_f(plist_id,ierr)
      call h5dclose_f(dset,ierr)
      call h5sclose_f(slabspace,ierr)
      call h5sclose_f(memspace,ierr)
      call h5fclose_f(file_id,ierr)
      !
      if(myid_loc == 0 .and. (present(time) .or. present(istep) .or. present(x_g) .or. present(y_g) .or. present(z_g))) then
        call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,ierr)
        if(present(time)) then
          call h5dopen_f(file_id,'meta/time',dset,ierr)
          call h5dread_f(dset,dtype_rp,meta_time,[1_HSIZE_T],ierr)
          call h5dclose_f(dset,ierr)
          time = meta_time(1)
        end if
        if(present(istep)) then
          call h5dopen_f(file_id,'meta/istep',dset,ierr)
          call h5dread_f(dset,dtype_int,meta_istep,[1_HSIZE_T],ierr)
          call h5dclose_f(dset,ierr)
          istep = meta_istep(1)
        end if
        if(present(x_g)) then
          call h5dopen_f(file_id,'grid/x',dset,ierr)
          call h5dread_f(dset,dtype_rp,x_g(1:ng(1)),[int(ng(1),HSIZE_T)],ierr)
          call h5dclose_f(dset,ierr)
        end if
        if(present(y_g)) then
          call h5dopen_f(file_id,'grid/y',dset,ierr)
          call h5dread_f(dset,dtype_rp,y_g(1:ng(2)),[int(ng(2),HSIZE_T)],ierr)
          call h5dclose_f(dset,ierr)
        end if
        if(present(z_g)) then
          call h5dopen_f(file_id,'grid/z',dset,ierr)
          call h5dread_f(dset,dtype_rp,z_g(1:ng(3)),[int(ng(3),HSIZE_T)],ierr)
          call h5dclose_f(dset,ierr)
        end if
        call h5fclose_f(file_id,ierr)
      end if
      if(present(time)) call MPI_BCAST(time,1,MPI_REAL_RP,0,comm,ierr)
      if(present(istep)) call MPI_BCAST(istep,1,MPI_INTEGER,0,comm,ierr)
      if(present(x_g)) call MPI_BCAST(x_g(1:ng(1)),ng(1),MPI_REAL_RP,0,comm,ierr)
      if(present(y_g)) call MPI_BCAST(y_g(1:ng(2)),ng(2),MPI_REAL_RP,0,comm,ierr)
      if(present(z_g)) call MPI_BCAST(z_g(1:ng(3)),ng(3),MPI_REAL_RP,0,comm,ierr)
    case('w')
      call h5screate_simple_f(ndims,dims,filespace,ierr)
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,ierr)
      call h5pset_fapl_mpio_f(plist_id,comm,MPI_INFO_NULL,ierr)
      if(first_write_) then
        ! guarantees overwriting if file exists
        call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,ierr,access_prp=plist_id)
        call h5gcreate_f(file_id,'fields',group_id,ierr)
        call h5gclose_f(group_id,ierr)
      else
        call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,ierr,access_prp=plist_id)
      end if
      call h5pclose_f(plist_id,ierr)
      !
      use_dcpl = .false.
      if(is_compress_) then
        !
        ! determine chunk dimensions to target about 2 MiB per chunk
        !
        call hdf5_checkpoint_chunk_dims(comm,n,ng,chunk_dims)
        select case(HDF5_COMPRESSION)
        case('gzip')
          call h5pcreate_f(H5P_DATASET_CREATE_F,dcpl_id,ierr)
          call h5pset_chunk_f(dcpl_id,ndims,chunk_dims,ierr)
          call h5pset_shuffle_f(dcpl_id,ierr) ! comment to disable shuffle filter
          call h5pset_deflate_f(dcpl_id,HDF5_COMPRESSION_LEVEL,ierr)
          use_dcpl = .true.
        case default
          if(first_write_ .and. myid_loc == 0) then
            print*, 'Warning: unknown HDF5 checkpoint compression preset.'
            print*, 'Writing uncompressed fields.'
          end if
        end select
      end if
      if(use_dcpl) then
        call h5gopen_f(file_id,'fields',group_id,ierr)
        call h5dcreate_f(group_id,trim(varname),dtype_rp,filespace,dset,ierr,dcpl_id=dcpl_id)
        call h5gclose_f(group_id,ierr)
        call h5pclose_f(dcpl_id,ierr)
      else
        call h5gopen_f(file_id,'fields',group_id,ierr)
        call h5dcreate_f(group_id,trim(varname),dtype_rp,filespace,dset,ierr)
        call h5gclose_f(group_id,ierr)
      end if
      call h5screate_simple_f(3,mem_dims,memspace,ierr)
      call h5dget_space_f(dset,slabspace,ierr)
      call h5sselect_hyperslab_f(slabspace,H5S_SELECT_SET_F,data_offset,data_count,ierr)
      call h5sselect_hyperslab_f(memspace,H5S_SELECT_SET_F,halo_offset,data_count,ierr)
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,ierr)
      call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,ierr)
      call h5dwrite_f(dset,dtype_rp,var,mem_dims,ierr,file_space_id=slabspace,mem_space_id=memspace,xfer_prp=plist_id)
      !
      call h5pclose_f(plist_id,ierr)
      call h5dclose_f(dset,ierr)
      call h5sclose_f(memspace,ierr)
      call h5sclose_f(slabspace,ierr)
      call h5sclose_f(filespace,ierr)
      call h5fclose_f(file_id,ierr)
      !
      ! write metadata
      !
      if(myid_loc == 0) then
        call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,ierr)
        !
        if(first_write_ .and. present(x_g) .and. present(y_g) .and. present(z_g)) then
          call h5gcreate_f(file_id,'grid',group_id,ierr)
          call h5screate_simple_f(1,[int(size(x_g),HSIZE_T)],filespace,ierr)
          call h5dcreate_f(group_id,'x',dtype_rp,filespace,dset,ierr)
          call h5dwrite_f(dset,dtype_rp,x_g,[int(size(x_g),HSIZE_T)],ierr)
          call h5dclose_f(dset,ierr)
          call h5sclose_f(filespace,ierr)
          call h5screate_simple_f(1,[int(size(y_g),HSIZE_T)],filespace,ierr)
          call h5dcreate_f(group_id,'y',dtype_rp,filespace,dset,ierr)
          call h5dwrite_f(dset,dtype_rp,y_g,[int(size(y_g),HSIZE_T)],ierr)
          call h5dclose_f(dset,ierr)
          call h5sclose_f(filespace,ierr)
          call h5screate_simple_f(1,[int(size(z_g),HSIZE_T)],filespace,ierr)
          call h5dcreate_f(group_id,'z',dtype_rp,filespace,dset,ierr)
          call h5dwrite_f(dset,dtype_rp,z_g,[int(size(z_g),HSIZE_T)],ierr)
          call h5dclose_f(dset,ierr)
          call h5sclose_f(filespace,ierr)
          call h5gclose_f(group_id,ierr)
        end if
        !
        if(first_write_ .and. (present(time) .or. present(istep))) then
          call h5gcreate_f(file_id,'meta',group_id,ierr)
          call h5gclose_f(group_id,ierr)
        end if
        if(first_write_ .and. present(time)) then
          call h5gopen_f(file_id,'meta',group_id,ierr)
          call h5screate_simple_f(1,[1_HSIZE_T],filespace,ierr)
          meta_time(1) = time
          call h5dcreate_f(group_id,'time',dtype_rp,filespace,dset,ierr)
          call h5dwrite_f(dset,dtype_rp,meta_time,[1_HSIZE_T],ierr)
          call h5dclose_f(dset,ierr)
          call h5sclose_f(filespace,ierr)
          call h5gclose_f(group_id,ierr)
        end if
        if(first_write_ .and. present(istep)) then
          call h5gopen_f(file_id,'meta',group_id,ierr)
          call h5screate_simple_f(1,[1_HSIZE_T],filespace,ierr)
          meta_istep(1) = istep
          call h5dcreate_f(group_id,'istep',dtype_int,filespace,dset,ierr)
          call h5dwrite_f(dset,dtype_int,meta_istep,[1_HSIZE_T],ierr)
          call h5dclose_f(dset,ierr)
          call h5sclose_f(filespace,ierr)
          call h5gclose_f(group_id,ierr)
        end if
        call h5fclose_f(file_id,ierr)
      end if
    end select
    call h5close_f(ierr)
  end subroutine io_field_hdf5
  !
  integer(HID_T) function HDF5_REAL_RP()
    !
    ! returns the HDF5 type matching rp
    !
    use hdf5
    implicit none
    if(rp == dp) then
      HDF5_REAL_RP = H5T_NATIVE_DOUBLE
    else
      HDF5_REAL_RP = H5T_NATIVE_REAL
    end if
  end function HDF5_REAL_RP
#endif
#if defined(_USE_ADIOS2)
  subroutine load_all_adios2(io,filename,comm,ng,nh,lo,hi,nscal,u,v,w,p,s,time,istep,x_g,y_g,z_g)
    use mod_param, only: is_use_compression
    !
    ! reads/writes a restart file for all fields using ADIOS2
    !
    implicit none
    character(len=1), intent(in) :: io
    character(len=*), intent(in) :: filename
    integer         , intent(in) :: comm
    integer , intent(in), dimension(3) :: ng,nh,lo,hi
    integer , intent(in)               :: nscal
    real(rp), intent(inout), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: u,v,w,p
    type(scalar), intent(inout), dimension(:) :: s
    real(rp), intent(inout) :: time
    integer , intent(inout) :: istep
    real(rp), intent(inout), dimension(1-nh(1):), optional :: x_g
    real(rp), intent(inout), dimension(1-nh(2):), optional :: y_g
    real(rp), intent(inout), dimension(1-nh(3):), optional :: z_g
    type(adios2_adios) :: adios
    type(adios2_io) :: io_handle
    type(adios2_engine) :: engine
    type(adios2_operator) :: compress
    type(adios2_compression_params) :: compress_params
    character(len=5) :: scalnum
    real(rp) :: meta_rp(1)
    integer :: meta_i4(1)
    integer :: iscal
    logical :: is_compress,is_pack
    !
    call adios2_open_engine(io,adios,io_handle,engine,filename,comm)
    is_compress = (io == 'w' .and. is_use_compression)
    !
    ! BP5 in ADIOS2 2.11.0 is unable to combine compression with SetMemorySelection to exclude halo regions when writing.
    ! this has been fixed and should be available in an upcoming release (seettps://github.com/ornladios/ADIOS2/issues/4965)
    ! as a temporary workaround, we pack the data into a contiguous buffer array for writing.
    !
    !is_pack = is_use_compression
    !
    ! ADIOS2 2.9.x on Ubuntu can mis-handle SetMemorySelection on haloed restart arrays even without compression.
    ! use the packed path for restart fields to keep read/write behavior consistent across versions.
    !
    is_pack = .true.
    if(is_compress) then
      compress_params = adios_get_default_compression_params('blosc')
      call adios2_define_compress(adios,compress,compress_params)
    end if
    !
    select case(io)
    case('r')
      call io_field_adios2(io,engine,io_handle,'u',ng,nh,lo,hi,u,is_pack)
      call io_field_adios2(io,engine,io_handle,'v',ng,nh,lo,hi,v,is_pack)
      call io_field_adios2(io,engine,io_handle,'w',ng,nh,lo,hi,w,is_pack)
      call io_field_adios2(io,engine,io_handle,'p',ng,nh,lo,hi,p,is_pack)
      do iscal=1,nscal
        write(scalnum,'(i3.3)') iscal
        call io_field_adios2(io,engine,io_handle,'s_'//scalnum,ng,nh,lo,hi,s(iscal)%val,is_pack)
      end do
      call io_field_adios2_1d(io,engine,io_handle,'time',meta_rp,comm)
      call io_field_adios2_1d(io,engine,io_handle,'istep',meta_i4,comm)
      if(present(x_g)) call io_field_adios2_1d(io,engine,io_handle,'x',x_g(1:ng(1)),comm)
      if(present(y_g)) call io_field_adios2_1d(io,engine,io_handle,'y',y_g(1:ng(2)),comm)
      if(present(z_g)) call io_field_adios2_1d(io,engine,io_handle,'z',z_g(1:ng(3)),comm)
      time = meta_rp(1)
      istep = meta_i4(1)
    case('w')
      if(is_compress) then
        call io_field_adios2(io,engine,io_handle,'u',ng,nh,lo,hi,u,is_pack,compress,compress_params)
        call io_field_adios2(io,engine,io_handle,'v',ng,nh,lo,hi,v,is_pack,compress,compress_params)
        call io_field_adios2(io,engine,io_handle,'w',ng,nh,lo,hi,w,is_pack,compress,compress_params)
        call io_field_adios2(io,engine,io_handle,'p',ng,nh,lo,hi,p,is_pack,compress,compress_params)
      else
        call io_field_adios2(io,engine,io_handle,'u',ng,nh,lo,hi,u,is_pack)
        call io_field_adios2(io,engine,io_handle,'v',ng,nh,lo,hi,v,is_pack)
        call io_field_adios2(io,engine,io_handle,'w',ng,nh,lo,hi,w,is_pack)
        call io_field_adios2(io,engine,io_handle,'p',ng,nh,lo,hi,p,is_pack)
      end if
      do iscal=1,nscal
        write(scalnum,'(i3.3)') iscal
        if(is_compress) then
          call io_field_adios2(io,engine,io_handle,'s_'//scalnum,ng,nh,lo,hi,s(iscal)%val,is_pack,compress, &
                               compress_params)
        else
          call io_field_adios2(io,engine,io_handle,'s_'//scalnum,ng,nh,lo,hi,s(iscal)%val,is_pack)
        end if
      end do
      meta_rp(1) = time
      meta_i4(1) = istep
      call io_field_adios2_1d(io,engine,io_handle,'time',meta_rp,comm)
      call io_field_adios2_1d(io,engine,io_handle,'istep',meta_i4,comm)
      if(present(x_g)) call io_field_adios2_1d(io,engine,io_handle,'x',x_g(1:ng(1)),comm)
      if(present(y_g)) call io_field_adios2_1d(io,engine,io_handle,'y',y_g(1:ng(2)),comm)
      if(present(z_g)) call io_field_adios2_1d(io,engine,io_handle,'z',z_g(1:ng(3)),comm)
    end select
    !
    call adios2_close_engine(adios,engine)
  end subroutine load_all_adios2
  !
  subroutine load_one_adios2(io,filename,comm,ng,nh,lo,hi,p,varname,time,istep,x_g,y_g,z_g)
    use mod_param, only: is_use_compression
    !
    ! reads/writes a restart file for a single field using ADIOS2
    !
    implicit none
    integer, parameter :: NAME_LEN_MAX = 256
    character(len=1), intent(in) :: io
    character(len=*), intent(in) :: filename
    integer         , intent(in) :: comm
    integer , intent(in), dimension(3) :: ng,nh,lo,hi
    real(rp), intent(inout), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: p
    character(len=*), intent(in) :: varname
    real(rp), intent(inout), optional :: time
    integer , intent(inout), optional :: istep
    real(rp), intent(inout), dimension(1-nh(1):), optional :: x_g
    real(rp), intent(inout), dimension(1-nh(2):), optional :: y_g
    real(rp), intent(inout), dimension(1-nh(3):), optional :: z_g
    type(adios2_adios) :: adios
    type(adios2_io) :: io_handle
    type(adios2_engine) :: engine
    type(adios2_operator) :: compress
    type(adios2_compression_params) :: compress_params
    real(rp) :: meta_rp(1)
    integer :: meta_i4(1)
    character(len=NAME_LEN_MAX) :: field_name
    logical :: is_compress,is_pack
    !
    field_name = trim(varname)
    call adios2_open_engine(io,adios,io_handle,engine,filename,comm)
    is_compress = io == 'w' .and. is_use_compression
    !is_pack = is_use_compression
    !
    ! see load_all_adios2: pack restart fields for ADIOS2 to avoid version-dependent SetMemorySelection issues on haloed arrays.
    !
    is_pack = .true.
    if(is_compress) then
      compress_params = adios_get_default_compression_params('blosc')
      call adios2_define_compress(adios,compress,compress_params)
    end if
    !
    select case(io)
    case('r')
      call io_field_adios2(io,engine,io_handle,field_name,ng,nh,lo,hi,p,is_pack)
      if(present(time)) then
        call io_field_adios2_1d(io,engine,io_handle,'time',meta_rp,comm)
        time = meta_rp(1)
      end if
      if(present(istep)) then
        call io_field_adios2_1d(io,engine,io_handle,'istep',meta_i4,comm)
        istep = meta_i4(1)
      end if
      if(present(x_g)) call io_field_adios2_1d(io,engine,io_handle,'x',x_g(1:ng(1)),comm)
      if(present(y_g)) call io_field_adios2_1d(io,engine,io_handle,'y',y_g(1:ng(2)),comm)
      if(present(z_g)) call io_field_adios2_1d(io,engine,io_handle,'z',z_g(1:ng(3)),comm)
    case('w')
      if(is_compress) then
        call io_field_adios2(io,engine,io_handle,field_name,ng,nh,lo,hi,p,is_pack,compress,compress_params)
      else
        call io_field_adios2(io,engine,io_handle,field_name,ng,nh,lo,hi,p,is_pack)
      end if
      if(present(time)) then
        meta_rp(1) = time
        call io_field_adios2_1d(io,engine,io_handle,'time',meta_rp,comm)
      end if
      if(present(istep)) then
        meta_i4(1) = istep
        call io_field_adios2_1d(io,engine,io_handle,'istep',meta_i4,comm)
      end if
      if(present(x_g)) call io_field_adios2_1d(io,engine,io_handle,'x',x_g(1:ng(1)),comm)
      if(present(y_g)) call io_field_adios2_1d(io,engine,io_handle,'y',y_g(1:ng(2)),comm)
      if(present(z_g)) call io_field_adios2_1d(io,engine,io_handle,'z',z_g(1:ng(3)),comm)
    end select
    !
    call adios2_close_engine(adios,engine)
  end subroutine load_one_adios2
  !
  subroutine adios2_open_engine(io,adios,io_handle,engine,filename,comm)
    !
    ! opens an ADIOS2 engine for restart I/O
    !
    implicit none
    character(len=1), intent(in) :: io
    type(adios2_adios), intent(out) :: adios
    type(adios2_io)   , intent(out) :: io_handle
    type(adios2_engine), intent(out) :: engine
    character(len=*), intent(in) :: filename
    integer         , intent(in) :: comm
    integer :: adios2_mode
    !
    select case(io)
    case('r')
      adios2_mode = adios2_mode_read
    case('w')
      adios2_mode = adios2_mode_write
    end select
    !
    call adios2_init(adios, comm, ierr)
    call adios2_declare_io(io_handle, adios, 'restart', ierr)
    call adios2_set_parameter(io_handle, 'Engine', 'BP5', ierr)
    call adios2_open(engine, io_handle, filename, adios2_mode, ierr)
    call adios2_begin_step(engine, ierr)
  end subroutine adios2_open_engine
  !
  subroutine adios2_close_engine(adios,engine)
    !
    ! closes an ADIOS2 engine for restart I/O
    !
    implicit none
    type(adios2_adios), intent(inout) :: adios
    type(adios2_engine), intent(inout) :: engine
    call adios2_end_step(engine, ierr)
    call adios2_close(engine, ierr)
    call adios2_finalize(adios, ierr)
  end subroutine adios2_close_engine
  !
  function adios_get_default_compression_params(name) result(params)
    !
    ! returns one of the supported ADIOS2 compression presets
    !
    implicit none
    character(len=*), intent(in), optional :: name
    type(adios2_compression_params) :: params
    character(len=16) :: name_
    name_ = 'blosc'
    if(present(name)) name_ = trim(name)
    params%keys(:) = ''
    params%values(:) = ''
    select case(trim(name_))
    case('blosc')
      params%name = 'blosc'
      params%keys(1) = 'clevel'
      params%keys(2) = 'doshuffle'
      params%values(1) = '5'
      params%values(2) = 'BLOSC_SHUFFLE'
    case('zfp')
      params%name = 'zfp'
      params%keys(1) = 'accuracy'
      params%values(1) = '1e-6'
    case('sz3')
      params%name = 'sz3'
      params%keys(1) = 'accuracy'
      params%values(1) = '1e-6'
    case default
      params%name = 'blosc'
      params%keys(1) = 'clevel'
      params%keys(2) = 'doshuffle'
      params%values(1) = '5'
      params%values(2) = 'BLOSC_SHUFFLE'
    end select
  end function adios_get_default_compression_params
  !
  subroutine adios2_define_compress(adios,compress,params)
    !
    ! defines the checkpoint compression operator requested by the user
    !
    implicit none
    type(adios2_adios), intent(in) :: adios
    type(adios2_operator), intent(out) :: compress
    type(adios2_compression_params), intent(in), optional :: params
    type(adios2_compression_params) :: params_
    params_ = adios_get_default_compression_params()
    if(present(params)) params_ = params
    call adios2_define_operator(compress,adios,'compress',trim(params_%name),ierr)
  end subroutine adios2_define_compress
  !
  subroutine adios2_add_compress(var_handle,compress,params)
    !
    ! attaches the configured operator parameters to a variable
    !
    implicit none
    type(adios2_variable), intent(in) :: var_handle
    type(adios2_operator), intent(in) :: compress
    type(adios2_compression_params), intent(in), optional :: params
    type(adios2_compression_params) :: params_
    integer :: operation_index
    integer :: first_param,iparam
    !
    params_ = adios_get_default_compression_params()
    if(present(params)) params_ = params
    first_param = 0
    do iparam=1,size(params_%keys)
      if(len_trim(params_%keys(iparam)) > 0) then
        first_param = iparam
        exit
      end if
    end do
    if(first_param > 0) then
      call adios2_add_operation(operation_index,var_handle,compress,trim(params_%keys(first_param)), &
                                trim(params_%values(first_param)),ierr)
      do iparam=first_param+1,size(params_%keys)
        if(len_trim(params_%keys(iparam)) > 0) then
          call adios2_set_operation_parameter(var_handle,operation_index,trim(params_%keys(iparam)), &
                                              trim(params_%values(iparam)),ierr)
        end if
      end do
    else
      call adios2_add_operation(operation_index,var_handle,compress,'','',ierr)
    end if
  end subroutine adios2_add_compress
  !
  subroutine io_field_adios2(io,engine,io_handle,varname,ng,nh,lo,hi,var,is_pack,compress,params)
    !
    ! reads/writes a halo-free 3D field using ADIOS2 selections
    !
    implicit none
    character(len=1), intent(in) :: io
    type(adios2_engine), intent(in) :: engine
    type(adios2_io)    , intent(in) :: io_handle
    character(len=*), intent(in) :: varname
    integer , intent(in), dimension(3) :: ng,nh,lo,hi
    real(rp), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: var
    logical , intent(in), optional :: is_pack
    type(adios2_operator), intent(in), optional :: compress
    type(adios2_compression_params), intent(in), optional :: params
    type(adios2_variable) :: var_handle
    integer(i8), dimension(3) :: shape,start,count
    integer(i8), dimension(3) :: mem_shape,mem_start
    integer , dimension(3) :: n
    logical, parameter :: adios2_constant_dims = .true.
    integer, parameter :: ndims = 3
    logical :: is_pack_
    real(rp), allocatable, dimension(:,:,:) :: tmp
    !
    n(:) = hi(:)-lo(:)+1
    shape(:) = ng(:)
    start(:) = lo(:)-1
    count(:) = n(:)
    mem_shape(:) = n(:)+2*nh(:)
    mem_start(:) = nh(:)
    is_pack_ = .false.
    if(present(is_pack)) is_pack_ = is_pack
    if(is_pack_) allocate(tmp(n(1),n(2),n(3)))
    !
    select case(io)
    case('r')
      call adios2_inquire_variable(var_handle, io_handle, varname, ierr)
      call adios2_set_selection(var_handle, ndims, start, count, ierr)
      if(is_pack_) then
        call adios2_get(engine, var_handle, tmp, adios2_mode_sync, ierr)
        var(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = tmp(:,:,:)
      else
        call adios2_set_memory_selection(var_handle, ndims, mem_start, mem_shape, ierr)
        call adios2_get(engine, var_handle, var, adios2_mode_sync, ierr)
      end if
    case('w')
      call adios2_define_variable(var_handle, io_handle, varname, ADIOS2_REAL_RP(), &
                                  ndims, shape, start, count, adios2_constant_dims, ierr)
      if(present(compress)) call adios2_add_compress(var_handle,compress,params)
      if(is_pack_) then
        tmp(:,:,:) = var(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
        call adios2_put(engine, var_handle, tmp, adios2_mode_sync, ierr)
      else
        call adios2_set_memory_selection(var_handle, ndims, mem_start, mem_shape, ierr)
        call adios2_put(engine, var_handle, var, adios2_mode_sync, ierr)
      end if
    end select
    if(is_pack_) deallocate(tmp)
  end subroutine io_field_adios2
  !
  subroutine io_field_adios2_1d_real(io,engine,io_handle,varname,var,comm)
    !
    ! reads/writes a 1D real metadata array using rank 0 ownership
    !
    implicit none
    character(len=1), intent(in) :: io
    type(adios2_engine), intent(in) :: engine
    type(adios2_io)    , intent(in) :: io_handle
    character(len=*), intent(in) :: varname
    real(rp), intent(inout), dimension(:) :: var
    integer         , intent(in) :: comm
    type(adios2_variable) :: var_handle
    integer(i8), dimension(1) :: shape,start,count
    real(rp), dimension(1) :: dummy
    logical, parameter :: adios2_constant_dims = .true.
    integer, parameter :: ndims = 1
    integer :: myid_loc
    !
    dummy = 0.
    shape(1) = size(var)
    start(1) = 0
    count(1) = 0
    call MPI_COMM_RANK(comm,myid_loc,ierr)
    if(myid_loc == 0) count(1) = shape(1)
    !
    select case(io)
    case('r')
      count(1) = shape(1)
      if(myid_loc == 0) then
        call adios2_inquire_variable(var_handle, io_handle, varname, ierr)
        call adios2_set_selection(var_handle, ndims, start, count, ierr)
        call adios2_get(engine, var_handle, var, adios2_mode_sync, ierr)
      end if
      call MPI_BCAST(var,size(var),MPI_REAL_RP,0,comm,ierr)
    case('w')
      call adios2_define_variable(var_handle, io_handle, varname, ADIOS2_REAL_RP(), &
                                  ndims, shape, start, count, adios2_constant_dims, ierr)
      if(myid_loc == 0) then
        call adios2_put(engine, var_handle, var, adios2_mode_sync, ierr)
      else
        call adios2_put(engine, var_handle, dummy, adios2_mode_sync, ierr)
      end if
    end select
  end subroutine io_field_adios2_1d_real
  !
  subroutine io_field_adios2_1d_int(io,engine,io_handle,varname,var,comm)
    !
    ! reads/writes a 1D integer metadata array using rank 0 ownership
    !
    implicit none
    character(len=1), intent(in) :: io
    type(adios2_engine), intent(in) :: engine
    type(adios2_io)    , intent(in) :: io_handle
    character(len=*), intent(in) :: varname
    integer , intent(inout), dimension(:) :: var
    integer         , intent(in) :: comm
    type(adios2_variable) :: var_handle
    integer(i8), dimension(1) :: shape,start,count
    integer, dimension(1) :: dummy
    logical, parameter :: adios2_constant_dims = .true.
    integer, parameter :: ndims = 1
    integer :: myid_loc
    !
    dummy = 0
    shape(1) = size(var)
    start(1) = 0
    count(1) = 0
    call MPI_COMM_RANK(comm,myid_loc,ierr)
    if(myid_loc == 0) count(1) = shape(1)
    !
    select case(io)
    case('r')
      count(1) = shape(1)
      if(myid_loc == 0) then
        call adios2_inquire_variable(var_handle, io_handle, varname, ierr)
        call adios2_set_selection(var_handle, ndims, start, count, ierr)
        call adios2_get(engine, var_handle, var, adios2_mode_sync, ierr)
      end if
      call MPI_BCAST(var,size(var),MPI_INTEGER,0,comm,ierr)
    case('w')
      call adios2_define_variable(var_handle, io_handle, varname, adios2_type_integer4, &
                                  ndims, shape, start, count, adios2_constant_dims, ierr)
      if(myid_loc == 0) then
        call adios2_put(engine, var_handle, var, adios2_mode_sync, ierr)
      else
        call adios2_put(engine, var_handle, dummy, adios2_mode_sync, ierr)
      end if
    end select
  end subroutine io_field_adios2_1d_int
  !
  integer function ADIOS2_REAL_RP()
    !
    ! returns the ADIOS2 type matching rp
    !
    implicit none
    if(rp == dp) then
      ADIOS2_REAL_RP = adios2_type_dp
    else
      ADIOS2_REAL_RP = adios2_type_real
    end if
  end function ADIOS2_REAL_RP
#endif
end module mod_load
