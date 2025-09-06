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
  use mod_common_mpi, only: myid,ierr
  use mod_types
  use mod_utils, only: f_sizeof
  use mod_scal, only: scalar
  implicit none
  private
  public load_one,io_field
  contains
  subroutine load_all(io,filename,comm,ng,nh,lo,hi,nscal,u,v,w,p,s,time,istep)
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
        if(myid == 0) print*, '*** Simulation aborted due a checkpoint file with incorrect size ***'
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
  end subroutine load_all
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
    ! n.b.: the Poisson solver buffers are being recycled here, meaning
    ! that I/O should use same precision as these buffers, in case this
    ! routine is used in the future;
    ! alternatively one could *temporarily* (i.e., during I/O) offload
    ! device memory and allocate larger buffers (and while at it, use a
    ! dedicated cuDecomp grid descriptor without axis-contiguous layout)
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
  subroutine load_one(io,filename,comm,ng,nh,lo,hi,p,time,istep)
    !
    ! reads/writes a restart file for a single field
    !
    implicit none
    character(len=1), intent(in) :: io
    character(len=*), intent(in) :: filename
    integer         , intent(in) :: comm
    integer , intent(in), dimension(3) :: ng,nh,lo,hi
    real(rp), intent(inout), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: p
    real(rp), intent(inout), optional :: time
    integer , intent(inout), optional :: istep
    real(rp), dimension(2) :: fldinfo
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
      good = (product(int(ng(:),MPI_OFFSET_KIND))*1+2)*f_sizeof(1._rp)
      if(filesize /= good) then
        if(myid == 0) print*, ''
        if(myid == 0) print*, '*** Simulation aborted due a checkpoint file with incorrect size ***'
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
      if(present(time) .and. present(istep)) then
        call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,MPI_REAL_RP,'native',MPI_INFO_NULL,ierr)
        nreals_myid = 0
        if(myid == 0) nreals_myid = 2
        call MPI_FILE_READ(fh,fldinfo,nreals_myid,MPI_REAL_RP,MPI_STATUS_IGNORE,ierr)
        call MPI_FILE_CLOSE(fh,ierr)
        call MPI_BCAST(fldinfo,2,MPI_REAL_RP,0,comm,ierr)
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
        call MPI_FILE_CLOSE(fh,ierr)
      end if
    end select
  end subroutine load_one
  !
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
                  p(1:n(1),1:n(2),1:n(3)), &
                  fldinfo(1:2)
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
                   p(1:n(1),1:n(2),1:n(3)), &
                   fldinfo(1:2)
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
    real(rp), intent(inout), optional :: time
    integer , intent(inout), optional :: istep
    real(rp), dimension(2) :: fldinfo
    integer :: iunit
    !
    select case(io)
    case('r')
      open(newunit=iunit,file=filename,action='read' ,access='stream',form='unformatted',status='old'    )
      read(iunit) p(1:n(1),1:n(2),1:n(3))
      if(present(time) .and. present(istep)) read(iunit)  fldinfo(1:2)
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
      if(present(time) .and. present(istep)) write(iunit) fldinfo(1:2)
      close(iunit)
    end select
  end subroutine load_one_local
  !
#if defined(_USE_HDF5)
  subroutine io_field_hdf5(io,filename,varname,ng,nh,lo,hi,var,meta,x_g,y_g,z_g)
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
    integer         , intent(in), dimension(3)   :: ng,nh,lo,hi
    real(rp), intent(inout), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: var
    real(rp), intent(inout), dimension(2), optional :: meta
    real(rp), intent(inout), dimension(0:), optional :: x_g,y_g,z_g
    integer , dimension(3) :: n
    integer , dimension(3) :: sizes,subsizes,starts
    !
    ! HDF5 variables
    !
    integer :: ndims, ierr
    integer(HID_T) :: file_id,group_id
    integer(HID_T) :: filespace
    integer(HID_T) :: slabspace
    integer(HID_T) :: memspace
    !
    integer(HID_T) :: dset
    !
    integer(HSIZE_T) :: dims(3)
    !
    integer(HID_T) :: plist_id
    integer(HSIZE_T) , dimension(3) :: data_count
    integer(HSSIZE_T), dimension(3) :: data_offset
    integer(HSSIZE_T), dimension(3) :: halo_offset
    !
    n(:)        = hi(:)-lo(:)+1
    sizes(:)    = ng(:)
    subsizes(:) = n(:)
    starts(:)   = lo(:) - 1 ! starts from 0
    !
    ndims = 3
    dims(:) = ng(:)
    data_count(:) = subsizes(:)
    data_offset(:) = starts(:)
    halo_offset(:) = nh(:)
    !
    select case(io)
    case('r')
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,ierr)
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,ierr)
      call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,ierr,access_prp=plist_id)
      call h5pclose_f(plist_id,ierr)
      !
      call h5dopen_f(file_id,'fields/'//varname,dset,ierr)
      call h5screate_simple_f(ndims,data_count+2*nh(:),memspace,ierr)
      !
      call h5dget_space_f(dset,slabspace,ierr)
      call h5sselect_hyperslab_f(slabspace,H5S_SELECT_SET_F,data_offset,data_count,ierr)
      call h5sselect_hyperslab_f(memspace,H5S_SELECT_SET_F,halo_offset,data_count,ierr)
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,ierr)
      call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,ierr)
      !
      call h5dread_f(dset,H5T_NATIVE_DOUBLE,var,dims,ierr,file_space_id=slabspace,mem_space_id=memspace,xfer_prp=plist_id)
      !
      call h5pclose_f(plist_id,ierr)
      call h5dclose_f(dset,ierr)
      call h5sclose_f(memspace,ierr)
      call h5fclose_f(file_id,ierr)
      !
      if(myid == 0) then
        call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,ierr)
        call h5dopen_f(file_id,'meta/time',dset,ierr)
        call h5dread_f(dset,H5T_NATIVE_DOUBLE,meta,[int(2,HSIZE_T)],ierr)
        call h5dclose_f(dset,ierr)
        call h5fclose_f(file_id,ierr)
      end if
      call MPI_Bcast(meta,2,MPI_REAL_RP,0,MPI_COMM_WORLD,ierr)
    case('w')
      call h5screate_simple_f(ndims,dims,filespace,ierr)
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,ierr)
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,ierr)
      call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,ierr,access_prp=plist_id)
      call h5pclose_f(plist_id,ierr)
      !
      call h5gcreate_f(file_id,'fields',group_id,ierr)
      call h5dcreate_f(group_id,varname,H5T_NATIVE_DOUBLE,filespace,dset,ierr)
      call h5screate_simple_f(ndims,data_count+2*nh(:),memspace,ierr)
      call h5dget_space_f(dset,slabspace,ierr)
      call h5sselect_hyperslab_f(slabspace,H5S_SELECT_SET_F,data_offset,data_count,ierr)
      call h5sselect_hyperslab_f(memspace,H5S_SELECT_SET_F,halo_offset,data_count,ierr)
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,ierr)
      call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,ierr)
      call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,var,dims,ierr,file_space_id=slabspace,mem_space_id=memspace,xfer_prp=plist_id)
      !
      call h5pclose_f(plist_id,ierr)
      call h5dclose_f(dset,ierr)
      call h5sclose_f(memspace,ierr)
      call h5sclose_f(slabspace,ierr)
      call h5sclose_f(filespace,ierr)
      call h5gclose_f(group_id,ierr)
      call h5fclose_f(file_id,ierr)
      !
      ! write metadata
      !
      if(myid == 0) then
        call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,ierr)
        !
        if(present(x_g) .and. present(y_g) .and. present(z_g)) then
          call h5gcreate_f(file_id,'grid',group_id,ierr)
          call h5screate_simple_f(1,[int(ng(1),hsize_t)],filespace,ierr)
          call h5dcreate_f(group_id,'x',h5t_native_double,filespace,dset,ierr)
          call h5dwrite_f(dset,h5t_native_double,x_g(1:ng(1)),[int(ng(1),hsize_t)],ierr)
          call h5screate_simple_f(1,[int(ng(2),hsize_t)],filespace,ierr)
          call h5dcreate_f(group_id,'y',h5t_native_double,filespace,dset,ierr)
          call h5dwrite_f(dset,h5t_native_double,y_g(1:ng(2)),[int(ng(2),hsize_t)],ierr)
          call h5screate_simple_f(1,[int(ng(3),hsize_t)],filespace,ierr)
          call h5dcreate_f(group_id,'z',h5t_native_double,filespace,dset,ierr)
          call h5dwrite_f(dset,h5t_native_double,z_g(1:ng(3)),[int(ng(3),hsize_t)],ierr)
          call h5dclose_f(dset,ierr)
          call h5gclose_f(group_id,ierr)
          call h5sclose_f(filespace,ierr)
        end if
        !
        if(present(meta)) then
          call h5gcreate_f(file_id,'meta',group_id,ierr)
          call h5screate_simple_f(1,[int(2,hsize_t)],filespace,ierr)
          call h5dcreate_f(group_id,'time',h5t_native_double,filespace,dset,ierr)
          call h5dwrite_f(dset,h5t_native_double,meta,[int(2,hsize_t)],ierr)
          call h5dclose_f(dset,ierr)
          call h5gclose_f(group_id,ierr)
          call h5sclose_f(filespace,ierr)
        end if
        call h5fclose_f(file_id,ierr)
      end if
    end select
  end subroutine io_field_hdf5
#endif
end module mod_load
