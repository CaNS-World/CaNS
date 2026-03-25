! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_output
  use mpi
  use decomp_2d_io
  use mod_common_mpi, only:ierr,myid
  use mod_load      , only: io_write_subset
  use mod_types
  implicit none
  private
  public out0d,gen_alias,out1d,out1d_chan,out2d,out3d,write_log_output,write_visu_2d,write_visu_3d,write_visu_subset
  character(len=*), parameter :: fmt_dp = '(*(es24.16e3,1x))', &
                                 fmt_sp = '(*(es15.8e2,1x))'
#if !defined(_SINGLE_PRECISION)
  character(len=*), parameter :: fmt_rp = fmt_dp
#else
  character(len=*), parameter :: fmt_rp = fmt_sp
#endif
  contains
  subroutine out0d(fname,n,var)
    !
    ! appends the first n entries of an array
    ! var to a file
    ! fname -> name of the file
    ! n     -> number of entries
    ! var   -> input array of real values
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in) :: n
    real(rp), intent(in), dimension(:) :: var
    integer :: iunit
    !
    if(myid == 0) then
      open(newunit=iunit,file=fname,position='append')
      write(iunit,fmt_rp) var(1:n)
      close(iunit)
    end if
  end subroutine out0d
  !
  subroutine gen_alias(myid,datadir,fname,fname_alias)
    !
    ! this subroutine generates a symlink with name `fname_alias`, pointing to
    ! file `datadir//fname` using the `execute_command_line` intrinsic;
    ! it is called by task `myid`
    !
    integer, intent(in) :: myid
    character(len=*), intent(in) :: datadir,fname,fname_alias
    if(myid == 0) call execute_command_line('ln -sf '//trim(fname)//' '//trim(datadir)//fname_alias)
  end subroutine gen_alias
  !
  subroutine out1d(fname,ng,lo,hi,idir,l,dl,z_g,dz,p)
    !
    ! writes the profile of a variable averaged
    ! over two domain directions
    !
    ! fname -> name of the file
    ! ng    -> global domain sizes
    ! lo,hi -> upper and lower extents of the input array
    ! idir  -> direction of the profile
    ! dl,l  -> uniform grid spacing and length arrays
    ! z_g   -> global z coordinate array (grid is non-uniform in z)
    ! dz    -> local z grid spacing array (should work also with the global one)
    ! p     -> 3D input scalar field
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: ng,lo,hi
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:       ) :: z_g
    real(rp), intent(in), dimension(lo(3)-1:) :: dz
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    real(rp), allocatable, dimension(:) :: p1d
    integer :: i,j,k
    integer :: iunit
    real(rp) :: grid_area_ratio,p1d_s
    !
    allocate(p1d(ng(idir)))
    !$acc enter data create(p1d)
    !$acc parallel loop default(present)
    do k=1,size(p1d)
      p1d(k) = 0._rp
    end do
    select case(idir)
    case(3)
      grid_area_ratio = dl(1)*dl(2)/(l(1)*l(2))
      !$acc parallel loop gang default(present) private(p1d_s)
      do k=lo(3),hi(3)
        p1d_s = 0._rp
        !$acc loop collapse(2) reduction(+:p1d_s)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            p1d_s = p1d_s + p(i,j,k)*grid_area_ratio
          end do
        end do
        p1d(k) = p1d_s
      end do
      !$acc exit data copyout(p1d)
      call MPI_ALLREDUCE(MPI_IN_PLACE,p1d(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myid == 0) then
        open(newunit=iunit,file=fname)
        do k=1,ng(3)
          write(iunit,fmt_rp) z_g(k),p1d(k)
        end do
        close(iunit)
      end if
    case(2)
      grid_area_ratio = dl(1)/(l(1)*l(3))
      !$acc parallel loop gang default(present) private(p1d_s)
      do j=lo(2),hi(2)
        p1d_s = 0._rp
        !$acc loop collapse(2) reduction(+:p1d_s)
        do k=lo(3),hi(3)
          do i=lo(1),hi(1)
            p1d_s = p1d_s + p(i,j,k)*dz(k)*grid_area_ratio
          end do
        end do
        p1d(j) = p1d_s
      end do
      !$acc exit data copyout(p1d)
      call MPI_ALLREDUCE(MPI_IN_PLACE,p1d(1),ng(2),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myid == 0) then
        open(newunit=iunit,file=fname)
        do j=1,ng(2)
          write(iunit,fmt_rp) (j-.5)*dl(2),p1d(j)
        end do
        close(iunit)
      end if
    case(1)
      grid_area_ratio = dl(2)/(l(2)*l(3))
      !$acc parallel loop gang default(present) private(p1d_s)
      do i=lo(1),hi(1)
        p1d_s = 0._rp
        !$acc loop collapse(2) reduction(+:p1d_s)
        do k=lo(3),hi(3)
          do j=lo(2),hi(2)
            p1d_s = p1d_s + p(i,j,k)*dz(k)*grid_area_ratio
          end do
        end do
        p1d(i) = p1d_s
      end do
      !$acc exit data copyout(p1d)
      call MPI_ALLREDUCE(MPI_IN_PLACE,p1d(1),ng(1),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myid == 0) then
        open(newunit=iunit,file=fname)
        do i=1,ng(1)
          write(iunit,fmt_rp) (i-.5)*dl(1),p1d(i)
        end do
        close(iunit)
      end if
    end select
  end subroutine out1d
  !
  subroutine out2d(fname,ng,nh,lo,hi,inorm,islice,p,varname,time,istep,x_g,y_g,z_g,is_pack)
    !
    ! saves a planar slice of a scalar field using the structured subset engine
    !
    ! fname  -> name of the output file
    ! ng     -> global domain sizes
    ! nh     -> local halo extents
    ! lo,hi  -> local interior ownership in global indices
    ! inorm  -> plane is perpendicular to direction inorm (1,2,3)
    ! islice -> plane is of constant index islice in direction inorm
    ! p      -> 3D input scalar field
    ! varname,time,istep,x_g,y_g,z_g -> subset metadata
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: ng,nh,lo,hi
    integer , intent(in) :: inorm,islice
    real(rp), intent(in), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: p
    character(len=*), intent(in) :: varname
    real(rp), intent(in) :: time
    integer , intent(in) :: istep
    real(rp), intent(in), dimension(1-nh(1):) :: x_g
    real(rp), intent(in), dimension(1-nh(2):) :: y_g
    real(rp), intent(in), dimension(1-nh(3):) :: z_g
    logical , intent(in), optional :: is_pack
    integer, dimension(3) :: lo_out,hi_out
    !
    select case(inorm)
    case(1)
      lo_out(:) = [islice,1     ,1   ]
      hi_out(:) = [islice,ng(2),ng(3)]
    case(2)
      lo_out(:) = [1    ,islice,1   ]
      hi_out(:) = [ng(1),islice,ng(3)]
    case(3)
      lo_out(:) = [1    ,1    ,islice]
      hi_out(:) = [ng(1),ng(2),islice]
    end select
    call io_write_subset(fname,varname,MPI_COMM_WORLD,ng,nh,lo,hi,lo_out,hi_out,[1,1,1],p,time,istep,x_g,y_g,z_g,is_pack)
  end subroutine out2d
  !
  subroutine out3d(fname,ng,nh,lo,hi,lo_out,hi_out,nskip,p,varname,time,istep,x_g,y_g,z_g,is_pack)
    !
    ! saves a structured subset of a scalar field using the subset engine
    !
    ! fname  -> name of the output file
    ! ng     -> global domain sizes
    ! nh     -> local halo extents
    ! lo,hi  -> local interior ownership in global indices
    ! lo_out -> first global point written in each direction
    ! hi_out -> last  global point written in each direction
    ! nskip  -> output stride in each direction
    ! p      -> 3D input scalar field
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: ng,nh,lo,hi,lo_out,hi_out
    integer , intent(in), dimension(3) :: nskip
    real(rp), intent(in), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: p
    character(len=*), intent(in) :: varname
    real(rp), intent(in) :: time
    integer , intent(in) :: istep
    real(rp), intent(in), dimension(1-nh(1):) :: x_g
    real(rp), intent(in), dimension(1-nh(2):) :: y_g
    real(rp), intent(in), dimension(1-nh(3):) :: z_g
    logical , intent(in), optional :: is_pack
    !
    call io_write_subset(fname,varname,MPI_COMM_WORLD,ng,nh,lo,hi,lo_out,hi_out,nskip,p,time,istep,x_g,y_g,z_g,is_pack)
  end subroutine out3d
  !
  subroutine write_log_output(fname,fname_fld,varname,lo_out,hi_out,nskip,time,istep)
    !
    ! appends information about a saved binary file to a log file
    ! this file is used to generate a xdmf file for visualization of field data
    !
    ! fname     -> name of the output log file
    ! fname_fld -> name of the saved binary file (excluding the directory)
    ! varname   -> name of the variable that is saved
    ! lo_out    -> first element of the field that is saved in each direction, e.g. [1,1,1]
    ! hi_out    -> last  element of the field that is saved in each direction, e.g. [ng(1),ng(2),ng(3)]
    ! nskip     -> step size between nmin and nmax, e.g. [1,1,1] if the whole array is saved
    ! time      -> physical time
    ! istep     -> time step number
    !
    implicit none
    character(len=*), intent(in) :: fname,fname_fld,varname
    integer , intent(in), dimension(3) :: lo_out,hi_out,nskip
    real(rp), intent(in)               :: time
    integer , intent(in)               :: istep
    character(len=100) :: cfmt
    integer :: iunit
    !
    write(cfmt, '(A)') '(A,A,A,9i5,E16.7e3,i7)'
    if(myid  ==  0) then
      open(newunit=iunit,file=fname,position='append')
      write(iunit,trim(cfmt)) trim(fname_fld),' ',trim(varname),lo_out,hi_out,nskip,time,istep
      close(iunit)
    end if
  end subroutine write_log_output
  !
  subroutine write_visu_subset(datadir,fname_bin,fname_log,varname,ng,nh,lo,hi,lo_out,hi_out,nskip,time,istep,p,x_g,y_g,z_g,is_pack)
    !
    ! wraps subset data output and subset metadata logging
    !
    implicit none
    character(len=*), intent(in)          :: datadir,fname_bin,fname_log,varname
    integer , intent(in), dimension(3)    :: ng,nh,lo,hi,lo_out,hi_out,nskip
    real(rp), intent(in)                  :: time
    integer , intent(in)                  :: istep
    real(rp), intent(in), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: p
    real(rp), intent(in), dimension(1-nh(1):)   :: x_g
    real(rp), intent(in), dimension(1-nh(2):)   :: y_g
    real(rp), intent(in), dimension(1-nh(3):)   :: z_g
    logical , intent(in), optional        :: is_pack
    !
    call out3d(trim(datadir)//trim(fname_bin),ng,nh,lo,hi,lo_out,hi_out,nskip,p,varname,time,istep,x_g,y_g,z_g,is_pack)
    call write_log_output(trim(datadir)//trim(fname_log),trim(fname_bin),trim(varname),lo_out,hi_out,nskip,time,istep)
  end subroutine write_visu_subset
  !
  subroutine write_visu_3d(datadir,fname_bin,fname_log,varname,ng,nh,lo,hi,lo_out,hi_out,nskip,time,istep,p,x_g,y_g,z_g,is_pack)
    !
    ! wraps the calls of out3d and write_log_output into the same subroutine
    !
    implicit none
    character(len=*), intent(in)          :: datadir,fname_bin,fname_log,varname
    integer , intent(in), dimension(3)    :: ng,nh,lo,hi,lo_out,hi_out,nskip
    real(rp), intent(in)                  :: time
    integer , intent(in)                  :: istep
    real(rp), intent(in), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: p
    real(rp), intent(in), dimension(1-nh(1):)   :: x_g
    real(rp), intent(in), dimension(1-nh(2):)   :: y_g
    real(rp), intent(in), dimension(1-nh(3):)   :: z_g
    logical , intent(in), optional        :: is_pack
    !
    call write_visu_subset(datadir,fname_bin,fname_log,varname,ng,nh,lo,hi,lo_out,hi_out,nskip,time,istep,p,x_g,y_g,z_g,is_pack)
  end subroutine write_visu_3d
  !
  subroutine write_visu_2d(datadir,fname_bin,fname_log,varname,ng,nh,lo,hi,inorm,nslice,time,istep,p,x_g,y_g,z_g,is_pack)
    !
    ! wraps the calls of out2d and write-log_output into the same subroutine
    !
    implicit none
    character(len=*), intent(in)          :: datadir,fname_bin,fname_log,varname
    integer , intent(in), dimension(3)    :: ng,nh,lo,hi
    integer , intent(in)                  :: inorm,nslice
    real(rp), intent(in)                  :: time
    integer , intent(in)                  :: istep
    real(rp), intent(in), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: p
    real(rp), intent(in), dimension(1-nh(1):)   :: x_g
    real(rp), intent(in), dimension(1-nh(2):)   :: y_g
    real(rp), intent(in), dimension(1-nh(3):)   :: z_g
    logical , intent(in), optional        :: is_pack
    integer , dimension(3) :: lo_out,hi_out
    !
    select case(inorm)
    case(1)
      lo_out(:) = [nslice,1     ,1   ]
      hi_out(:) = [nslice,ng(2),ng(3)]
    case(2)
      lo_out(:) = [1    ,nslice,1    ]
      hi_out(:) = [ng(1),nslice,ng(3)]
    case(3)
      lo_out(:) = [1    ,1    ,nslice]
      hi_out(:) = [ng(1),ng(2),nslice]
    end select
    call write_visu_subset(datadir,fname_bin,fname_log,varname,ng,nh,lo,hi,lo_out,hi_out,[1,1,1],time,istep,p,x_g,y_g,z_g,is_pack)
  end subroutine write_visu_2d
  !
  subroutine out1d_chan(fname,ng,lo,hi,idir,l,dl,z_g,u,v,w) ! e.g. for a channel with streamwise dir in x
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: ng,lo,hi
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: z_g
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: u,v,w
    real(rp), allocatable, dimension(:) :: um,vm,wm,u2,v2,w2,uw
    integer :: i,j,k
    integer :: iunit
    integer :: q
    real(rp) :: grid_area_ratio
    !
    q = ng(idir)
    select case(idir)
    case(3)
      grid_area_ratio = dl(1)*dl(2)/(l(1)*l(2))
      allocate(um(0:q+1),vm(0:q+1),wm(0:q+1),u2(0:q+1),v2(0:q+1),w2(0:q+1),uw(0:q+1))
      um(:) = 0.
      vm(:) = 0.
      wm(:) = 0.
      u2(:) = 0.
      v2(:) = 0.
      w2(:) = 0.
      uw(:) = 0.
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            um(k) = um(k) + u(i,j,k)
            vm(k) = vm(k) + v(i,j,k)
            wm(k) = wm(k) + 0.50*(w(i,j,k-1) + w(i,j,k))
            u2(k) = u2(k) + u(i,j,k)**2
            v2(k) = v2(k) + v(i,j,k)**2
            w2(k) = w2(k) + 0.50*(w(i,j,k)**2+w(i,j,k-1)**2)
            uw(k) = uw(k) + 0.25*(u(i-1,j,k) + u(i,j,k))* &
                                 (w(i,j,k-1) + w(i,j,k))
          end do
        end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,um(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vm(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,wm(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,u2(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,v2(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,w2(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,uw(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      um(:) = um(:)*grid_area_ratio
      vm(:) = vm(:)*grid_area_ratio
      wm(:) = wm(:)*grid_area_ratio
      u2(:) = u2(:)*grid_area_ratio - um(:)**2
      v2(:) = v2(:)*grid_area_ratio - vm(:)**2
      w2(:) = w2(:)*grid_area_ratio - wm(:)**2
      uw(:) = uw(:)*grid_area_ratio - um(:)*wm(:)
      if(myid == 0) then
        open(newunit=iunit,file=fname)
        do k=1,ng(3)
          write(iunit,fmt_rp) z_g(k),um(k),vm(k),wm(k), &
                                     u2(k),v2(k),w2(k), &
                                     uw(k)
        end do
        close(iunit)
      end if
    case(2)
    case(1)
    end select
  end subroutine out1d_chan
  !
  subroutine out2d_duct(fname,ng,lo,hi,idir,l,dl,z_g,u,v,w) ! e.g. for a duct
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: ng,lo,hi
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: z_g
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: u,v,w
    real(rp), allocatable, dimension(:,:) :: um,vm,wm,u2,v2,w2,uv,uw,vw
    integer :: i,j,k
    integer :: iunit
    integer :: p,q
    real(rp) :: x_g,y_g,grid_area_ratio
    !
    select case(idir) ! streamwise direction
    case(3)
    case(2)
      grid_area_ratio = dl(2)/l(2)
      p = ng(1)
      q = ng(3)
      allocate(um(p,q),vm(p,q),wm(p,q),u2(p,q),v2(p,q),w2(p,q),uv(p,q),vw(p,q))
      !
      um(:,:) = 0.
      vm(:,:) = 0.
      wm(:,:) = 0.
      u2(:,:) = 0.
      v2(:,:) = 0.
      w2(:,:) = 0.
      uv(:,:) = 0.
      vw(:,:) = 0.
      do k=lo(3),hi(3)
        do i=lo(1),hi(1)
          do j=lo(2),hi(2)
            um(i,k) = um(i,k) + 0.5*(u(i-1,j,k)+u(i,j,k))
            vm(i,k) = vm(i,k) + v(i,j,k)
            wm(i,k) = wm(i,k) + 0.5*(w(i,j,k-1)+w(i,j,k))
            u2(i,k) = u2(i,k) + 0.5*(u(i-1,j,k)**2+u(i,j,k)**2)
            v2(i,k) = v2(i,k) + v(i,j,k)**2
            w2(i,k) = w2(i,k) + 0.5*(w(i,j,k-1)**2+w(i,j,k)**2)
            vw(i,k) = vw(i,k) + 0.25*(v(i,j-1,k) + v(i,j,k))* &
                                     (w(i,j,k-1) + w(i,j,k))
            uv(i,k) = uv(i,k) + 0.25*(u(i-1,j,k) + u(i,j,k))* &
                                     (v(i,j-1,k) + v(i,j,k))
          end do
        end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,um(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vm(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,wm(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,u2(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,v2(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,w2(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vw(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,uv(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      um(:,:) = um(:,:)*grid_area_ratio
      vm(:,:) = vm(:,:)*grid_area_ratio
      wm(:,:) = wm(:,:)*grid_area_ratio
      u2(:,:) = u2(:,:)*grid_area_ratio - um(:,:)**2
      v2(:,:) = v2(:,:)*grid_area_ratio - vm(:,:)**2
      w2(:,:) = w2(:,:)*grid_area_ratio - wm(:,:)**2
      vw(:,:) = vw(:,:)*grid_area_ratio - vm(:,:)*wm(:,:)
      uv(:,:) = uv(:,:)*grid_area_ratio - um(:,:)*vm(:,:)
      if(myid == 0) then
        open(newunit=iunit,file=fname)
        do k=1,ng(3)
          do i=1,ng(1)
            x_g = (i-.5)*dl(1)
            write(iunit,fmt_rp) x_g,z_g(k),um(i,k),vm(i,k),wm(i,k), &
                                           u2(i,k),v2(i,k),w2(i,k), &
                                           vw(i,k),uv(i,k)
          end do
        end do
        close(iunit)
      end if
    case(1)
      grid_area_ratio = dl(1)/l(1)
      p = ng(2)
      q = ng(3)
      allocate(um(p,q),vm(p,q),wm(p,q),u2(p,q),v2(p,q),w2(p,q),uv(p,q),uw(p,q))
      !
      um(:,:) = 0.
      vm(:,:) = 0.
      wm(:,:) = 0.
      u2(:,:) = 0.
      v2(:,:) = 0.
      w2(:,:) = 0.
      uv(:,:) = 0.
      uw(:,:) = 0.
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            um(j,k) = um(j,k) + u(i,j,k)
            vm(j,k) = vm(j,k) + 0.5*(v(i,j-1,k)+v(i,j,k))
            wm(j,k) = wm(j,k) + 0.5*(w(i,j,k-1)+w(i,j,k))
            u2(j,k) = u2(j,k) + u(i,j,k)**2
            v2(j,k) = v2(j,k) + 0.5*(v(i,j-1,k)**2+v(i,j,k)**2)
            w2(j,k) = w2(j,k) + 0.5*(w(i,j,k-1)**2+w(i,j,k)**2)
            uv(j,k) = uv(j,k) + 0.25*(u(i-1,j,k) + u(i,j,k))* &
                                     (v(i,j-1,k) + v(i,j,k))
            uw(j,k) = uw(j,k) + 0.25*(u(i-1,j,k) + u(i,j,k))* &
                                     (w(i,j,k-1) + w(i,j,k))
          end do
        end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,um(1,1),ng(2)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vm(1,1),ng(2)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,wm(1,1),ng(2)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,u2(1,1),ng(2)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,v2(1,1),ng(2)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,w2(1,1),ng(2)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,uv(1,1),ng(2)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,uw(1,1),ng(2)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      um(:,:) = um(:,:)*grid_area_ratio
      vm(:,:) = vm(:,:)*grid_area_ratio
      wm(:,:) = wm(:,:)*grid_area_ratio
      u2(:,:) = u2(:,:)*grid_area_ratio - um(:,:)**2
      v2(:,:) = v2(:,:)*grid_area_ratio - vm(:,:)**2
      w2(:,:) = w2(:,:)*grid_area_ratio - wm(:,:)**2
      uv(:,:) = uv(:,:)*grid_area_ratio - um(:,:)*vm(:,:)
      uw(:,:) = uw(:,:)*grid_area_ratio - um(:,:)*wm(:,:)
      if(myid == 0) then
        open(newunit=iunit,file=fname)
        do k=1,ng(3)
          do j=1,ng(2)
            y_g = (j-.5)*dl(2)
            write(iunit,fmt_rp) y_g,z_g(k),um(j,k),vm(j,k),wm(j,k), &
                                           u2(j,k),v2(j,k),w2(j,k), &
                                           uv(j,k),uw(j,k)
          end do
        end do
        close(iunit)
      end if
    end select
  end subroutine out2d_duct
end module mod_output
