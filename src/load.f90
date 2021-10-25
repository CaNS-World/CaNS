module mod_load
  use mpi_f08
  use mod_common_mpi, only:myid_block
  use mod_types
  implicit none
  private
  public load,io_field
  contains
  subroutine load(io,filename,comm,ng,nh,lo,hi,u,v,w,p,time,istep)
    !
    ! reads/writes a restart file
    !
    implicit none
    character(len=1), intent(in) :: io
    character(len=*), intent(in) :: filename
    type(MPI_COMM) , intent(in)         :: comm
    integer , intent(in), dimension(3) :: ng,nh,lo,hi
    real(rp), intent(inout), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: u,v,w,p
    real(rp), intent(inout) :: time
    integer , intent(inout) :: istep
    real(rp), dimension(2) :: fldinfo
    type(MPI_FILE) :: fh
    integer :: nreals_myid
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp,good
    !
    select case(io)
    case('r')
      call MPI_FILE_OPEN(comm, filename, &
           MPI_MODE_RDONLY, MPI_INFO_NULL,fh)
      !
      ! check file size first
      !
      call MPI_FILE_GET_SIZE(fh,filesize)
      good = (product(ng)*4+2)*int(storage_size(1._rp)/8,MPI_OFFSET_KIND)
      if(filesize /= good) then
        if(myid_block == 0) write(stderr,*) ''
        if(myid_block == 0) write(stderr,*) '*** Simulation aborted due a checkpoint file with incorrect size ***'
        if(myid_block == 0) write(stderr,*) '    file: ', filename, ' | expected size: ', good, '| actual size: ', filesize
        call MPI_FINALIZE()
        error stop
      endif
      !
      ! read
      !
      disp = 0_MPI_OFFSET_KIND
      call io_field('r',fh,ng,nh,lo,hi,disp,u)
      call io_field('r',fh,ng,nh,lo,hi,disp,v)
      call io_field('r',fh,ng,nh,lo,hi,disp,w)
      call io_field('r',fh,ng,nh,lo,hi,disp,p)
      call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,MPI_REAL_RP,'native',MPI_INFO_NULL)
      nreals_myid = 0
      if(myid_block == 0) nreals_myid = 2
      call MPI_FILE_READ(fh,fldinfo,nreals_myid,MPI_REAL_RP,MPI_STATUS_IGNORE)
      call MPI_FILE_CLOSE(fh)
      call MPI_BCAST(fldinfo,2,MPI_REAL_RP,0,comm)
      time  =      fldinfo(1)
      istep = nint(fldinfo(2))
    case('w')
      !
      ! write
      !
      call MPI_FILE_OPEN(comm, filename                 , &
           MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh)
      filesize = 0_MPI_OFFSET_KIND
      call MPI_FILE_SET_SIZE(fh,filesize)
      disp = 0_MPI_OFFSET_KIND
      call io_field('w',fh,ng,nh,lo,hi,disp,u)
      call io_field('w',fh,ng,nh,lo,hi,disp,v)
      call io_field('w',fh,ng,nh,lo,hi,disp,w)
      call io_field('w',fh,ng,nh,lo,hi,disp,p)
      call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,MPI_REAL_RP,'native',MPI_INFO_NULL)
      fldinfo = [time,1._rp*istep]
      nreals_myid = 0
      if(myid_block == 0) nreals_myid = 2
      call MPI_FILE_WRITE(fh,fldinfo,nreals_myid,MPI_REAL_RP,MPI_STATUS_IGNORE)
      call MPI_FILE_CLOSE(fh)
    end select
  end subroutine load
  subroutine io_field(io,fh,ng,nh,lo,hi,disp,var)
    implicit none
    character(len=1), intent(in)                 :: io
    type(MPI_FILE)  , intent(in)                 :: fh
    integer , intent(in), dimension(3)           :: ng,nh,lo,hi
    integer(kind=MPI_OFFSET_KIND), intent(inout) :: disp
    real(rp), intent(inout), dimension(lo(1)-nh(1):,lo(2)-nh(2):,lo(3)-nh(3):) :: var
    integer , dimension(3) :: n
    integer , dimension(3) :: sizes,subsizes,starts
    type(MPI_DATATYPE) :: type_glob,type_loc
    n(:)        = hi(:)-lo(:)+1
    sizes(:)    = ng(:)
    subsizes(:) = n(:)
    starts(:)   = lo(:) - 1 ! starts from 0
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,MPI_REAL_RP,type_glob)
    call MPI_TYPE_COMMIT(type_glob)
    sizes(:)    = n(:) + 2*nh(:)
    subsizes(:) = n(:)
    starts(:)   = 0 + nh(:)
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,MPI_REAL_RP,type_loc )
    call MPI_TYPE_COMMIT(type_loc)
    select case(io)
    case('r')
      call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,type_glob,'native',MPI_INFO_NULL)
      call MPI_FILE_READ_ALL(fh,var,1,type_loc,MPI_STATUS_IGNORE)
    case('w')
      call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,type_glob,'native',MPI_INFO_NULL)
      call MPI_FILE_WRITE_ALL(fh,var,1,type_loc,MPI_STATUS_IGNORE)
    end select
    disp = disp+product(ng)*storage_size(1._rp)/8
    call MPI_TYPE_FREE(type_glob)
    call MPI_TYPE_FREE(type_loc )
  end subroutine io_field
end module mod_load
