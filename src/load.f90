module mod_load
  use mpi
  use mod_common_mpi, only: ierr,ipencil,dims,myid
  use decomp_2d
  use decomp_2d_io
  use mod_types
  implicit none
  private
  public load
  contains
  subroutine load(io,filename,n,u,v,w,p,time,istep)
    !
    ! reads/writes a restart file
    !
    implicit none
    character(len=1)  , intent(in) :: io
    character(len=*), intent(in) :: filename
    integer , intent(in), dimension(3) :: n
    real(rp), intent(inout), dimension(n(1),n(2),n(3)) :: u,v,w,p
    real(rp), intent(inout) :: time,istep
    real(rp), dimension(2) :: fldinfo
    integer :: fh
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp,good
    integer(8), dimension(3) :: ng
    integer :: lenr
    !
    select case(io)
    case('r')
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename                 , &
           MPI_MODE_RDONLY, MPI_INFO_NULL,fh, ierr)
      !
      ! check file size first
      !
      call MPI_FILE_GET_SIZE(fh,filesize,ierr)
      ng(1:3) = n(1:3)*dims(1:3)
      lenr = sizeof(time)
      good = (product(ng)*4+2)*lenr
      if(filesize.ne.good) then
        if(myid.eq.0) print*, ''
        if(myid.eq.0) print*, '*** Simulation aborted due a checkpoint file with incorrect size ***'
        if(myid.eq.0) print*, '    file: ', filename, ' | expected size: ', good, '| actual size: ', filesize
        call decomp_2d_finalize
        call MPI_FINALIZE(ierr)
        call exit
      endif
      !
      ! read
      !
      disp = 0_MPI_OFFSET_KIND
      call decomp_2d_read_var(fh,disp,ipencil,u)
      call decomp_2d_read_var(fh,disp,ipencil,v)
      call decomp_2d_read_var(fh,disp,ipencil,w)
      call decomp_2d_read_var(fh,disp,ipencil,p)
      call decomp_2d_read_scalar(fh,disp,2,fldinfo)
      time  = fldinfo(1)
      istep = fldinfo(2)
      call MPI_FILE_CLOSE(fh,ierr)
    case('w')
      !
      ! write
      !
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename                 , &
           MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
      filesize = 0_MPI_OFFSET_KIND
      call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
      disp = 0_MPI_OFFSET_KIND
      call decomp_2d_write_var(fh,disp,ipencil,u)
      call decomp_2d_write_var(fh,disp,ipencil,v)
      call decomp_2d_write_var(fh,disp,ipencil,w)
      call decomp_2d_write_var(fh,disp,ipencil,p)
      fldinfo = (/time,istep/)
      call decomp_2d_write_scalar(fh,disp,2,fldinfo)
      call MPI_FILE_CLOSE(fh,ierr)
    end select
    return
  end subroutine load
end module mod_load
