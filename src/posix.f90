! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
!
! file-system operations through the C standard library (zero on success)
!
module mod_posix
  use, intrinsic :: iso_c_binding, only: c_char,c_int,c_long,c_null_char,c_size_t
  implicit none
  private
  public link_file,rename_file,remove_file,read_link
  interface
    function c_symlink(path_target,path_link) bind(C,name='symlink') result(istat)
      import :: c_char,c_int
      character(kind=c_char), intent(in) :: path_target(*),path_link(*)
      integer(c_int) :: istat
    end function c_symlink
    function c_rename(path_old,path_new) bind(C,name='rename') result(istat)
      import :: c_char,c_int
      character(kind=c_char), intent(in) :: path_old(*),path_new(*)
      integer(c_int) :: istat
    end function c_rename
    function c_unlink(path) bind(C,name='unlink') result(istat)
      import :: c_char,c_int
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int) :: istat
    end function c_unlink
    function c_readlink(path,buf,bufsiz) bind(C,name='readlink') result(nbytes)
      import :: c_char,c_long,c_size_t
      character(kind=c_char), intent(in) :: path(*)
      character(kind=c_char) :: buf(*)
      integer(c_size_t), value :: bufsiz
      integer(c_long) :: nbytes
    end function c_readlink
  end interface
contains
  integer function link_file(ftarget,flink) result(istat)
    implicit none
    character(len=*), intent(in) :: ftarget,flink
    istat = c_symlink(trim(ftarget)//c_null_char,trim(flink)//c_null_char)
  end function link_file
  integer function rename_file(fold,fnew) result(istat)
    implicit none
    character(len=*), intent(in) :: fold,fnew
    istat = c_rename(trim(fold)//c_null_char,trim(fnew)//c_null_char)
  end function rename_file
  integer function remove_file(fname) result(istat)
    implicit none
    character(len=*), intent(in) :: fname
    istat = c_unlink(trim(fname)//c_null_char)
  end function remove_file
  integer function read_link(flink,ftarget) result(istat)
    implicit none
    character(len=*), intent(in ) :: flink
    character(len=*), intent(out) :: ftarget
    character(kind=c_char) :: buf(len(ftarget))
    integer :: nbytes,i
    nbytes = int(c_readlink(trim(flink)//c_null_char,buf,int(len(ftarget),c_size_t)))
    ftarget = ''
    do i=1,nbytes
      ftarget(i:i) = buf(i)
    end do
    istat = merge(0,1,nbytes >= 0)
  end function read_link
end module mod_posix
