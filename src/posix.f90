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
  use, intrinsic :: iso_c_binding, only: c_char,c_int,c_null_char
  implicit none
  private
  public link_file,rename_file,remove_file
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
end module mod_posix
