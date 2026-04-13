! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017 Matthias Redies
! SPDX-License-Identifier: MIT
!
! -
!
! adapted from https://github.com/MRedies/NPY-for-Fortran (copyright in the header)
! Pedro Costa
!

module mod_npy
  use, intrinsic :: iso_fortran_env, only: r4 => real32, &
                                           r8 => real64, &
                                           i4 => int32
  implicit none
  character, parameter                :: magic_num = achar(147) ! x93
  character, parameter                :: major = achar(2)   !major *.npy version
  character, parameter                :: minor = achar(0)   !minor *.npy version
  character(len=*), parameter         :: magic_str = "NUMPY"
  integer(i4)                         :: iunit
  interface save_npy
     module procedure write_npy_1d_r4,write_npy_2d_r4,write_npy_3d_r4,write_npy_4d_r4, &
                      write_npy_1d_r8,write_npy_2d_r8,write_npy_3d_r8,write_npy_4d_r8, &
                      write_npy_1d_c4,write_npy_1d_c8,write_npy_2d_c4,write_npy_2d_c8, &
                      write_npy_3d_c4,write_npy_3d_c8,write_npy_4d_c4,write_npy_4d_c8
  end interface save_npy
  private
  public save_npy
  contains
  !
  subroutine write_npy_1d_r4(filename,arr)
    implicit none
    character(len=*), intent(in) :: filename
    real(r4), intent(in)         :: arr(:)
    character(len=*), parameter  :: var_type = "<f4"
    integer(i4)                  :: header_len
    !
    header_len = len(dict_str(var_type,shape(arr)))
    open(newunit=iunit, file=filename, access="stream", form="unformatted", status="replace")
    write(iunit) magic_num,magic_str,major,minor
    write(iunit) header_len
    write(iunit) dict_str(var_type,shape(arr))
    write(iunit) arr
    close(unit=iunit)
  end subroutine write_npy_1d_r4
  !
  subroutine write_npy_1d_r8(filename,arr)
    implicit none
    character(len=*), intent(in) :: filename
    real(r8), intent(in)         :: arr(:)
    character(len=*), parameter  :: var_type = "<f8"
    integer(i4)                  :: header_len
    !
    header_len = len(dict_str(var_type,shape(arr)))
    open(newunit=iunit, file=filename, access="stream", form="unformatted", status="replace")
    write(iunit) magic_num,magic_str,major,minor
    write(iunit) header_len
    write(iunit) dict_str(var_type,shape(arr))
    write(iunit) arr
    close(unit=iunit)
  end subroutine write_npy_1d_r8
  !
  subroutine write_npy_2d_r4(filename,arr)
    implicit none
    character(len=*), intent(in) :: filename
    real(r4), intent(in)         :: arr(:,:)
    character(len=*), parameter  :: var_type = "<f4"
    integer(i4)                  :: header_len
    !
    header_len = len(dict_str(var_type,shape(arr)))
    open(newunit=iunit, file=filename, access="stream", form="unformatted", status="replace")
    write(iunit) magic_num,magic_str,major,minor
    write(iunit) header_len
    write(iunit) dict_str(var_type,shape(arr))
    write(iunit) arr
    close(unit=iunit)
  end subroutine write_npy_2d_r4
  !
  subroutine write_npy_2d_r8(filename,arr)
    implicit none
    character(len=*), intent(in) :: filename
    real(r8), intent(in)         :: arr(:,:)
    character(len=*), parameter  :: var_type = "<f8"
    integer(i4)                  :: header_len
    !
    header_len = len(dict_str(var_type,shape(arr)))
    open(newunit=iunit, file=filename, access="stream", form="unformatted", status="replace")
    write(iunit) magic_num,magic_str,major,minor
    write(iunit) header_len
    write(iunit) dict_str(var_type,shape(arr))
    write(iunit) arr
    close(unit=iunit)
  end subroutine write_npy_2d_r8
  !
  subroutine write_npy_3d_r4(filename,arr)
    implicit none
    character(len=*), intent(in) :: filename
    real(r4), intent(in)         :: arr(:,:,:)
    character(len=*), parameter  :: var_type = "<f4"
    integer(i4)                  :: header_len
    !
    header_len = len(dict_str(var_type,shape(arr)))
    open(newunit=iunit, file=filename, access="stream", form="unformatted", status="replace")
    write(iunit) magic_num,magic_str,major,minor
    write(iunit) header_len
    write(iunit) dict_str(var_type,shape(arr))
    write(iunit) arr
    close(unit=iunit)
  end subroutine write_npy_3d_r4
  !
  subroutine write_npy_3d_r8(filename,arr)
    implicit none
    character(len=*), intent(in) :: filename
    real(r8), intent(in)         :: arr(:,:,:)
    character(len=*), parameter  :: var_type = "<f8"
    integer(i4)                  :: header_len
    !
    header_len = len(dict_str(var_type,shape(arr)))
    open(newunit=iunit, file=filename, access="stream", form="unformatted", status="replace")
    write(iunit) magic_num,magic_str,major,minor
    write(iunit) header_len
    write(iunit) dict_str(var_type,shape(arr))
    write(iunit) arr
    close(unit=iunit)
  end subroutine write_npy_3d_r8
  !
  subroutine write_npy_4d_r4(filename,arr)
    implicit none
    character(len=*), intent(in) :: filename
    real(r4), intent(in)         :: arr(:,:,:,:)
    character(len=*), parameter  :: var_type = "<f4"
    integer(i4)                  :: header_len
    !
    header_len = len(dict_str(var_type,shape(arr)))
    open(newunit=iunit, file=filename, access="stream", form="unformatted", status="replace")
    write(iunit) magic_num,magic_str,major,minor
    write(iunit) header_len
    write(iunit) dict_str(var_type,shape(arr))
    write(iunit) arr
    close(unit=iunit)
  end subroutine write_npy_4d_r4
  !
  subroutine write_npy_4d_r8(filename,arr)
    implicit none
    character(len=*), intent(in) :: filename
    real(r8), intent(in)         :: arr(:,:,:,:)
    character(len=*), parameter  :: var_type = "<f8"
    integer(i4)                  :: header_len
    !
    header_len = len(dict_str(var_type,shape(arr)))
    open(newunit=iunit, file=filename, access="stream", form="unformatted", status="replace")
    write(iunit) magic_num,magic_str,major,minor
    write(iunit) header_len
    write(iunit) dict_str(var_type,shape(arr))
    write(iunit) arr
    close(unit=iunit)
  end subroutine write_npy_4d_r8
  !
  subroutine write_npy_1d_c4(filename,arr)
    implicit none
    character(len=*), intent(in) :: filename
    complex(r4), intent(in)      :: arr(:)
    character(len=*), parameter  :: var_type = "<c8"
    integer(i4)                  :: header_len
    !
    header_len = len(dict_str(var_type,shape(arr)))
    open(newunit=iunit, file=filename, access="stream", form="unformatted", status="replace")
    write(iunit) magic_num,magic_str,major,minor
    write(iunit) header_len
    write(iunit) dict_str(var_type,shape(arr))
    write(iunit) arr
    close(unit=iunit)
  end subroutine write_npy_1d_c4
  !
  subroutine write_npy_1d_c8(filename,arr)
    implicit none
    character(len=*), intent(in) :: filename
    complex(r8), intent(in)      :: arr(:)
    character(len=*), parameter  :: var_type = "<c16"
    integer(i4)                  :: header_len
    !
    header_len = len(dict_str(var_type,shape(arr)))
    open(newunit=iunit, file=filename, access="stream", form="unformatted", status="replace")
    write(iunit) magic_num,magic_str,major,minor
    write(iunit) header_len
    write(iunit) dict_str(var_type,shape(arr))
    write(iunit) arr
    close(unit=iunit)
  end subroutine write_npy_1d_c8
  !
  subroutine write_npy_2d_c4(filename,arr)
    implicit none
    character(len=*), intent(in)      :: filename
    complex(r4), intent(in)           :: arr(:,:)
    character(len=*), parameter       :: var_type = "<c8"
    integer(i4)                       :: header_len
    !
    header_len = len(dict_str(var_type,shape(arr)))
    open(newunit=iunit, file=filename, access="stream", form="unformatted", status="replace")
    write(iunit) magic_num,magic_str,major,minor
    write(iunit) header_len
    write(iunit) dict_str(var_type,shape(arr))
    write(iunit) arr
    close(unit=iunit)
  end subroutine write_npy_2d_c4
  !
  subroutine write_npy_2d_c8(filename,arr)
    implicit none
    character(len=*), intent(in) :: filename
    complex(r8), intent(in)      :: arr(:,:)
    character(len=*), parameter  :: var_type = "<c16"
    integer(i4)                  :: header_len
    !
    header_len = len(dict_str(var_type,shape(arr)))
    open(newunit=iunit, file=filename, access="stream", form="unformatted", status="replace")
    write(iunit) magic_num,magic_str,major,minor
    write(iunit) header_len
    write(iunit) dict_str(var_type,shape(arr))
    write(iunit) arr
    close(unit=iunit)
  end subroutine write_npy_2d_c8
  !
  subroutine write_npy_3d_c4(filename,arr)
    implicit none
    character(len=*), intent(in) :: filename
    complex(r4), intent(in)      :: arr(:,:,:)
    character(len=*), parameter  :: var_type = "<c8"
    integer(i4)                  :: header_len
    !
    header_len = len(dict_str(var_type,shape(arr)))
    open(newunit=iunit, file=filename, access="stream", form="unformatted", status="replace")
    write(iunit) magic_num,magic_str,major,minor
    write(iunit) header_len
    write(iunit) dict_str(var_type,shape(arr))
    write(iunit) arr
    close(unit=iunit)
  end subroutine write_npy_3d_c4
  !
  subroutine write_npy_3d_c8(filename,arr)
    implicit none
    character(len=*), intent(in) :: filename
    complex(r8), intent(in)      :: arr(:,:,:)
    character(len=*), parameter  :: var_type = "<c16"
    integer(i4)                  :: header_len
    !
    header_len = len(dict_str(var_type,shape(arr)))
    open(newunit=iunit, file=filename, access="stream", form="unformatted", status="replace")
    write(iunit) magic_num,magic_str,major,minor
    write(iunit) header_len
    write(iunit) dict_str(var_type,shape(arr))
    write(iunit) arr
    close(unit=iunit)
  end subroutine write_npy_3d_c8
  !
  subroutine write_npy_4d_c4(filename,arr)
    implicit none
    character(len=*), intent(in) :: filename
    complex(r4), intent(in)      :: arr(:,:,:,:)
    character(len=*), parameter  :: var_type = "<c8"
    integer(i4)                  :: header_len
    !
    header_len = len(dict_str(var_type,shape(arr)))
    open(newunit=iunit, file=filename, access="stream", form="unformatted", status="replace")
    write(iunit) magic_num,magic_str,major,minor
    write(iunit) header_len
    write(iunit) dict_str(var_type,shape(arr))
    write(iunit) arr
    close(unit=iunit)
  end subroutine write_npy_4d_c4
  !
  subroutine write_npy_4d_c8(filename,arr)
    implicit none
    character(len=*), intent(in) :: filename
    complex(r8), intent(in)      :: arr(:,:,:,:)
    character(len=*), parameter  :: var_type = "<c16"
    integer(i4)                  :: header_len
    !
    header_len = len(dict_str(var_type,shape(arr)))
    open(newunit=iunit, file=filename, access="stream", form="unformatted", status="replace")
    write(iunit) magic_num,magic_str,major,minor
    write(iunit) header_len
    write(iunit) dict_str(var_type,shape(arr))
    write(iunit) arr
    close(unit=iunit)
  end subroutine write_npy_4d_c8
  !
  function dict_str(var_type,var_shape) result(str)
    implicit none
    character(len=*), intent(in)  :: var_type
    integer(i4), intent(in)       :: var_shape(:)
    character(len=:), allocatable :: str
    integer(i4)                   :: cnt
    !
    cnt = len("{'descr': '")
    cnt = cnt + len(var_type)
    cnt = cnt + len("', 'fortran_order': True, 'shape': (")
    cnt = cnt + len(shape_str(var_shape))
    cnt = cnt + len(",), }")
    do while (mod(cnt + 10, 16) /= 0)
      cnt = cnt + 1
    end do
    !
    allocate(character(cnt) :: str)
    !
    str = "{'descr': '"//var_type// &
          "', 'fortran_order': True, 'shape': ("// &
          shape_str(var_shape)//"), }"
    !
    do while (mod(len(str) + 11, 16) /= 0)
      str = str//" "
    end do
    !
    str = str//achar(10)
  end function dict_str
  !
  function shape_str(var_shape) result(fin_str)
    implicit none
    integer(i4), intent(in)       :: var_shape(:)
    character(len=:), allocatable :: str, small_str, fin_str
    integer(i4)                   :: i, length, start, halt
    !
    length = 14*size(var_shape)
    allocate(character(length) :: str)
    allocate(character(14)     :: small_str)
    str = " "
    !
    do i = 1,size(var_shape)
      start = (i - 1)*length + 1
      halt = i*length + 1
      write(small_str,"(i13,a)") var_shape(i), ","
      str = trim(str)//adjustl(small_str)
    end do
    !
    fin_str = trim(str)
  end function shape_str
end module mod_npy
