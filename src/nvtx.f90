! -
!
! SPDX-FileCopyrightText: Copyright (c) 2022 Pedro Costa. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
!
! NVTX Fortran Module, adapted from https://github.com/maxcuda/NVTX_example (MIT)
!
module mod_nvtx
#if defined(_USE_NVTX)
  use, intrinsic :: iso_c_binding
  implicit none
  private
  public nvtxStartRange,nvtxEndRange
  enum, bind(c)
    enumerator :: COLOR_G = 1 ! green
    enumerator :: COLOR_B = 2 ! blue
    enumerator :: COLOR_Y = 3 ! yellow
    enumerator :: COLOR_M = 4 ! magenta
    enumerator :: COLOR_C = 5 ! cyan
    enumerator :: COLOR_R = 6 ! red
    enumerator :: COLOR_W = 7 ! white
  end enum
  integer(kind=C_INT32_T), private :: col(7) = [ &
                                                 int(Z'0000ff00',kind=C_INT32_T), & ! 1 -> green
                                                 int(Z'000000ff',kind=C_INT32_T), & ! 2 -> blue
                                                 int(Z'00ffff00',kind=C_INT32_T), & ! 3 -> yellow
                                                 int(Z'00ff00ff',kind=C_INT32_T), & ! 4 -> magenta
                                                 int(Z'0000ffff',kind=C_INT32_T), & ! 5 -> cyan
                                                 int(Z'00ff0000',kind=C_INT32_T), & ! 6 -> red
                                                 int(Z'00ffffff',kind=C_INT32_T)  & ! 7 -> white
                                               ]
  integer  , parameter :: NAME_LEN_MAX = 256
  character, target :: tempName(NAME_LEN_MAX)
  type, bind(C) :: nvtxEventAttributes
    integer(C_INT16_T) :: version=1
    integer(C_INT16_T) :: size=48 !
    integer(C_INT) :: category=0
    integer(C_INT) :: colorType=1 ! NVTX_COLOR_ARGB = 1
    integer(C_INT) :: color
    integer(C_INT) :: payloadType=0 ! NVTX_PAYLOAD_UNKNOWN = 0
    integer(C_INT) :: reserved0
    integer(C_INT64_T) :: payload   ! union uint,int,double
    integer(C_INT) :: messageType=1  ! NVTX_MESSAGE_TYPE_ASCII     = 1
    type(C_PTR) :: message  ! ascii char
  end type
  interface nvtxRangePush
    ! push range with custom label and standard color
    subroutine nvtxRangePushA(name) bind(C, name='nvtxRangePushA')
    use, intrinsic :: iso_c_binding
    character(kind=C_CHAR) :: name(*)
    end subroutine
    ! push range with custom label and custom color
    subroutine nvtxRangePushEx(event) bind(C, name='nvtxRangePushEx')
    use, intrinsic :: iso_c_binding
    import :: nvtxEventAttributes
    type(nvtxEventAttributes) :: event
    end subroutine
  end interface
  interface nvtxRangePop
    subroutine nvtxRangePop() bind(C, name='nvtxRangePop')
    end subroutine
  end interface
contains
  subroutine nvtxStartRange(name,id,color)
    character(kind=c_char,len=*) :: name
    integer, optional :: id
    character(len=1), optional :: color ! g/b/y/m/c/r/w following matplotlib's convention
    type(nvtxEventAttributes) :: event
    character(kind=c_char,len=NAME_LEN_MAX) :: trimmed_name
    integer :: i,icolor
    trimmed_name=trim(name)//c_null_char
    ! move scalar trimmed_name into character array tempName
    do i=1,len(trim(name))+1
      tempName(i) = trimmed_name(i:i)
    end do
    if(present(color)) then
      select case(color)
      case('g')
        icolor = COLOR_G
      case('b')
        icolor = COLOR_B
      case('y')
        icolor = COLOR_Y
      case('m')
        icolor = COLOR_M
      case('c')
        icolor = COLOR_C
      case('r')
        icolor = COLOR_R
      case('w')
        icolor = COLOR_W
      case default
        icolor = COLOR_W
      end select
    else if(present(id)) then
      icolor = mod(id-1,size(col))+1
    end if
    if (present(id).or.present(color)) then
      event%color=col(icolor)
      event%message=c_loc(tempName)
      call nvtxRangePushEx(event)
    else
      call nvtxRangePush(tempName)
    end if
  end subroutine
  subroutine nvtxEndRange
    call nvtxRangePop
  end subroutine
#endif
end module mod_nvtx
