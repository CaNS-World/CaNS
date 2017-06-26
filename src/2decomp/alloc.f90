!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2012 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Utility routine to help allocate 3D arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! X-pencil real arrays
  subroutine alloc_x_real(var, opt_decomp, opt_global)

    implicit none

    real(mytype), allocatable, dimension(:,:,:) :: var
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    logical, intent(IN), optional :: opt_global

    TYPE(DECOMP_INFO) :: decomp
    logical :: global
    integer :: alloc_stat, errorcode

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    if (present(opt_global)) then
       global = opt_global
    else
       global = .false.
    end if

    if (global) then
       allocate(var(decomp%xst(1):decomp%xen(1), &
            decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)), &
            stat=alloc_stat)
    else
       allocate(var(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)), &
            stat=alloc_stat)
    end if
    
    if (alloc_stat /= 0) then
       errorcode = 8
       call decomp_2d_abort(errorcode, &
            'Memory allocation failed when creating new arrays')
    end if

    return
  end subroutine alloc_x_real


  ! X-pencil complex arrays
  subroutine alloc_x_complex(var, opt_decomp, opt_global)

    implicit none

    complex(mytype), allocatable, dimension(:,:,:) :: var
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    logical, intent(IN), optional :: opt_global

    TYPE(DECOMP_INFO) :: decomp
    logical :: global
    integer :: alloc_stat, errorcode

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    if (present(opt_global)) then
       global = opt_global
    else
       global = .false.
    end if

    if (global) then
       allocate(var(decomp%xst(1):decomp%xen(1), &
            decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)), &
            stat=alloc_stat)
    else
       allocate(var(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)), &
            stat=alloc_stat)
    end if
    
    if (alloc_stat /= 0) then
       errorcode = 8
       call decomp_2d_abort(errorcode, &
            'Memory allocation failed when creating new arrays')
    end if

    return
  end subroutine alloc_x_complex


  ! Y-pencil real arrays
  subroutine alloc_y_real(var, opt_decomp, opt_global)

    implicit none

    real(mytype), allocatable, dimension(:,:,:) :: var
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    logical, intent(IN), optional :: opt_global

    TYPE(DECOMP_INFO) :: decomp
    logical :: global
    integer :: alloc_stat, errorcode

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    if (present(opt_global)) then
       global = opt_global
    else
       global = .false.
    end if

    if (global) then
       allocate(var(decomp%yst(1):decomp%yen(1), &
            decomp%yst(2):decomp%yen(2), decomp%yst(3):decomp%yen(3)), &
            stat=alloc_stat)
    else
       allocate(var(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)), &
            stat=alloc_stat)
    end if
    
    if (alloc_stat /= 0) then
       errorcode = 8
       call decomp_2d_abort(errorcode, &
            'Memory allocation failed when creating new arrays')
    end if

    return
  end subroutine alloc_y_real


  ! Y-pencil complex arrays
  subroutine alloc_y_complex(var, opt_decomp, opt_global)

    implicit none

    complex(mytype), allocatable, dimension(:,:,:) :: var
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    logical, intent(IN), optional :: opt_global

    TYPE(DECOMP_INFO) :: decomp
    logical :: global
    integer :: alloc_stat, errorcode

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    if (present(opt_global)) then
       global = opt_global
    else
       global = .false.
    end if

    if (global) then
       allocate(var(decomp%yst(1):decomp%yen(1), &
            decomp%yst(2):decomp%yen(2), decomp%yst(3):decomp%yen(3)), &
            stat=alloc_stat)
    else
       allocate(var(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)), &
            stat=alloc_stat)
    end if
    
    if (alloc_stat /= 0) then
       errorcode = 8
       call decomp_2d_abort(errorcode, &
            'Memory allocation failed when creating new arrays')
    end if

    return
  end subroutine alloc_y_complex


  ! Z-pencil real arrays
  subroutine alloc_z_real(var, opt_decomp, opt_global)

    implicit none

    real(mytype), allocatable, dimension(:,:,:) :: var
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    logical, intent(IN), optional :: opt_global

    TYPE(DECOMP_INFO) :: decomp
    logical :: global
    integer :: alloc_stat, errorcode

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    if (present(opt_global)) then
       global = opt_global
    else
       global = .false.
    end if

    if (global) then
       allocate(var(decomp%zst(1):decomp%zen(1), &
            decomp%zst(2):decomp%zen(2), decomp%zst(3):decomp%zen(3)), &
            stat=alloc_stat)
    else
       allocate(var(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3)), &
            stat=alloc_stat)
    end if
    
    if (alloc_stat /= 0) then
       errorcode = 8
       call decomp_2d_abort(errorcode, &
            'Memory allocation failed when creating new arrays')
    end if

    return
  end subroutine alloc_z_real


  ! Z-pencil complex arrays
  subroutine alloc_z_complex(var, opt_decomp, opt_global)

    implicit none

    complex(mytype), allocatable, dimension(:,:,:) :: var
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    logical, intent(IN), optional :: opt_global

    TYPE(DECOMP_INFO) :: decomp
    logical :: global
    integer :: alloc_stat, errorcode

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    if (present(opt_global)) then
       global = opt_global
    else
       global = .false.
    end if

    if (global) then
       allocate(var(decomp%zst(1):decomp%zen(1), &
            decomp%zst(2):decomp%zen(2), decomp%zst(3):decomp%zen(3)), &
            stat=alloc_stat)
    else
       allocate(var(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3)), &
            stat=alloc_stat)
    end if
    
    if (alloc_stat /= 0) then
       errorcode = 8
       call decomp_2d_abort(errorcode, &
            'Memory allocation failed when creating new arrays')
    end if

    return
  end subroutine alloc_z_complex
