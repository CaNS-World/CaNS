!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2021 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This file contains subroutines that generate ACML plans
! for several types of 1D multiple FFTs.

! Note most ACML plans can be shared by forward/backward transforms

  ! Return an ACML plan for multiple 1D c2c FFTs in X direction
  subroutine c2c_1m_x_plan(comm, decomp)

    implicit none
    complex(mytype), allocatable, dimension(:), intent(OUT) :: comm
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    complex(mytype), allocatable, dimension(:,:,:) :: dummy
    integer :: info

#ifdef DOUBLE_PREC
    allocate(comm(3*decomp%xsz(1)+100))
#else
    allocate(comm(5*decomp%xsz(1)+100))
#endif

    allocate(dummy(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)))

#ifdef DOUBLE_PREC
    call zfft1mx(plan_type,scale,.true.,decomp%xsz(2)*decomp%xsz(3), &
         decomp%xsz(1), dummy,1,decomp%xsz(1),dummy,1,decomp%xsz(1), &
         comm,info)
#else
    call cfft1mx(plan_type,scale,.true.,decomp%xsz(2)*decomp%xsz(3), &
         decomp%xsz(1), dummy,1,decomp%xsz(1),dummy,1,decomp%xsz(1), &
         comm,info)
#endif

    deallocate(dummy)

    return
  end subroutine c2c_1m_x_plan

  ! Return an ACML plan for multiple 1D c2c FFTs in Y direction
  subroutine c2c_1m_y_plan(comm, decomp)

    implicit none
    complex(mytype), allocatable, dimension(:), intent(OUT) :: comm
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    complex(mytype), allocatable, dimension(:,:,:) :: dummy
    integer :: info

#ifdef DOUBLE_PREC
    allocate(comm(3*decomp%ysz(2)+100))
#else
    allocate(comm(5*decomp%ysz(2)+100))
#endif

    allocate(dummy(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))

#ifdef DOUBLE_PREC
    call zfft1mx(plan_type,scale,.true.,decomp%ysz(1),decomp%ysz(2), &
            dummy,decomp%ysz(1),1,dummy,decomp%ysz(1),1,comm,info)
#else
    call cfft1mx(plan_type,scale,.true.,decomp%ysz(1),decomp%ysz(2), &
            dummy,decomp%ysz(1),1,dummy,decomp%ysz(1),1,comm,info)
#endif

    deallocate(dummy)

    return
  end subroutine c2c_1m_y_plan

  ! Return an ACML plan for multiple 1D c2c FFTs in Z direction  
  subroutine c2c_1m_z_plan(comm, decomp)

    implicit none
    complex(mytype), allocatable, dimension(:), intent(OUT) :: comm
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    complex(mytype), allocatable, dimension(:,:,:) :: dummy
    integer :: info

#ifdef DOUBLE_PREC
    allocate(comm(3*decomp%zsz(3)+100))
#else
    allocate(comm(5*decomp%zsz(3)+100))
#endif

    allocate(dummy(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3)))

#ifdef DOUBLE_PREC
    call zfft1mx(plan_type,scale,.true.,decomp%zsz(1)*decomp%zsz(2), &
         decomp%zsz(3),dummy,decomp%zsz(1)*decomp%zsz(2),1,dummy, &
         decomp%zsz(1)*decomp%zsz(2),1,comm,info)
#else
    call cfft1mx(plan_type,scale,.true.,decomp%zsz(1)*decomp%zsz(2), &
         decomp%zsz(3),dummy,decomp%zsz(1)*decomp%zsz(2),1,dummy, &
         decomp%zsz(1)*decomp%zsz(2),1,comm,info)
#endif

    deallocate(dummy)

    return
  end subroutine c2c_1m_z_plan


  ! Return an ACML plan for multiple 1D r2c FFTs in X direction
  subroutine r2c_1m_x_plan(comm, decomp)

    implicit none
    real(mytype), allocatable, dimension(:), intent(OUT) :: comm
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    real(mytype), allocatable, dimension(:) :: dummy
    integer :: info

    allocate(comm(3*decomp%xsz(1)+100))

    allocate(dummy(decomp%xsz(1)))

#ifdef DOUBLE_PREC
    call dzfft(plan_type,decomp%xsz(1),dummy,comm,info)
#else
    call scfft(plan_type,decomp%xsz(1),dummy,comm,info)
#endif

    deallocate(dummy)

    return
  end subroutine r2c_1m_x_plan

  ! Return an ACML plan for multiple 1D c2r FFTs in X direction
  subroutine c2r_1m_x_plan(comm, decomp)

    implicit none
    real(mytype), allocatable, dimension(:), intent(OUT) :: comm
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    real(mytype), allocatable, dimension(:) :: dummy
    integer :: info

    allocate(comm(3*decomp%xsz(1)+100))

    allocate(dummy(decomp%xsz(1)))

#ifdef DOUBLE_PREC
    call zdfft(plan_type,decomp%xsz(1),dummy,comm,info)
#else
    call csfft(plan_type,decomp%xsz(1),dummy,comm,info)
#endif

    deallocate(dummy)

    return
  end subroutine c2r_1m_x_plan

  ! Return an ACML plan for multiple 1D r2c FFTs in Z direction
  subroutine r2c_1m_z_plan(comm, decomp)

    implicit none
    real(mytype), allocatable, dimension(:), intent(OUT) :: comm
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    real(mytype), allocatable, dimension(:) :: dummy
    integer :: info

    allocate(comm(3*decomp%zsz(3)+100))

    allocate(dummy(decomp%zsz(3)))

#ifdef DOUBLE_PREC
    call dzfft(plan_type,decomp%zsz(3),dummy,comm,info)
#else
    call scfft(plan_type,decomp%zsz(3),dummy,comm,info)
#endif

    deallocate(dummy)

    return
  end subroutine r2c_1m_z_plan

  ! Return an ACML plan for multiple 1D c2r FFTs in Z direction
  subroutine c2r_1m_z_plan(comm, decomp)

    implicit none
    real(mytype), allocatable, dimension(:), intent(OUT) :: comm
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    real(mytype), allocatable, dimension(:) :: dummy
    integer :: info

    allocate(comm(3*decomp%zsz(3)+100))

    allocate(dummy(decomp%zsz(3)))

#ifdef DOUBLE_PREC
    call zdfft(plan_type,decomp%zsz(3),dummy,comm,info)
#else
    call csfft(plan_type,decomp%zsz(3),dummy,comm,info)
#endif

    deallocate(dummy)

    return
  end subroutine c2r_1m_z_plan
