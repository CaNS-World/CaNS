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

! This is the ESSL implementation of the FFT library

module decomp_2d_fft
  
  use decomp_2d  ! 2D decomposition module
  use glassman
  
  implicit none
  
  private        ! Make everything private unless declared public

  ! engine-specific global variables

  ! always perform unscaled FFTs
  real(mytype), parameter :: scale=1.0_mytype

  ! number of supported transform lengths
  integer, parameter :: MAGIC=451

  ! work space
  integer, save :: naux1, naux2, naux3
  double precision, allocatable, dimension(:) :: aux1, aux2, aux3
  complex(mytype), allocatable, dimension(:) :: buf, scratch

  ! ESSL has a list of acceptable lengths for transforms
  !  - if winthin this list, use ESSL calls directly
  !  - if not, use the less efficient generic implementation
  logical, save :: x_goodsize, y_goodsize, z_goodsize

  ! common code used for all engines, including global variables, 
  ! generic interface definitions and several subroutines
#include "fft_common.f90"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time initialisations for the FFT engine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_fft_engine

    implicit none

    integer, dimension(MAGIC) :: lengths
    integer :: bufsize, tmp1, tmp2

    if (nrank==0) then
       write(*,*) ' '
       write(*,*) '***** Using the ESSL engine *****'
       write(*,*) ' '
    end if

    ! Most ESSL FFT routines work on a limited set of sample size
    call get_acceptable_length(lengths)

    if (binarySearch(lengths, nx_fft)==0) then
       x_goodsize = .false.
    else
       x_goodsize = .true.
    end if
    if (binarySearch(lengths, ny_fft)==0) then
       y_goodsize = .false.
    else
       y_goodsize = .true.
    end if
    if (binarySearch(lengths, nz_fft)==0) then
       z_goodsize = .false.
    else
       z_goodsize = .true.
    end if

    if (nrank==0) then
       if ( (.not.x_goodsize) .or. (.not.y_goodsize) &
            .or. (.not.z_goodsize) ) then
          write(*,*) '***WARNING***: sample size not supported by ESSL.'
          write(*,*) 'Please consult the ESSL documentation for detail.'
          write(*,*) '2DECOMP&FFT will use alternative FFT algorithms.'
          write(*,*) 'You will still get the correct results.'
          write(*,*) 'But performance of your application may suffer.'
       end if
    end if

    ! allocate ESSL work space
    call cftd_wk(naux1, naux2)
    call cft_wk(tmp1, tmp2)
    naux1 = max(naux1, tmp1)
    naux2 = max(naux2, tmp2)
    call rcft_wk(tmp1, tmp2)
    naux1 = max(naux1, tmp1)
    naux2 = max(naux2, tmp2)
    naux3 = 1 ! not in use by ESSL but exists for migration purpose
    !if (nrank==0) then
    !   write(*,*) '3D problem size:', ph%xsz(1), ph%ysz(2), ph%zsz(3)
    !   write(*,*) 'Work array size:', naux1, naux2
    !end if
    allocate(aux1(naux1))
    allocate(aux2(naux2))
    allocate(aux3(naux3))

    ! work space by generic algorithm for FFT sizes unsupported by ESSL
    bufsize = max(ph%xsz(1),ph%zsz(3))
    allocate(buf(bufsize))
    allocate(scratch(bufsize))

    return
  end subroutine init_fft_engine
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time finalisations for the FFT engine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine finalize_fft_engine

    implicit none

    deallocate(aux1, aux2, aux3)
    deallocate(buf, scratch)

    return
  end subroutine finalize_fft_engine


  ! Following routines calculate multiple one-dimensional FFTs to form 
  ! the basis of three-dimensional FFTs.

  ! *** NOTE, in ESSL forward transforms have isign>0 and backward ones
  ! have isign<0. This is different from other FFT libraries. 

  ! c2c transform, multiple 1D FFTs in x direction
  subroutine c2c_1m_x(inout, isign, decomp)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    if (x_goodsize) then
#ifdef DOUBLE_PREC
       call dcft(1,inout,1,decomp%xsz(1),inout,1,decomp%xsz(1), &
            decomp%xsz(1),decomp%xsz(2)*decomp%xsz(3) &
            ,-isign,scale,aux1,naux1,aux2,naux2)
       call dcft(0,inout,1,decomp%xsz(1),inout,1,decomp%xsz(1), &
            decomp%xsz(1),decomp%xsz(2)*decomp%xsz(3) &
            ,-isign,scale,aux1,naux1,aux2,naux2)
#else
       call scft(1,inout,1,decomp%xsz(1),inout,1,decomp%xsz(1), &
            decomp%xsz(1),decomp%xsz(2)*decomp%xsz(3) &
            ,-isign,scale,aux1,naux1,aux2,naux2)
       call scft(0,inout,1,decomp%xsz(1),inout,1,decomp%xsz(1), &
            decomp%xsz(1),decomp%xsz(2)*decomp%xsz(3) &
            ,-isign,scale,aux1,naux1,aux2,naux2)
#endif
    else
#ifdef DOUBLE_PREC
       call dcftd(1,1,inout,1,decomp%xsz(1),inout,1,decomp%xsz(1), &
            decomp%xsz(1),decomp%xsz(2)*decomp%xsz(3) &
            ,-isign,scale,aux1,naux1,aux2,naux2)
       call dcftd(0,1,inout,1,decomp%xsz(1),inout,1,decomp%xsz(1), &
            decomp%xsz(1),decomp%xsz(2)*decomp%xsz(3) &
            ,-isign,scale,aux1,naux1,aux2,naux2)
#else
       call scftd(1,1,inout,1,decomp%xsz(1),inout,1,decomp%xsz(1), &
            decomp%xsz(1),decomp%xsz(2)*decomp%xsz(3) &
            ,-isign,scale,aux1,naux1,aux2,naux2)
       call scftd(0,1,inout,1,decomp%xsz(1),inout,1,decomp%xsz(1), &
            decomp%xsz(1),decomp%xsz(2)*decomp%xsz(3) &
            ,-isign,scale,aux1,naux1,aux2,naux2)
#endif          
    end if
    
    return
    
  end subroutine c2c_1m_x


  ! c2c transform, multiple 1D FFTs in y direction
  subroutine c2c_1m_y(inout, isign, decomp)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: k

    if (y_goodsize) then
#ifdef DOUBLE_PREC
       call dcft(1,inout,decomp%ysz(1),1,inout,decomp%ysz(1),1, &
            decomp%ysz(2),decomp%ysz(1), &
            -isign,scale,aux1,naux1,aux2,naux2)
       do k=1,decomp%ysz(3)
          call dcft(0,inout(:,:,k),decomp%ysz(1),1,inout(:,:,k), &
               decomp%ysz(1),1,decomp%ysz(2),decomp%ysz(1), &
               -isign,scale,aux1,naux1,aux2,naux2)
       end do
#else
       call scft(1,inout,decomp%ysz(1),1,inout,decomp%ysz(1),1, &
            decomp%ysz(2),decomp%ysz(1), &
            -isign,scale,aux1,naux1,aux2,naux2)
       do k=1,decomp%ysz(3)
          call scft(0,inout(:,:,k),decomp%ysz(1),1,inout(:,:,k), &
               decomp%ysz(1),1,decomp%ysz(2),decomp%ysz(1), &
               -isign,scale,aux1,naux1,aux2,naux2)
       end do
#endif    
    else
#ifdef DOUBLE_PREC
       call dcftd(1,1,inout,decomp%ysz(1),1,inout,decomp%ysz(1),1, &
            decomp%ysz(2),decomp%ysz(1), &
            -isign,scale,aux1,naux1,aux2,naux2)
       do k=1,decomp%ysz(3)
          call dcftd(0,1,inout(:,:,k),decomp%ysz(1),1,inout(:,:,k), &
               decomp%ysz(1),1,decomp%ysz(2),decomp%ysz(1), &
               -isign,scale,aux1,naux1,aux2,naux2)
       enddo
#else
       call scftd(1,1,inout,decomp%ysz(1),1,inout,decomp%ysz(1),1, &
            decomp%ysz(2),decomp%ysz(1), &
            -isign,scale,aux1,naux1,aux2,naux2)
       do k=1,decomp%ysz(3)
          call scftd(0,1,inout(:,:,k),decomp%ysz(1),1,inout(:,:,k), &
               decomp%ysz(1),1,decomp%ysz(2),decomp%ysz(1), &
               -isign,scale,aux1,naux1,aux2,naux2)
       enddo
#endif
    end if

    return

  end subroutine c2c_1m_y


  ! c2c transform, multiple 1D FFTs in z direction
  subroutine c2c_1m_z(inout, isign, decomp)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT)  :: inout
    integer, intent(IN) :: isign
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    if (z_goodsize) then
#ifdef DOUBLE_PREC
       call dcft(1,inout,decomp%zsz(1)*decomp%zsz(2),1,inout, &
            decomp%zsz(1)*decomp%zsz(2),1,decomp%zsz(3), &
            decomp%zsz(1)*decomp%zsz(2), &
            -isign,scale,aux1,naux1,aux2,naux2)
       call dcft(0,inout,decomp%zsz(1)*decomp%zsz(2),1,inout, &
            decomp%zsz(1)*decomp%zsz(2),1,decomp%zsz(3), &
            decomp%zsz(1)*decomp%zsz(2), &
            -isign,scale,aux1,naux1,aux2,naux2)
#else
       call scft(1,inout,decomp%zsz(1)*decomp%zsz(2),1,inout, &
            decomp%zsz(1)*decomp%zsz(2),1,decomp%zsz(3), &
            decomp%zsz(1)*decomp%zsz(2), &
            -isign,scale,aux1,naux1,aux2,naux2)
       call scft(0,inout,decomp%zsz(1)*decomp%zsz(2),1,inout, &
            decomp%zsz(1)*decomp%zsz(2),1,decomp%zsz(3), &
            decomp%zsz(1)*decomp%zsz(2), &
            -isign,scale,aux1,naux1,aux2,naux2)
#endif
    else
#ifdef DOUBLE_PREC
       call dcftd(1,1,inout,decomp%zsz(1)*decomp%zsz(2),1,inout, &
            decomp%zsz(1)*decomp%zsz(2),1,decomp%zsz(3), &
            decomp%zsz(1)*decomp%zsz(2), &
            -isign,scale,aux1,naux1,aux2,naux2)
       call dcftd(0,1,inout,decomp%zsz(1)*decomp%zsz(2),1,inout, &
            decomp%zsz(1)*decomp%zsz(2),1,decomp%zsz(3), &
            decomp%zsz(1)*decomp%zsz(2), &
            -isign,scale,aux1,naux1,aux2,naux2)
#else
       call scftd(1,1,inout,decomp%zsz(1)*decomp%zsz(2),1,inout, &
            decomp%zsz(1)*decomp%zsz(2),1,decomp%zsz(3), &
            decomp%zsz(1)*decomp%zsz(2), &
            -isign,scale,aux1,naux1,aux2,naux2)
       call scftd(0,1,inout,decomp%zsz(1)*decomp%zsz(2),1,inout, &
            decomp%zsz(1)*decomp%zsz(2),1,decomp%zsz(3), &
            decomp%zsz(1)*decomp%zsz(2), &
            -isign,scale,aux1,naux1,aux2,naux2)
#endif
    end if

    return

  end subroutine c2c_1m_z


  ! r2c transform, multiple 1D FFTs in x direction
  subroutine r2c_1m_x(input, output)

    implicit none
    
    real(mytype), dimension(:,:,:), intent(IN)  ::  input
    complex(mytype), dimension(:,:,:), intent(OUT) :: output

    integer :: i,j,k
    
    if (x_goodsize) then
#ifdef DOUBLE_PREC
       call drcft(1,input,ph%xsz(1),output,sp%xsz(1), &
            ph%xsz(1),ph%xsz(2)*ph%xsz(3),1,scale, &
            aux1,naux1,aux2,naux2,aux3,naux3)
       call drcft(0,input,ph%xsz(1),output,sp%xsz(1), &
            ph%xsz(1),ph%xsz(2)*ph%xsz(3),1,scale, &
            aux1,naux1,aux2,naux2,aux3,naux3)
#else
       call srcft(1,input,ph%xsz(1),output,sp%xsz(1), &
            ph%xsz(1),ph%xsz(2)*ph%xsz(3),1,scale, &
            aux1,naux1,aux2,naux2,aux3,naux3)
       call srcft(0,input,ph%xsz(1),output,sp%xsz(1), &
            ph%xsz(1),ph%xsz(2)*ph%xsz(3),1,scale, &
            aux1,naux1,aux2,naux2,aux3,naux3)
#endif
    else
       ! For tranfer lengths not supported by ESSL, using generic
       do k=1,ph%xsz(3)
          do j=1,ph%xsz(2)
             do i=1,ph%xsz(1)
                buf(i) = cmplx(input(i,j,k),0._mytype, kind=mytype)
             end do
             call spcfft(buf,ph%xsz(1),-1,scratch)
             do i=1,sp%xsz(1)
                output(i,j,k) = buf(i)
             end do
          end do
       end do
    end if

    return

  end subroutine r2c_1m_x


  ! r2c transform, multiple 1D FFTs in z direction
  subroutine r2c_1m_z(input, output)

    implicit none

    real(mytype), dimension(:,:,:), intent(IN)  ::  input
    complex(mytype), dimension(:,:,:), intent(OUT) :: output

    integer :: i,j,k

    ! ESSL doesn't support non-stride-1 r2c transform
    ! use the 'generic' implementation here
    
    do j=1,ph%zsz(2)
       do i=1,ph%zsz(1)
          do k=1,ph%zsz(3)
             buf(k) = cmplx(input(i,j,k),0._mytype, kind=mytype)
          end do
          call spcfft(buf,ph%zsz(3),-1,scratch)
          do k=1,sp%zsz(3)
             output(i,j,k) = buf(k)
          end do
       end do
    end do

    return

  end subroutine r2c_1m_z


  ! c2r transform, multiple 1D FFTs in x direction
  subroutine c2r_1m_x(input, output)

    implicit none

    complex(mytype), dimension(:,:,:), intent(IN)  ::  input
    real(mytype), dimension(:,:,:), intent(OUT) :: output

    integer :: i,j,k

    if (x_goodsize) then
#ifdef DOUBLE_PREC
       call dcrft(1,input,sp%xsz(1),output,ph%xsz(1), &
            ph%xsz(1),ph%xsz(2)*ph%xsz(3),-1,scale, &
            aux1,naux1,aux2,naux2,aux3,naux3)
       call dcrft(0,input,sp%xsz(1),output,ph%xsz(1), &
            ph%xsz(1),ph%xsz(2)*ph%xsz(3),-1,scale, &
            aux1,naux1,aux2,naux2,aux3,naux3)
#else
       call scrft(1,input,sp%xsz(1),output,ph%xsz(1), &
            ph%xsz(1),ph%xsz(2)*ph%xsz(3),-1,scale, &
            aux1,naux1,aux2,naux2,aux3,naux3)
       call scrft(0,input,sp%xsz(1),output,ph%xsz(1), &
            ph%xsz(1),ph%xsz(2)*ph%xsz(3),-1,scale, &
            aux1,naux1,aux2,naux2,aux3,naux3)
#endif
    else
       ! For tranfer lengths not supported by ESSL, using generic
       do k=1,ph%xsz(3)
          do j=1,ph%xsz(2)
             do i=1,ph%xsz(1)/2+1
                buf(i) = input(i,j,k)
             end do
             do i=ph%xsz(1)/2+2,ph%xsz(1)
                buf(i) =  conjg(buf(ph%xsz(1)+2-i))
             end do
             call spcfft(buf,ph%xsz(1),1,scratch)
             do i=1,ph%xsz(1)
                output(i,j,k) = real(buf(i), kind=mytype)
             end do
          end do
       end do
    end if

    return

  end subroutine c2r_1m_x


  ! c2r transform, multiple 1D FFTs in z direction
  subroutine c2r_1m_z(input, output)

    implicit none

    complex(mytype), dimension(:,:,:), intent(IN)  ::  input
    real(mytype), dimension(:,:,:), intent(OUT) :: output

    integer :: i,j,k

    ! ESSL doesn't support non-stride-1 c2r transform
    ! use the 'generic' implementation here
    do j=1,ph%zsz(2)
       do i=1,ph%zsz(1)
          do k=1,ph%zsz(3)/2+1
             buf(k) = input(i,j,k)
          end do
          do k=ph%zsz(3)/2+2,ph%zsz(3)
             buf(k) =  conjg(buf(ph%zsz(3)+2-k))
          end do
          call spcfft(buf,ph%zsz(3),1,scratch)
          do k=1,ph%zsz(3)
             output(i,j,k) = real(buf(k), kind=mytype)
          end do
       end do
    end do

    return

  end subroutine c2r_1m_z


#include "fft_common_3d.f90"


  ! Acceptable length for transforms - read ESSL document for details
  subroutine get_acceptable_length(values)

    implicit none

    integer, dimension(MAGIC), intent(OUT) :: values
    integer, parameter :: nh=25,ni=2,nj=1,nk=1,nm=1
    integer, parameter :: nmax=37748736
    integer, parameter :: dbl=kind(0.0D0)
    
    integer :: h,i,j,k,m, counter
    double precision :: tmp
    
    counter = 1
    do h=1,nh
       do i=0,ni
          do j=0,nj
             do k=0,nk
                do m=0,nm
                   ! use a double precision number to avoid overflow
                   tmp = real(2**h, kind=dbl)*3**i*5**j*7**k*11**m
                   if (tmp>real(nmax+1, kind=dbl)) then
                      exit
                   else
                      values(counter) = 2**h*3**i*5**j*7**k*11**m
                      counter = counter + 1
                   end if
                end do
             end do
          end do
       end do
    end do
    
    call hpsort(counter-1,values) ! sort acceptable lengths

    return
  end subroutine get_acceptable_length


  ! Sorts an array in ascending order by the Heapsort method
  ! The Heapsort method is a N Log2 N algorithm
  subroutine hpsort(N,A)
    
    implicit none
    
    integer, intent(IN) :: N
    integer, dimension(N), intent(INOUT) :: A

    integer :: L,IR,RA,I,J
    
    L=N/2+1
    IR=N
    
    do
       if (L>1) then
          L=L-1
          RA=A(L)
       else
          RA=A(IR)
          A(IR)=A(1)
          IR=IR-1
          if (IR==1) then
             A(1)=RA
             exit
          end if
       end if
       I=L
       J=L+L
       
       do
          if (J>IR) exit
          if (J<IR) then
             if (A(J)<A(J+1)) J=J+1
          end if
          if (RA<A(J)) then
             A(I)=A(J)
             I=J 
             J=J+J
          else
             J=IR+1
          end if
       end do
       
       A(I)=RA
    end do
    
    return
  end subroutine hpsort
  
  ! binary serach
  function binarySearch(a, value)
    integer                     :: binarySearch
    integer, intent(in), target :: a(:)
    integer, intent(in)         :: value
    integer, pointer            :: p(:)
    integer                     :: mid, offset
    
    p => a
    binarySearch = 0
    offset = 0
    do while (size(p) > 0)
       mid = size(p)/2 + 1
       if (p(mid) > value) then
          p => p(:mid-1)
       else if (p(mid) < value) then
          offset = offset + mid
          p => p(mid+1:)
       else
          binarySearch = offset + mid    ! SUCCESS!!
          return
       end if
    end do
  end function binarySearch
  

  ! work array sizes for SCFTD/DCFTD
  subroutine cftd_wk(naux_1, naux_2)

    implicit none
    
    integer, intent(OUT) :: naux_1, naux_2
    integer :: maxdim, mindim

    ! enough work space for the largest dimension
    maxdim = max(ph%xsz(1), ph%ysz(2))
    maxdim = max(maxdim, ph%zsz(3))
    mindim = min(ph%xsz(1), ph%ysz(2))
    mindim = min(mindim, ph%zsz(3))

#ifdef DOUBLE_PREC
    if (maxdim>1024) then
       naux_1 = 60000 + int(28.24*maxdim)
    else
       naux_1 = 30000
    end if
    if (maxdim<252) then
       naux_2 = 20000
    else
       naux_2 = 20000 + (2*maxdim+256)*int(min(64,mindim)+17.12)
    end if
#else
    if (maxdim>2048) then
       naux_1 = 60000 + int(14.12*maxdim)
    else
       naux_1 = 30000
    end if
    if (maxdim<252) then
       naux_2 = 20000
    else
       naux_2 = 20000 + (maxdim+256)*int(min(64,mindim)+8.56)
    end if
#endif

    return
  end subroutine cftd_wk

  ! work array sizes for SCFT/DCFT
  subroutine cft_wk(naux_1, naux_2)

    implicit none
    
    integer, intent(OUT) :: naux_1, naux_2
    integer :: maxdim, tmp

    ! enough work space for the largest dimension
    maxdim = max(ph%xsz(1), ph%ysz(2))
    maxdim = max(maxdim, ph%zsz(3))

#ifdef DOUBLE_PREC
    if (maxdim>2048) then
       naux_1 = 20000 + int(2.28*maxdim)
    else
       naux_1 = 20000
    end if
    ! naux2 also dependent on input stride, treat each directions
    ! X (always stride 1)
    if (ph%xsz(1)>2048) then
       naux_2 = 20000 + int(2.28*ph%xsz(1))
    else
       naux_2 = 20000
    end if
    ! Y
    if (ph%ysz(2)>2048) then
       tmp = 20000 + int(2.28*ph%ysz(2))
    else
       tmp = 20000
    end if
    if (ph%ysz(2)>=252) tmp = tmp + (2*ph%ysz(2)+256)*min(64,ph%ysz(1))
    naux_2 = max(naux_2, tmp)
    ! Z
    if (ph%zsz(3)>2048) then
       tmp = 20000 + int(2.28*ph%zsz(3))
    else
       tmp = 20000
    end if
    if (ph%zsz(3)>=252) tmp = tmp &
         + (2*ph%zsz(3)+256)*min(64,ph%zsz(1)*ph%zsz(2))
    naux_2 = max(naux_2, tmp)
#else
    if (maxdim>8192) then
       naux_1 = 20000 + int(1.14*maxdim)
    else
       naux_1 = 20000
    end if
    ! naux2 also dependent on input stride, treat each directions
    ! X (always stride 1)
    if (ph%xsz(1)>8192) then
       naux_2 = 20000 + int(1.14*ph%xsz(1))
    else
       naux_2 = 20000
    end if
    ! Y
    if (ph%ysz(2)>8192) then
       tmp = 20000 + int(1.14*ph%ysz(2))
    else
       tmp = 20000
    end if
    if (ph%ysz(2)>=252) tmp = tmp + (ph%ysz(2)+256)*min(64,ph%ysz(1))
    naux_2 = max(naux_2, tmp)
    ! Z
    if (ph%zsz(3)>8192) then
       tmp = 20000 + int(1.14*ph%zsz(3))
    else
       tmp = 20000
    end if
    if (ph%zsz(3)>=252) tmp = tmp &
         + (ph%zsz(3)+256)*min(64,ph%zsz(1)*ph%zsz(2))
    naux_2 = max(naux_2, tmp)
#endif

    ! for r2c/c2r the complex arrays are smaller so enough work space

    return
  end subroutine cft_wk

  ! work array sizes for SRCFT/DRCFT
  subroutine rcft_wk(naux_1, naux_2)

    implicit none
    
    integer, intent(OUT) :: naux_1, naux_2

    ! only do r2c for 1D FFTs in X
#ifdef DOUBLE_PREC
    if (ph%xsz(1)>4096) then
       naux_1 = 20000 + int(1.64*ph%xsz(1))
    else
       naux_1 = 22000
    end if
    if (ph%xsz(1)>4096) then
       naux_2 = 20000 + int(1.14*ph%xsz(1))
    else
       naux_1 = 20000
    end if
#else    
    if (ph%xsz(1)>16384) then
       naux_1 = 20000 + int(0.82*ph%xsz(1))
    else
       naux_1 = 25000
    end if
    if (ph%xsz(1)>16384) then
       naux_2 = 20000 + int(0.57*ph%xsz(1))
    else
       naux_1 = 20000
    end if
#endif
    return
  end subroutine rcft_wk


end module decomp_2d_fft
