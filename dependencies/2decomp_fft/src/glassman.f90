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

! This module contains a few 'generic' FFT routines, making the 
! 2DECOMP&FFT library not dependent on any external libraries

module glassman

  use decomp_2d, only : mytype

  implicit none

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Following is a FFT implementation based on algorithm proposed by
  ! Glassman, a general FFT algorithm supporting arbitrary input length.
  !
  ! W. E. Ferguson, Jr., "A simple derivation of Glassman general-n fast
  ! Fourier transform," Comput. and Math. with Appls., vol. 8, no. 6, pp.
  ! 401-411, 1982.
  !
  ! Original implemtation online at http://www.jjj.de/fft/fftpage.html
  !
  ! Updated  
  !  -  to handle double-precision as well
  !  -  unnecessary scaling code removed
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE SPCFFT(U,N,ISIGN,WORK)
    
    IMPLICIT NONE
    
    LOGICAL :: INU
    INTEGER :: A,B,C,N,I,ISIGN
    COMPLEX(mytype) :: U(*),WORK(*)
    
    A = 1
    B = N
    C = 1
    INU = .TRUE.
    
    DO WHILE ( B .GT. 1 )
       A = C * A
       C = 2
       DO WHILE ( MOD(B,C) .NE. 0 )
          C = C + 1
       END DO
       B = B / C
       IF ( INU ) THEN
          CALL SPCPFT (A,B,C,U,WORK,ISIGN)
       ELSE
          CALL SPCPFT (A,B,C,WORK,U,ISIGN)
       END IF
       INU = ( .NOT. INU )
    END DO
    
    IF ( .NOT. INU ) THEN
       DO I = 1, N
          U(I) = WORK(I)
       END DO
    END IF
    
    RETURN
  END SUBROUTINE SPCFFT
  
  
  SUBROUTINE SPCPFT( A, B, C, UIN, UOUT, ISIGN )
    
    IMPLICIT NONE
    
    INTEGER :: ISIGN,A,B,C,IA,IB,IC,JCR,JC
    
    DOUBLE PRECISION :: ANGLE
    
    COMPLEX(mytype) :: UIN(B,C,A),UOUT(B,A,C),DELTA,OMEGA,SUM
    
    ANGLE = 8.D0*DATAN(1.D0) / REAL( A * C, kind=mytype )
    OMEGA = CMPLX( 1.0, 0.0, kind=mytype )
    
    IF( ISIGN .EQ. 1 ) THEN
       DELTA = CMPLX( DCOS(ANGLE), DSIN(ANGLE), kind=mytype )
    ELSE
       DELTA = CMPLX( DCOS(ANGLE), -DSIN(ANGLE), kind=mytype )
    END IF
    
    DO IC = 1, C
       DO IA = 1, A
          DO IB = 1, B
             SUM = UIN( IB, C, IA )
             DO JCR = 2, C
                JC = C + 1 - JCR
                SUM = UIN( IB, JC, IA ) + OMEGA * SUM
             END DO
             UOUT( IB, IA, IC ) = SUM
          END DO
          OMEGA = DELTA * OMEGA
       END DO
    END DO
    
    RETURN
  END SUBROUTINE SPCPFT


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! A 3D real-to-complex routine implemented using the 1D FFT above
  !   Input:   nx*ny*nz real numbers
  !   Output:  (nx/2+1)*ny*nz complex numbers
  ! Just like big FFT libraries (such as FFTW) do
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine glassman_3d_r2c(in_r,nx,ny,nz,out_c)

    implicit none
    
    integer, intent(IN) :: nx,ny,nz
    real(mytype), dimension(nx,ny,nz) :: in_r
    complex(mytype), dimension(nx/2+1,ny,nz) :: out_c
    
    complex(mytype), allocatable, dimension(:) :: buf, scratch
    integer :: maxsize, i,j,k
    
    maxsize = max(nx, max(ny,nz))
    allocate(buf(maxsize))
    allocate(scratch(maxsize))
    
    ! ===== 1D FFTs in X =====
    do k=1,nz
       do j=1,ny
          ! Glassman's 1D FFT is c2c only, 
          ! needing some pre- and post-processing for r2c
          ! pack real input in complex storage
          do i=1,nx
             buf(i) = cmplx(in_r(i,j,k),0._mytype, kind=mytype)
          end do
          call spcfft(buf,nx,-1,scratch)
          ! simply drop the redundant part of the complex output
          do i=1,nx/2+1
             out_c(i,j,k) = buf(i)
          end do
       end do
    end do
    
    ! ===== 1D FFTs in Y =====
    do k=1,nz
       do i=1,nx/2+1
          do j=1,ny
             buf(j) = out_c(i,j,k)
          end do
          call spcfft(buf,ny,-1,scratch)
          do j=1,ny
             out_c(i,j,k) = buf(j)
          end do
       end do
    end do
    
    ! ===== 1D FFTs in Z =====
    do j=1,ny
       do i=1,nx/2+1
          do k=1,nz
             buf(k) = out_c(i,j,k)
          end do
          call spcfft(buf,nz,-1,scratch)
          do k=1,nz
             out_c(i,j,k) = buf(k)
          end do
       end do
    end do
    
    deallocate(buf,scratch)
    
    return
  end subroutine glassman_3d_r2c


end module glassman

