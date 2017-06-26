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

! A few utility routines to find factors of integer numbers

  subroutine findfactor(num, factors, nfact)
    
    implicit none

    integer, intent(IN) :: num
    integer, intent(OUT), dimension(*) :: factors
    integer, intent(OUT) :: nfact
    integer :: i, m

    ! find the factors <= sqrt(num)
    m = int(sqrt(real(num)))
    nfact = 1
    do i=1,m
       if (num/i*i == num) then
          factors(nfact) = i
          nfact = nfact + 1
       end if
    end do
    nfact = nfact - 1

    ! derive those > sqrt(num)
    if (factors(nfact)**2/=num) then
       do i=nfact+1, 2*nfact
          factors(i) = num / factors(2*nfact-i+1)
       end do
       nfact = nfact * 2
    else
       do i=nfact+1, 2*nfact-1
          factors(i) = num / factors(2*nfact-i)
       end do
       nfact = nfact * 2 - 1
    endif
       
    return

  end subroutine findfactor


  subroutine primefactors(num, factors, nfact)

    implicit none
  
    integer, intent(IN) :: num
    integer, intent(OUT), dimension(*) :: factors
    integer, intent(INOUT) :: nfact

    integer :: i, n
    
    i = 2  
    nfact = 1
    n = num 
    do
       if (mod(n,i) == 0) then
          factors(nfact) = i
          nfact = nfact + 1
          n = n / i
       else
          i = i + 1
       end if
       if (n == 1) then
          nfact = nfact - 1
          exit
       end if
    end do
    
    return

  end subroutine primefactors
  
