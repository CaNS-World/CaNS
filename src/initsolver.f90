! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_initsolver
  use, intrinsic :: iso_c_binding, only: C_PTR
  use mod_fft  , only: fftini
  use mod_types
  implicit none
  private
  public initsolver
  contains
  subroutine initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci,dzfi,cbc,bc,lambdaxy,c_or_f,a,b,c,arrplan,normfft, &
                        rhsbx,rhsby,rhsbz)
    !
    ! initializes the Poisson/Helmholtz solver
    !
    implicit none
    integer , intent(in), dimension(3) :: ng,n_x_fft,n_y_fft,lo_z,hi_z
    real(rp), intent(in), dimension(3 ) :: dli
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    real(rp)        , intent(in), dimension(0:1,3) :: bc
    real(rp), intent(out), dimension(lo_z(1):,lo_z(2):) :: lambdaxy
    character(len=1), intent(in), dimension(3) :: c_or_f
    real(rp), intent(out), dimension(:) :: a,b,c
#if !defined(_OPENACC)
    type(C_PTR), intent(out), dimension(2,2) :: arrplan
#else
    integer    , intent(out), dimension(2,2) :: arrplan
#endif
    real(rp), intent(out), dimension(:,:,0:) :: rhsbx
    real(rp), intent(out), dimension(:,:,0:) :: rhsby
    real(rp), intent(out), dimension(:,:,0:) :: rhsbz
    real(rp), intent(out) :: normfft
    real(rp), dimension(3)        :: dl
    real(rp), dimension(0:ng(3)+1) :: dzc,dzf
    integer :: i,j
    real(rp), dimension(ng(1))      :: lambdax
    real(rp), dimension(ng(2))      :: lambday
    !
    ! generating eigenvalues consistent with the BCs
    !
    call eigenvalues(ng(1),cbc(:,1),c_or_f(1),lambdax)
    lambdax(:) = lambdax(:)*dli(1)**2
    call eigenvalues(ng(2),cbc(:,2),c_or_f(2),lambday)
    lambday(:) = lambday(:)*dli(2)**2
    !
    ! add eigenvalues
    !
    do j=lo_z(2),hi_z(2)
      do i=lo_z(1),hi_z(1)
        lambdaxy(i,j) = lambdax(i)+lambday(j)
      end do
    end do
    !
    ! compute coefficients for tridiagonal solver
    !
    call tridmatrix(cbc(:,3),ng(3),dli(3),dzci,dzfi,c_or_f(3),a,b,c)
    !
    ! compute values to be added to the right hand side
    !
    dl(:)  = dli( :)**(-1)
    dzc(:) = dzci(:)**(-1)
    dzf(:) = dzfi(:)**(-1)
    call bc_rhs(cbc(:,1),bc(:,1),[dl(1) ,dl(1)      ],[dl(1) ,dl(1)    ],c_or_f(1),rhsbx)
    call bc_rhs(cbc(:,2),bc(:,2),[dl(2) ,dl(2)      ],[dl(2) ,dl(2)    ],c_or_f(2),rhsby)
    if(     c_or_f(3) == 'c') then
      call bc_rhs(cbc(:,3),bc(:,3),[dzc(0),dzc(ng(3)  )],[dzf(1),dzf(ng(3))],c_or_f(3),rhsbz)
    else if(c_or_f(3) == 'f') then
      call bc_rhs(cbc(:,3),bc(:,3),[dzc(1),dzc(ng(3)-1)],[dzf(1),dzf(ng(3))],c_or_f(3),rhsbz)
    end if
    !
    ! prepare ffts
    !
    call fftini(ng,n_x_fft,n_y_fft,cbc(:,1:2),c_or_f(1:2),arrplan,normfft)
  end subroutine initsolver
  !
  subroutine eigenvalues(n,bc,c_or_f,lambda)
    use mod_param, only: pi
    implicit none
    integer , intent(in ) :: n
    character(len=1), intent(in), dimension(0:1) :: bc
    character(len=1), intent(in) :: c_or_f ! c -> cell-centered; f -> face-centered
    real(rp), intent(out), dimension(n) :: lambda
    integer :: l
    select case(bc(0)//bc(1))
    case('PP')
      do l=1,n
        lambda(l  )   = -2.*(1.-cos((2*(l-1))*pi/(1.*n)))
      end do
#if defined(_OPENACC)
      block
        !
        ! new format: (r[0],r[n],r[1],i[1],...,r[n-1],i[n-1])
        ! note that i[0] = i[n] = 0 in a R2C DFT
        !
        integer :: nh,iswap(n)
        nh = (n+1)/2
        iswap(1) = 1
        iswap(2) = nh+(1-mod(n,2))
        do l=2,n-1
          if(l <= nh) then ! real eigenvalue
            iswap(2*l-1                  ) = l
          else             ! imaginary eigenvalue
            iswap(n-2*(l-(nh+1))-mod(n,2)) = l+1
          end if
        end do
        lambda(:) = lambda(iswap(:))
      end block
#endif
    case('NN')
      if(     c_or_f == 'c') then
        do l=1,n
          lambda(l)   = -2.*(1.-cos((l-1  )*pi/(1.*n)))
        end do
      else if(c_or_f == 'f') then
        do l=1,n
          lambda(l)   = -2.*(1.-cos((l-1  )*pi/(1.*(n-1+1))))
        end do
      end if
    case('DD')
      if(     c_or_f == 'c') then
        do l=1,n
          lambda(l)   = -2.*(1.-cos((l    )*pi/(1.*n)))
        end do
      else if(c_or_f == 'f') then
        do l=1,n-1 ! point at n is a boundary and is excluded here
          lambda(l)   = -2.*(1.-cos((l    )*pi/(1.*(n+1-1))))
        end do
      end if
    case('ND','DN')
      do l=1,n
        lambda(l)     = -2.*(1.-cos((2*l-1)*pi/(2.*n)))
      end do
    end select
  end subroutine eigenvalues
  !
  subroutine tridmatrix(bc,n,dzi,dzci,dzfi,c_or_f,a,b,c)
    implicit none
    character(len=1), intent(in), dimension(0:1) :: bc
    integer , intent(in) :: n
    real(rp), intent(in) :: dzi
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    character(len=1), intent(in) :: c_or_f ! c -> cell-centered; f-face-centered
    real(rp), intent(out), dimension(n) :: a,b,c
    integer :: k
    integer :: ibound
    real(rp), dimension(0:1) :: factor
    select case(c_or_f)
    case('c')
      do k=1,n
        a(k) = dzfi(k)*dzci(k-1)
        c(k) = dzfi(k)*dzci(k)
      end do
    case('f')
      do k = 1,n
        a(k) = dzfi(k)*dzci(k)
        c(k) = dzfi(k+1)*dzci(k)
      end do
    end select
    b(:) = -(a(:)+c(:))
    do ibound = 0,1
      select case(bc(ibound))
      case('P')
        factor(ibound) = 0.
      case('D')
        factor(ibound) = -1.
      case('N')
        factor(ibound) = 1.
      end select
    end do
    select case(c_or_f)
    case('c')
      b(1) = b(1) + factor(0)*a(1)
      b(n) = b(n) + factor(1)*c(n)
    case('f')
      if(bc(0) == 'N') b(1) = b(1) + factor(0)*a(1)
      if(bc(1) == 'N') b(n) = b(n) + factor(1)*c(n)
    end select
  end subroutine tridmatrix
  !
  subroutine bc_rhs(cbc,bc,dlc,dlf,c_or_f,rhs)
    implicit none
    character(len=1), intent(in), dimension(0:1) :: cbc
    real(rp), intent(in), dimension(0:1) :: bc
    real(rp), intent(in), dimension(0:1) :: dlc,dlf
    real(rp), intent(out), dimension(:,:,0:) :: rhs
    character(len=1), intent(in) :: c_or_f ! c -> cell-centered; f -> face-centered
    real(rp), dimension(0:1) :: factor
    real(rp) :: sgn
    integer :: ibound
    !
    select case(c_or_f)
    case('c')
      do ibound = 0,1
        select case(cbc(ibound))
        case('P')
          factor(ibound) = 0.
        case('D')
          factor(ibound) = -2.*bc(ibound)
        case('N')
          if(ibound == 0) sgn =  1.
          if(ibound == 1) sgn = -1.
          factor(ibound) = sgn*dlc(ibound)*bc(ibound)
        end select
      end do
    case('f')
      do ibound = 0,1
        select case(cbc(ibound))
        case('P')
          factor(ibound) = 0.
        case('D')
          factor(ibound) = -bc(ibound)
        case('N')
          if(ibound == 0) sgn =  1.
          if(ibound == 1) sgn = -1.
          factor(ibound) = sgn*dlf(ibound)*bc(ibound)
        end select
      end do
    end select
    do concurrent(ibound=0:1)
      rhs(:,:,ibound) = factor(ibound)/dlc(ibound)/dlf(ibound)
      rhs(:,:,ibound) = factor(ibound)/dlc(ibound)/dlf(ibound)
    end do
  end subroutine bc_rhs
end module mod_initsolver
