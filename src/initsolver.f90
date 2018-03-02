module mod_initsolver
  use iso_c_binding, only: C_PTR
  use mod_common_mpi, only: coord
  use mod_fft       , only: fftini
  use mod_param     , only: pi,dims
  implicit none
  private
  public initsolver
  contains
  subroutine initsolver(n,dli,dzci,dzfi,bc,lambdaxy,c_or_f,a,b,c,arrplan,normfft)
    implicit none
    integer, intent(in), dimension(3) :: n
    real(8), intent(in), dimension(3) :: dli
    real(8), intent(in), dimension(0:) :: dzci,dzfi
    character(len=1), intent(in), dimension(0:1,3) :: bc
    real(8), intent(out), dimension(n(1),n(2)) :: lambdaxy
    character(len=1), intent(in), dimension(3) :: c_or_f
    real(8), intent(out), dimension(n(3)) :: a,b,c
    type(C_PTR), intent(out), dimension(2,2) :: arrplan
    real(8), intent(out) :: normfft
    integer :: i,j
    real(8), dimension(n(1)*dims(1))      :: lambdax
    real(8), dimension(n(2)*dims(2))      :: lambday
    integer, dimension(3) :: ng
    integer :: ii,jj
    !
    ! generating eigenvalues consistent with the BCs
    !
    ! note, for periodic this is redundant and the tri-diagonal solver could be optimized, lets leave this for later
    ! call eigenvalues(bctype)
    !
    ! add eigenvalues
    !
    ng(:) = n(:)
    ng(1:2) = ng(1:2)*dims(1:2)
    call eigenvalues(ng(1),bc(:,1),c_or_f(1),lambdax)
    lambdax(:) = lambdax(:)*dli(1)**2
    call eigenvalues(ng(2),bc(:,2),c_or_f(2),lambday)
    lambday(:) = lambday(:)*dli(2)**2
    do j=1,n(2)
      jj = coord(2)*n(2)+j
      do i=1,n(1)
        ii = coord(1)*n(1)+i
        lambdaxy(i,j) = lambdax(ii)+lambday(jj)
      enddo
    enddo
    !
    ! compute coefficients for tridiagonal solver
    !
    call tridmatrix(bc(:,3),n(3),dli(3),dzci,dzfi,c_or_f(3),a,b,c)
    !
    ! prepare ffts
    !
    call fftini(ng(1),ng(2),ng(3),bc(:,1:2),c_or_f(1:2),arrplan,normfft)
    return
  end subroutine initsolver
  !
  subroutine eigenvalues(n,bc,c_or_f,lambda)
    implicit none
    integer, intent(in ) :: n
    character(len=1), intent(in), dimension(0:1) :: bc
    character(len=1), intent(in) :: c_or_f ! c -> cell-centered; f-face-centered
    real(8), intent(out), dimension(n) :: lambda
    integer :: l 
    select case(bc(0)//bc(1))
    case('PP')
      l = 1
      lambda(l)     = -4.*sin((1.*(l-1))*pi/(1.*n))**2
      do l=2,n/2
        lambda(l  )   = -4.*sin((1.*(l-1))*pi/(1.*n))**2
        lambda(n-l+2) = lambda(l) ! according to the half-complex format of fftw
      enddo
      l = n/2+1
      lambda(l)    = -4.*sin((1.*(l-1))*pi/(1.*n))**2
    case('NN')
      if(    c_or_f.eq.'c') then
        do l=1,n
          lambda(l)   = -4.*sin((1.*(l-1))*pi/(2.*n))**2
        enddo
      elseif(c_or_f.eq.'f') then
        do l=1,n
          lambda(l)   = -4.*sin((1.*(l-1))*pi/(2.*(n-1)))**2
        enddo
      endif
    case('DD')
      if(    c_or_f.eq.'c') then
        do l=1,n
          lambda(l)   = -4.*sin((1.*(l-0))*pi/(2.*n))**2
        enddo
      elseif(c_or_f.eq.'f') then
        do l=1,n
          lambda(l)   = -4.*sin((1.*(l-0))*pi/(2.*(n+1)))**2
        enddo
      endif
    case('ND')
      do l=1,n
        lambda(l)   = -4.*sin((1.*(2*l-1))*pi/(4.*n))**2
      enddo
    end select   
    return
  end subroutine eigenvalues
  !
  subroutine tridmatrix(bc,n,dzi,dzci,dzfi,c_or_f,a,b,c)
    implicit none
    real(8), parameter :: eps = 1.e-10
    character(len=1), intent(in), dimension(0:1) :: bc
    integer, intent(in) :: n
    real(8), intent(in) :: dzi
    real(8), intent(in), dimension(0:) :: dzci,dzfi
    character(len=1), intent(in) :: c_or_f ! c -> cell-centered; f-face-centered
    real(8), intent(out), dimension(n) :: a,b,c
    integer :: k
    real(8) :: factor1,factor2
    ! bc(l) =  1 -> Neumann   BC in the lth direction
    ! bc(l) =  0 -> Periodic  BC in the lth direction
    ! bc(l) = -1 -> Dirichlet BC in the lth direction
    select case(c_or_f)
    case('c')
      do k=1,n
        a(k) = dzfi(k)*dzci(k-1)
        c(k) = dzfi(k)*dzci(k)
      enddo
    case('f')
      do k = 1,n-1 ! needs to be changed; factor should go to the RHS for solving the system of eqs right!
        a(k) = dzfi(k)*dzci(k)
        c(k) = dzfi(k+1)*dzci(k)
      enddo
    end select
      b(:) = -(a(:)+c(:))
    select case(bc(-))
    case('P')
      factor1 = 0.
    case('D')
      factor1 = -1.
    case('N')
      factor1 = 1.
    end select
    select case(bc(1))
    case('P')
      factor2 = 0.
    case('D')
      factor2 = -1.
    case('N')
      factor2 = 1.
    end select
    select case(c_or_f)
    case('c')
      b(1) = b(1) + factor1*a(1)
      b(n) = b(n) + factor2*c(n)
    case('f')
      ! needs to be changed; factor should go to the RHS for solving the system of eqs right!
    end select
    a(:) = a(:) + eps
    b(:) = b(:) + eps
    c(:) = c(:) + eps
    return
  end subroutine tridmatrix
end module mod_initsolver
