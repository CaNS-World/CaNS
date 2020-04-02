module mod_initsolver
  use decomp_2d
  use iso_c_binding, only: C_PTR
  use mod_common_mpi, only: dims,n_x,n_y,n_z,ijk_start_z
  use mod_fft       , only: fftini
  use mod_param     , only: pi
  use mod_types
  implicit none
  private
  public initsolver
  contains
  subroutine initsolver(n,dli,dzci,dzfi,cbc,bc,lambdaxy,c_or_f,a,b,c,arrplan,normfft,rhsbx,rhsby,rhsbz)
    !
    ! initializes the Poisson/Helmholtz solver
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3) :: dli
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    real(rp)         , intent(in), dimension(0:1,3) :: bc
    real(rp), intent(out), dimension(:,:) :: lambdaxy
    character(len=1), intent(in), dimension(3) :: c_or_f
    real(rp), intent(out), dimension(n_z(3)) :: a,b,c
    type(C_PTR), intent(out), dimension(2,2) :: arrplan
    real(rp), intent(out), dimension(n(2),n(3),0:1) :: rhsbx
    real(rp), intent(out), dimension(n(1),n(3),0:1) :: rhsby
    real(rp), intent(out), dimension(n(1),n(2),0:1) :: rhsbz
    real(rp), intent(out) :: normfft
    real(rp), dimension(3)        :: dl
    real(rp), dimension(0:n(3)*dims(3)+1) :: dzc,dzf
    integer :: i,j
    real(rp), dimension(n(1)*dims(1))      :: lambdax
    real(rp), dimension(n(2)*dims(2))      :: lambday
    integer, dimension(3) :: ng
    integer :: ii,jj
    !
    ! generating eigenvalues consistent with the BCs
    !
    ng(:) = n(:)*dims(:)
    call eigenvalues(ng(1),cbc(:,1),c_or_f(1),lambdax)
    lambdax(:) = lambdax(:)*dli(1)**2
    call eigenvalues(ng(2),cbc(:,2),c_or_f(2),lambday)
    lambday(:) = lambday(:)*dli(2)**2
    !
    ! add eigenvalues
    !
    do j=1,n_z(2)
      jj = ijk_start_z(2) + j
      do i=1,n_z(1)
        ii = ijk_start_z(1) + i
        lambdaxy(i,j) = lambdax(ii)+lambday(jj)
      enddo
    enddo
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
    call bc_rhs(cbc(:,1),n(1),bc(:,1),(/dl(1) ,dl(1)    /),(/dl(1) ,dl(1)    /),c_or_f(1),rhsbx)
    call bc_rhs(cbc(:,2),n(2),bc(:,2),(/dl(2) ,dl(2)    /),(/dl(2) ,dl(2)    /),c_or_f(2),rhsby)
    if(    c_or_f(3).eq.'c') then
      call bc_rhs(cbc(:,3),n(3),bc(:,3),(/dzc(0),dzc(ng(3)  )/),(/dzf(1),dzf(ng(3))/),c_or_f(3),rhsbz)
    elseif(c_or_f(3).eq.'f') then
      call bc_rhs(cbc(:,3),n(3),bc(:,3),(/dzc(1),dzc(ng(3)-1)/),(/dzf(1),dzf(ng(3))/),c_or_f(3),rhsbz)
    endif
    !
    ! prepare ffts
    !
    call fftini(n_x,n_y,cbc(:,1:2),c_or_f(1:2),arrplan,normfft)
    return
  end subroutine initsolver
  !
  subroutine eigenvalues(n,bc,c_or_f,lambda)
    implicit none
    integer , intent(in ) :: n
    character(len=1), intent(in), dimension(0:1) :: bc
    character(len=1), intent(in) :: c_or_f ! c -> cell-centered; f -> face-centered
    real(rp), intent(out), dimension(n) :: lambda
    integer :: l 
    select case(bc(0)//bc(1))
    case('PP')
      do l=1,n
        lambda(l  )   = -4.*sin((1.*(l-1))*pi/(1.*n))**2
      enddo
    case('NN')
      if(    c_or_f.eq.'c') then
        do l=1,n
          lambda(l)   = -4.*sin((1.*(l-1))*pi/(2.*n))**2
        enddo
      elseif(c_or_f.eq.'f') then
        do l=1,n
          lambda(l)   = -4.*sin((1.*(l-1))*pi/(2.*(n-1+1)))**2
        enddo
      endif
    case('DD')
      if(    c_or_f.eq.'c') then
        do l=1,n
          lambda(l)   = -4.*sin((1.*(l-0))*pi/(2.*n))**2
        enddo
      elseif(c_or_f.eq.'f') then
        do l=1,n-1 ! point at n is a boundary and is excluded here
          lambda(l)   = -4.*sin((1.*(l-0))*pi/(2.*(n+1-1)))**2
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
    real(rp), parameter :: eps = epsilon(1._rp)*10.
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
      enddo
    case('f')
      do k = 1,n
        a(k) = dzfi(k)*dzci(k)
        c(k) = dzfi(k+1)*dzci(k)
      enddo
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
    enddo
    select case(c_or_f)
    case('c')
      b(1) = b(1) + factor(0)*a(1)
      b(n) = b(n) + factor(1)*c(n)
    case('f')
      if(bc(0).eq.'N') b(1) = b(1) + factor(0)*a(1)
      if(bc(1).eq.'N') b(n) = b(n) + factor(1)*c(n)
    end select
    ! n.b.: a(1) and c(n) not set to zero here;
    !       the values are not used in the solver unless
    !       the direction is periodic
    a(:) = a(:)! + eps
    b(:) = b(:)! + eps
    c(:) = c(:)! + eps
    return
  end subroutine tridmatrix
  !
  subroutine bc_rhs(cbc,n,bc,dlc,dlf,c_or_f,rhs)
    implicit none
    character(len=1), intent(in), dimension(0:1) :: cbc
    integer , intent(in) :: n
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
          if(ibound.eq.0) sgn =  1.
          if(ibound.eq.1) sgn = -1.
          factor(ibound) = sgn*dlc(ibound)*bc(ibound)
        end select
      enddo
    case('f')
      do ibound = 0,1
        select case(cbc(ibound))
        case('P')
          factor(ibound) = 0.
        case('D')
          factor(ibound) = -bc(ibound)
        case('N')
          if(ibound.eq.0) sgn =  1.
          if(ibound.eq.1) sgn = -1.
          factor(ibound) = sgn*dlf(ibound)*bc(ibound)
        end select
      enddo
    end select
    forall(ibound=0:1)
      rhs(:,:,ibound) = factor(ibound)/dlc(ibound)/dlf(ibound)
      rhs(:,:,ibound) = factor(ibound)/dlc(ibound)/dlf(ibound)
    end forall
    return
  end subroutine bc_rhs
end module mod_initsolver
