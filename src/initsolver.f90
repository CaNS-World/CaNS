module mod_initsolver
  use iso_c_binding, only: C_PTR
  use decomp_2d    , only: xsize,ysize
  use mod_fft      , only: fftini
  use mod_param    , only: pi
  use mod_types
  implicit none
  private
  public initsolver
  interface eigenvalues
    module procedure eigenvalues_sp,eigenvalues_dp
  end interface
  interface initsolver
    module procedure initsolver_sp,initsolver_dp
  end interface
  interface tridmatrix
    module procedure tridmatrix_sp,tridmatrix_dp
  end interface
  contains
  !
#define MYREAL real(sp)
  subroutine initsolver_sp(ng,lo_z,hi_z,dli,dzci,dzfi,cbc,bc,lambdaxy,c_or_f,a,b,c,arrplan,normfft, &
                           rhsbx,rhsby,rhsbz)
#include "initsolver_initsolver-inc.f90"
  end subroutine initsolver_sp
#undef MYREAL
#define MYREAL real(dp)
  subroutine initsolver_dp(ng,lo_z,hi_z,dli,dzci,dzfi,cbc,bc,lambdaxy,c_or_f,a,b,c,arrplan,normfft, &
                           rhsbx,rhsby,rhsbz)
#include "initsolver_initsolver-inc.f90"
  end subroutine initsolver_dp
#undef MYREAL
  !
#define MYREAL real(sp)
  subroutine eigenvalues_sp(n,bc,c_or_f,lambda)
#include "initsolver_eigenvalues-inc.f90"
  end subroutine eigenvalues_sp
#undef MYREAL
#define MYREAL real(dp)
  subroutine eigenvalues_dp(n,bc,c_or_f,lambda)
#include "initsolver_eigenvalues-inc.f90"
  end subroutine eigenvalues_dp
#undef MYREAL
  !
#define MYREAL real(sp)
  subroutine tridmatrix_sp(bc,n,dzi,dzci,dzfi,c_or_f,a,b,c)
#include "initsolver_tridmatrix-inc.f90"
  end subroutine tridmatrix_sp
#undef MYREAL
#define MYREAL real(dp)
  subroutine tridmatrix_dp(bc,n,dzi,dzci,dzfi,c_or_f,a,b,c)
#include "initsolver_tridmatrix-inc.f90"
  end subroutine tridmatrix_dp
#undef MYREAL
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
