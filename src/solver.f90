! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_solver
  use, intrinsic :: iso_c_binding, only: C_PTR
  use decomp_2d
  use mod_fft   , only: fft
  use mod_types
  implicit none
  private
  public solver
#if defined(_IMPDIFF_1D)
  public solver_gaussel_z
#endif
  interface gaussel
    module procedure gaussel_sp,gaussel_dp
  end interface
  interface gaussel_periodic
    module procedure gaussel_periodic_sp,gaussel_periodic_dp
  end interface
  interface dgtsv_homebrewed
    module procedure dgtsv_homebrewed_sp,dgtsv_homebrewed_dp
  end interface
  contains
  subroutine solver(n,ng,arrplan,normfft,lambdaxy,a,b,c,bc,c_or_f,p)
    implicit none
    integer , intent(in), dimension(3) :: n,ng ! ng is here just for compatibility with the argument list of solver_gpu
    type(C_PTR), intent(in), dimension(2,2) :: arrplan
    real(rp), intent(in) :: normfft
    real(gp), intent(in), dimension(:,:) :: lambdaxy
    real(gp), intent(in), dimension(:) :: a,b,c
    character(len=1), dimension(0:1,3), intent(in) :: bc
    character(len=1), intent(in), dimension(3) :: c_or_f
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    real(gp), dimension(xsize(1),xsize(2),xsize(3)) :: px
    real(gp), dimension(ysize(1),ysize(2),ysize(3)) :: py
    real(gp), dimension(zsize(1),zsize(2),zsize(3)) :: pz
    integer :: q
    integer, dimension(3) :: n_z
    !
    n_z(:) = zsize(:)
#if !defined(_DECOMP_Y) && !defined(_DECOMP_Z)
    !$OMP PARALLEL WORKSHARE
    px(:,:,:) = p(1:n(1),1:n(2),1:n(3))
    !$OMP END PARALLEL WORKSHARE
#elif defined(_DECOMP_Y)
    !$OMP PARALLEL WORKSHARE
    py(:,:,:) = p(1:n(1),1:n(2),1:n(3))
    !$OMP END PARALLEL WORKSHARE
    call transpose_y_to_x(py,px)
#elif defined(_DECOMP_Z)
    !$OMP PARALLEL WORKSHARE
    pz(:,:,:) = p(1:n(1),1:n(2),1:n(3))
    !$OMP END PARALLEL WORKSHARE
    !call transpose_z_to_x(pz,px)
    call transpose_z_to_y(pz,py)
    call transpose_y_to_x(py,px)
#endif
    call fft(arrplan(1,1),px) ! fwd transform in x
    !
    call transpose_x_to_y(px,py)
    call fft(arrplan(1,2),py) ! fwd transform in y
    !
    call transpose_y_to_z(py,pz)
    q = 0
    if(c_or_f(3) == 'f'.and.bc(1,3) == 'D') q = 1
    if(bc(0,3)//bc(1,3) == 'PP') then
      call gaussel_periodic(n_z(1),n_z(2),n_z(3)-q,0,a,b,c,pz,lambdaxy)
    else
      call gaussel(         n_z(1),n_z(2),n_z(3)-q,0,a,b,c,pz,lambdaxy)
    end if
    !
    call transpose_z_to_y(pz,py)
    call fft(arrplan(2,2),py) ! bwd transform in y
    !
    call transpose_y_to_x(py,px)
    call fft(arrplan(2,1),px) ! bwd transform in x
    !
#if !defined(_DECOMP_Y) && !defined(_DECOMP_Z)
    !$OMP PARALLEL WORKSHARE
    p(1:n(1),1:n(2),1:n(3)) = px(:,:,:)*normfft
    !$OMP END PARALLEL WORKSHARE
#elif defined(_DECOMP_Y)
    call transpose_x_to_y(px,py)
    !$OMP PARALLEL WORKSHARE
    p(1:n(1),1:n(2),1:n(3)) = py(:,:,:)*normfft
    !$OMP END PARALLEL WORKSHARE
#elif defined(_DECOMP_Z)
    !call transpose_x_to_z(px,pz)
    call transpose_x_to_y(px,py)
    call transpose_y_to_z(py,pz)
    !$OMP PARALLEL WORKSHARE
    p(1:n(1),1:n(2),1:n(3)) = pz(:,:,:)*normfft
    !$OMP END PARALLEL WORKSHARE
#endif
  end subroutine solver
  !
  subroutine gaussel_sp(nx,ny,n,nh,a,b,c,p,lambdaxy)
    use mod_types, only: wp => sp
    use mod_param, only: eps => eps_sp
#include "solver_gaussel-inc.f90"
  end subroutine gaussel_sp
  subroutine gaussel_dp(nx,ny,n,nh,a,b,c,p,lambdaxy)
    use mod_types, only: wp => dp
    use mod_param, only: eps => eps_dp
#include "solver_gaussel-inc.f90"
  end subroutine gaussel_dp
  !
  subroutine gaussel_periodic_sp(nx,ny,n,nh,a,b,c,p,lambdaxy)
    use mod_types, only: wp => sp
    use mod_param, only: eps => eps_sp
#include "solver_gaussel_periodic-inc.f90"
  end subroutine gaussel_periodic_sp
  subroutine gaussel_periodic_dp(nx,ny,n,nh,a,b,c,p,lambdaxy)
    use mod_types, only: wp => dp
    use mod_param, only: eps => eps_dp
#include "solver_gaussel_periodic-inc.f90"
  end subroutine gaussel_periodic_dp
  !
  subroutine dgtsv_homebrewed_sp(n,a,b,c,p)
    use mod_types, only: wp => sp
    use mod_param, only: eps => eps_sp
#include "solver_dgtsv_homebrewed-inc.f90"
  end subroutine dgtsv_homebrewed_sp
  subroutine dgtsv_homebrewed_dp(n,a,b,c,p)
    use mod_types, only: wp => dp
    use mod_param, only: eps => eps_dp
#include "solver_dgtsv_homebrewed-inc.f90"
  end subroutine dgtsv_homebrewed_dp
  !
#if defined(_IMPDIFF_1D)
  subroutine solver_gaussel_z(n,a,b,c,bcz,c_or_f,p)
    use mod_param, only: eps => eps_rp
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(:) :: a,b,c
    character(len=1), dimension(0:1), intent(in) :: bcz
    character(len=1), intent(in), dimension(3) :: c_or_f
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
#if !defined(_DECOMP_Z)
    real(rp), dimension(xsize(1),xsize(2),xsize(3)) :: px
    real(rp), dimension(ysize(1),ysize(2),ysize(3)) :: py
    real(rp), dimension(zsize(1),zsize(2),zsize(3)) :: pz
#endif
    integer :: q
    integer, dimension(3) :: n_z
    !
    n_z(:) = zsize(:)
#if !defined(_DECOMP_Y) && !defined(_DECOMP_Z)
    !$OMP PARALLEL WORKSHARE
    px(:,:,:) = p(1:n(1),1:n(2),1:n(3))
    !$OMP END PARALLEL WORKSHARE
    !call transpose_x_to_z(px,pz)
    call transpose_x_to_y(px,py)
    call transpose_y_to_z(py,pz)
#elif defined(_DECOMP_Y)
    !$OMP PARALLEL WORKSHARE
    py(:,:,:) = p(1:n(1),1:n(2),1:n(3))
    !$OMP END PARALLEL WORKSHARE
    call transpose_y_to_z(py,pz)
#endif
    q = 0
    if(c_or_f(3) == 'f'.and.bcz(1) == 'D') q = 1
#if !defined(_DECOMP_Z)
    if(bcz(0)//bcz(1) == 'PP') then
      call gaussel_periodic(n_z(1),n_z(2),n_z(3)-q,0,a,b,c,pz)
    else
      call gaussel(         n_z(1),n_z(2),n_z(3)-q,0,a,b,c,pz)
    end if
#else
    if(bcz(0)//bcz(1) == 'PP') then
      call gaussel_periodic(n_z(1),n_z(2),n_z(3)-q,1,a,b,c,p)
    else
      call gaussel(         n_z(1),n_z(2),n_z(3)-q,1,a,b,c,p)
    end if
#endif
    !
#if !defined(_DECOMP_Y) && !defined(_DECOMP_Z)
    !call transpose_z_to_x(pz,px)
    call transpose_z_to_y(pz,py)
    call transpose_y_to_x(py,px)
    !$OMP PARALLEL WORKSHARE
    p(1:n(1),1:n(2),1:n(3)) = px(:,:,:)
    !$OMP END PARALLEL WORKSHARE
#elif defined(_DECOMP_Y)
    call transpose_z_to_y(pz,py)
    !$OMP PARALLEL WORKSHARE
    p(1:n(1),1:n(2),1:n(3)) = py(:,:,:)
    !$OMP END PARALLEL WORKSHARE
#endif
  end subroutine solver_gaussel_z
  !
#if 0
  subroutine gaussel_lapack(nx,ny,n,a,b,c,p)
    implicit none
#if !defined(_SINGLE_PRECISION)
    external :: dgttrf,dgttrs
    procedure(), pointer :: gttrf => dgttrf, gttrs => dgttrs
#else
    external :: sgttrf,sgttrs
    procedure(), pointer :: gttrf => sgttrf, gttrs => sgttrs
#endif
    integer , intent(in) :: nx,ny,n
    real(rp), intent(in), dimension(:) :: a,b,c
    real(rp), intent(inout), dimension(:,:,:) :: p
    real(rp), allocatable, dimension(:) :: aa,bb,cc,ccc
    integer , allocatable, dimension(:) :: ipiv
    integer :: i,j,info
    !real(rp), dimension(n,nx,ny) :: p_t
    !
    allocate(aa,source=a(2:n  ))
    allocate(bb,source=b(1:n  ))
    allocate(cc,source=c(1:n-1))
    allocate(ccc(n-2),ipiv(n))
    call gttrf(n,aa,bb,cc,ccc,ipiv,info)
    do j=1,ny
      do i=1,nx
        call gttrs('N',n,1,aa,bb,cc,ccc,ipiv,p(i,j,1:n),n,info)
      end do
    end do
    !p_t = reshape(p(1:nx,1:ny,1:n),shape(p_t),order=[2,3,1])
    !call gttrs('N',n,nx*ny,aa,bb,cc,ccc,ipiv,p_t(1:n,:,:),n,info)
    !p(1:nx,1:ny,1:n) = reshape(p_t,shape(p(1:nx,1:ny,1:n)),order=[3,1,2])
  end subroutine gaussel_lapack
#endif
#endif
end module mod_solver
