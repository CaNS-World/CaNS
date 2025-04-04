! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_solve_helmholtz
  use, intrinsic :: iso_c_binding, only: C_PTR
  use mod_types
  use mod_bound, only: updt_rhs_b
#if !defined(_OPENACC)
  use mod_solver    , only: solver
#if defined(_IMPDIFF_1D)
  use mod_solver    , only: solver_gaussel_z
#endif
#else
  use mod_solver_gpu, only: solver => solver_gpu
#if defined(_IMPDIFF_1D)
  use mod_solver_gpu, only: solver_gaussel_z => solver_gaussel_z_gpu
#endif
#endif
  implicit none
  private
  type rhs_bound
    real(rp), allocatable, dimension(:,:,:) :: x
    real(rp), allocatable, dimension(:,:,:) :: y
    real(rp), allocatable, dimension(:,:,:) :: z
  end type rhs_bound
  public solve_helmholtz,rhs_bound
  contains
  subroutine solve_helmholtz(n,ng,arrplan,normfft,alpha, &
                             lambdaxyi,ai,bi,ci,rhsbxi,rhsbyi,rhsbzi,is_bound,cbc,c_or_f,p)
    !
    ! this is a wrapper subroutine to solve 1D/3D helmholtz problems: p/alpha + lap(p) = rhs
    !
    integer ,    intent(in   ), dimension(3)                :: n,ng
#if !defined(_OPENACC)
    type(C_PTR), intent(in   ), dimension(2,2),    optional :: arrplan
#else
    integer    , intent(in   ), dimension(2,2),    optional :: arrplan
#endif
    real(rp),    intent(in   ),                    optional :: normfft
    real(rp),    intent(in   )                              :: alpha
    real(rp),    intent(in   ), dimension(:,:),    optional :: lambdaxyi
    real(rp),    intent(in   ), dimension(:)                :: ai,bi,ci
    real(rp),    intent(in   ), dimension(:,:,0:), optional :: rhsbxi,rhsbyi,rhsbzi
    logical ,    intent(in   ), dimension(2,3)              :: is_bound
    character(len=1), intent(in), dimension(0:1,3)          :: cbc
    character(len=1), intent(in), dimension(3)              :: c_or_f
    real(rp),    intent(inout), dimension(:,:,:)            :: p
    real(rp), allocatable, dimension(:,:)  , save :: lambdaxy
    real(rp), allocatable, dimension(:)    , save :: a,b,c
    real(rp), allocatable, dimension(:,:,:), save :: rhsbx,rhsby,rhsbz
    !
    logical, save :: is_first = .true.
    !
    ! initialization
    !
    if(is_first) then ! leverage save attribute to allocate these arrays on the device only once
      is_first = .false.
      if(present(lambdaxyi)) allocate(lambdaxy,mold=lambdaxyi)
      allocate(a,mold=ai)
      allocate(b,mold=bi)
      allocate(c,mold=ci)
      if(present(rhsbxi)) allocate(rhsbx(n(2),n(3),0:1)) ! allocate(rhsbx,mold=rhsbxi) ! gfortran 11.4.0 bug
      if(present(rhsbyi)) allocate(rhsby(n(1),n(3),0:1)) ! allocate(rhsby,mold=rhsbyi) ! gfortran 11.4.0 bug
      if(present(rhsbzi)) allocate(rhsbz(n(1),n(2),0:1)) ! allocate(rhsbz,mold=rhsbzi) ! gfortran 11.4.0 bug
      !$acc enter data create(lambdaxy,a,b,c,rhsbx,rhsby,rhsbz) async(1)
    end if
    !
    !$acc kernels default(present) async(1)
    !$OMP PARALLEL WORKSHARE
#if !defined(_IMPDIFF_1D)
    rhsbx(:,:,0:1) = rhsbxi(:,:,0:1)*alpha
    rhsby(:,:,0:1) = rhsbyi(:,:,0:1)*alpha
#endif
    rhsbz(:,:,0:1) = rhsbzi(:,:,0:1)*alpha
    !$OMP END PARALLEL WORKSHARE
    !$acc end kernels
    call updt_rhs_b(c_or_f,cbc,n,is_bound,rhsbx,rhsby,rhsbz,p)
    !$acc kernels default(present) async(1)
    !$OMP PARALLEL WORKSHARE
    a(:) = ai(:)*alpha
    b(:) = bi(:)*alpha + 1.
    c(:) = ci(:)*alpha
#if !defined(_IMPDIFF_1D)
    lambdaxy(:,:) = lambdaxyi(:,:)*alpha
#endif
    !$OMP END PARALLEL WORKSHARE
    !$acc end kernels
#if !defined(_IMPDIFF_1D)
    call solver(          n,ng,arrplan,normfft,lambdaxy,a,b,c,cbc     ,c_or_f,p)
#else
    call solver_gaussel_z(n                            ,a,b,c,cbc(:,3),c_or_f,p)
#endif
  end subroutine solve_helmholtz
end module mod_solve_helmholtz
