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
  use mod_param, only: is_impdiff_1d
#if !defined(_OPENACC)
  use mod_solver    , only: solver
  use mod_solver    , only: solver_gaussel_z
#else
  use mod_solver_gpu, only: solver => solver_gpu
  use mod_solver_gpu, only: solver_gaussel_z => solver_gaussel_z_gpu
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
  subroutine solve_helmholtz(n,ng,hi,arrplan,normfft,alpha,lambdaxy,a,b,c,rhsbx,rhsby,rhsbz,is_bound,cbc,c_or_f,p)
    !
    ! this is a wrapper subroutine to solve 1D/3D helmholtz problems: p/alpha + lap(p) = rhs
    !
    integer ,    intent(in   ), dimension(3)                :: n,ng,hi
#if !defined(_OPENACC)
    type(C_PTR), intent(in   ), dimension(2,2),    optional :: arrplan
#else
    integer    , intent(in   ), dimension(2,2),    optional :: arrplan
#endif
    real(rp),    intent(in   ),                    optional :: normfft
    real(rp),    intent(in   )                              :: alpha
    real(rp),    intent(in   ), dimension(:,:),    optional :: lambdaxy
    real(rp),    intent(in   ), dimension(:)                :: a,b,c
    real(rp),    intent(in   ), dimension(:,:,0:), optional :: rhsbx,rhsby,rhsbz
    logical ,    intent(in   ), dimension(2,3)              :: is_bound
    character(len=1), intent(in), dimension(0:1,3)          :: cbc
    character(len=1), intent(in), dimension(3)              :: c_or_f
    real(rp),    intent(inout), dimension(:,:,:)            :: p
    real(rp), allocatable, dimension(:), save :: bb
    real(rp) :: alphai
    integer :: k
    !
    logical, save :: is_first = .true.
    !
    ! initialization
    !
    if(is_first) then ! leverage save attribute to allocate these arrays on the device only once
      is_first = .false.
      allocate(bb,mold=b)
      !$acc enter data create(bb) async(1)
    end if
    !
    call updt_rhs_b(c_or_f,cbc,n,is_bound,rhsbx,rhsby,rhsbz,p,alpha)
    !
    alphai = alpha**(-1)
    !$acc parallel loop default(present) async(1)
    !$OMP PARALLEL DO   DEFAULT(shared)
    do k=1,size(b)
      bb(k) = b(k) + alphai
    end do
    !
    if(.not.is_impdiff_1d) then
      call solver(n,ng,arrplan,normfft*alphai,lambdaxy,a,bb,c,cbc,c_or_f,p)
    else
      call solver_gaussel_z(n,ng,hi,a,bb,c,cbc(:,3),c_or_f,alphai,p)
    end if
  end subroutine solve_helmholtz
end module mod_solve_helmholtz
