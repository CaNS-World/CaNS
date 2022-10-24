! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_debug
  use mpi
  use mod_common_mpi, only: myid,ierr
  use mod_param     , only: dims
  use mod_types
  implicit none
  private
  public chk_helmholtz
  contains
  !
  subroutine chk_helmholtz(lo,hi,dli,dzci,dzfi,alpha,fp,fpp,bc,is_bound,c_or_f,diffmax)
    !
    ! this subroutine checks if the implementation of implicit diffusion is
    ! correct under sanity.f90
    !
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(2) :: dli
    real(rp), intent(in) :: alpha
    real(rp), intent(in), dimension(lo(3)-1:) :: dzci,dzfi
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: fp,fpp
    character(len=1), intent(in), dimension(0:1,3) :: bc
    character(len=1), intent(in), dimension(3) :: c_or_f
    logical         , intent(in), dimension(0:1,3) :: is_bound
    real(rp), intent(out) :: diffmax
    real(rp) :: val
    integer :: i,j,k
    integer :: idir
    integer, dimension(3) :: q
    real(rp) :: dxi,dyi
    dxi = dli(1); dyi = dli(2)
    q(:) = 0
    do idir = 1,3
      if(bc(1,idir) /= 'P'.and.c_or_f(idir) == 'f'.and.is_bound(1,idir)) q(idir) = 1
    end do
    diffmax = 0.
    !$acc wait
    !$acc data copy(diffmax,q)
    select case(c_or_f(3))
    case('c')
      !$acc parallel loop collapse(3) default(present) private(val) reduction(max:diffmax)
      do k=lo(3),hi(3)-q(3)
        do j=lo(2),hi(2)-q(2)
          do i=lo(1),hi(1)-q(1)
            val =  fpp(i,j,k)+(1./alpha)*( &
                  (fpp(i+1,j,k)-2.*fpp(i,j,k)+fpp(i-1,j,k))*(dxi**2) + &
                  (fpp(i,j+1,k)-2.*fpp(i,j,k)+fpp(i,j-1,k))*(dyi**2) + &
                 ((fpp(i,j,k+1)-fpp(i,j,k  ))*dzci(k ) - &
                  (fpp(i,j,k  )-fpp(i,j,k-1))*dzci(k-1))*dzfi(k) )
            val = val*alpha
            diffmax = max(diffmax,abs(val-fp(i,j,k)))
            !if(abs(val-fp(i,j,k)) > 1.e-8) print*, 'Large difference : ', val-fp(i,j,k),i,j,k
          end do
        end do
      end do
    case('f')
      !$acc parallel loop collapse(3) default(present) private(val) reduction(max:diffmax)
      do k=lo(3),hi(3)-q(3)
        do j=lo(2),hi(2)-q(2)
          do i=lo(1),hi(1)-q(1)
            val =  fpp(i,j,k)+(1./alpha)*( &
                  (fpp(i+1,j,k)-2.*fpp(i,j,k)+fpp(i-1,j,k))*(dxi**2) + &
                  (fpp(i,j+1,k)-2.*fpp(i,j,k)+fpp(i,j-1,k))*(dyi**2) + &
                 ((fpp(i,j,k+1)-fpp(i,j,k  ))*dzfi(k+1) - &
                  (fpp(i,j,k  )-fpp(i,j,k-1))*dzfi(k ))*dzci(k) )
            val = val*alpha
            diffmax = max(diffmax,abs(val-fp(i,j,k)))
            !if(abs(val-fp(i,j,k)) > 1.e-8) print*, 'Large difference : ', val,fp(i,j,k),i,j,k
          end do
        end do
      end do
    end select
    !$acc end data
    call MPI_ALLREDUCE(MPI_IN_PLACE,diffmax,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
  end subroutine chk_helmholtz
  subroutine chk_poisson(lo,hi,dli,dzci,dzfi,fp,fpp,diffmax)
    !
    ! this subroutine checks if the Poisson equation is correctly solved;
    ! the inputs should be downcasted to precision `gp` to ensure proper
    ! testing when single precision is used
    !
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(2) :: dli
    real(rp), intent(in), dimension(lo(3)-1:) :: dzci,dzfi
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: fp,fpp
    real(rp), intent(out) :: diffmax
    real(rp) :: val
    integer :: i,j,k
    real(rp) :: dxi,dyi
    dxi = dli(1); dyi = dli(2)
    diffmax = 0.
    !$acc wait
    !$acc data copy(diffmax)
    !$acc parallel loop collapse(3) default(present) private(val) reduction(max:diffmax)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          val = (fpp(i+1,j,k)-2.*fpp(i,j,k)+fpp(i-1,j,k))*(dxi**2) + &
                (fpp(i,j+1,k)-2.*fpp(i,j,k)+fpp(i,j-1,k))*(dyi**2) + &
               ((fpp(i,j,k+1)-fpp(i,j,k  ))*dzci(k ) - &
                (fpp(i,j,k  )-fpp(i,j,k-1))*dzci(k-1))*dzfi(k)
          diffmax = max(diffmax,abs(val-fp(i,j,k)))
          !if(abs(val-fp(i,j,k)) > 1.e-8) print*, 'Large difference : ', val-fp(i,j,k),i,j,k
        end do
      end do
    end do
    !$acc end data
    call MPI_ALLREDUCE(MPI_IN_PLACE,diffmax,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
  end subroutine chk_poisson
end module mod_debug
