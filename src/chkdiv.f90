! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_chkdiv
  use mpi
  use mod_common_mpi, only: ierr
  use mod_types
  implicit none
  private
  public chkdiv
  contains
  subroutine chkdiv(lo,hi,dli,dzfi,u,v,w,divtot,divmax)
    !
    ! checks the divergence of the velocity field
    !
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(3) :: dli
    real(rp), intent(in), dimension(lo(3)-1:) :: dzfi
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: u,v,w
    real(rp), intent(out) :: divtot,divmax
    real(rp) :: dxi,dyi,div!,dzi
    integer :: i,j,k
    !
    dxi = dli(1)
    dyi = dli(2)
    !dzi = dli(3)
    divtot = 0.
    divmax = 0.
    !$acc data copy(divtot,divmax) async(1)
    !$acc parallel loop collapse(3) default(present) private(div) reduction(+:divtot) reduction(max:divmax) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  PRIVATE(div) REDUCTION(+:divtot) REDUCTION(max:divmax)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          div = (w(i,j,k)-w(i,j,k-1))*dzfi(k) + &
                (v(i,j,k)-v(i,j-1,k))*dyi     + &
                (u(i,j,k)-u(i-1,j,k))*dxi
          divmax = max(divmax,abs(div))
          divtot = divtot + div
          !if(abs(div) >= 1.e-12) print*,div,'Large divergence at grid cell: ',i,j,k,div
        end do
      end do
    end do
    !$acc end data
    !$acc wait(1)
    call MPI_ALLREDUCE(MPI_IN_PLACE,divtot,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,divmax,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
  end subroutine chkdiv
end module mod_chkdiv
