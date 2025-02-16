! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_updatep
  use mod_types
  use mod_param, only: is_impdiff,is_impdiff_1d
  implicit none
  private
  public updatep
  contains
  subroutine updatep(n,dli,dzci,dzfi,alpha,pp,p)
    !
    ! updates the final pressure
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    real(rp), intent(in   ), dimension(3 ) :: dli
    real(rp), intent(in   ), dimension(0:) :: dzci,dzfi
    real(rp), intent(in   ) :: alpha
    real(rp), intent(in   ), dimension(0:,0:,0:) :: pp
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    real(rp) :: dxi,dyi
    integer :: i,j,k
    real(rp) :: lap_pp
    !
    if(is_impdiff) then
      dxi = dli(1); dyi = dli(2)
#if !defined(_LOOP_UNSWITCHING)
      !$acc parallel loop collapse(3) default(present) private(lap_pp) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared ) private(lap_pp)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            lap_pp = 0.
            if(.not.is_impdiff_1d) then
              lap_pp = lap_pp + (pp(i+1,j,k)-2.*pp(i,j,k)+pp(i-1,j,k))*(dxi**2) + &
                                (pp(i,j+1,k)-2.*pp(i,j,k)+pp(i,j-1,k))*(dyi**2)
            end if
            lap_pp = lap_pp + ((pp(i,j,k+1)-pp(i,j,k  ))*dzci(k  ) - &
                               (pp(i,j,k  )-pp(i,j,k-1))*dzci(k-1))*dzfi(k)
            p(i,j,k) = p(i,j,k) + pp(i,j,k) + alpha*lap_pp
          end do
        end do
      end do
#else
      if(is_impdiff_1d) then
        !$acc parallel loop collapse(3) default(present) private(lap_pp) async(1)
        !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared ) private(lap_pp)
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              lap_pp = ((pp(i,j,k+1)-pp(i,j,k  ))*dzci(k) - &
                        (pp(i,j,k  )-pp(i,j,k-1))*dzci(k-1))*dzfi(k)
              p(i,j,k) = p(i,j,k) + pp(i,j,k) + alpha*lap_pp
            end do
          end do
        end do
      else
        !$acc parallel loop collapse(3) default(present) private(lap_pp) async(1)
        !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared ) private(lap_pp)
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              lap_pp = (pp(i+1,j,k)-2.*pp(i,j,k)+pp(i-1,j,k))*(dxi**2) + &
                       (pp(i,j+1,k)-2.*pp(i,j,k)+pp(i,j-1,k))*(dyi**2) + &
                      ((pp(i,j,k+1)-pp(i,j,k  ))*dzci(k) - &
                       (pp(i,j,k  )-pp(i,j,k-1))*dzci(k-1))*dzfi(k)
              p(i,j,k) = p(i,j,k) + pp(i,j,k) + alpha*lap_pp
            end do
          end do
        end do
      end if
#endif
    else
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            p(i,j,k) = p(i,j,k) + pp(i,j,k)
          end do
        end do
      end do
    end if
  end subroutine updatep
end module mod_updatep
