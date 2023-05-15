! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_post
  use mod_types
  implicit none
  private
  public vorticity,rotation_rate,strain_rate,q_criterion
  contains
  subroutine vorticity(n,dli,dzci,ux,uy,uz,vox,voy,voz)
    !
    ! computes the vorticity field
    !
    implicit none
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(3)        :: dli
    real(rp), intent(in ), dimension(0:)       :: dzci
    real(rp), intent(in ), dimension(0:,0:,0:) :: ux ,uy ,uz
    real(rp), intent(out), dimension( :, :, :) :: vox,voy,voz
    real(rp) :: dxi,dyi
    integer :: i,j,k
    dxi = dli(1)
    dyi = dli(2)
    !$acc parallel loop collapse(3) default(present)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          ! x component of the vorticity at cell center
          !
          vox(i,j,k) = 0.25_rp*( &
                                (uz(i,j+1,k  )-uz(i,j  ,k  ))*dyi - (uy(i,j  ,k+1)-uy(i,j  ,k  ))*dzci(k  ) + &
                                (uz(i,j+1,k-1)-uz(i,j  ,k-1))*dyi - (uy(i,j  ,k  )-uy(i,j  ,k-1))*dzci(k-1) + &
                                (uz(i,j  ,k  )-uz(i,j-1,k  ))*dyi - (uy(i,j-1,k+1)-uy(i,j-1,k  ))*dzci(k  ) + &
                                (uz(i,j  ,k-1)-uz(i,j-1,k-1))*dyi - (uy(i,j-1,k  )-uy(i,j-1,k-1))*dzci(k-1) &
                               )
          !
          ! y component of the vorticity at cell center
          !
          voy(i,j,k) = 0.25_rp*( &
                                (ux(i  ,j,k+1)-ux(i  ,j,k  ))*dzci(k  ) - (uz(i+1,j,k  )-uz(i  ,j,k  ))*dxi + &
                                (ux(i  ,j,k  )-ux(i  ,j,k-1))*dzci(k-1) - (uz(i+1,j,k-1)-uz(i  ,j,k-1))*dxi + &
                                (ux(i-1,j,k+1)-ux(i-1,j,k  ))*dzci(k  ) - (uz(i  ,j,k  )-uz(i-1,j,k  ))*dxi + &
                                (ux(i-1,j,k  )-ux(i-1,j,k-1))*dzci(k-1) - (uz(i  ,j,k-1)-uz(i-1,j,k-1))*dxi &
                               )
          !
          ! z component of the vorticity at cell center
          !
          voz(i,j,k) = 0.25_rp*( &
                                (uy(i+1,j  ,k)-uy(i  ,j  ,k))*dxi - (ux(i  ,j+1,k)-ux(i  ,j  ,k))*dyi + &
                                (uy(i+1,j-1,k)-uy(i  ,j-1,k))*dxi - (ux(i  ,j  ,k)-ux(i  ,j-1,k))*dyi + &
                                (uy(i  ,j  ,k)-uy(i-1,j  ,k))*dxi - (ux(i-1,j+1,k)-ux(i-1,j  ,k))*dyi + &
                                (uy(i  ,j-1,k)-uy(i-1,j-1,k))*dxi - (ux(i-1,j  ,k)-ux(i-1,j-1,k))*dyi &
                               )
        end do
      end do
    end do
  end subroutine vorticity
  !
  subroutine strain_rate(n,dli,dzci,dzfi,ux,uy,uz,str)
    implicit none
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(3)        :: dli
    real(rp), intent(in ), dimension(0:)       :: dzci,dzfi
    real(rp), intent(in ), dimension(0:,0:,0:) :: ux,uy,uz
    real(rp), intent(out), dimension(1:,1:,1:) :: str
    real(rp) :: s11,s22,s33,s12,s13,s23
    real(rp) :: dxi,dyi
    integer :: i,j,k
    dxi = dli(1)
    dyi = dli(2)
    !
    ! compute sijsij, where sij = (1/2)(du_i/dx_j + du_j/dx_i)
    !
    !$acc parallel loop collapse(3) default(present) private(s11,s12,s13,s22,s23,s33)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  PRIVATE(s11,s12,s13,s22,s23,s33)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          s11 = ((ux(i,j,k)-ux(i-1,j,k))*dxi    )**2
          s22 = ((uy(i,j,k)-uy(i,j-1,k))*dyi    )**2
          s33 = ((uz(i,j,k)-uz(i,j,k-1))*dzfi(k))**2
          s12 = .25_rp*( &
                        ((ux(i  ,j+1,k)-ux(i  ,j  ,k))*dyi + (uy(i+1,j  ,k)-uy(i  ,j  ,k))*dxi)**2 + &
                        ((ux(i  ,j  ,k)-ux(i  ,j-1,k))*dyi + (uy(i+1,j-1,k)-uy(i  ,j-1,k))*dxi)**2 + &
                        ((ux(i-1,j+1,k)-ux(i-1,j  ,k))*dyi + (uy(i  ,j  ,k)-uy(i-1,j  ,k))*dxi)**2 + &
                        ((ux(i-1,j  ,k)-ux(i-1,j-1,k))*dyi + (uy(i  ,j-1,k)-uy(i-1,j-1,k))*dxi)**2 &
                       )*.25_rp
          s13 = .25_rp*( &
                         ((ux(i  ,j,k+1)-ux(i  ,j,k  ))*dzci(k  ) + (uz(i+1,j,k  )-uz(i  ,j,k  ))*dxi)**2 + &
                         ((ux(i  ,j,k  )-ux(i  ,j,k-1))*dzci(k-1) + (uz(i+1,j,k-1)-uz(i  ,j,k-1))*dxi)**2 + &
                         ((ux(i-1,j,k+1)-ux(i-1,j,k  ))*dzci(k  ) + (uz(i  ,j,k  )-uz(i-1,j,k  ))*dxi)**2 + &
                         ((ux(i-1,j,k  )-ux(i-1,j,k-1))*dzci(k-1) + (uz(i  ,j,k-1)-uz(i-1,j,k-1))*dxi)**2 &
                       )*.25_rp
          s23 = .25_rp*( &
                        ((uy(i,j  ,k+1)-uy(i,j  ,k  ))*dzci(k  ) + (uz(i,j+1,k  )-uz(i,j  ,k  ))*dyi)**2 + &
                        ((uy(i,j  ,k  )-uy(i,j  ,k-1))*dzci(k-1) + (uz(i,j+1,k-1)-uz(i,j  ,k-1))*dyi)**2 + &
                        ((uy(i,j-1,k+1)-uy(i,j-1,k  ))*dzci(k  ) + (uz(i,j  ,k  )-uz(i,j-1,k  ))*dyi)**2 + &
                        ((uy(i,j-1,k  )-uy(i,j-1,k-1))*dzci(k-1) + (uz(i,j  ,k-1)-uz(i,j-1,k-1))*dyi)**2 &
                       )*.25_rp
          str(i,j,k) = s11+s22+s33 + 2*(s12+s13+s23)
        end do
      end do
    end do
  end subroutine strain_rate
  !
  subroutine rotation_rate(n,dli,dzci,ux,uy,uz,ens)
    implicit none
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(3)        :: dli
    real(rp), intent(in ), dimension(0:)       :: dzci
    real(rp), intent(in ), dimension(0:,0:,0:) :: ux,uy,uz
    real(rp), intent(out), dimension(1:,1:,1:) :: ens
    real(rp) :: e12,e13,e23
    real(rp) :: dxi,dyi
    integer :: i,j,k
    !
    ! compute wijwij, where wij = (1/2)(du_i/dx_j - du_j/dx_i)
    !
    dxi = dli(1)
    dyi = dli(2)
    !$acc parallel loop collapse(3) default(present) private(e12,e13,e23)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  PRIVATE(e12,e13,e23)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          e12 = .25_rp*( &
                        ((ux(i  ,j+1,k)-ux(i  ,j  ,k))*dyi - (uy(i+1,j  ,k)-uy(i  ,j  ,k))*dxi)**2 + &
                        ((ux(i  ,j  ,k)-ux(i  ,j-1,k))*dyi - (uy(i+1,j-1,k)-uy(i  ,j-1,k))*dxi)**2 + &
                        ((ux(i-1,j+1,k)-ux(i-1,j  ,k))*dyi - (uy(i  ,j  ,k)-uy(i-1,j  ,k))*dxi)**2 + &
                        ((ux(i-1,j  ,k)-ux(i-1,j-1,k))*dyi - (uy(i  ,j-1,k)-uy(i-1,j-1,k))*dxi)**2 &
                       )*.25_rp
          e13 = .25_rp*( &
                        ((ux(i  ,j,k+1)-ux(i  ,j,k  ))*dzci(k  ) - (uz(i+1,j,k  )-uz(i  ,j,k  ))*dxi)**2 + &
                        ((ux(i  ,j,k  )-ux(i  ,j,k-1))*dzci(k-1) - (uz(i+1,j,k-1)-uz(i  ,j,k-1))*dxi)**2 + &
                        ((ux(i-1,j,k+1)-ux(i-1,j,k  ))*dzci(k  ) - (uz(i  ,j,k  )-uz(i-1,j,k  ))*dxi)**2 + &
                        ((ux(i-1,j,k  )-ux(i-1,j,k-1))*dzci(k-1) - (uz(i  ,j,k-1)-uz(i-1,j,k-1))*dxi)**2 &
                       )*.25_rp
          e23 = .25_rp*( &
                        ((uy(i,j  ,k+1)-uy(i,j  ,k  ))*dzci(k  ) - (uz(i,j+1,k  )-uz(i,j  ,k  ))*dyi)**2 + &
                        ((uy(i,j  ,k  )-uy(i,j  ,k-1))*dzci(k-1) - (uz(i,j+1,k-1)-uz(i,j  ,k-1))*dyi)**2 + &
                        ((uy(i,j-1,k+1)-uy(i,j-1,k  ))*dzci(k  ) - (uz(i,j  ,k  )-uz(i,j-1,k  ))*dyi)**2 + &
                        ((uy(i,j-1,k  )-uy(i,j-1,k-1))*dzci(k-1) - (uz(i,j  ,k-1)-uz(i,j-1,k-1))*dyi)**2 &
                       )*.25_rp
          ens(i,j,k) =  2._rp*(e12+e13+e23)
        end do
      end do
    end do
  end subroutine rotation_rate
  !
  subroutine q_criterion(n,ens,str,qcr)
    implicit none
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(1:,1:,1:) :: ens,str
    real(rp), intent(out), dimension(0:,0:,0:) :: qcr
    integer  :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          qcr(i,j,k) = .5_rp*(ens(i,j,k)-str(i,j,k))
        end do
      end do
    end do
  end subroutine q_criterion
end module mod_post
