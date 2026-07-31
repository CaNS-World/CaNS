! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
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
    ! computes the vorticity components at their natural edge centers:
    !   vox: (cell,face,face), voy: (face,cell,face), voz: (face,face,cell)
    ! valid ranges are vox(1:nx,0:ny,0:nz), voy(0:nx,1:ny,0:nz),
    ! and voz(0:nx,0:ny,1:nz); velocity halos must be current.
    !
    implicit none
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(3)        :: dli
    real(rp), intent(in ), dimension(0:)       :: dzci
    real(rp), intent(in ), dimension(0:,0:,0:) :: ux ,uy ,uz
    real(rp), intent(out), dimension(0:,0:,0:) :: vox,voy,voz
    real(rp) :: dxi,dyi
    integer :: i,j,k
    dxi = dli(1)
    dyi = dli(2)
    !
    ! x component at x-directed edge centers
    !
    !$acc parallel loop collapse(3) default(present)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=0,n(3)
      do j=0,n(2)
        do i=1,n(1)
          vox(i,j,k) = (uz(i,j+1,k)-uz(i,j,k))*dyi - &
                       (uy(i,j,k+1)-uy(i,j,k))*dzci(k)
        end do
      end do
    end do
    !
    ! y component at y-directed edge centers
    !
    !$acc parallel loop collapse(3) default(present)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=0,n(3)
      do j=1,n(2)
        do i=0,n(1)
          voy(i,j,k) = (ux(i,j,k+1)-ux(i,j,k))*dzci(k) - &
                       (uz(i+1,j,k)-uz(i,j,k))*dxi
        end do
      end do
    end do
    !
    ! z component at z-directed edge centers
    !
    !$acc parallel loop collapse(3) default(present)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=1,n(3)
      do j=0,n(2)
        do i=0,n(1)
          voz(i,j,k) = (uy(i+1,j,k)-uy(i,j,k))*dxi - &
                       (ux(i,j+1,k)-ux(i,j,k))*dyi
        end do
      end do
    end do
    !
    ! Previous cell-centered formulation, retained for reference:
    !
    ! vox(i,j,k) = 0.25_rp*( &
    !   (uz(i,j+1,k  )-uz(i,j  ,k  ))*dyi - (uy(i,j  ,k+1)-uy(i,j  ,k  ))*dzci(k  ) + &
    !   (uz(i,j+1,k-1)-uz(i,j  ,k-1))*dyi - (uy(i,j  ,k  )-uy(i,j  ,k-1))*dzci(k-1) + &
    !   (uz(i,j  ,k  )-uz(i,j-1,k  ))*dyi - (uy(i,j-1,k+1)-uy(i,j-1,k  ))*dzci(k  ) + &
    !   (uz(i,j  ,k-1)-uz(i,j-1,k-1))*dyi - (uy(i,j-1,k  )-uy(i,j-1,k-1))*dzci(k-1))
    ! voy(i,j,k) = 0.25_rp*( &
    !   (ux(i  ,j,k+1)-ux(i  ,j,k  ))*dzci(k  ) - (uz(i+1,j,k  )-uz(i  ,j,k  ))*dxi + &
    !   (ux(i  ,j,k  )-ux(i  ,j,k-1))*dzci(k-1) - (uz(i+1,j,k-1)-uz(i  ,j,k-1))*dxi + &
    !   (ux(i-1,j,k+1)-ux(i-1,j,k  ))*dzci(k  ) - (uz(i  ,j,k  )-uz(i-1,j,k  ))*dxi + &
    !   (ux(i-1,j,k  )-ux(i-1,j,k-1))*dzci(k-1) - (uz(i  ,j,k-1)-uz(i-1,j,k-1))*dxi)
    ! voz(i,j,k) = 0.25_rp*( &
    !   (uy(i+1,j  ,k)-uy(i  ,j  ,k))*dxi - (ux(i  ,j+1,k)-ux(i  ,j  ,k))*dyi + &
    !   (uy(i+1,j-1,k)-uy(i  ,j-1,k))*dxi - (ux(i  ,j  ,k)-ux(i  ,j-1,k))*dyi + &
    !   (uy(i  ,j  ,k)-uy(i-1,j  ,k))*dxi - (ux(i-1,j+1,k)-ux(i-1,j  ,k))*dyi + &
    !   (uy(i  ,j-1,k)-uy(i-1,j-1,k))*dxi - (ux(i-1,j  ,k)-ux(i-1,j-1,k))*dyi)
  end subroutine vorticity
  !
  subroutine strain_rate(n,dli,dzci,dzfi,ux,uy,uz,str)
    !
    ! computes Sij*Sij at cell centers, evaluating each shear strain once
    ! at its natural edge center before interpolation
    !
    implicit none
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(3)        :: dli
    real(rp), intent(in ), dimension(0:)       :: dzci,dzfi
    real(rp), intent(in ), dimension(0:,0:,0:) :: ux,uy,uz
    real(rp), intent(out), dimension(1:,1:,1:) :: str
    real(rp), allocatable, dimension(:,:,:) :: s12e,s13e,s23e
    real(rp) :: s11,s22,s33,s12c,s13c,s23c
    real(rp) :: dxi,dyi
    integer :: i,j,k
    dxi = dli(1)
    dyi = dli(2)
    !
    ! compute the off-diagonal strains at edge centers
    !
    allocate(s12e(0:n(1),0:n(2),1:n(3)), &
             s13e(0:n(1),1:n(2),0:n(3)), &
             s23e(1:n(1),0:n(2),0:n(3)))
    !$acc enter data create(s12e,s13e,s23e)
    !$acc parallel loop collapse(3) default(present)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=1,n(3)
      do j=0,n(2)
        do i=0,n(1)
          s12e(i,j,k) = .5_rp*((ux(i,j+1,k)-ux(i,j,k))*dyi + &
                               (uy(i+1,j,k)-uy(i,j,k))*dxi)
        end do
      end do
    end do
    !$acc parallel loop collapse(3) default(present)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=0,n(3)
      do j=1,n(2)
        do i=0,n(1)
          s13e(i,j,k) = .5_rp*((ux(i,j,k+1)-ux(i,j,k))*dzci(k) + &
                               (uz(i+1,j,k)-uz(i,j,k))*dxi)
        end do
      end do
    end do
    !$acc parallel loop collapse(3) default(present)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=0,n(3)
      do j=0,n(2)
        do i=1,n(1)
          s23e(i,j,k) = .5_rp*((uy(i,j,k+1)-uy(i,j,k))*dzci(k) + &
                               (uz(i,j+1,k)-uz(i,j,k))*dyi)
        end do
      end do
    end do
    !
    ! interpolate squared edge strains only when forming the cell-centered invariant
    !
    !$acc parallel loop collapse(3) default(present) private(s11,s12c,s13c,s22,s23c,s33)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  PRIVATE(s11,s12c,s13c,s22,s23c,s33)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          s11  = ((ux(i,j,k)-ux(i-1,j,k))*dxi    )**2
          s22  = ((uy(i,j,k)-uy(i,j-1,k))*dyi    )**2
          s33  = ((uz(i,j,k)-uz(i,j,k-1))*dzfi(k))**2
          s12c = .25_rp*(s12e(i,j,k)**2+s12e(i-1,j,k)**2+s12e(i,j-1,k)**2+s12e(i-1,j-1,k)**2)
          s13c = .25_rp*(s13e(i,j,k)**2+s13e(i-1,j,k)**2+s13e(i,j,k-1)**2+s13e(i-1,j,k-1)**2)
          s23c = .25_rp*(s23e(i,j,k)**2+s23e(i,j-1,k)**2+s23e(i,j,k-1)**2+s23e(i,j-1,k-1)**2)
          str(i,j,k) = s11+s22+s33 + 2._rp*(s12c+s13c+s23c)
        end do
      end do
    end do
    !$acc exit data delete(s12e,s13e,s23e)
    deallocate(s12e,s13e,s23e)
    !
    ! Previous direct cell-centered shear-strain formulation, retained for reference:
    !
    ! s12 = .25_rp*( &
    !   ((ux(i  ,j+1,k)-ux(i  ,j  ,k))*dyi + (uy(i+1,j  ,k)-uy(i  ,j  ,k))*dxi)**2 + &
    !   ((ux(i  ,j  ,k)-ux(i  ,j-1,k))*dyi + (uy(i+1,j-1,k)-uy(i  ,j-1,k))*dxi)**2 + &
    !   ((ux(i-1,j+1,k)-ux(i-1,j  ,k))*dyi + (uy(i  ,j  ,k)-uy(i-1,j  ,k))*dxi)**2 + &
    !   ((ux(i-1,j  ,k)-ux(i-1,j-1,k))*dyi + (uy(i  ,j-1,k)-uy(i-1,j-1,k))*dxi)**2)*.25_rp
    ! s13 = .25_rp*( &
    !   ((ux(i  ,j,k+1)-ux(i  ,j,k  ))*dzci(k  ) + (uz(i+1,j,k  )-uz(i  ,j,k  ))*dxi)**2 + &
    !   ((ux(i  ,j,k  )-ux(i  ,j,k-1))*dzci(k-1) + (uz(i+1,j,k-1)-uz(i  ,j,k-1))*dxi)**2 + &
    !   ((ux(i-1,j,k+1)-ux(i-1,j,k  ))*dzci(k  ) + (uz(i  ,j,k  )-uz(i-1,j,k  ))*dxi)**2 + &
    !   ((ux(i-1,j,k  )-ux(i-1,j,k-1))*dzci(k-1) + (uz(i  ,j,k-1)-uz(i-1,j,k-1))*dxi)**2)*.25_rp
    ! s23 = .25_rp*( &
    !   ((uy(i,j  ,k+1)-uy(i,j  ,k  ))*dzci(k  ) + (uz(i,j+1,k  )-uz(i,j  ,k  ))*dyi)**2 + &
    !   ((uy(i,j  ,k  )-uy(i,j  ,k-1))*dzci(k-1) + (uz(i,j+1,k-1)-uz(i,j  ,k-1))*dyi)**2 + &
    !   ((uy(i,j-1,k+1)-uy(i,j-1,k  ))*dzci(k  ) + (uz(i,j  ,k  )-uz(i,j-1,k  ))*dyi)**2 + &
    !   ((uy(i,j-1,k  )-uy(i,j-1,k-1))*dzci(k-1) + (uz(i,j  ,k-1)-uz(i,j-1,k-1))*dyi)**2)*.25_rp
  end subroutine strain_rate
  !
  subroutine rotation_rate(n,dli,dzci,ux,uy,uz,ens)
    !
    ! computes Wij*Wij at cell centers from edge-centered vorticity
    !
    implicit none
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(3)        :: dli
    real(rp), intent(in ), dimension(0:)       :: dzci
    real(rp), intent(in ), dimension(0:,0:,0:) :: ux,uy,uz
    real(rp), intent(out), dimension(1:,1:,1:) :: ens
    real(rp), allocatable, dimension(:,:,:) :: vox,voy,voz
    integer :: i,j,k
    !
    ! compute each vorticity component once at its natural edge center
    !
    allocate(vox(0:n(1),0:n(2),0:n(3)), &
             voy(0:n(1),0:n(2),0:n(3)), &
             voz(0:n(1),0:n(2),0:n(3)))
    !$acc enter data create(vox,voy,voz)
    call vorticity(n,dli,dzci,ux,uy,uz,vox,voy,voz)
    !
    ! Wij*Wij = |vorticity|^2/2; interpolate squared edge values to cells
    !
    !$acc parallel loop collapse(3) default(present)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          ens(i,j,k) = .125_rp*( &
            vox(i,j,k  )**2+vox(i,j-1,k  )**2+vox(i,j,k-1)**2+vox(i,j-1,k-1)**2 + &
            voy(i,j,k  )**2+voy(i-1,j,k  )**2+voy(i,j,k-1)**2+voy(i-1,j,k-1)**2 + &
            voz(i,j,k  )**2+voz(i-1,j,k  )**2+voz(i,j-1,k)**2+voz(i-1,j-1,k)**2)
        end do
      end do
    end do
    !$acc exit data delete(vox,voy,voz)
    deallocate(vox,voy,voz)
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
