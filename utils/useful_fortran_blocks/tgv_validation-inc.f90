! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
      block
        real(rp) :: dx,dy
        real(rp) :: ycl,yfl,xcl,xfl
        real(rp) :: err1u,err2u,err1v,err2v,err1p,err2p
        real(rp) :: utru,vtru,ptru,erru,errv,errp
        integer :: i,j,k
        dx = dl(1)
        dy = dl(2)
        err1u = 0._rp
        err2u = 0._rp
        err1v = 0._rp
        err2v = 0._rp
        err1p = 0._rp
        err2p = 0._rp
        do k=lo(3),hi(3)
          do j=lo(2),hi(2)
            ycl = (j+lo(2)-1-.5)*dy
            yfl = (j+lo(2)-1-.0)*dy
            do i=lo(1),hi(1)
              xcl = (i+lo(1)-1-.5)*dx
              xfl = (i+lo(1)-1-.0)*dx
              utru =  cos(xfl)*sin(ycl)*uref*exp(-2.*visc*time)
              vtru = -sin(xcl)*cos(yfl)*uref*exp(-2.*visc*time)
              ptru = -(cos(2.*xcl)+cos(2.*ycl))/4.*exp(-4.*visc*time)*uref**2
              !
              erru  = abs(u(i,j,k) - utru)
              errv  = abs(v(i,j,k) - vtru)
              errp  = abs(p(i,j,k) - ptru)
              err1u = err1u + erru**1*dx*dy*dzf(k)/product(l)
              err2u = err2u + erru**2*dx*dy*dzf(k)/product(l)
              err1v = err1v + errv**1*dx*dy*dzc(k)/product(l)
              err2v = err2v + errv**2*dx*dy*dzc(k)/product(l)
              err1p = err1p + errp**1*dx*dy*dzf(k)/product(l)
              err2p = err2p + errp**2*dx*dy*dzf(k)/product(l)
            end do
          end do
        end do
        call MPI_ALLREDUCE(MPI_IN_PLACE,err1u,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,err2u,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,err1v,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,err2v,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,err1p,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,err2p,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
        var(1) = 1._rp*istep
        var(2) = time
        var(3:4) = [err1u,sqrt(err2u)]
        var(5:6) = [err1v,sqrt(err2v)]
        var(7:8) = [err1p,sqrt(err2p)]
        call out0d(trim(datadir)//'err_tgv.out',8,var)
      end block
