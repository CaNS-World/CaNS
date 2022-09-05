! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
      block
        real(rp), allocatable, dimension(:,:,:) :: str,rot,qcr
        allocate(str(n(1),n(2),n(3)),rot(n(1),n(2),n(3)),qcr(n(1),n(2),n(3)))
        !$acc enter data copyin(str,ens,qcr)
        call strain_rate(n,dli,dzci,dzfi,u,v,w,str)
        call rotation_rate(n,dli,dzci,u,v,w,rot)
        call q_criterion(n,rot,str,qcr)
        !$acc exit data copyout(str,ens,qcr)
        call write_visu_3d(datadir,'qcr_fld_'//fldnum//'.bin','log_visu_3d.out','Q-criterion', &
                   [1,1,1],[ng(1),ng(2),ng(3)],[1,1,1],time,istep,qcr(1:n(1),1:n(2),1:n(3)))
        deallocate(str,rot,qcr)
      end block
