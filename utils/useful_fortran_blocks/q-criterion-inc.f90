! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
      block
        real(rp), allocatable, dimension(:,:,:) :: str,rot,qcr
        allocate(str(n(1),n(2),n(3)),rot(n(1),n(2),n(3)), &
                 qcr(0:n(1)+1,0:n(2)+1,0:n(3)+1))
        !$acc enter data create(str,rot,qcr)
        call strain_rate(n,dli,dzci,dzfi,u,v,w,str)
        call rotation_rate(n,dli,dzci,u,v,w,rot)
        call q_criterion(n,rot,str,qcr)
        !$acc exit data copyout(str,rot,qcr)
        call write_visu_3d(datadir,'qcr_fld_'//fldnum//trim(io_ext),'log_visu_3d.out','Q-criterion', &
                           ng,[1,1,1],lo,hi,[1,1,1],ng,[1,1,1],time,istep,qcr,xc_g,yc_g,zc_g,.true.)
        deallocate(str,rot,qcr)
      end block
