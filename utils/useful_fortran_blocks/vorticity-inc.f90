! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
      block
        real(rp), allocatable, dimension(:,:,:) :: vox,voy,voz
        allocate(vox(n(1),n(2),n(3)),voy(n(1),n(2),n(3)),voz(n(1),n(2),n(3)))
        !$acc enter data copyin(vox,voy,voz)
        call vorticity(n,dli,dzci,u,v,w,vox,voy,voz)
        !$acc exit data copyout(vox,voy,voz)
        call write_visu_3d(datadir,'vox_fld_'//fldnum//'.bin','log_visu_3d.out','Vorticity_X', &
                   [1,1,1],[ng(1),ng(2),ng(3)],[1,1,1],time,istep,vox(1:n(1),1:n(2),1:n(3)))
        call write_visu_3d(datadir,'vox_fld_'//fldnum//'.bin','log_visu_3d.out','Vorticity_Y', &
                   [1,1,1],[ng(1),ng(2),ng(3)],[1,1,1],time,istep,voy(1:n(1),1:n(2),1:n(3)))
        call write_visu_3d(datadir,'vox_fld_'//fldnum//'.bin','log_visu_3d.out','Vorticity_Z', &
                   [1,1,1],[ng(1),ng(2),ng(3)],[1,1,1],time,istep,voy(1:n(1),1:n(2),1:n(3)))
        deallocate(vox,voy,voz)
      end block
