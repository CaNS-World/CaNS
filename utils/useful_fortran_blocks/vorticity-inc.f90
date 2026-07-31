! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
      block
        real(rp), allocatable, dimension(:,:,:) :: vox,voy,voz
        allocate(vox(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
                 voy(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
                 voz(0:n(1)+1,0:n(2)+1,0:n(3)+1))
        !$acc enter data create(vox,voy,voz)
        call vorticity(n,dli,dzci,u,v,w,vox,voy,voz)
        !$acc exit data copyout(vox,voy,voz)
        call write_visu_3d(datadir,'vox_fld_'//fldnum//trim(io_ext),'log_visu_3d.out','Vorticity_X', &
                           ng,[1,1,1],lo,hi,[1,1,1],ng,[1,1,1],time,istep,vox,xc_g,yf_g,zf_g,.true.)
        call write_visu_3d(datadir,'voy_fld_'//fldnum//trim(io_ext),'log_visu_3d.out','Vorticity_Y', &
                           ng,[1,1,1],lo,hi,[1,1,1],ng,[1,1,1],time,istep,voy,xf_g,yc_g,zf_g,.true.)
        call write_visu_3d(datadir,'voz_fld_'//fldnum//trim(io_ext),'log_visu_3d.out','Vorticity_Z', &
                           ng,[1,1,1],lo,hi,[1,1,1],ng,[1,1,1],time,istep,voz,xf_g,yf_g,zc_g,.true.)
        deallocate(vox,voy,voz)
      end block
