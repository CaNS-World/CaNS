! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
program gen_grid
!
! this program generates the grid files to be read by the xdmf file
! the output file from CaNS 'grid.bin' must be in this folder
!
! Pedro Costa (p.simoes.costa@gmail.com)
!
implicit none
include 'param.h90'
!
integer :: iunit,i,j,lenr
real(8), dimension(nz) :: dummy,zc
!
iunit = 99
inquire (iolength=lenr) dummy(1)
!
! generate cell-centered grid in x (uniform)
!
open(iunit,file='x.bin',access='direct',recl=nx*lenr)
write(iunit,rec=1) ((i-0.5)*dx,i=1,nx)
close(iunit)
!
! generate cell-centered grid in y (uniform)
!
open(iunit,file='y.bin',access='direct',recl=ny*lenr)
write(iunit,rec=1) ((j-0.5)*dy,j=1,ny)
close(iunit)
!
! generate cell-centered grid in z (non-uniform) from file 'grid.bin'
!
open(iunit,file='grid.bin',action='read',form='unformatted',access='stream',status='old')
read(iunit) dummy(1:nz),dummy(1:nz),zc(1:nz),dummy(1:nz)
close(iunit)
open(iunit,file='z.bin',action='write',form='unformatted',access='stream',status='replace')
write(iunit) zc(1:nz)
close(iunit)
end program gen_grid
