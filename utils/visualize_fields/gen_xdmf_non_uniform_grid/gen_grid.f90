program gen_grid
!
! this program generates a the grid files to be read by the xdmf file
! the output file from CaNS 'grid.bin' must be on this folder
!
! Pedro Costa (p.simoes.costa@gmail.com)
!
implicit none
include 'param.h90'
!
integer :: iunit,i,j,k,lenr
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
open(iunit,file='grid.bin',access='direct',recl=4*nz*lenr)
read(iunit,rec=1) dummy(1:nz),dummy(1:nz),zc(1:nz),dummy(1:nz)
close(iunit)
open(iunit,file='z.bin',access='direct',recl=nz*lenr)
write(iunit,rec=1) zc(1:nz)
close(iunit)
end program gen_grid
