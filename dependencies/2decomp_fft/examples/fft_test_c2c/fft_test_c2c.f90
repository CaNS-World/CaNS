!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main test program for the FFT interface
!  - use input data from a NAG FFT library for validation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program fft_test_c2c

use decomp_2d
use decomp_2d_fft

implicit none

integer, parameter :: nx=2, ny=3, nz=4
integer, parameter :: p_row=2, p_col=2

complex(mytype), allocatable, dimension(:,:,:) :: in, out

complex(mytype), dimension(nx,ny,nz) :: in1, out1
integer :: ierror, i,j,k

interface
   subroutine assemble_global(ndir,local,global,nx,ny,nz)
     use decomp_2d
     integer, intent(IN) :: ndir
     integer, intent(IN) :: nx,ny,nz
     complex(mytype), dimension(:,:,:), intent(IN) :: local
     complex(mytype), dimension(nx,ny,nz), intent(OUT) :: global
   end subroutine assemble_global
end interface

call MPI_INIT(ierror)
call decomp_2d_init(nx,ny,nz,p_row,p_col)
call decomp_2d_fft_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! (1) Testing the complex-to-complex interface (c2c) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  input is X-pencil data
! output is Z-pencil data
allocate (in(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
allocate (out(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Following is the testing input for NAG library C06FXF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
in1(1,1,1) = (1.000,  0.000)
in1(1,1,2) = (0.999, -0.040)
in1(1,1,3) = (0.987, -0.159)
in1(1,1,4) = (0.936, -0.352)
in1(1,2,1) = (0.994, -0.111)
in1(1,2,2) = (0.989, -0.151)
in1(1,2,3) = (0.963, -0.268)
in1(1,2,4) = (0.891, -0.454)
in1(1,3,1) = (0.903, -0.430)
in1(1,3,2) = (0.885, -0.466)
in1(1,3,3) = (0.823, -0.568)
in1(1,3,4) = (0.694, -0.720)
in1(2,1,1) = (0.500,  0.500)
in1(2,1,2) = (0.499,  0.040)
in1(2,1,3) = (0.487,  0.159)
in1(2,1,4) = (0.436,  0.352)
in1(2,2,1) = (0.494,  0.111)
in1(2,2,2) = (0.489,  0.151)
in1(2,2,3) = (0.463,  0.268)
in1(2,2,4) = (0.391,  0.454)
in1(2,3,1) = (0.403,  0.430)
in1(2,3,2) = (0.385,  0.466)
in1(2,3,3) = (0.323,  0.568)
in1(2,3,4) = (0.194,  0.720)

! each processor gets its local portion of global data
do k=xstart(3),xend(3)
   do j=xstart(2),xend(2)
      do i=xstart(1),xend(1)
         in(i,j,k) = in1(i,j,k)
      end do
   end do
end do

! write out input, to match the format of NAG example result file
if (nrank==0) then
   write(*,*) 'C06FXF Example Program Results'
   write(*,*) ''
   write(*,*) 'Original data values'
   write(*,*) ''
   call print_global(in1,nx,ny,nz)
end if

! ===== 3D forward FFT =====
call decomp_2d_fft_3d(in, out, DECOMP_2D_FFT_FORWARD)

! normalisation - note FFTW doesn't normalise 
do k=zstart(3),zend(3)
   do j=zstart(2),zend(2)
      do i=zstart(1),zend(1)
         out(i,j,k) = out(i,j,k) / sqrt(real(nx*ny*nz))
      end do
   end do
end do

call assemble_global(3,out,out1,nx,ny,nz)

! write out forward FFT result
if (nrank==0) then
   write(*,*) 'Components of discrete Fourier transform'
   write(*,*) ''
   call print_global(out1,nx,ny,nz)
end if

! ===== 3D inverse FFT =====
call decomp_2d_fft_3d(out, in, DECOMP_2D_FFT_BACKWARD)

! normalisation - note FFTW doesn't normalise
do k=xstart(3),xend(3)
   do j=xstart(2),xend(2)
      do i=xstart(1),xend(1)
         in(i,j,k) = in(i,j,k) / sqrt(real(nx*ny*nz))
      end do
   end do
end do

call assemble_global(1,in,in1,nx,ny,nz)

! write out inverse FFT result
if (nrank==0) then
   write(*,*) 'Original sequence as restored by inverse transform'
   write(*,*) ''
   call print_global(in1,nx,ny,nz)
end if

call decomp_2d_fft_finalize
call decomp_2d_finalize
call MPI_FINALIZE(ierror)

end program fft_test_c2c


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Collect data from each processor and assemble into a global array
! at the master rank
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine assemble_global(ndir,local,global,nx,ny,nz)
  
  use decomp_2d
  use MPI
  
  implicit none
  
  integer, intent(IN) :: ndir  ! 1 = X-pencil; 3 = Z-pencil
  integer, intent(IN) :: nx,ny,nz
  complex(mytype), dimension(:,:,:), intent(IN) :: local
  complex(mytype), dimension(nx,ny,nz), intent(OUT) :: global
  
  complex(mytype), allocatable, dimension(:,:,:) :: rbuf
  integer, dimension(9) :: sbuf1, rbuf1
  
  integer :: ierror, i,j,k,m, i1,i2,j1,j2,k1,k2, count
  integer, dimension(MPI_STATUS_SIZE) :: status
  
  if (nrank==0) then
     ! master writes its own data to a global array
     if (ndir==3) then  ! Z-pencil 
        i1 = zstart(1)
        i2 = zend(1)
        j1 = zstart(2)
        j2 = zend(2)
        k1 = zstart(3)
        k2 = zend(3)
     else if (ndir==1) then  ! X-pencil
        i1 = xstart(1)
        i2 = xend(1)
        j1 = xstart(2)
        j2 = xend(2)
        k1 = xstart(3)
        k2 = xend(3)
     end if
     do k=k1,k2
        do j=j1,j2
           do i=i1,i2
              ! 'local' is assumbed shape array
              ! but it is OK as starting index for rank 0 always 1
              global(i,j,k)=local(i,j,k)
           end do
        end do
     end do
     ! then loop through all other ranks to collect data
     do m=1,nproc-1
        CALL MPI_RECV(rbuf1,9,MPI_INTEGER,m,m,MPI_COMM_WORLD, &
             status,ierror)
        allocate(rbuf(rbuf1(1):rbuf1(2),rbuf1(4):rbuf1(5), &
             rbuf1(7):rbuf1(8)))
        CALL MPI_RECV(rbuf,rbuf1(3)*rbuf1(6)*rbuf1(9),complex_type,m, &
             m+nproc,MPI_COMM_WORLD,status,ierror)
        do k=rbuf1(7),rbuf1(8)
           do j=rbuf1(4),rbuf1(5)
              do i=rbuf1(1),rbuf1(2)
                 global(i,j,k)=rbuf(i,j,k)
              end do
           end do
        end do
        deallocate(rbuf)
     end do
  else
     ! slaves send data to mater
     if (ndir==3) then  ! Z-pencil
        sbuf1(1) = zstart(1)
        sbuf1(2) = zend(1)
        sbuf1(3) = zsize(1)
        sbuf1(4) = zstart(2)
        sbuf1(5) = zend(2)
        sbuf1(6) = zsize(2)
        sbuf1(7) = zstart(3)
        sbuf1(8) = zend(3)
        sbuf1(9) = zsize(3)
        count = zsize(1)*zsize(2)*zsize(3)
     else if (ndir==1) then  ! X-pencil
        sbuf1(1) = xstart(1)
        sbuf1(2) = xend(1)
        sbuf1(3) = xsize(1)
        sbuf1(4) = xstart(2)
        sbuf1(5) = xend(2)
        sbuf1(6) = xsize(2)
        sbuf1(7) = xstart(3)
        sbuf1(8) = xend(3)
        sbuf1(9) = xsize(3)
        count = xsize(1)*xsize(2)*xsize(3)
     end if
     ! send partition information
     CALL MPI_SEND(sbuf1,9,MPI_INTEGER,0,nrank,MPI_COMM_WORLD,ierror)
     ! send data array
     CALL MPI_SEND(local,count,complex_type,0, &
          nrank+nproc,MPI_COMM_WORLD,ierror)
  end if
  
  return
end subroutine assemble_global


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Print out a global data array using special format that matches
! NAG library C06FXF Example Program Results for validation purpose
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine print_global(data,nx,ny,nz)

use decomp_2d

implicit none

integer, intent(IN) :: nx,ny,nz
complex(mytype), dimension(nx,ny,nz), intent(IN) :: data

integer :: i,j,k

do i=1,nx
   write(*,10) i
   write(*,*) ''
   do j=1,ny
      write(*,20) (real(data(i,j,k)),k=1,nz)
      write(*,21) (aimag(data(i,j,k)),k=1,nz)
      write(*,*) ''
   end do
end do
10 format(1x,'z(i,j,k) for i =', I6)
20 format(1x,'Real ', 4F10.3)
21 format(1x,'Imag ', 4F10.3)

return
end subroutine print_global

