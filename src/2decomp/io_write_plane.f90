!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2012 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This file contain common code to be included by subroutines 
! 'mpiio_write_plane_3d_...' in io.f90

    ! It is much easier to implement if all mpi ranks participate I/O.
    ! Transpose the 3D data if necessary.

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       call get_decomp_info(decomp)
    end if

    if (iplane==1) then
       allocate(wk(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)))
       if (ipencil==1) then
          wk = var
       else if (ipencil==2) then
          call transpose_y_to_x(var,wk,decomp)
       else if (ipencil==3) then
          allocate(wk2(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))
          call transpose_z_to_y(var,wk2,decomp)
          call transpose_y_to_x(wk2,wk,decomp)
          deallocate(wk2)
       end if
       allocate(wk2d(1,decomp%xsz(2),decomp%xsz(3)))
       do k=1,decomp%xsz(3)
          do j=1,decomp%xsz(2)
             wk2d(1,j,k)=wk(n,j,k)
          end do
       end do
       sizes(1) = 1
       sizes(2) = decomp%ysz(2)
       sizes(3) = decomp%zsz(3)
       subsizes(1) = 1
       subsizes(2) = decomp%xsz(2)
       subsizes(3) = decomp%xsz(3)
       starts(1) = 0
       starts(2) = decomp%xst(2)-1
       starts(3) = decomp%xst(3)-1

    else if (iplane==2) then
       allocate(wk(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))
       if (ipencil==1) then
          call transpose_x_to_y(var,wk,decomp)
       else if (ipencil==2) then
          wk = var
       else if (ipencil==3) then
          call transpose_z_to_y(var,wk,decomp)
       end if
       allocate(wk2d(decomp%ysz(1),1,decomp%ysz(3)))
       do k=1,decomp%ysz(3)
          do i=1,decomp%ysz(1)
             wk2d(i,1,k)=wk(i,n,k)
          end do
       end do
       sizes(1) = decomp%xsz(1)
       sizes(2) = 1
       sizes(3) = decomp%zsz(3)
       subsizes(1) = decomp%ysz(1)
       subsizes(2) = 1
       subsizes(3) = decomp%ysz(3)
       starts(1) = decomp%yst(1)-1
       starts(2) = 0
       starts(3) = decomp%yst(3)-1

    else if (iplane==3) then
       allocate(wk(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3)))
       if (ipencil==1) then
          allocate(wk2(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))
          call transpose_x_to_y(var,wk2,decomp)
          call transpose_y_to_z(wk2,wk,decomp)
          deallocate(wk2)
       else if (ipencil==2) then
          call transpose_y_to_z(var,wk,decomp)
       else if (ipencil==3) then
          wk = var
       end if
       allocate(wk2d(decomp%zsz(1),decomp%zsz(2),1))
       do j=1,decomp%zsz(2)
          do i=1,decomp%zsz(1) 
             wk2d(i,j,1)=wk(i,j,n)
          end do
       end do
       sizes(1) = decomp%xsz(1)
       sizes(2) = decomp%ysz(2)
       sizes(3) = 1
       subsizes(1) = decomp%zsz(1)
       subsizes(2) = decomp%zsz(2)
       subsizes(3) = 1
       starts(1) = decomp%zst(1)-1
       starts(2) = decomp%zst(2)-1
       starts(3) = 0
    end if

    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, data_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
         fh, ierror)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_VIEW(fh,disp,data_type, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, wk2d, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         data_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    
    deallocate(wk,wk2d)
