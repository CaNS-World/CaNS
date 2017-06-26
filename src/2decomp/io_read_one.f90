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
! 'mpiio_read_one_...' in io.f90

    ! Using MPI-IO to write a distributed 3D array into a file

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       call get_decomp_info(decomp)
    end if
    
    ! determine subarray parameters
    sizes(1) = decomp%xsz(1)
    sizes(2) = decomp%ysz(2)
    sizes(3) = decomp%zsz(3)
    
    if (ipencil == 1) then
       subsizes(1) = decomp%xsz(1)
       subsizes(2) = decomp%xsz(2)
       subsizes(3) = decomp%xsz(3)
       starts(1) = decomp%xst(1)-1  ! 0-based index
       starts(2) = decomp%xst(2)-1
       starts(3) = decomp%xst(3)-1
    else if (ipencil == 2) then
       subsizes(1) = decomp%ysz(1)
       subsizes(2) = decomp%ysz(2)
       subsizes(3) = decomp%ysz(3)
       starts(1) = decomp%yst(1)-1
       starts(2) = decomp%yst(2)-1
       starts(3) = decomp%yst(3)-1
    else if (ipencil == 3) then
       subsizes(1) = decomp%zsz(1)
       subsizes(2) = decomp%zsz(2)
       subsizes(3) = decomp%zsz(3)
       starts(1) = decomp%zst(1)-1
       starts(2) = decomp%zst(2)-1
       starts(3) = decomp%zst(3)-1
    endif
    
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, data_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, &
         fh, ierror)
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_VIEW(fh,disp,data_type, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         data_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
