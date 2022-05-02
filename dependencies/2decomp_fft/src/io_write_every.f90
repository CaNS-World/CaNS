!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2021 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This file contain common code to be included by subroutines 
! 'write_every_...' in io.f90

  ! To write every few points of a 3D array to a file

    ! work out the distribution parameters, which may be different from 
    ! the default distribution used by the decomposition library
    !  For exmample if nx=17 and p_row=4
    !    distribution is: 4 4 4 5

    ! If writing from the 1st element
    !  If saving every 3 points, then 5 points to be saved (17/3)
    !    default distribution would be 1 1 1 2
    !    However, 1st block (1-4) contains the 3rd point
    !             2nd block (5-8) contains the 6th point
    !             3rd block (9-12) contains the 9th and 12th point
    !             4th block (13-17) contains then 15th point
    !    giving a 1 1 2 1 distribution
    !    So cannot use the base decomposition library for such IO

    ! If writing from the n-th element (n=?skip)
    !  If saving every 3 points, then 6 points to be saved
    !    However, 1st block (1-4) contains the 1st & 4th point
    !             2nd block (5-8) contains the 7th point
    !             3rd block (9-12) contains the 10th point
    !             4th block (13-17) contains then 12th & 15th point
    !    giving a 1 2 2 1 distribution

    skip(1)=iskip
    skip(2)=jskip
    skip(3)=kskip

    do i=1,3
       if (from1) then
          xst(i) = (xstart(i)+skip(i)-1)/skip(i)
          if (mod(xstart(i)+skip(i)-1,skip(i))/=0) xst(i)=xst(i)+1
          xen(i) = (xend(i)+skip(i)-1)/skip(i)
       else
          xst(i) = xstart(i)/skip(i)
          if (mod(xstart(i),skip(i))/=0) xst(i)=xst(i)+1
          xen(i) = xend(i)/skip(i)
       end if
       xsz(i) = xen(i)-xst(i)+1
    end do
       
    do i=1,3
       if (from1) then
          yst(i) = (ystart(i)+skip(i)-1)/skip(i)
          if (mod(ystart(i)+skip(i)-1,skip(i))/=0) yst(i)=yst(i)+1
          yen(i) = (yend(i)+skip(i)-1)/skip(i)
       else
          yst(i) = ystart(i)/skip(i)
          if (mod(ystart(i),skip(i))/=0) yst(i)=yst(i)+1
          yen(i) = yend(i)/skip(i)
       end if
       ysz(i) = yen(i)-yst(i)+1
    end do

    do i=1,3
       if (from1) then
          zst(i) = (zstart(i)+skip(i)-1)/skip(i)
          if (mod(zstart(i)+skip(i)-1,skip(i))/=0) zst(i)=zst(i)+1
          zen(i) = (zend(i)+skip(i)-1)/skip(i)
       else
          zst(i) = zstart(i)/skip(i)
          if (mod(zstart(i),skip(i))/=0) zst(i)=zst(i)+1
          zen(i) = zend(i)/skip(i)
       end if
       zsz(i) = zen(i)-zst(i)+1
    end do

    ! if 'skip' value is large it is possible that some ranks do not 
    ! contain any points to be written. Subarray constructor requires 
    ! nonzero size so it is not possible to use MPI_COMM_WORLD for IO.
    ! Create a sub communicator for this...
    color = 1
    key = 0  ! rank order doesn't matter
    if (ipencil==1) then
       if (xsz(1)==0 .or. xsz(2)==0 .or. xsz(3)==0) then
          color = 2
       end if
    else if (ipencil==2) then
       if (ysz(1)==0 .or. ysz(2)==0 .or. ysz(3)==0) then
          color = 2
       end if
    else if (ipencil==3) then
       if (zsz(1)==0 .or. zsz(2)==0 .or. zsz(3)==0) then
          color = 2
       end if
    end if
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,newcomm,ierror)

    if (color==1) then ! only ranks in this group do IO collectively
       
       ! generate subarray information
       sizes(1) = xsz(1)
       sizes(2) = ysz(2)
       sizes(3) = zsz(3)
       if (ipencil==1) then
          subsizes(1) = xsz(1)
          subsizes(2) = xsz(2)
          subsizes(3) = xsz(3)
          starts(1) = xst(1)-1
          starts(2) = xst(2)-1
          starts(3) = xst(3)-1
       else if (ipencil==2) then
          subsizes(1) = ysz(1)
          subsizes(2) = ysz(2)
          subsizes(3) = ysz(3)
          starts(1) = yst(1)-1
          starts(2) = yst(2)-1
          starts(3) = yst(3)-1
       else if (ipencil==3) then
          subsizes(1) = zsz(1)
          subsizes(2) = zsz(2)
          subsizes(3) = zsz(3)
          starts(1) = zst(1)-1
          starts(2) = zst(2)-1
          starts(3) = zst(3)-1
       end if
       
       ! copy data from original array
       ! needs a copy of original array in global coordinate 
       if (ipencil==1) then
          allocate(wk(xst(1):xen(1),xst(2):xen(2),xst(3):xen(3)))
          allocate(wk2(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
          wk2=var
          if (from1) then
             do k=xst(3),xen(3)
                do j=xst(2),xen(2)
                   do i=xst(1),xen(1)
                      wk(i,j,k) = wk2((i-1)*iskip+1,(j-1)*jskip+1,(k-1)*kskip+1)
                   end do
                end do
             end do
          else
             do k=xst(3),xen(3)
                do j=xst(2),xen(2)
                   do i=xst(1),xen(1)
                      wk(i,j,k) = wk2(i*iskip,j*jskip,k*kskip)
                   end do
                end do
             end do
          end if   
       else if (ipencil==2) then
          allocate(wk(yst(1):yen(1),yst(2):yen(2),yst(3):yen(3)))
          allocate(wk2(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
          wk2=var
          if (from1) then
             do k=yst(3),yen(3)
                do j=yst(2),yen(2)
                   do i=yst(1),yen(1)
                      wk(i,j,k) = wk2((i-1)*iskip+1,(j-1)*jskip+1,(k-1)*kskip+1)
                   end do
                end do
             end do
          else
             do k=yst(3),yen(3)
                do j=yst(2),yen(2)
                   do i=yst(1),yen(1)
                      wk(i,j,k) = wk2(i*iskip,j*jskip,k*kskip)
                   end do
                end do
             end do
          end if
       else if (ipencil==3) then
          allocate(wk(zst(1):zen(1),zst(2):zen(2),zst(3):zen(3)))
          allocate(wk2(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
          wk2=var
          if (from1) then
             do k=zst(3),zen(3)
                do j=zst(2),zen(2)
                   do i=zst(1),zen(1)
                      wk(i,j,k) = wk2((i-1)*iskip+1,(j-1)*jskip+1,(k-1)*kskip+1)
                   end do
                end do
             end do
          else
             do k=zst(3),zen(3)
                do j=zst(2),zen(2)
                   do i=zst(1),zen(1)
                      wk(i,j,k) = wk2(i*iskip,j*jskip,k*kskip)
                   end do
                end do
             end do
          end if
       end if
       deallocate(wk2)

       ! MPI-IO
       call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
            MPI_ORDER_FORTRAN, data_type, newtype, ierror)
       call MPI_TYPE_COMMIT(newtype,ierror)
       call MPI_FILE_OPEN(newcomm, filename, &
            MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
            fh, ierror)
       filesize = 0_MPI_OFFSET_KIND
       call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
       disp = 0_MPI_OFFSET_KIND
       call MPI_FILE_SET_VIEW(fh,disp,data_type, &
            newtype,'native',MPI_INFO_NULL,ierror)
       call MPI_FILE_WRITE_ALL(fh, wk, &
            subsizes(1)*subsizes(2)*subsizes(3), &
            data_type, MPI_STATUS_IGNORE, ierror)
       call MPI_FILE_CLOSE(fh,ierror)
       call MPI_TYPE_FREE(newtype,ierror)

       deallocate(wk)

    end if ! color==1

    call MPI_BARRIER(MPI_COMM_WORLD, ierror)
