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
! 'update_halo_...' in halo.f90

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    if (present(opt_global)) then
       global = opt_global
    else
       global = .false.
    end if

    s1 = size(in,1)
    s2 = size(in,2)
    s3 = size(in,3)

    ! Calculate the starting index and ending index of output
    if (s1==decomp%xsz(1)) then  ! X-pencil input
       if (global) then
          xs = decomp%xst(1) 
          xe = decomp%xen(1)
          ys = decomp%xst(2) - level
          ye = decomp%xen(2) + level
          zs = decomp%xst(3) - level
          ze = decomp%xen(3) + level
       else
          xs = 1 
          xe = s1
          ys = 1 - level
          ye = s2 + level 
          zs = 1 - level
          ze = s3 + level
       end if
    else if (s2==decomp%ysz(2)) then  ! Y-pencil input
       if (global) then
          xs = decomp%yst(1) - level 
          xe = decomp%yen(1) + level
          ys = decomp%yst(2)
          ye = decomp%yen(2)
          zs = decomp%yst(3) - level
          ze = decomp%yen(3) + level
       else
          xs = 1 - level
          xe = s1 + level
          ys = 1
          ye = s2
          zs = 1 - level
          ze = s3 + level
       end if
    else if (s3==decomp%zsz(3)) then  ! Z-pencil input
       if (global) then
          xs = decomp%zst(1) - level 
          xe = decomp%zen(1) + level
          ys = decomp%zst(2) - level
          ye = decomp%zen(2) + level
          zs = decomp%zst(3)
          ze = decomp%zen(3)
       else
          xs = 1 - level
          xe = s1 + level
          ys = 1 - level
          ye = s2 + level
          zs = 1
          ze = s3
       end if
    else
       ! invalid input
       call decomp_2d_abort(10, &
            'Invalid data passed to update_halo')
    end if


    allocate(out(xs:xe, ys:ye, zs:ze))
!    out = -1.0_mytype ! fill the halo for debugging

    ! copy input data to output
    if (global) then
       ! using global coordinate
       ! note the input array passed in always has index starting from 1 
       ! need to work out the corresponding global index
       if (s1==decomp%xsz(1)) then
          do k=decomp%xst(3),decomp%xen(3)
             do j=decomp%xst(2),decomp%xen(2)
                do i=1,s1  ! x all local
                   out(i,j,k) = in(i,j-decomp%xst(2)+1,k-decomp%xst(3)+1)
                end do
             end do
          end do
       else if (s2==decomp%ysz(2)) then
          do k=decomp%yst(3),decomp%yen(3)
             do j=1,s2  ! y all local
                do i=decomp%yst(1),decomp%yen(1)
                   out(i,j,k) = in(i-decomp%yst(1)+1,j,k-decomp%yst(3)+1)
                end do
             end do
          end do
       else if (s3==decomp%zsz(3)) then
          do k=1,s3  ! z all local
             do j=decomp%zst(2),decomp%zen(2)
                do i=decomp%zst(1),decomp%zen(1)
                   out(i,j,k) = in(i-decomp%zst(1)+1,j-decomp%zst(2)+1,k)
                end do
             end do
          end do
       end if
    else
       ! not using global coordinate
       do k=1,s3
          do j=1,s2
             do i=1,s1
                out(i,j,k) = in(i,j,k)
             end do
          end do
       end do
    end if

    ! If needed, define MPI derived data type to pack halo data,
    ! then call MPI send/receive to exchange halo data
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! X-pencil
    if (s1==decomp%xsz(1)) then

#ifdef HALO_DEBUG
       if (nrank==4) then
          write(*,*) 'X-pencil input'
          write(*,*) '=============='
          write(*,*) 'Data on a y-z plane is shown'
          write(*,*) 'Before halo exchange'
          do j=ye,ys,-1
             write(*,'(10F4.0)') (out(1,j,k),k=zs,ze)
          end do
       end if
#endif

       ! *** east/west ***
       ! all data in local memory already, no halo exchange

       ! *** north/south *** 
       tag_s = coord(1)
       if (coord(1)==dims(1)-1 .AND. periodic_y) then
          tag_n = 0
       else
          tag_n = coord(1) + 1
       end if
       icount = s3 + 2*level
       ilength = level * s1
       ijump = s1*(s2+2*level)
       call MPI_TYPE_VECTOR(icount, ilength, ijump, &
            data_type, halo12, ierror)
       call MPI_TYPE_COMMIT(halo12, ierror)
       ! receive from south
       call MPI_IRECV(out(xs,ys,zs), 1, halo12, &
            neighbour(1,4), tag_s, DECOMP_2D_COMM_CART_X, &
            requests(1), ierror)
       ! receive from north
       call MPI_IRECV(out(xs,ye-level+1,zs), 1, halo12, &
            neighbour(1,3), tag_n, DECOMP_2D_COMM_CART_X, &
            requests(2), ierror)
       ! send to south
       call MPI_ISSEND(out(xs,ys+level,zs), 1, halo12, &
            neighbour(1,4), tag_s, DECOMP_2D_COMM_CART_X, &
            requests(3), ierror)
       ! send to north
       call MPI_ISSEND(out(xs,ye-level-level+1,zs), 1, halo12, &
            neighbour(1,3), tag_n, DECOMP_2D_COMM_CART_X, &
            requests(4), ierror)
       call MPI_WAITALL(4, requests, status, ierror)
       call MPI_TYPE_FREE(halo12, ierror)
#ifdef HALO_DEBUG
       if (nrank==4) then
          write(*,*) 'After exchange in Y'
          do j=ye,ys,-1
             write(*,'(10F4.0)') (out(1,j,k),k=zs,ze)
          end do
       end if
#endif

       ! *** top/bottom ***
       ! no need to define derived data type as data on xy-planes
       ! all contiguous in memory, which can be sent/received using
       ! MPI directly
       tag_b = coord(2)
       if (coord(2)==dims(2)-1 .AND. periodic_z) then
          tag_t = 0
       else
          tag_t = coord(2) + 1
       end if
       icount = (s1 * (s2+2*level)) * level
       ! receive from bottom
       call MPI_IRECV(out(xs,ys,zs), icount, data_type, &
            neighbour(1,6), tag_b, DECOMP_2D_COMM_CART_X, &
            requests(1), ierror)
       ! receive from top
       call MPI_IRECV(out(xs,ys,ze-level+1), icount, data_type, &
            neighbour(1,5), tag_t, DECOMP_2D_COMM_CART_X, &
            requests(2), ierror)
       ! send to bottom
       call MPI_ISSEND(out(xs,ys,zs+level), icount, data_type, &
            neighbour(1,6), tag_b, DECOMP_2D_COMM_CART_X, &
            requests(3), ierror)
       ! send to top
       call MPI_ISSEND(out(xs,ys,ze-level-level+1), icount, data_type, &
            neighbour(1,5), tag_t, DECOMP_2D_COMM_CART_X, &
            requests(4), ierror)
       call MPI_WAITALL(4, requests, status, ierror)
#ifdef HALO_DEBUG       
       if (nrank==4) then
          write(*,*) 'After exchange in Z'
          do j=ye,ys,-1
             write(*,'(10F4.0)') (out(1,j,k),k=zs,ze)
          end do
       end if
#endif       

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    ! Y-pencil   
    else if (s2==decomp%ysz(2)) then

#ifdef HALO_DEBUG
       if (nrank==4) then
          write(*,*) 'Y-pencil input'
          write(*,*) '=============='
          write(*,*) 'Data on a x-z plane is shown'
          write(*,*) 'Before halo exchange'
          do i=xe,xs,-1
             write(*,'(10F4.0)') (out(i,1,k),k=zs,ze)
          end do
       end if
#endif

       ! *** east/west ***
       tag_w = coord(1)
       if (coord(1)==dims(1)-1 .AND. periodic_x) then
          tag_e = 0
       else
          tag_e = coord(1) + 1
       end if
       icount = s2*(s3+2*level)
       ilength = level
       ijump = s1+2*level
       call MPI_TYPE_VECTOR(icount, ilength, ijump, &
            data_type, halo21, ierror)
       call MPI_TYPE_COMMIT(halo21, ierror)
       ! receive from west
       call MPI_IRECV(out(xs,ys,zs), 1, halo21, &
            neighbour(2,2), tag_w, DECOMP_2D_COMM_CART_Y, &
            requests(1), ierror)
       ! receive from east
       call MPI_IRECV(out(xe-level+1,ys,zs), 1, halo21, &
            neighbour(2,1), tag_e, DECOMP_2D_COMM_CART_Y, &
            requests(2), ierror)
       ! send to west
       call MPI_ISSEND(out(xs+level,ys,zs), 1, halo21, &
            neighbour(2,2), tag_w, DECOMP_2D_COMM_CART_Y, &
            requests(3), ierror)
       ! send to east
       call MPI_ISSEND(out(xe-level-level+1,ys,zs), 1, halo21, &
            neighbour(2,1), tag_e, DECOMP_2D_COMM_CART_Y, &
            requests(4), ierror)
       call MPI_WAITALL(4, requests, status, ierror)
       call MPI_TYPE_FREE(halo21, ierror)
#ifdef HALO_DEBUG       
       if (nrank==4) then
          write(*,*) 'After exchange in X'
          do i=xe,xs,-1
             write(*,'(10F4.0)') (out(i,1,k),k=zs,ze)
          end do
       end if
#endif

       ! *** north/south ***
       ! all data in local memory already, no halo exchange

       ! *** top/bottom ***
       ! no need to define derived data type as data on xy-planes
       ! all contiguous in memory, which can be sent/received using
       ! MPI directly
       tag_b = coord(2)
       if (coord(2)==dims(2)-1 .AND. periodic_z) then
          tag_t = 0
       else
          tag_t = coord(2) + 1
       end if
       icount = (s2 * (s1+2*level)) * level
       ! receive from bottom
       call MPI_IRECV(out(xs,ys,zs), icount, data_type, &
            neighbour(2,6), tag_b, DECOMP_2D_COMM_CART_Y, &
            requests(1), ierror)
       ! receive from top
       call MPI_IRECV(out(xs,ys,ze-level+1), icount, data_type, &
            neighbour(2,5), tag_t, DECOMP_2D_COMM_CART_Y, &
            requests(2), ierror)
       ! send to bottom
       call MPI_ISSEND(out(xs,ys,zs+level), icount, data_type, &
            neighbour(2,6), tag_b, DECOMP_2D_COMM_CART_Y, &
            requests(3), ierror)
       ! send to top
       call MPI_ISSEND(out(xs,ys,ze-level-level+1), icount, data_type, &
            neighbour(2,5), tag_t, DECOMP_2D_COMM_CART_Y, &
            requests(4), ierror)
       call MPI_WAITALL(4, requests, status, ierror)
#ifdef HALO_DEBUG
       if (nrank==4) then
          write(*,*) 'After exchange in Z'
          do i=xe,xs,-1
             write(*,'(10F4.0)') (out(i,1,k),k=zs,ze)
          end do
       end if
#endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    ! Z-pencil
    else if (s3==decomp%zsz(3)) then   

#ifdef HALO_DEBUG
       if (nrank==4) then
          write(*,*) 'Z-pencil input'
          write(*,*) '=============='
          write(*,*) 'Data on a x-y plane is shown'
          write(*,*) 'Before halo exchange'
          do i=xe,xs,-1
             write(*,'(10F4.0)') (out(i,j,1),j=ys,ye)
          end do
       end if
#endif

       ! *** east/west ***
       tag_w = coord(1)
       if (coord(1)==dims(1)-1 .AND. periodic_x) then
          tag_e = 0
       else
          tag_e = coord(1) + 1
       end if
       icount = (s2+2*level)*s3
       ilength = level
       ijump = s1+2*level
       call MPI_TYPE_VECTOR(icount, ilength, ijump, &
            data_type, halo31, ierror)
       call MPI_TYPE_COMMIT(halo31, ierror)
       ! receive from west
       call MPI_IRECV(out(xs,ys,zs), 1, halo31, &
            neighbour(3,2), tag_w, DECOMP_2D_COMM_CART_Z, &
            requests(1), ierror)
       ! receive from east
       call MPI_IRECV(out(xe-level+1,ys,zs), 1, halo31, &
            neighbour(3,1), tag_e, DECOMP_2D_COMM_CART_Z, &
            requests(2), ierror)
       ! send to west
       call MPI_ISSEND(out(xs+level,ys,zs), 1, halo31, &
            neighbour(3,2), tag_w, DECOMP_2D_COMM_CART_Z, &
            requests(3), ierror)
       ! send to east
       call MPI_ISSEND(out(xe-level-level+1,ys,zs), 1, halo31, &
            neighbour(3,1), tag_e, DECOMP_2D_COMM_CART_Z, &
            requests(4), ierror)
       call MPI_WAITALL(4, requests, status, ierror)
       call MPI_TYPE_FREE(halo31, ierror)       
#ifdef HALO_DEBUG
       if (nrank==4) then
          write(*,*) 'After exchange in X'
          do i=xe,xs,-1
             write(*,'(10F4.0)') (out(i,j,1),j=ys,ye)
          end do
       end if
#endif

       ! *** north/south *** 
       tag_s = coord(2)
       if (coord(2)==dims(2)-1 .AND. periodic_y) then
          tag_n = 0
       else
          tag_n = coord(2) + 1
       end if
       icount = s3
       ilength = level * (s1+2*level)
       ijump = (s1+2*level) * (s2+2*level)
       call MPI_TYPE_VECTOR(icount, ilength, ijump, &
            data_type, halo32, ierror)
       call MPI_TYPE_COMMIT(halo32, ierror)
       ! receive from south
       call MPI_IRECV(out(xs,ys,zs), 1, halo32, &
            neighbour(3,4), tag_s, DECOMP_2D_COMM_CART_Z, &
            requests(1), ierror)
       ! receive from north
       call MPI_IRECV(out(xs,ye-level+1,zs), 1, halo32, &
            neighbour(3,3), tag_n, DECOMP_2D_COMM_CART_Z, &
            requests(2), ierror)
       ! send to south
       call MPI_ISSEND(out(xs,ys+level,zs), 1, halo32, &
            neighbour(3,4), tag_s, DECOMP_2D_COMM_CART_Z, &
            requests(3), ierror)
       ! send to north
       call MPI_ISSEND(out(xs,ye-level-level+1,zs), 1, halo32, &
            neighbour(3,3), tag_n, DECOMP_2D_COMM_CART_Z, &
            requests(4), ierror)
       call MPI_WAITALL(4, requests, status, ierror)
       call MPI_TYPE_FREE(halo32, ierror)
#ifdef HALO_DEBUG
       if (nrank==4) then
          write(*,*) 'After exchange in Y'
          do i=xe,xs,-1
             write(*,'(10F4.0)') (out(i,j,1),j=ys,ye)
          end do
       end if
#endif

       ! *** top/bottom ***
       ! all data in local memory already, no halo exchange

    end if  ! pencil
