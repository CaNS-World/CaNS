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

! This file contains the routines that transpose data from X to Y pencil

  subroutine transpose_x_to_y_real(src, dst, opt_decomp)

    implicit none
    
    real(mytype), dimension(:,:,:), intent(IN) :: src
    real(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

#ifdef SHM
    real(mytype) :: work1(*), work2(*)
    POINTER  (work1_p, work1), (work2_p, work2)  ! Cray pointers
#endif
    
    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)

    ! rearrange source array as send buffer
#ifdef SHM
    work1_p = decomp%COL_INFO%SND_P
    call mem_split_xy_real(src, s1, s2, s3, work1, dims(1), &
         decomp%x1dist, decomp)
#else
    call mem_split_xy_real(src, s1, s2, s3, work1_r, dims(1), &
         decomp%x1dist, decomp)
#endif

    ! define receive buffer
#ifdef SHM
    work2_p = decomp%COL_INFO%RCV_P
    call MPI_BARRIER(decomp%COL_INFO%CORE_COMM, ierror)
#endif
    
    ! transpose using MPI_ALLTOALL(V)
#ifdef SHM
    if (decomp%COL_INFO%CORE_ME==1) THEN
       call MPI_ALLTOALLV(work1, decomp%x1cnts_s, decomp%x1disp_s, &
            real_type, work2, decomp%y1cnts_s, decomp%y1disp_s, &
            real_type, decomp%COL_INFO%SMP_COMM, ierror)
    end if
#else
#ifdef EVEN
    call MPI_ALLTOALL(work1_r, decomp%x1count, &
         real_type, work2_r, decomp%y1count, &
         real_type, DECOMP_2D_COMM_COL, ierror)
#else
    call MPI_ALLTOALLV(work1_r, decomp%x1cnts, decomp%x1disp, &
         real_type, work2_r, decomp%y1cnts, decomp%y1disp, &
         real_type, DECOMP_2D_COMM_COL, ierror)
#endif
#endif

    ! rearrange receive buffer
#ifdef SHM
    call MPI_BARRIER(decomp%COL_INFO%CORE_COMM, ierror)
    call mem_merge_xy_real(work2, d1, d2, d3, dst, dims(1), &
         decomp%y1dist, decomp)
#else
    call mem_merge_xy_real(work2_r, d1, d2, d3, dst, dims(1), &
         decomp%y1dist, decomp)
#endif
    
    return
  end subroutine transpose_x_to_y_real


#ifdef OCC
  subroutine transpose_x_to_y_real_start(handle, src, dst, sbuf, rbuf, &
       opt_decomp)

    implicit none
    
    integer :: handle
    real(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

    integer :: s1,s2,s3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)

    ! rearrange source array as send buffer
    call mem_split_xy_real(src, s1, s2, s3, sbuf, dims(1), &
         decomp%x1dist, decomp)
    
#ifdef EVEN
    call NBC_IALLTOALL(sbuf, decomp%x1count, real_type, &
         rbuf, decomp%y1count, real_type, &
         DECOMP_2D_COMM_COL, handle, ierror)
#else
    call NBC_IALLTOALLV(sbuf, decomp%x1cnts, decomp%x1disp, real_type, &
         rbuf, decomp%y1cnts, decomp%y1disp, real_type, &
         DECOMP_2D_COMM_COL, handle, ierror)
#endif

    return
  end subroutine transpose_x_to_y_real_start


  subroutine transpose_x_to_y_real_wait(handle, src, dst, sbuf, rbuf, &
       opt_decomp)

    implicit none
    
    integer :: handle
    real(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

    integer :: d1,d2,d3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)

    call NBC_WAIT(handle, ierror)

    ! rearrange receive buffer
    call mem_merge_xy_real(rbuf, d1, d2, d3, dst, dims(1), &
         decomp%y1dist, decomp)

    return
  end subroutine transpose_x_to_y_real_wait
#endif


  subroutine transpose_x_to_y_complex(src, dst, opt_decomp)

    implicit none
    
    complex(mytype), dimension(:,:,:), intent(IN) :: src
    complex(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

#ifdef SHM
    complex(mytype) :: work1(*), work2(*)
    POINTER  (work1_p, work1), (work2_p, work2)  ! Cray pointers
#endif
    
    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)
    
    ! rearrange source array as send buffer
#ifdef SHM
    work1_p = decomp%COL_INFO%SND_P_c
    call mem_split_xy_complex(src, s1, s2, s3, work1, dims(1), &
         decomp%x1dist, decomp)
#else
    call mem_split_xy_complex(src, s1, s2, s3, work1_c, dims(1), &
         decomp%x1dist, decomp)
#endif
    
    ! define receive buffer
#ifdef SHM
    work2_p = decomp%COL_INFO%RCV_P_c
    call MPI_BARRIER(decomp%COL_INFO%CORE_COMM, ierror)
#endif
    
    ! transpose using MPI_ALLTOALL(V)
#ifdef SHM
    if (decomp%COL_INFO%CORE_ME==1) THEN
       call MPI_ALLTOALLV(work1, decomp%x1cnts_s, decomp%x1disp_s, &
            complex_type, work2, decomp%y1cnts_s, decomp%y1disp_s, &
            complex_type, decomp%COL_INFO%SMP_COMM, ierror)
    end if
#else
#ifdef EVEN
    call MPI_ALLTOALL(work1_c, decomp%x1count, &
         complex_type, work2_c, decomp%y1count, &
         complex_type, DECOMP_2D_COMM_COL, ierror)
#else
    call MPI_ALLTOALLV(work1_c, decomp%x1cnts, decomp%x1disp, &
         complex_type, work2_c, decomp%y1cnts, decomp%y1disp, &
         complex_type, DECOMP_2D_COMM_COL, ierror)
#endif
#endif

    ! rearrange receive buffer
#ifdef SHM
    call MPI_BARRIER(decomp%COL_INFO%CORE_COMM, ierror)
    call mem_merge_xy_complex(work2, d1, d2, d3, dst, dims(1), &
         decomp%y1dist, decomp)
#else
    call mem_merge_xy_complex(work2_c, d1, d2, d3, dst, dims(1), &
         decomp%y1dist, decomp)
#endif

    return
  end subroutine transpose_x_to_y_complex


#ifdef OCC
  subroutine transpose_x_to_y_complex_start(handle, src, dst, sbuf, &
       rbuf, opt_decomp)

    implicit none
    
    integer :: handle
    complex(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

    integer :: s1,s2,s3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    
    ! rearrange source array as send buffer
    call mem_split_xy_complex(src, s1, s2, s3, sbuf, dims(1), &
         decomp%x1dist, decomp)

#ifdef EVEN
    call NBC_IALLTOALL(sbuf, decomp%x1count, &
         complex_type, rbuf, decomp%y1count, &
         complex_type, DECOMP_2D_COMM_COL, handle, ierror)
#else
    call NBC_IALLTOALLV(sbuf, decomp%x1cnts, decomp%x1disp, &
         complex_type, rbuf, decomp%y1cnts, decomp%y1disp, &
         complex_type, DECOMP_2D_COMM_COL, handle, ierror)
#endif

    return
  end subroutine transpose_x_to_y_complex_start


  subroutine transpose_x_to_y_complex_wait(handle, src, dst, sbuf, &
       rbuf, opt_decomp)

    implicit none
    
    integer :: handle
    complex(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

    integer :: d1,d2,d3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)

    call NBC_WAIT(handle, ierror)

    ! rearrange receive buffer
    call mem_merge_xy_complex(rbuf, d1, d2, d3, dst, dims(1), &
         decomp%y1dist, decomp)

    return
  end subroutine transpose_x_to_y_complex_wait
#endif


  ! pack/unpack ALLTOALL(V) buffers

  subroutine mem_split_xy_real(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none

    integer, intent(IN) :: n1,n2,n3
    real(mytype), dimension(n1,n2,n3), intent(IN) :: in
    real(mytype), dimension(*), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2,pos

    do m=0,iproc-1
       if (m==0) then 
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

#ifdef SHM
       pos = decomp%x1disp_o(m) + 1
#else
#ifdef EVEN
       pos = m * decomp%x1count + 1
#else
       pos = decomp%x1disp(m) + 1
#endif
#endif

       do k=1,n3
          do j=1,n2
             do i=i1,i2
                out(pos) = in(i,j,k)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_split_xy_real


  subroutine mem_split_xy_complex(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none

    integer, intent(IN) :: n1,n2,n3
    complex(mytype), dimension(n1,n2,n3), intent(IN) :: in
    complex(mytype), dimension(*), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2,pos

    do m=0,iproc-1
       if (m==0) then 
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

#ifdef SHM
       pos = decomp%x1disp_o(m) + 1
#else
#ifdef EVEN
       pos = m * decomp%x1count + 1
#else
       pos = decomp%x1disp(m) + 1
#endif
#endif

       do k=1,n3
          do j=1,n2
             do i=i1,i2
                out(pos) = in(i,j,k)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_split_xy_complex


  subroutine mem_merge_xy_real(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none
    
    integer, intent(IN) :: n1,n2,n3
    real(mytype), dimension(*), intent(IN) :: in
    real(mytype), dimension(n1,n2,n3), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2, pos

    do m=0,iproc-1
       if (m==0) then
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

#ifdef SHM
       pos = decomp%y1disp_o(m) + 1
#else
#ifdef EVEN
       pos = m * decomp%y1count + 1
#else
       pos = decomp%y1disp(m) + 1
#endif
#endif

       do k=1,n3
          do j=i1,i2
             do i=1,n1
                out(i,j,k) = in(pos)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_merge_xy_real


  subroutine mem_merge_xy_complex(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none
    
    integer, intent(IN) :: n1,n2,n3
    complex(mytype), dimension(*), intent(IN) :: in
    complex(mytype), dimension(n1,n2,n3), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2, pos

    do m=0,iproc-1
       if (m==0) then
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

#ifdef SHM
       pos = decomp%y1disp_o(m) + 1
#else
#ifdef EVEN
       pos = m * decomp%y1count + 1
#else
       pos = decomp%y1disp(m) + 1
#endif
#endif

       do k=1,n3
          do j=i1,i2
             do i=1,n1
                out(i,j,k) = in(pos)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_merge_xy_complex
