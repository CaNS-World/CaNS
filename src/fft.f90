! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_fft
  use, intrinsic :: iso_c_binding , only: C_INT
  use mod_common_mpi, only: ierr
  use mod_fftw_param
  use mod_types
#if defined(_OPENACC)
  use mod_utils     , only: f_sizeof
#endif
  !$ use omp_lib
  private
  public fftini,fftend,fft
  !@acc public signal_processing,fftf_gpu,fftb_gpu
  !@acc integer(i8), public :: wsize_fft
  contains
  subroutine fftini(ng,n_x,n_y,bcxy,c_or_f,arrplan,normfft)
    implicit none
    integer, intent(in), dimension(3) :: ng,n_x,n_y
    character(len=1), intent(in), dimension(0:1,2) :: bcxy
    character(len=1), intent(in), dimension(2) :: c_or_f
#if !defined(_OPENACC)
    type(C_PTR), intent(out), dimension(2,2) :: arrplan
#else
    integer    , intent(out), dimension(2,2) :: arrplan
#endif
    real(rp), intent(out) :: normfft
    real(gp), dimension(n_x(1),n_x(2),n_x(3))  :: arrx
    real(gp), dimension(n_y(1),n_y(2),n_y(3))  :: arry
#if !defined(_OPENACC)
    type(C_PTR) :: plan_fwd_x,plan_bwd_x, &
                   plan_fwd_y,plan_bwd_y
    type(fftw_iodim), dimension(1) :: iodim
    type(fftw_iodim), dimension(2) :: iodim_howmany
#else
    integer     :: plan_fwd_x,plan_bwd_x, &
                   plan_fwd_y,plan_bwd_y
#endif
    integer :: kind_fwd,kind_bwd
    real(rp), dimension(2) :: norm
    integer(C_INT) :: nx_x,ny_x,nz_x, &
                      nx_y,ny_y,nz_y
    integer :: ix,iy
#if defined(_OPENACC)
    integer :: istat,batch
    integer(int_ptr_kind()) :: wsize,max_wsize
#endif
#if defined(_SINGLE_PRECISION) || defined(_SINGLE_PRECISION_POISSON)
    !$ call sfftw_init_threads(ierr)
    !$ call sfftw_plan_with_nthreads(omp_get_max_threads())
#else
    !$ call dfftw_init_threads(ierr)
    !$ call dfftw_plan_with_nthreads(omp_get_max_threads())
#endif
#if !defined(_OPENACC)
    nx_x = n_x(1)
    ny_x = n_x(2)
    nz_x = n_x(3)
    nx_y = n_y(1)
    ny_y = n_y(2)
    nz_y = n_y(3)
#endif
    normfft = 1.
    !
    ! fft in x
    !
    call find_fft(bcxy(:,1),c_or_f(1),kind_fwd,kind_bwd,norm)
    ix = 0
    ! size of transform reduced by 1 point with Dirichlet BC in face
    if(bcxy(0,1)//bcxy(1,1) == 'DD'.and.c_or_f(1) == 'f') ix = 1
#if !defined(_OPENACC)
    !
    ! prepare plans with guru interface
    !
    iodim(1)%n  = nx_x-ix
    iodim(1)%is = 1
    iodim(1)%os = 1
    iodim_howmany(1)%n  = ny_x
    iodim_howmany(1)%is = nx_x
    iodim_howmany(1)%os = nx_x
    iodim_howmany(2)%n  = nz_x
    iodim_howmany(2)%is = nx_x*ny_x
    iodim_howmany(2)%os = nx_x*ny_x
    plan_fwd_x=fftw_plan_guru_r2r(1,iodim,2,iodim_howmany,arrx,arrx,kind_fwd,FFTW_ESTIMATE)
    plan_bwd_x=fftw_plan_guru_r2r(1,iodim,2,iodim_howmany,arrx,arrx,kind_bwd,FFTW_ESTIMATE)
#else
    max_wsize = -1
    batch = product(n_x(2:3)) ! padded & axis-contiguous layout
    nx_x = n_x(1)             ! padded & axis-contiguous layout
    istat = cufftCreate(plan_fwd_x)
    istat = cufftSetAutoAllocation(plan_fwd_x,0)
    istat = cufftMakePlanMany(plan_fwd_x,1,ng(1),null(),1,nx_x,null(),1,ng(1)/2+1,CUFFT_FWD_TYPE,batch,wsize)
    max_wsize = max(wsize,max_wsize)
    !
    istat = cufftCreate(plan_bwd_x)
    istat = cufftSetAutoAllocation(plan_bwd_x,0)
    istat = cufftMakePlanMany(plan_bwd_x,1,ng(1),null(),1,ng(1)/2+1,null(),1,nx_x,CUFFT_BWD_TYPE,batch,wsize)
    max_wsize = max(wsize,max_wsize)
#endif
    normfft = normfft*norm(1)*(ng(1)+norm(2)-ix)
    !
    ! fft in y
    !
    call find_fft(bcxy(:,2),c_or_f(2),kind_fwd,kind_bwd,norm)
    iy = 0
    ! size of transform reduced by 1 point with Dirichlet BC in face
    if(bcxy(0,2)//bcxy(1,2) == 'DD'.and.c_or_f(2) == 'f') iy = 1
#if !defined(_OPENACC)
    !
    ! prepare plans with guru interface
    !
    iodim(1)%n  = ny_y-iy
    iodim(1)%is = nx_y
    iodim(1)%os = nx_y
    iodim_howmany(1)%n  = nx_y
    iodim_howmany(1)%is = 1
    iodim_howmany(1)%os = 1
    iodim_howmany(2)%n  = nz_y
    iodim_howmany(2)%is = nx_y*ny_y
    iodim_howmany(2)%os = nx_y*ny_y
    plan_fwd_y=fftw_plan_guru_r2r(1,iodim,2,iodim_howmany,arry,arry,kind_fwd,FFTW_ESTIMATE)
    plan_bwd_y=fftw_plan_guru_r2r(1,iodim,2,iodim_howmany,arry,arry,kind_bwd,FFTW_ESTIMATE)
#else
    batch = product(n_y(2:3)) ! padded & axis-contiguous layout
    ny_y = n_y(1)             ! padded & axis-contiguous layout
    istat = cufftCreate(plan_fwd_y)
    istat = cufftSetAutoAllocation(plan_fwd_y,0)
    istat = cufftMakePlanMany(plan_fwd_y,1,ng(2),null(),1,ny_y,null(),1,ng(2)/2+1,CUFFT_FWD_TYPE,batch,wsize)
    max_wsize = max(wsize,max_wsize)
    !
    istat = cufftCreate(plan_bwd_y)
    istat = cufftSetAutoAllocation(plan_bwd_y,0)
    istat = cufftMakePlanMany(plan_bwd_y,1,ng(2),null(),1,ng(2)/2+1,null(),1,ny_y,CUFFT_BWD_TYPE,batch,wsize)
    max_wsize = max(wsize,max_wsize)
    wsize_fft = max_wsize/f_sizeof(1._gp)
#endif
    normfft = normfft*norm(1)*(ng(2)+norm(2)-iy)
    !
    arrplan(1,1) = plan_fwd_x
    arrplan(2,1) = plan_bwd_x
    arrplan(1,2) = plan_fwd_y
    arrplan(2,2) = plan_bwd_y
    normfft = normfft**(-1)
  end subroutine fftini
  !
  subroutine fftend(arrplan)
    implicit none
#if !defined(_OPENACC)
    type(C_PTR), intent(in), dimension(:,:) :: arrplan
#else
    integer    , intent(in), dimension(:,:) :: arrplan
#endif
    integer :: i,j
    !@acc integer :: istat
#if !defined(_OPENACC)
#if defined(_SINGLE_PRECISION) || defined(_SINGLE_PRECISION_POISSON)
    do j=1,size(arrplan,2)
      do i=1,size(arrplan,1)
        call sfftw_destroy_plan(arrplan(i,j))
      end do
    end do
    !$ call sfftw_cleanup_threads(ierr)
#else
    do j=1,size(arrplan,2)
      do i=1,size(arrplan,1)
        call dfftw_destroy_plan(arrplan(i,j))
      end do
    end do
    !$ call dfftw_cleanup_threads(ierr)
#endif
#else
    do j=1,size(arrplan,2)
      do i=1,size(arrplan,1)
        istat = cufftDestroy(arrplan(i,j))
      end do
    end do
#endif
  end subroutine fftend
  !
  subroutine fft(plan,arr)
    implicit none
    type(C_PTR), intent(in) :: plan
    real(gp), intent(inout), dimension(:,:,:) :: arr
#if !defined(_OPENACC)
#if defined(_SINGLE_PRECISION) || defined(_SINGLE_PRECISION_POISSON)
    call sfftw_execute_r2r(plan,arr,arr)
#else
    call dfftw_execute_r2r(plan,arr,arr)
#endif
#endif
  end subroutine fft
  !
  subroutine find_fft(bc,c_or_f,kind_fwd,kind_bwd,norm)
  implicit none
  character(len=1), intent(in), dimension(0:1) :: bc
  character(len=1), intent(in) :: c_or_f
  integer , intent(out) :: kind_fwd,kind_bwd
  real(rp), intent(out), dimension(2) :: norm
  if(c_or_f == 'c') then
    select case(bc(0)//bc(1))
    case('PP')
      kind_fwd = FFTW_R2HC
      kind_bwd = FFTW_HC2R
      norm = [1.,0.]
    case('NN')
      kind_fwd = FFTW_REDFT10
      kind_bwd = FFTW_REDFT01
      norm = [2.,0.]
    case('DD')
      kind_fwd = FFTW_RODFT10
      kind_bwd = FFTW_RODFT01
      norm = [2.,0.]
    case('ND')
      kind_fwd = FFTW_REDFT11
      kind_bwd = FFTW_REDFT11
      norm = [2.,0.]
    case('DN')
      kind_fwd = FFTW_RODFT11
      kind_bwd = FFTW_RODFT11
      norm = [2.,0.]
    end select
  else if(c_or_f == 'f') then
    select case(bc(0)//bc(1))
    case('PP')
      kind_fwd = FFTW_R2HC
      kind_bwd = FFTW_HC2R
      norm = [1.,0.]
    case('NN')
      kind_fwd = FFTW_REDFT00
      kind_bwd = FFTW_REDFT00
      norm = [2.,-1.]
    case('DD')
      kind_fwd = FFTW_RODFT00
      kind_bwd = FFTW_RODFT00
      norm = [2.,1.]
    case('ND')
      kind_fwd = FFTW_REDFT10
      kind_bwd = FFTW_REDFT01
      norm = [2.,0.]
    case('DN')
      kind_fwd = FFTW_RODFT01
      kind_bwd = FFTW_RODFT10
      norm = [2.,0.]
    end select
  end if
  end subroutine find_fft
#if defined(_OPENACC)
  subroutine fftf_gpu(plan,arr)
    implicit none
    integer , intent(in) :: plan
    real(gp), intent(inout), dimension(:,:,:) :: arr
    integer :: istat
    !$acc host_data use_device(arr)
#if defined(_SINGLE_PRECISION) || defined(_SINGLE_PRECISION_POISSON)
    istat = cufftExecR2C(plan,arr,arr)
#else
    istat = cufftExecD2Z(plan,arr,arr)
#endif
    !$acc end host_data
  end subroutine fftf_gpu
  subroutine fftb_gpu(plan,arr)
    implicit none
    integer , intent(in) :: plan
    real(gp), intent(inout), dimension(:,:,:) :: arr
    integer :: istat
    !$acc host_data use_device(arr)
#if defined(_SINGLE_PRECISION) || defined(_SINGLE_PRECISION_POISSON)
    istat = cufftExecC2R(plan,arr,arr)
#else
    istat = cufftExecZ2D(plan,arr,arr)
#endif
    !$acc end host_data
  end subroutine fftb_gpu
  subroutine posp_fftf(nn,n,idir,arr)
    !
    ! post-processing of a signal following a forward FFT
    ! to order the data as follows:
    ! (r[0],r[n],r[1],i[1],...,r[n-1],i[n-1])
    !
    implicit none
    integer , intent(in   ) :: nn
    integer , intent(in   ), dimension(3) :: n       ! dimensions of input/output array
    integer , intent(in   ) :: idir                  ! direction where the transform is taken
    real(gp), intent(inout), dimension(:,:,:) :: arr ! input/output array
    integer :: j,k,n_2,n_3
    !
    select case(idir)
    case(1)
      n_2 = n(2); n_3 = n(3)
      !$acc parallel loop collapse(2) default(present) async(1)
      do k=1,n_3
        do j=1,n_2
          arr(2,j,k) = arr(nn+1,j,k)
        end do
      end do
    end select
  end subroutine posp_fftf
  subroutine prep_fftb(nn,n,idir,arr)
    !
    ! pre-processing of a signal preciding a backward FFT
    ! to order the data as follows:
    ! (r[0],i[0],r[1],i[1],...,r[n-1],i[n-1],r[n],i[n])
    ! note that i[0] = i[n] = 0
    !
    implicit none
    integer , intent(in   ) :: nn
    integer , intent(in   ), dimension(3) :: n       ! dimensions of input/output array
    integer , intent(in   ) :: idir                  ! direction where the transform is taken
    real(gp), intent(inout), dimension(:,:,:) :: arr ! input/output array
    integer :: j,k,n_2,n_3
    !
    select case(idir)
    case(1)
      n_2 = n(2); n_3 = n(3)
      !$acc parallel loop collapse(2) default(present) async(1)
      do k=1,n_3
        do j=1,n_2
          arr(nn+1,j,k) = arr(2,j,k)
          arr(2   ,j,k) = 0.
        end do
      end do
    end select
  end subroutine prep_fftb
  subroutine prep_dctiif(nn,n,idir,arr,is_swap_order,is_negate_even)
    !
    ! pre-processing of a signal to perform a fast forward
    ! discrete cosine transform (DCT) with FFTs (see Makhoul 1980)
    !
    ! the input signal x(n) is pre-processed into a signal v(n)
    ! as follows (now done using the subroutine remap with ib=0):
    !
    ! v(n) = x(2n       ),              0 <= n <= floor((N-1)/2)
    !      = x(2N -2n -1), floor((N+1)/2) <= n <= N-1
    ! with n = 0,...,N-1 and N is the total number of elements
    ! of the signal.
    !
    ! pre-processing required for computing the corresponding
    ! discrete sine transform (DST) may also be performed
    ! if one of the last boolean input variables is .true.
    !
    implicit none
    integer , intent(in   ) :: nn
    integer , intent(in   ), dimension(3) :: n       ! dimensions of input/output array
    integer , intent(in   ) :: idir                  ! array direction where the transform is taken
    real(gp), intent(inout), dimension(:,:,:) :: arr ! input/output array
    logical, intent(in) :: is_swap_order  ! swap order of the elements of the input array? (for DST)
    logical, intent(in) :: is_negate_even ! negate every other element of the input array?
    integer :: i,j,k
    !
    select case(idir)
    case(1)
      if(is_swap_order ) call swap_order( nn,n(2),n(3),arr)
      if(is_negate_even) call negate_even(nn,n(2),n(3),arr)
      call remap(0,nn,n(2),n(3),arr)
    end select
  end subroutine prep_dctiif
  subroutine posp_dctiif(nn,n,idir,arr,is_swap_order,is_negate_even)
    !
    ! post-processing of a signal to perform a fast forward discrete
    ! cosine transform with FFTs (see Makhoul 1980)
    !
    ! post-processing required for computing the corresponding
    ! discrete sine transform (DST) may also be performed
    ! if one of the last boolean input variables is .true.
    !
    use mod_param          , only: pi_rp => pi
    use mod_common_cudecomp, only: buf => work
    implicit none
    integer , intent(in   ) :: nn
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   ) :: idir
    real(gp), intent(inout), dimension(:,:,:) :: arr
    logical, intent(in) :: is_swap_order  ! swap order of the elements of the input array? (for DST)
    logical, intent(in) :: is_negate_even ! negate every other element of the input array?
    real(gp), pointer, contiguous, dimension(:,:,:) :: arr_tmp
    integer :: i,j,k,ii
    real(gp) :: arg
    integer :: n_2,n_3
    real(gp) :: pi
    pi = pi_rp ! converts double to single if needed
    !
    select case(idir)
    case(1)
      arr_tmp(0:n(idir)-1,1:n(2),1:n(3)) => buf(1:product(n(:)))
      n_2 = n(2); n_3 = n(3)
      !$acc parallel loop collapse(3) default(present) private(i,arg) async(1)
      do k=1,n_3
        do j=1,n_2
          do ii=0,nn/2
            i = 2*ii+1
            !arr_tmp(ii   ,j,k) =    real( &
            !                         2.*exp(-ri_unit*pi*ii/(2.*nn))*cmplx(arr(i,j,k),arr(i+1,j,k),gp) &
            !                        )
            !arr_tmp(nn-ii,j,k) = - aimag( &
            !                         2.*exp(-ri_unit*pi*ii/(2.*nn))*cmplx(arr(i,j,k),arr(i+1,j,k),gp) &
            !                        ) ! = 0 for ii=0
            arg = -pi*ii/(2.*nn)
            arr_tmp(ii   ,j,k) =  2.*(cos(arg)*arr(i,j,k) - sin(arg)*arr(i+1,j,k))
            arr_tmp(nn-ii,j,k) = -2.*(sin(arg)*arr(i,j,k) + cos(arg)*arr(i+1,j,k))
          end do
        end do
      end do
      !$acc kernels default(present) async(1)
      arr(:,:,:) = arr_tmp(:,:,:)
      !$acc end kernels
      if(is_swap_order ) call swap_order( nn,n(2),n(3),arr)
      if(is_negate_even) call negate_even(nn,n(2),n(3),arr)
    end select
  end subroutine posp_dctiif
  subroutine prep_dctiib(nn,n,idir,arr,is_swap_order,is_negate_even)
    !
    ! pre-processing of a signal to perform a fast backward
    ! discrete cosine transform (DST) with FFTs (see Makhoul 1980)
    !
    ! pre-processing required for computing the corresponding
    ! discrete sine transform (DST) may also be performed
    ! if one of the last boolean input variables is .true.
    !
    use mod_param          , only: pi_rp => pi
    use mod_common_cudecomp, only: buf => work
    implicit none
    integer , intent(in   ) :: nn
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   ) :: idir
    real(gp), intent(inout), dimension(:,:,:) :: arr
    logical, intent(in) :: is_swap_order  ! swap order of the elements of the input array? (for DST)
    logical, intent(in) :: is_negate_even ! negate every other element of the input array?
    real(gp), pointer, contiguous, dimension(:,:,:) :: arr_tmp
    integer :: i,j,k,ii
    real(gp) :: arg
    integer :: n_2,n_3
    real(gp) :: pi
    pi = pi_rp ! converts double to single if needed
    !
    select case(idir)
    case(1)
      if(is_swap_order ) call swap_order( nn,n(2),n(3),arr)
      if(is_negate_even) call negate_even(nn,n(2),n(3),arr)
      !$acc parallel loop collapse(2) default(present) async(1)
      do k=1,n_3
        do j=1,n_2
          do ii = 1,2-mod(nn,2)
            arr(nn+ii,j,k) = 0.
          end do
        end do
      end do
      arr_tmp(0:n(idir)-1,1:n(2),1:n(3)) => buf(1:product(n(:)))
      n_2 = n(2); n_3 = n(3)
      !$acc parallel loop collapse(3) default(present) private(ii,arg) async(1)
      do k=1,n_3
        do j=1,n_2
          do ii=0,nn/2
            !arr_tmp(2*ii  ,j,k)  = real( 1.*exp(ri_unit*pi*ii/(2.*nn))*(arr(ii+1,j,k)-ri_unit*arr(nn-ii+1,j,k)),gp)
            !arr_tmp(2*ii+1,j,k)  = aimag(1.*exp(ri_unit*pi*ii/(2.*nn))*(arr(ii+1,j,k)-ri_unit*arr(nn-ii+1,j,k)),gp)
            arg = pi*ii/(2.*nn)
            arr_tmp(2*ii  ,j,k) = 1.*(cos(arg)*arr(ii+1,j,k) + sin(arg)*arr(nn-ii+1,j,k))
            arr_tmp(2*ii+1,j,k) = 1.*(sin(arg)*arr(ii+1,j,k) - cos(arg)*arr(nn-ii+1,j,k))
          end do
        end do
      end do
      !$acc kernels default(present) async(1)
      arr(:,:,:) = arr_tmp(:,:,:)
      !$acc end kernels
    end select
  end subroutine prep_dctiib
  subroutine posp_dctiib(nn,n,idir,arr,is_swap_order,is_negate_even)
    !
    ! post-processing of a signal to perform a fast forward
    ! discrete cosine transform (DCT) with FFTs (see Makhoul 1980)
    !
    ! the input signal v(n) is post-processed into a signal x(n)
    ! as follows (now done using the subroutine remap with ib=1):
    !
    ! v(n) = x(2n       ),              0 <= n <= floor((N-1)/2)
    !      = x(2N -2n -1), floor((N+1)/2) <= n <= N-1
    ! with n = 0,...,N-1 and N is the total number of elements
    ! of the signal.
    !
    ! post-processing required for computing the corresponding
    ! discrete sine transform (DST) may also be performed
    ! if one of the last boolean input variables is .true.
    !
    implicit none
    integer , intent(in   ) :: nn
    integer , intent(in   ), dimension(3) :: n       ! dimensions of input/output array
    integer , intent(in   ) :: idir                  ! array direction where the transform is taken
    real(gp), intent(inout), dimension(:,:,:) :: arr ! input/output array
    logical, intent(in) :: is_swap_order  ! swap order of the elements of the input array? (for DST)
    logical, intent(in) :: is_negate_even ! negate every other element of the input array?
    integer :: i,j,k
    !
    select case(idir)
    case(1)
      call remap(1,nn,n(2),n(3),arr)
      if(is_swap_order ) call swap_order( nn,n(2),n(3),arr)
      if(is_negate_even) call negate_even(nn,n(2),n(3),arr)
    end select
  end subroutine posp_dctiib
  subroutine signal_processing(pre_or_pos,f_or_b,cbc,c_or_f,nn,n,idir,arr)
    implicit none
    !
    ! wrapper subroutine for signal processing to compute FFT-based transforms
    ! (can also be done with pointers to a subroutine like in initgrid.f90)
    !
    integer,          intent(in) :: pre_or_pos ! prior (0) or after (1) fft
    character(len=1), intent(in) :: f_or_b     ! forward or backward transform
    character(len=2), intent(in) :: cbc        ! type of boundary condition
    character(len=1), intent(in) :: c_or_f     ! cell- or face-centred BC?
    integer, intent(in)                       :: nn ! number of points in the signal
    integer, intent(in), dimension(3)         :: n
    integer, intent(in)                       :: idir
    real(gp), intent(inout), dimension(:,:,:) :: arr
    integer :: istat
    select case(cbc)
    case('PP')
      select case(f_or_b)
      case('F')
        if(     pre_or_pos == 0) then
          return
        else if(pre_or_pos == 1) then
          call posp_fftf(nn,n,idir,arr)
        else
        end if
      case('B')
        if(     pre_or_pos == 0) then
          call prep_fftb(nn,n,idir,arr)
        else if(pre_or_pos == 1) then
          return
        else
        end if
      end select
    case('NN')
      if(c_or_f == 'c') then
        select case(f_or_b)
        case('F')
          if(     pre_or_pos == 0) then
            call prep_dctiif(nn,n,idir,arr,.false.,.false.)
          else if(pre_or_pos == 1) then
            call posp_dctiif(nn,n,idir,arr,.false.,.false.)
          else
          end if
        case('B')
          if(    pre_or_pos == 0) then
            call prep_dctiib(nn,n,idir,arr,.false.,.false.)
          else if(pre_or_pos == 1) then
            call posp_dctiib(nn,n,idir,arr,.false.,.false.)
          else
          end if
        case default
        end select
      end if
    case('DD')
      if(c_or_f == 'c') then
        select case(f_or_b)
        case('F')
          if(     pre_or_pos == 0) then
            call prep_dctiif(nn,n,idir,arr,.false.,.true. )
          else if(pre_or_pos == 1) then
            call posp_dctiif(nn,n,idir,arr,.true. ,.false.)
          else
          end if
        case('B')
          if(     pre_or_pos == 0) then
            call prep_dctiib(nn,n,idir,arr,.true. ,.false.)
          else if(pre_or_pos == 1) then
            call posp_dctiib(nn,n,idir,arr,.false.,.true. )
          else
          end if
        case default
        end select
      end if
    case default
      ! ERROR  (trap this in sanity check in sanity.f90)
    end select
  end subroutine signal_processing
  subroutine negate_even(n,n2,n3,arr)
    implicit none
    integer , intent(in   ) :: n,n2,n3
    real(gp), intent(inout) :: arr(:,:,:)
    integer :: i,j,k
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,n3
      do j=1,n2
        do i=1,n/2
          arr(2*i,j,k) = - arr(2*i,j,k)
        end do
      end do
    end do
  end subroutine negate_even
  subroutine swap_order(n,n2,n3,arr)
    implicit none
    integer , intent(in   ) :: n,n2,n3
    real(gp), intent(inout) :: arr(:,:,:)
    real(gp) :: tmp
    integer  :: i,j,k
    !$acc parallel loop collapse(3) default(present) private(tmp) async(1)
    do k=1,n3
      do j=1,n2
        do i=1,n/2
          tmp            = arr(i    ,j,k)
          arr(i    ,j,k) = arr(n-i+1,j,k)
          arr(n-i+1,j,k) = tmp
        end do
      end do
    end do
  end subroutine swap_order
  subroutine remap(ib,n,n2,n3,arr)
    !
    ! maps a signal x to v (ib = 0), or v to x (ib = 1)
    ! where:
    ! v(n) = x(2n       ),              0 <= n <= floor((N-1)/2)
    !      = x(2N -2n -1), floor((N+1)/2) <= n <= N-1
    ! with n = 0,...,N-1; N = size(v) = size(x)
    !
    use mod_common_cudecomp, only: buf => work
    implicit none
    integer , intent(in   ) :: ib,n,n2,n3
    real(gp), intent(inout) :: arr(:,:,:)
    real(gp), pointer, contiguous :: arr_tmp(:,:,:)
    integer :: i,j,k
    integer :: nh
    arr_tmp(1:n,1:n2,1:n3) => buf(1:n*n2*n3)
    nh = (n+1)/2
    select case(ib)
    case(0)
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1,n3
        do j=1,n2
          do i=1,nh
            arr_tmp(i   ,j,k) = arr(2*i-1        ,j,k)
          end do
        end do
      end do
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1,n3
        do j=1,n2
          do i=1,n/2
            arr_tmp(i+nh,j,k) = arr(2*(n-(i+nh)+1),j,k)
          end do
        end do
      end do
      !$acc kernels default(present) async(1)
      arr(:,:,:) = arr_tmp(:,:,:)
      !$acc end kernels
    case(1)
      !$acc kernels default(present) async(1)
      arr_tmp(:,:,:) = arr(:,:,:)
      !$acc end kernels
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1,n3
        do j=1,n2
          do i=1,nh
            arr(2*i-1         ,j,k) = arr_tmp(i   ,j,k)
          end do
        end do
      end do
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1,n3
        do j=1,n2
          do i=1,n/2
            arr(2*(n-(i+nh)+1),j,k) = arr_tmp(i+nh,j,k)
          end do
        end do
      end do
    end select
  end subroutine remap
#endif
end module mod_fft
