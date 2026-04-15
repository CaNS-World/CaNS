! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_fft
  use, intrinsic :: iso_c_binding, only: C_INT,c_intptr_t,C_PTR,c_loc,c_null_ptr
  use mod_common_mpi, only: ierr
  use mod_fftw_param
  use mod_types
#if defined(_OPENACC)
  use mod_utils     , only: f_sizeof
#endif
  !$ use omp_lib
  private
  public fftini,fftend,fft
#if defined(_OPENACC)
  public fft_gpu
  integer(i8), public :: wsize_fft,wsize_tmp
  real(rp), allocatable, target :: sincos_theta_x(:,:),sincos_theta_y(:,:)
#endif
  contains
  subroutine fftini(ng,n_x,n_y,bcxy,c_or_f,arrplan,normfft)
    use mod_param, only: pi
    implicit none
    integer , target, intent(in), dimension(3) :: ng,n_x,n_y
    character(len=1), intent(in), dimension(0:1,2) :: bcxy
    character(len=1), intent(in), dimension(2) :: c_or_f
#if !defined(_OPENACC) || defined(_USE_HIP)
    type(C_PTR), intent(out), dimension(2,2) :: arrplan
#else
    integer    , intent(out), dimension(2,2) :: arrplan
#endif
    real(rp), intent(out) :: normfft
    real(rp), dimension(n_x(1),n_x(2),n_x(3))  :: arrx
    real(rp), dimension(n_y(1),n_y(2),n_y(3))  :: arry
#if !defined(_OPENACC)
    type(fftw_iodim), dimension(1) :: iodim
    type(fftw_iodim), dimension(2) :: iodim_howmany
#endif
#if !defined(_OPENACC) || defined(_USE_HIP)
    type(C_PTR) :: plan_fwd_x,plan_bwd_x, &
                   plan_fwd_y,plan_bwd_y
#else
    integer     :: plan_fwd_x,plan_bwd_x, &
                   plan_fwd_y,plan_bwd_y
#endif
    integer :: kind_fwd,kind_bwd
    real(rp), dimension(2) :: norm
    real(rp) :: theta
    integer(C_INT) :: nx_x,ny_x,nz_x, &
                      nx_y,ny_y,nz_y
    integer :: ix,iy
#if defined(_OPENACC)
    integer :: istat,batch,ii
    integer(c_intptr_t) :: wsize,max_wsize
#endif
#if defined(_SINGLE_PRECISION)
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
    !
    ! store sine/cosine values for real-to-real transforms on GPUs
    !
    ! this assumes that all variables have the same number of points along each direction;
    ! in the future, if real-to-real with (semi-) collocated boundary conditions are implemented
    ! (e.g., for implicit diffusion on GPUs), we need to initialize arrays for each solved variable
    ! and store it.
    !
    if(bcxy(0,1)//bcxy(1,1) /= 'PP') then
      if(.not.allocated(sincos_theta_x)) then
        allocate(sincos_theta_x(0:ng(1)/2,1:2))
        do ii=0,ng(1)/2
          theta = pi*ii/(2._rp*ng(1))
          sincos_theta_x(ii,1) = sin(theta)
          sincos_theta_x(ii,2) = cos(theta)
        end do
        !$acc enter data copyin(sincos_theta_x)
      end if
    end if
    batch = product(n_x(2:3)) ! padded & axis-contiguous layout
    nx_x = n_x(1)             ! padded & axis-contiguous layout
    istat = cufftCreate(plan_fwd_x)
    istat = cufftSetAutoAllocation(plan_fwd_x,0)
#if !defined(_USE_HIP)
    istat = cufftMakePlanMany(plan_fwd_x,1,ng(1),null(),1,nx_x,null(),1,ng(1)/2+1,CUFFT_FWD_TYPE,batch,wsize)
#else
    istat = cufftMakePlanMany(plan_fwd_x,1,c_loc(ng(1)),c_null_ptr,1,nx_x,c_null_ptr,1,ng(1)/2+1,CUFFT_FWD_TYPE,batch,c_loc(wsize))
#endif
    max_wsize = max(wsize,max_wsize)
    !
    istat = cufftCreate(plan_bwd_x)
    istat = cufftSetAutoAllocation(plan_bwd_x,0)
#if !defined(_USE_HIP)
    istat = cufftMakePlanMany(plan_bwd_x,1,ng(1),null(),1,ng(1)/2+1,null(),1,nx_x,CUFFT_BWD_TYPE,batch,wsize)
#else
    istat = cufftMakePlanMany(plan_bwd_x,1,c_loc(ng(1)),c_null_ptr,1,ng(1)/2+1,c_null_ptr,1,nx_x,CUFFT_BWD_TYPE,batch,c_loc(wsize))
#endif
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
    !
    ! store sine/cosine values for real-to-real transforms on GPUs
    !
    ! see the comment above on prospective extensions of real-to-real transforms for implicit diffusion
    !
    if(bcxy(0,2)//bcxy(1,2) /= 'PP') then
      if(.not.allocated(sincos_theta_y)) then
        allocate(sincos_theta_y(0:ng(2)/2,1:2))
        do ii=0,ng(2)/2
          theta = pi*ii/(2._rp*ng(2))
          sincos_theta_y(ii,1) = sin(theta)
          sincos_theta_y(ii,2) = cos(theta)
        end do
        !$acc enter data copyin(sincos_theta_y)
      end if
    end if
    batch = product(n_y(2:3)) ! padded & axis-contiguous layout
    ny_y = n_y(1)             ! padded & axis-contiguous layout
    istat = cufftCreate(plan_fwd_y)
    istat = cufftSetAutoAllocation(plan_fwd_y,0)
#if !defined(_USE_HIP)
    istat = cufftMakePlanMany(plan_fwd_y,1,ng(2),null(),1,ny_y,null(),1,ng(2)/2+1,CUFFT_FWD_TYPE,batch,wsize)
#else
    istat = cufftMakePlanMany(plan_fwd_y,1,c_loc(ng(2)),c_null_ptr,1,ny_y,c_null_ptr,1,ng(2)/2+1,CUFFT_FWD_TYPE,batch,c_loc(wsize))
#endif
    max_wsize = max(wsize,max_wsize)
    !
    istat = cufftCreate(plan_bwd_y)
    istat = cufftSetAutoAllocation(plan_bwd_y,0)
#if !defined(_USE_HIP)
    istat = cufftMakePlanMany(plan_bwd_y,1,ng(2),null(),1,ng(2)/2+1,null(),1,ny_y,CUFFT_BWD_TYPE,batch,wsize)
#else
    istat = cufftMakePlanMany(plan_bwd_y,1,c_loc(ng(2)),c_null_ptr,1,ng(2)/2+1,c_null_ptr,1,ny_y,CUFFT_BWD_TYPE,batch,c_loc(wsize))
#endif
    max_wsize = max(wsize,max_wsize)
    wsize_fft = max_wsize/f_sizeof(1._rp)
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
#if !defined(_OPENACC) || defined(_USE_HIP)
    type(C_PTR), intent(in), dimension(:,:) :: arrplan
#else
    integer    , intent(in), dimension(:,:) :: arrplan
#endif
    integer :: i,j
#if defined(_OPENACC)
    integer :: istat
#endif
#if !defined(_OPENACC)
#if defined(_SINGLE_PRECISION)
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
    real(rp), intent(inout), dimension(:,:,:) :: arr
#if !defined(_OPENACC)
#if defined(_SINGLE_PRECISION)
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
#if !defined(_USE_HIP)
    integer    , intent(in) :: plan
#else
    type(C_PTR), intent(in) :: plan
#endif
    real(rp), target, intent(inout), dimension(:,:,:) :: arr
    integer :: istat
    !$acc host_data use_device(arr)
#if !defined(_USE_HIP)
#if defined(_SINGLE_PRECISION)
    istat = cufftExecR2C(plan,arr,arr)
#else
    istat = cufftExecD2Z(plan,arr,arr)
#endif
#else
    !$acc wait(1)
#if defined(_SINGLE_PRECISION)
    istat = cufftExecR2C(plan,c_loc(arr),c_loc(arr))
#else
    istat = cufftExecD2Z(plan,c_loc(arr),c_loc(arr))
#endif
    istat = hipDeviceSynchronize()
#endif
    !$acc end host_data
  end subroutine fftf_gpu
  subroutine fftb_gpu(plan,arr)
    implicit none
#if !defined(_USE_HIP)
    integer    , intent(in) :: plan
#else
    type(C_PTR), intent(in) :: plan
#endif
    real(rp), target, intent(inout), dimension(:,:,:) :: arr
    integer :: istat
    !$acc host_data use_device(arr)
#if !defined(_USE_HIP)
#if defined(_SINGLE_PRECISION)
    istat = cufftExecC2R(plan,arr,arr)
#else
    istat = cufftExecZ2D(plan,arr,arr)
#endif
#else
    !$acc wait(1)
#if defined(_SINGLE_PRECISION)
    istat = cufftExecC2R(plan,c_loc(arr),c_loc(arr))
#else
    istat = cufftExecZ2D(plan,c_loc(arr),c_loc(arr))
#endif
    istat = hipDeviceSynchronize()
#endif
    !$acc end host_data
  end subroutine fftb_gpu
  subroutine fft_gpu(f_or_b,cbc,c_or_f,nn,n,plan,arr,arr_tmp)
    implicit none
    character(len=1), intent(in) :: f_or_b,c_or_f
    character(len=2), intent(in) :: cbc
    integer         , intent(in) :: nn
    integer         , intent(in), dimension(3) :: n
#if !defined(_USE_HIP)
    integer    , intent(in) :: plan
#else
    type(C_PTR), intent(in) :: plan
#endif
    real(rp), intent(inout), target, dimension(:,:,:) :: arr
    real(rp), intent(inout), target, dimension(:,:,:) :: arr_tmp
    if(c_or_f == 'c'.and.cbc /= 'PP') then
      call signal_processing(0,f_or_b,cbc,c_or_f,nn,n,1,arr,arr_tmp)
      select case(f_or_b)
      case('F')
        call fftf_gpu(plan,arr_tmp)
      case('B')
        call fftb_gpu(plan,arr_tmp)
      end select
      call signal_processing(1,f_or_b,cbc,c_or_f,nn,n,1,arr_tmp,arr)
    else
      call signal_processing(0,f_or_b,cbc,c_or_f,nn,n,1,arr)
      select case(f_or_b)
      case('F')
        call fftf_gpu(plan,arr)
      case('B')
        call fftb_gpu(plan,arr)
      end select
      call signal_processing(1,f_or_b,cbc,c_or_f,nn,n,1,arr)
    end if
  end subroutine fft_gpu
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
    real(rp), intent(inout), dimension(:,:,:) :: arr ! input/output array
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
    real(rp), intent(inout), dimension(:,:,:) :: arr ! input/output array
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
  subroutine prep_dctiif(nn,n,idir,arr,arr_out,is_swap_order,is_negate_even)
    !
    ! pre-processing of a signal to perform a fast forward
    ! discrete cosine transform (DCT) with FFTs (see Makhoul 1980)
    !
    ! the input signal x(n) is pre-processed into a signal v(n)
    ! as follows (now done using the subroutine remap with ib=0):
    !
    ! v(n) = x(2n       ),              0 <= n <= floor((N-1)/2)
    !      = x(2N -2n -1), floor((N+1)/2) <= n <= N-1
    ! with n = 0,...,N-1 and N being the total number of elements of the
    ! signal.
    !
    ! pre-processing required for computing the corresponding
    ! discrete sine transform (DST) may also be performed
    ! if one of the last boolean input variables is .true.
    !
    implicit none
    integer , intent(in   ) :: nn
    integer , intent(in   ), dimension(3) :: n            ! dimensions of input/output array
    integer , intent(in   ) :: idir                       ! array direction where the transform is taken
    real(rp), intent(inout), dimension(: ,:,:) :: arr     ! input/output array
    real(rp), intent(out  ), dimension(0:,:,:) :: arr_out ! output array
    logical, intent(in) :: is_swap_order  ! swap order of the elements of the input array? (for DST)
    logical, intent(in) :: is_negate_even ! negate every other element of the input array?
    integer :: i,j,k
    !
    select case(idir)
    case(1)
      if(is_swap_order ) call swap_order( nn,n(2),n(3),arr)
      if(is_negate_even) call negate_even(nn,n(2),n(3),arr)
      call remap(0,nn,n(2),n(3),arr,arr_out)
    end select
  end subroutine prep_dctiif
  subroutine posp_dctiif(nn,n,idir,arr,arr_out,is_swap_order,is_negate_even)
    !
    ! post-processing of a signal to perform a fast forward discrete
    ! cosine transform with FFTs (see Makhoul 1980)
    !
    ! post-processing required for computing the corresponding
    ! discrete sine transform (DST) may also be performed
    ! if one of the last boolean input variables is .true.
    !
    implicit none
    integer , intent(in   ) :: nn
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   ) :: idir
    real(rp), intent(in   ), dimension(: ,:,:) :: arr
    real(rp), intent(out  ), dimension(0:,:,:) :: arr_out
    logical, intent(in) :: is_swap_order  ! swap order of the elements of the input array? (for DST)
    logical, intent(in) :: is_negate_even ! negate every other element of the input array?
    integer :: i,j,k,ii
    real(rp), pointer, contiguous :: sin_theta(:),cos_theta(:)
    integer :: n_2,n_3
    !
    select case(idir)
    case(1)
      if(allocated(sincos_theta_x) .or. allocated(sincos_theta_y)) then
        if(allocated(sincos_theta_x)) then
          if(ubound(sincos_theta_x,1) == nn/2) then
            sin_theta(0:nn/2) => sincos_theta_x(:,1)
            cos_theta(0:nn/2) => sincos_theta_x(:,2)
          end if
        end if
        if(.not.associated(sin_theta) .and. allocated(sincos_theta_y)) then
          if(ubound(sincos_theta_y,1) == nn/2) then
            sin_theta(0:nn/2) => sincos_theta_y(:,1)
            cos_theta(0:nn/2) => sincos_theta_y(:,2)
          end if
        end if
      else
        error stop 'ERROR: sincos arrays were not computed.'
      end if
      n_2 = n(2); n_3 = n(3)
      !$acc parallel loop collapse(3) default(present) private(i) async(1)
      do k=1,n_3
        do j=1,n_2
          do ii=0,nn/2
            i = 2*ii+1
            !arr_out(ii   ,j,k) =    real( &
            !                         2.*exp(-ri_unit*pi*ii/(2.*nn))*cmplx(arr(i,j,k),arr(i+1,j,k),rp) &
            !                        )
            !arr_out(nn-ii,j,k) = - aimag( &
            !                         2.*exp(-ri_unit*pi*ii/(2.*nn))*cmplx(arr(i,j,k),arr(i+1,j,k),rp) &
            !                        ) ! = 0 for ii=0
            arr_out(ii   ,j,k) =  2.*(cos_theta(ii)*arr(i,j,k) + sin_theta(ii)*arr(i+1,j,k))
            arr_out(nn-ii,j,k) =  2.*(sin_theta(ii)*arr(i,j,k) - cos_theta(ii)*arr(i+1,j,k))
          end do
        end do
      end do
      if(is_swap_order ) call swap_order( nn,n(2),n(3),arr_out)
      if(is_negate_even) call negate_even(nn,n(2),n(3),arr_out)
    end select
  end subroutine posp_dctiif
  subroutine prep_dctiib(nn,n,idir,arr,arr_out,is_swap_order,is_negate_even)
    !
    ! pre-processing of a signal to perform a fast backward
    ! discrete cosine transform (DST) with FFTs (see Makhoul 1980)
    !
    ! pre-processing required for computing the corresponding
    ! discrete sine transform (DST) may also be performed
    ! if one of the last boolean input variables is .true.
    !
    implicit none
    integer , intent(in   ) :: nn
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   ) :: idir
    real(rp), intent(inout), dimension(: ,:,:) :: arr
    real(rp), intent(out  ), dimension(0:,:,:) :: arr_out
    logical, intent(in) :: is_swap_order  ! swap order of the elements of the input array? (for DST)
    logical, intent(in) :: is_negate_even ! negate every other element of the input array?
    integer :: i,j,k,ii
    real(rp), pointer, contiguous :: sin_theta(:),cos_theta(:)
    integer :: n_2,n_3
    !
    select case(idir)
    case(1)
      if(allocated(sincos_theta_x) .or. allocated(sincos_theta_y)) then
        if(allocated(sincos_theta_x)) then
          if(ubound(sincos_theta_x,1) == nn/2) then
            sin_theta(0:nn/2) => sincos_theta_x(:,1)
            cos_theta(0:nn/2) => sincos_theta_x(:,2)
          end if
        end if
        if(.not.associated(sin_theta) .and. allocated(sincos_theta_y)) then
          if(ubound(sincos_theta_y,1) == nn/2) then
            sin_theta(0:nn/2) => sincos_theta_y(:,1)
            cos_theta(0:nn/2) => sincos_theta_y(:,2)
          end if
        end if
      else
        error stop 'ERROR: sincos arrays were not computed.'
      end if
      if(is_swap_order ) call swap_order( nn,n(2),n(3),arr)
      if(is_negate_even) call negate_even(nn,n(2),n(3),arr)
      n_2 = n(2); n_3 = n(3)
      !$acc parallel loop collapse(2) default(present) async(1)
      do k=1,n_3
        do j=1,n_2
          do ii = 1,2-mod(nn,2)
            arr(nn+ii,j,k) = 0.
          end do
        end do
      end do
      !$acc parallel loop collapse(3) default(present) private(ii) async(1)
      do k=1,n_3
        do j=1,n_2
          do ii=0,nn/2
            !arr_out(2*ii  ,j,k)  = real( &
            !                         1.*exp(ri_unit*pi*ii/(2.*nn))*(arr(ii+1,j,k)-ri_unit*arr(nn-ii+1,j,k)),rp &
            !                        )
            !arr_out(2*ii+1,j,k)  = aimag( &
            !                         1.*exp(ri_unit*pi*ii/(2.*nn))*(arr(ii+1,j,k)-ri_unit*arr(nn-ii+1,j,k)),rp &
            !                        )
            arr_out(2*ii  ,j,k) = cos_theta(ii)*arr(ii+1,j,k) + sin_theta(ii)*arr(nn-ii+1,j,k)
            arr_out(2*ii+1,j,k) = sin_theta(ii)*arr(ii+1,j,k) - cos_theta(ii)*arr(nn-ii+1,j,k)
          end do
        end do
      end do
    end select
  end subroutine prep_dctiib
  subroutine posp_dctiib(nn,n,idir,arr,arr_out,is_swap_order,is_negate_even)
    !
    ! post-processing of a signal to perform a fast forward
    ! discrete cosine transform (DCT) with FFTs (see Makhoul 1980)
    !
    ! the input signal v(n) is post-processed into a signal x(n)
    ! as follows (now done using the subroutine remap with ib=1):
    !
    ! v(n) = x(2n       ),              0 <= n <= floor((N-1)/2)
    !      = x(2N -2n -1), floor((N+1)/2) <= n <= N-1
    ! with n = 0,...,N-1 and N being the total number of elements of the
    ! signal.
    !
    ! post-processing required for computing the corresponding
    ! discrete sine transform (DST) may also be performed
    ! if one of the last boolean input variables is .true.
    !
    implicit none
    integer , intent(in   ) :: nn
    integer , intent(in   ), dimension(3) :: n            ! dimensions of input/output array
    integer , intent(in   ) :: idir                       ! array direction where the transform is taken
    real(rp), intent(in   ), dimension(:,:,:) :: arr     ! input/output array
    real(rp), intent(out  ), dimension(:,:,:) :: arr_out ! output array
    logical, intent(in) :: is_swap_order  ! swap order of the elements of the input array? (for DST)
    logical, intent(in) :: is_negate_even ! negate every other element of the input array?
    integer :: i,j,k
    !
    select case(idir)
    case(1)
      call remap(1,nn,n(2),n(3),arr,arr_out)
      if(is_swap_order ) call swap_order( nn,n(2),n(3),arr_out)
      if(is_negate_even) call negate_even(nn,n(2),n(3),arr_out)
    end select
  end subroutine posp_dctiib
  subroutine signal_processing(pre_or_pos,f_or_b,cbc,c_or_f,nn,n,idir,arr,arr_out)
    implicit none
    !
    ! wrapper subroutine for signal processing to compute FFT-based transforms
    !
    integer,          intent(in) :: pre_or_pos ! prior (0) or after (1) fft
    character(len=1), intent(in) :: f_or_b     ! forward or backward transform
    character(len=2), intent(in) :: cbc        ! type of boundary condition
    character(len=1), intent(in) :: c_or_f     ! cell- or face-centred BC?
    integer, intent(in)                       :: nn ! number of points in the signal
    integer, intent(in), dimension(3)         :: n
    integer, intent(in)                       :: idir
    real(rp), intent(inout), dimension(:,:,:) :: arr
    real(rp), intent(out  ), dimension(:,:,:), optional :: arr_out
    select case(cbc)
    case('PP')
      select case(f_or_b)
      case('F')
        if(pre_or_pos == 1) call posp_fftf(nn,n,idir,arr)
      case('B')
        if(pre_or_pos == 0) call prep_fftb(nn,n,idir,arr)
      end select
    case('NN')
      if(c_or_f == 'c') then
        select case(f_or_b)
        case('F')
          if(pre_or_pos == 0) call prep_dctiif(nn,n,idir,arr,arr_out,.false.,.false.)
          if(pre_or_pos == 1) call posp_dctiif(nn,n,idir,arr,arr_out,.false.,.false.)
        case('B')
          if(pre_or_pos == 0) call prep_dctiib(nn,n,idir,arr,arr_out,.false.,.false.)
          if(pre_or_pos == 1) call posp_dctiib(nn,n,idir,arr,arr_out,.false.,.false.)
        end select
      end if
    case('DD')
      if(c_or_f == 'c') then
        select case(f_or_b)
        case('F')
          if(pre_or_pos == 0) call prep_dctiif(nn,n,idir,arr,arr_out,.false.,.true. )
          if(pre_or_pos == 1) call posp_dctiif(nn,n,idir,arr,arr_out,.true. ,.false.)
        case('B')
          if(pre_or_pos == 0) call prep_dctiib(nn,n,idir,arr,arr_out,.true. ,.false.)
          if(pre_or_pos == 1) call posp_dctiib(nn,n,idir,arr,arr_out,.false.,.true. )
        end select
      end if
    case default
      error stop 'ERROR: unsupported boundary condition' ! should be trapped before under `sanity.f90`
    end select
  end subroutine signal_processing
  subroutine negate_even(n,n2,n3,arr)
    implicit none
    integer , intent(in   ) :: n,n2,n3
    real(rp), intent(inout) :: arr(:,:,:)
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
    real(rp), intent(inout) :: arr(:,:,:)
    real(rp) :: tmp
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
  subroutine remap(ib,n,n2,n3,arr,arr_out)
    !
    ! maps a signal x to v (ib = 0), or v to x (ib = 1)
    ! where:
    ! v(n) = x(2n       ),              0 <= n <= floor((N-1)/2)
    !      = x(2N -2n -1), floor((N+1)/2) <= n <= N-1
    ! with n = 0,...,N-1; N = size(v) = size(x)
    !
    implicit none
    integer , intent(in ) :: ib,n,n2,n3
    real(rp), intent(in ) :: arr(:,:,:)
    real(rp), intent(out) :: arr_out(:,:,:)
    integer :: i,j,k
    integer :: nh
    nh = (n+1)/2
    select case(ib)
    case(0)
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1,n3
        do j=1,n2
          do i=1,nh
            arr_out(i,j,k) = arr(2*i-1,j,k)
            if(i+nh <= n) arr_out(i+nh,j,k) = arr(2*(n-(i+nh)+1),j,k)
          end do
        end do
      end do
    case(1)
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1,n3
        do j=1,n2
          do i=1,nh
            arr_out(2*i-1,j,k) = arr(i,j,k)
            if(i+nh <= n) arr_out(2*(n-(i+nh)+1),j,k) = arr(i+nh,j,k)
          end do
        end do
      end do
    end select
  end subroutine remap
#endif
end module mod_fft
