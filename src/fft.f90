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
#if defined(_OPENACC) || defined(_OPENMP)
  use mod_utils     , only: f_sizeof
#endif
  private
  public fftini,fftend,fft,prep_dctviii,posp_dctviii
#if !(defined(_OPENACC) || defined(_OPENMP))
  real(rp), allocatable, target :: fft_work(:)
  integer :: nplans = 0
#endif
#if defined(_OPENACC) || defined(_OPENMP)
  public fft_gpu,fft_gpu_layout
  public signal_processing,fftf_gpu,fftb_gpu
  integer(i8), public :: wsize_fft = 0,wsize_tmp = 0
  real(rp), allocatable, target :: sincos_theta_x(:,:,:),sincos_theta_y(:,:,:)
  integer :: n_sincos_theta_x(0:3),n_sincos_theta_y(0:3)
#endif
  contains
  subroutine fftini(ng,n_x,n_y,bcxy,c_or_f,arrplan,normfft)
    implicit none
    integer , target, intent(in), dimension(3) :: ng,n_x,n_y
    character(len=1), intent(in), dimension(0:1,2) :: bcxy
    character(len=1), intent(in), dimension(2) :: c_or_f
#if !(defined(_OPENACC) || defined(_OPENMP)) || defined(_USE_HIP)
    type(C_PTR), intent(out), dimension(2,2) :: arrplan
#else
    integer    , intent(out), dimension(2,2) :: arrplan
#endif
    real(rp), intent(out) :: normfft
#if !(defined(_OPENACC) || defined(_OPENMP))
    integer(i8) :: nw_x(3),nw_y(3)
#endif
    !
#if !(defined(_OPENACC) || defined(_OPENMP))
    !
    ! one buffer serves both pencils on the fixed grid;
    ! planning and execution use the same allocation to preserve
    ! FFTW alignment across all plans
    !
    if((any((c_or_f(:) == 'f').and.(bcxy(0,:) /= bcxy(1,:)))).and.(.not.allocated(fft_work))) then
      nw_x = int(n_x,i8); nw_x(1) = 2*nw_x(1)-1
      nw_y = int(n_y,i8); nw_y(2) = 2*nw_y(2)-1
      allocate(fft_work(max(product(nw_x),product(nw_y))))
    end if
#endif
    normfft = 1.
    call fftini_axis(1,ng(1),n_x,bcxy(:,1),c_or_f(1),arrplan(:,1),normfft)
    call fftini_axis(2,ng(2),n_y,bcxy(:,2),c_or_f(2),arrplan(:,2),normfft)
#if !(defined(_OPENACC) || defined(_OPENMP))
    nplans = nplans+size(arrplan)
#endif
    normfft = normfft**(-1)
  end subroutine fftini
  !
  subroutine fftini_axis(idir,nn,n,bc,c_or_f,arrplan,normfft)
    !
    ! initialize one X/Y plan pair and accumulate its round-trip factor
    !
    implicit none
    integer , intent(in) :: idir,nn,n(3)
    character(len=1), intent(in) :: bc(0:1),c_or_f
#if !(defined(_OPENACC) || defined(_OPENMP)) || defined(_USE_HIP)
    type(C_PTR), intent(out), dimension(2) :: arrplan
#else
    integer    , intent(out), dimension(2) :: arrplan
#endif
    real(rp), intent(inout) :: normfft
    integer :: kind_fwd,kind_bwd,iexcl
    real(rp) :: norm(2)
#if !(defined(_OPENACC) || defined(_OPENMP))
    real(rp), target :: arr(n(1),n(2),n(3))
    real(rp), pointer, contiguous :: arrwork(:,:,:)
    type(fftw_iodim) :: iodim(1),iodim_howmany(2)
    integer(C_INT) :: nx,ny,nz
#else
    integer :: istat,batch,nw(3)
    integer(C_INT), target :: nfft
    integer(c_intptr_t), target :: wsize
    integer(c_intptr_t) :: max_wsize
#endif
    !
    call find_fft(bc,c_or_f,kind_fwd,kind_bwd,norm)
    iexcl = 0
    !
    ! exclude the dependent endpoint for face DD/NN
    !
    if((c_or_f == 'f').and.(any(bc(0)//bc(1) == ['DD','NN']))) iexcl = 1
#if !(defined(_OPENACC) || defined(_OPENMP))
    !
    ! prepare plans with guru interface
    !
    nx = n(1); ny = n(2); nz = n(3)
    arrwork(1:nx,1:ny,1:nz) => arr(:,:,:)
    if((c_or_f == 'f').and.(bc(0) /= bc(1))) then
      if(idir == 1) nx = 2*nn-1
      if(idir == 2) ny = 2*nn-1
      arrwork(1:nx,1:ny,1:nz) => fft_work(1:1_i8*nx*ny*nz)
    end if
    !
    ! X is unit-stride; Y strides over nx; keep full extents for batches
    !
    iodim(1)%n  = merge(nx,ny,idir == 1)-iexcl
    iodim(1)%is = merge(1 ,nx,idir == 1)
    iodim(1)%os = iodim(1)%is
    iodim_howmany(1)%n  = merge(ny,nx,idir == 1)
    iodim_howmany(1)%is = merge(nx,1 ,idir == 1)
    iodim_howmany(1)%os = iodim_howmany(1)%is
    iodim_howmany(2)%n  = nz
    iodim_howmany(2)%is = nx*ny
    iodim_howmany(2)%os = iodim_howmany(2)%is
    arrplan(1)=fftw_plan_guru_r2r(1,iodim,2,iodim_howmany,arrwork,arrwork,kind_fwd,FFTW_ESTIMATE)
    arrplan(2)=fftw_plan_guru_r2r(1,iodim,2,iodim_howmany,arrwork,arrwork,kind_bwd,FFTW_ESTIMATE)
#else
    !
    ! idir selects the physical axis; both GPU pencils are axis-contiguous
    !
    if(bc(0)//bc(1) /= 'PP') then
      select case(idir)
      case(1)
        call init_sincos_theta(nn,sincos_theta_x,n_sincos_theta_x)
      case(2)
        call init_sincos_theta(nn,sincos_theta_y,n_sincos_theta_y)
      end select
    end if
    batch = product(n(2:3))
    call fft_gpu_layout(nn,n,bc(0)//bc(1),c_or_f,nfft,nw)
    wsize_tmp = max(wsize_tmp,product(int(nw(:),i8)))
    max_wsize = -1
    istat = cufftCreate(arrplan(1))
    istat = cufftSetAutoAllocation(arrplan(1),0)
#if !defined(_USE_HIP)
    istat = cufftMakePlanMany(arrplan(1),1,nfft       ,null()    ,1,nw(1),null()    ,1,nfft/2+1,CUFFT_FWD_TYPE,batch,wsize)
#else
    istat = cufftMakePlanMany(arrplan(1),1,c_loc(nfft),c_null_ptr,1,nw(1),c_null_ptr,1,nfft/2+1,CUFFT_FWD_TYPE,batch,c_loc(wsize))
#endif
    max_wsize = max(wsize,max_wsize)
    !
    istat = cufftCreate(arrplan(2))
    istat = cufftSetAutoAllocation(arrplan(2),0)
#if !defined(_USE_HIP)
    istat = cufftMakePlanMany(arrplan(2),1,nfft       ,null()    ,1,nfft/2+1,null()    ,1,nw(1),CUFFT_BWD_TYPE,batch,wsize)
#else
    istat = cufftMakePlanMany(arrplan(2),1,c_loc(nfft),c_null_ptr,1,nfft/2+1,c_null_ptr,1,nw(1),CUFFT_BWD_TYPE,batch,c_loc(wsize))
#endif
    max_wsize = max(wsize,max_wsize)
    !
    ! all axes and fields share the workspace allocated after initialization
    !
    wsize_fft = max(wsize_fft,(max_wsize+f_sizeof(1._rp)-1)/f_sizeof(1._rp))
#endif
    normfft = normfft*norm(1)*(nn+norm(2)-iexcl)
  end subroutine fftini_axis
  !
#if defined(_OPENACC) || defined(_OPENMP)
  subroutine init_sincos_theta(nn,sincos_theta,n_sincos_theta)
    use mod_param, only: pi
    implicit none
    integer , intent(in) :: nn
    real(rp), allocatable, target, intent(inout) :: sincos_theta(:,:,:)
    integer , intent(inout) :: n_sincos_theta(0:3)
    integer :: ii,ip,n_theta
    real(rp) :: theta
    !
    if(allocated(sincos_theta)) return
    !
    ! keep phases for cell (N), face NN (N-1), and mixed (2N,2N-1) transforms
    !
    allocate(sincos_theta(0:nn,1:2,0:3))
    n_sincos_theta(:) = [nn,max(1,nn-1),2*nn,2*nn-1]
    sincos_theta(:,:,:) = 0.
    do ip=0,3
      n_theta = n_sincos_theta(ip)
      do ii=0,n_theta/2
        theta = pi*ii/(2._rp*n_theta)
        sincos_theta(ii,1,ip) = sin(theta)
        sincos_theta(ii,2,ip) = cos(theta)
      end do
    end do
    !$acc        enter data copyin(sincos_theta)
    !$omp target enter data map(to:sincos_theta)
  end subroutine init_sincos_theta
  !
#endif
  subroutine fftend(arrplan)
    implicit none
#if !(defined(_OPENACC) || defined(_OPENMP)) || defined(_USE_HIP)
    type(C_PTR), intent(in), dimension(:,:) :: arrplan
#else
    integer    , intent(in), dimension(:,:) :: arrplan
#endif
    integer :: i,j
#if defined(_OPENACC) || defined(_OPENMP)
    integer :: istat
#endif
#if !(defined(_OPENACC) || defined(_OPENMP))
#if defined(_SINGLE_PRECISION)
    do j=1,size(arrplan,2)
      do i=1,size(arrplan,1)
        call sfftw_destroy_plan(arrplan(i,j))
      end do
    end do
#else
    do j=1,size(arrplan,2)
      do i=1,size(arrplan,1)
        call dfftw_destroy_plan(arrplan(i,j))
      end do
    end do
#endif
    !
    ! temporary self-test plans share FFTW state with the main solver plans
    !
    nplans = nplans-size(arrplan)
    if(nplans == 0) then
      if(allocated(fft_work)) deallocate(fft_work)
    end if
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
#if !(defined(_OPENACC) || defined(_OPENMP))
#if defined(_SINGLE_PRECISION)
    call sfftw_execute_r2r(plan,arr,arr)
#else
    call dfftw_execute_r2r(plan,arr,arr)
#endif
#endif
  end subroutine fft
  !
  subroutine prep_dctviii(f_or_b,cbc,idir,arr,arr_out)
    !
    ! M=N-1 independent values; DCT-II of [x,0,-reverse(x)] has only
    ! odd modes; half those coefficients give the ND DCT-VIII; embed them
    ! at odd modes for DCT-III; the round-trip factor is 2M+1
    ! DN reverses physical values before/after the same construction
    !
    implicit none
    character(len=1), intent(in) :: f_or_b
    character(len=2), intent(in) :: cbc
    integer, intent(in) :: idir
    real(rp), intent(in) :: arr(:,:,:)
    real(rp), pointer, contiguous, intent(out) :: arr_out(:,:,:)
#if !(defined(_OPENACC) || defined(_OPENMP))
    integer :: n(3),nw(3),m,nfft,i,j,k,ii,jj
    logical :: is_reverse
    !
    n = shape(arr); m = n(idir)-1; nfft = 2*m+1
    nw = n; nw(idir) = nfft
    arr_out(1:nw(1),1:nw(2),1:nw(3)) => fft_work(1:product(int(nw,i8)))
    is_reverse = cbc == 'DN'
    arr_out(:,:,:) = 0.
    select case(idir)
    case(1)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,m
            if(f_or_b == 'F') then
              ii = i
              if(is_reverse) ii = m-i+1
              arr_out(i,j,k) = arr(ii,j,k)
              arr_out(nfft-i+1,j,k) = -arr(ii,j,k)
            else
              arr_out(2*i,j,k) = arr(i,j,k)
            end if
          end do
        end do
      end do
    case(2)
      do k=1,n(3)
        do j=1,m
          do i=1,n(1)
            if(f_or_b == 'F') then
              jj = j
              if(is_reverse) jj = m-j+1
              arr_out(i,j,k) = arr(i,jj,k)
              arr_out(i,nfft-j+1,k) = -arr(i,jj,k)
            else
              arr_out(i,2*j,k) = arr(i,j,k)
            end if
          end do
        end do
      end do
    end select
#endif
  end subroutine prep_dctviii
  !
  subroutine posp_dctviii(f_or_b,cbc,idir,arr,arr_out)
    !
    ! extract odd coefficients or restore physical values and the endpoint
    !
    implicit none
    character(len=1), intent(in) :: f_or_b
    character(len=2), intent(in) :: cbc
    integer, intent(in) :: idir
    real(rp), intent(in) :: arr(:,:,:)
    real(rp), intent(out) :: arr_out(:,:,:)
#if !(defined(_OPENACC) || defined(_OPENMP))
    integer :: n(3),m,i,j,k,ii,jj
    logical :: is_reverse
    !
    n = shape(arr_out); m = n(idir)-1
    is_reverse = cbc == 'DN'
    select case(idir)
    case(1)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,m
            if(f_or_b == 'F') then
              arr_out(i,j,k) = 0.5_rp*arr(2*i,j,k)
            else
              ii = i
              if(is_reverse) ii = m-i+1
              arr_out(ii,j,k) = arr(i,j,k)
            end if
          end do
        end do
      end do
      do k=1,n(3)
        do j=1,n(2)
          arr_out(m+1,j,k) = 0.
          if((f_or_b == 'B').and.(is_reverse)) arr_out(m+1,j,k) = arr_out(m,j,k)
        end do
      end do
    case(2)
      do k=1,n(3)
        do j=1,m
          do i=1,n(1)
            if(f_or_b == 'F') then
              arr_out(i,j,k) = 0.5_rp*arr(i,2*j,k)
            else
              jj = j
              if(is_reverse) jj = m-j+1
              arr_out(i,jj,k) = arr(i,j,k)
            end if
          end do
        end do
      end do
      do k=1,n(3)
        do i=1,n(1)
          arr_out(i,m+1,k) = 0.
          if((f_or_b == 'B').and.(is_reverse)) arr_out(i,m+1,k) = arr_out(i,m,k)
        end do
      end do
    end select
#endif
  end subroutine posp_dctviii
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
        kind_fwd = FFTW_REDFT10
        kind_bwd = FFTW_REDFT01
        norm = [2.,0.]
      case('DD')
        kind_fwd = FFTW_RODFT00
        kind_bwd = FFTW_RODFT00
        norm = [2.,1.]
      case('ND','DN')
        kind_fwd = FFTW_REDFT10
        kind_bwd = FFTW_REDFT01
        norm = [2.,-0.5]
      end select
    end if
  end subroutine find_fft
  !
#if defined(_OPENACC) || defined(_OPENMP)
  subroutine fft_gpu_layout(nn,n,cbc,c_or_f,nfft,nwork)
    !
    ! FFT length and local scratch shape; the distributed pencil shape stays unchanged
    !
    implicit none
    integer , intent(in) :: nn
    integer , intent(in), dimension(3) :: n
    character(len=2), intent(in) :: cbc
    character(len=1), intent(in) :: c_or_f
    integer , intent(out) :: nfft
    integer , intent(out), dimension(3) :: nwork
    !
    nfft = nn
    nwork(:) = n(:)
    if(    ((cbc == 'DD'                 ).and.(c_or_f == 'f')) .or. &
           ((cbc == 'ND' .or. cbc == 'DN').and.(c_or_f == 'c'))) then
      nfft = 2*nn
    else if((cbc == 'ND' .or. cbc == 'DN').and.(c_or_f == 'f')) then
      nfft = 2*nn-1
    else if((cbc == 'NN'                 ).and.(c_or_f == 'f')) then
      nfft = nn-1
    end if
    nwork(1) = 2*(nfft/2+1)
  end subroutine fft_gpu_layout
  !
  subroutine get_sincos_theta(nn,sin_theta,cos_theta)
    implicit none
    integer, intent(in) :: nn
    real(rp), pointer, contiguous, intent(out) :: sin_theta(:),cos_theta(:)
    integer :: ip
    !
    nullify(sin_theta,cos_theta)
    do ip=0,3
      if(allocated(sincos_theta_x)) then
        if(n_sincos_theta_x(ip) == nn) then
          sin_theta(0:nn/2) => sincos_theta_x(0:nn/2,1,ip)
          cos_theta(0:nn/2) => sincos_theta_x(0:nn/2,2,ip)
          return
        end if
      end if
      if(allocated(sincos_theta_y)) then
        if(n_sincos_theta_y(ip) == nn) then
          sin_theta(0:nn/2) => sincos_theta_y(0:nn/2,1,ip)
          cos_theta(0:nn/2) => sincos_theta_y(0:nn/2,2,ip)
          return
        end if
      end if
    end do
    error stop 'ERROR: sincos arrays were not computed for this transform length.'
  end subroutine get_sincos_theta
  !
  subroutine fftf_gpu(plan,arr)
    implicit none
#if !defined(_USE_HIP)
    integer    , intent(in) :: plan
#else
    type(C_PTR), intent(in) :: plan
#endif
    real(rp), target, intent(inout), dimension(:,:,:) :: arr
    integer :: istat
    !$acc   host_data use_device(     arr)
    !$omp target data use_device_addr(arr)
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
    !$omp end target data
    !$acc end   host_data
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
    !$acc   host_data use_device(     arr)
    !$omp target data use_device_addr(arr)
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
    !$omp end target data
    !$acc end   host_data
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
    !
    if(cbc /= 'PP') then
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
  !
  subroutine prep_dsti(f_or_b,nn,n,arr,arr_out)
    !
    ! face DD has N-1 independent values; for a forward DST-I, form the
    ! length-2N odd extension [0,x(1),...,x(N-1),0,-x(N-1),...,-x(1)]
    ! for the inverse, form a purely imaginary Hermitian half-spectrum
    !
    implicit none
    character(len=1), intent(in) :: f_or_b
    integer , intent(in) :: nn
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in ), dimension(:,:,:) :: arr
    real(rp), intent(out), dimension(:,:,:) :: arr_out
    integer :: i,j,k,n_2,n_3
    !
    n_2 = n(2); n_3 = n(3)
    select case(f_or_b)
    case('F')
      !$acc parallel     loop collapse(3) default(present) async(1)
      !$omp target teams loop collapse(3)
      do k=1,n_3
        do j=1,n_2
          do i=1,2*nn+2
            arr_out(i,j,k) = 0.
            if((i >= 2   ).and.(i <= nn  )) arr_out(i,j,k) =  arr(     i-1,j,k)
            if((i >= nn+2).and.(i <= 2*nn)) arr_out(i,j,k) = -arr(2*nn-i+1,j,k)
          end do
        end do
      end do
    case('B')
      !$acc parallel     loop collapse(3) default(present) async(1)
      !$omp target teams loop collapse(3)
      do k=1,n_3
        do j=1,n_2
          do i=1,nn+1
            arr_out(2*i-1,j,k) = 0.
            arr_out(2*i  ,j,k) = 0.
            if((i >= 2).and.(i <= nn)) arr_out(2*i,j,k) = -arr(i-1,j,k)
          end do
        end do
      end do
    end select
  end subroutine prep_dsti
  !
  subroutine posp_dsti(f_or_b,nn,n,arr,arr_out)
    !
    ! extract DST-I coefficients or interior values, with round-trip factor 2N;
    ! the excluded Dirichlet boundary slot and FFT padding are set to zero
    !
    implicit none
    character(len=1), intent(in) :: f_or_b
    integer , intent(in) :: nn
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in ), dimension(:,:,:) :: arr
    real(rp), intent(out), dimension(:,:,:) :: arr_out
    integer :: i,j,k,n_1,n_2,n_3
    !
    n_1 = n(1); n_2 = n(2); n_3 = n(3)
    select case(f_or_b)
    case('F')
      !$acc parallel     loop collapse(3) default(present) async(1)
      !$omp target teams loop collapse(3)
      do k=1,n_3
        do j=1,n_2
          do i=1,n_1
            if(i >= nn) then
              arr_out(i,j,k) = 0.
            else
              arr_out(i,j,k) = -arr(2*i+2,j,k)
            end if
          end do
        end do
      end do
    case('B')
      !$acc parallel     loop collapse(3) default(present) async(1)
      !$omp target teams loop collapse(3)
      do k=1,n_3
        do j=1,n_2
          do i=1,n_1
            if(i >= nn) then
              arr_out(i,j,k) = 0.
            else
              arr_out(i,j,k) =  arr(i+1  ,j,k)
            end if
          end do
        end do
      end do
    end select
  end subroutine posp_dsti
  !
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
      !$acc parallel     loop collapse(2) default(present) async(1)
      !$omp target teams loop collapse(2)
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
      !$acc parallel     loop collapse(2) default(present) async(1)
      !$omp target teams loop collapse(2)
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
      call get_sincos_theta(nn,sin_theta,cos_theta)
      n_2 = n(2); n_3 = n(3)
      !$acc parallel     loop collapse(3) default(present) private(i) async(1)
      !$omp target teams loop collapse(3)                  private(i)
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
      call get_sincos_theta(nn,sin_theta,cos_theta)
      if(is_swap_order ) call swap_order( nn,n(2),n(3),arr)
      if(is_negate_even) call negate_even(nn,n(2),n(3),arr)
      n_2 = n(2); n_3 = n(3)
      !$acc parallel     loop collapse(2) default(present) async(1)
      !$omp target teams loop collapse(2)
      do k=1,n_3
        do j=1,n_2
          do ii = 1,2-mod(nn,2)
            arr(nn+ii,j,k) = 0.
          end do
        end do
      end do
      !$acc parallel     loop collapse(3) default(present) private(ii) async(1)
      !$omp target teams loop collapse(3)                  private(ii)
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
    integer , intent(in   ), dimension(3) :: n           ! dimensions of input/output array
    integer , intent(in   ) :: idir                      ! array direction where the transform is taken
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
  !
  subroutine prep_dctiv_viii(f_or_b,cbc,c_or_f,nn,n,arr,arr_out)
    !
    ! take half the odd DCT-II coefficients of [x,-reverse(x)] for cell
    ! DCT-IV transforms, or [x,0,-reverse(x)] for DCT-VIII on M=N-1 face interiors
    ! fuse the extension with the FFT reorder; the inverse embeds input
    ! in odd DCT-II modes; cell DN negates even points and reverses modes;
    ! face DN reverses physical points before/after the ND transform
    !
    implicit none
    character(len=1), intent(in) :: f_or_b,c_or_f
    character(len=2), intent(in) :: cbc
    integer , intent(in) :: nn,n(3)
    real(rp), intent(in ) :: arr(:,:,:)
    real(rp), intent(out) :: arr_out(:,:,:)
    logical :: is_sine,is_reverse,is_face
    integer :: i,j,k,ii,ip,im,n_2,n_3,m,nfft
    real(rp) :: s,a,b
    real(rp), pointer, contiguous :: sin_theta(:),cos_theta(:)
    !
    is_face = c_or_f == 'f'
    is_sine    = (cbc == 'DN').and.(.not.is_face)
    is_reverse = (cbc == 'DN').and.(     is_face)
    m = nn; nfft = 2*nn
    if(is_face) then
      m = nn-1; nfft = 2*nn-1
    end if
    n_2 = n(2); n_3 = n(3)
    select case(f_or_b)
    case('F')
      !$acc parallel     loop collapse(3) default(present) private(ip,ii,s) async(1)
      !$omp target teams loop collapse(3)                  private(ip,ii,s)
      do k=1,n_3
        do j=1,n_2
          do i=1,2*(nfft/2+1)
            arr_out(i,j,k) = 0.
            if(i <= nfft) then
              ip = 2*i-1
              if(i > (nfft+1)/2) ip = 2*(nfft-i+1)
              ii = ip; s = 1.
              if(ip > m) then
                ii = nfft-ip+1; s = -1.
              end if
              if(ii <= m) then
                if(is_sine.and.(mod(ii,2) == 0)) s = -s
                if(is_reverse) ii = m-ii+1
                arr_out(i,j,k) = s*arr(ii,j,k)
              end if
            end if
          end do
        end do
      end do
    case('B')
      call get_sincos_theta(nfft,sin_theta,cos_theta)
      !$acc parallel     loop collapse(3) default(present) private(ip,im,a,b) async(1)
      !$omp target teams loop collapse(3)                  private(ip,im,a,b)
      do k=1,n_3
        do j=1,n_2
          do i=0,nfft/2
            a = 0.; b = 0.
            if(mod(i,2) == 1) then
              ip = (i+1)/2
              if(is_sine) ip = m-ip+1
              a = arr(ip,j,k)
            end if
            if((i > 0).and.(mod(nfft-i,2) == 1)) then
              im = (nfft-i+1)/2
              if(is_sine) im = m-im+1
              b = arr(im,j,k)
            end if
            arr_out(2*i+1,j,k) = cos_theta(i)*a + sin_theta(i)*b
            arr_out(2*i+2,j,k) = sin_theta(i)*a - cos_theta(i)*b
            if((mod(nfft,2) == 0).and.(i == nfft/2)) arr_out(2*i+2,j,k) = 0.
          end do
        end do
      end do
    end select
  end subroutine prep_dctiv_viii
  !
  subroutine posp_dctiv_viii(f_or_b,cbc,c_or_f,nn,n,arr,arr_out)
    implicit none
    character(len=1), intent(in) :: f_or_b,c_or_f
    character(len=2), intent(in) :: cbc
    integer , intent(in) :: nn,n(3)
    real(rp), intent(in ) :: arr(:,:,:)
    real(rp), intent(out) :: arr_out(:,:,:)
    logical :: is_sine,is_reverse,is_face
    integer :: i,j,k,ii,ip,n_1,n_2,n_3,m,nfft
    real(rp), pointer, contiguous :: sin_theta(:),cos_theta(:)
    !
    is_face = c_or_f == 'f'
    is_sine    = (cbc == 'DN').and.(.not.is_face)
    is_reverse = (cbc == 'DN').and.(     is_face)
    m = nn; nfft = 2*nn
    if(is_face) then
      m = nn-1; nfft = 2*nn-1
    end if
    n_1 = n(1); n_2 = n(2); n_3 = n(3)
    !
    ! clear padding before it enters a transform along the other axis
    !
    select case(f_or_b)
    case('F')
      !
      ! extract half the odd DCT-II coefficients from the half-spectrum
      ! the round-trip factor is the extension length (2N or 2N-1)
      !
      call get_sincos_theta(nfft,sin_theta,cos_theta)
      !$acc parallel     loop collapse(3) default(present) private(ii,ip) async(1)
      !$omp target teams loop collapse(3)                  private(ii,ip)
      do k=1,n_3
        do j=1,n_2
          do i=1,n_1
            if(i > m) then
              arr_out(i,j,k) = 0.
            else
              ii = 2*i-1
              if(is_sine) ii = 2*(m-i+1)-1
              if(ii <= nfft/2) then
                ip = 2*ii+1
                arr_out(i,j,k) = cos_theta(ii)*arr(ip,j,k) + sin_theta(ii)*arr(ip+1,j,k)
              else
                ii = nfft-ii; ip = 2*ii+1
                arr_out(i,j,k) = sin_theta(ii)*arr(ip,j,k) - cos_theta(ii)*arr(ip+1,j,k)
              end if
            end if
          end do
        end do
      end do
    case('B')
      !$acc parallel     loop collapse(3) default(present) private(ip,ii) async(1)
      !$omp target teams loop collapse(3)                  private(ip,ii)
      do k=1,n_3
        do j=1,n_2
          do i=1,n_1
            if(i > m) then
              arr_out(i,j,k) = 0.
              if(is_reverse.and.(i == nn)) arr_out(i,j,k) = arr(1,j,k)
            else
              ip = (i+1)/2
              if(mod(i,2) == 0) ip = nfft-i/2+1
              ii = i
              if(is_reverse) ii = m-i+1
              arr_out(ii,j,k) = arr(ip,j,k)
              if(is_sine.and.(mod(i,2) == 0)) arr_out(ii,j,k) = -arr_out(ii,j,k)
            end if
          end do
        end do
      end do
    end select
  end subroutine posp_dctiv_viii
  !
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
    integer :: m,j,k,n_2,n_3
    !
    select case(cbc)
    case('PP')
      select case(f_or_b)
      case('F')
        if(pre_or_pos == 1) call posp_fftf(nn,n,idir,arr)
      case('B')
        if(pre_or_pos == 0) call prep_fftb(nn,n,idir,arr)
      end select
    case('NN')
      m = nn
      if(c_or_f == 'f') m = nn-1
      select case(f_or_b)
      case('F')
        if(pre_or_pos == 0) call prep_dctiif(m,n,idir,arr,arr_out,.false.,.false.)
        if(pre_or_pos == 1) call posp_dctiif(m,n,idir,arr,arr_out,.false.,.false.)
      case('B')
        if(pre_or_pos == 0) call prep_dctiib(m,n,idir,arr,arr_out,.false.,.false.)
        if(pre_or_pos == 1) call posp_dctiib(m,n,idir,arr,arr_out,.false.,.false.)
      end select
      if((c_or_f == 'f').and.(pre_or_pos == 1).and.(f_or_b == 'B')) then
        !
        ! restore the dependent boundary value after the inverse transform
        !
        n_2 = n(2); n_3 = n(3)
        !$acc parallel     loop collapse(2) default(present) async(1)
        !$omp target teams loop collapse(2)
        do k=1,n_3
          do j=1,n_2
            arr_out(nn,j,k) = arr_out(m,j,k)
          end do
        end do
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
      else if(c_or_f == 'f') then
        if(pre_or_pos == 0) call prep_dsti(f_or_b,nn,n,arr,arr_out)
        if(pre_or_pos == 1) call posp_dsti(f_or_b,nn,n,arr,arr_out)
      end if
    case('ND','DN')
      if(pre_or_pos == 0) call prep_dctiv_viii(f_or_b,cbc,c_or_f,nn,n,arr,arr_out)
      if(pre_or_pos == 1) call posp_dctiv_viii(f_or_b,cbc,c_or_f,nn,n,arr,arr_out)
    case default
      error stop 'ERROR: unsupported boundary condition' ! should be trapped before under `sanity.f90`
    end select
  end subroutine signal_processing
  subroutine negate_even(n,n2,n3,arr)
    implicit none
    integer , intent(in   ) :: n,n2,n3
    real(rp), intent(inout) :: arr(:,:,:)
    integer :: i,j,k
    !$acc parallel     loop collapse(3) default(present) async(1)
    !$omp target teams loop collapse(3)
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
    !$acc parallel     loop collapse(3) default(present) private(tmp) async(1)
    !$omp target teams loop collapse(3)                  private(tmp)
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
      !$acc parallel     loop collapse(3) default(present) async(1)
      !$omp target teams loop collapse(3)
      do k=1,n3
        do j=1,n2
          do i=1,nh
            arr_out(i,j,k) = arr(2*i-1,j,k)
            if(i+nh <= n) arr_out(i+nh,j,k) = arr(2*(n-(i+nh)+1),j,k)
            !
            ! initialize real padding read by cuFFT
            !
            if(i == 1) then
              arr_out(n+1,j,k) = 0.
              if(mod(n,2) == 0) arr_out(n+2,j,k) = 0.
            end if
          end do
        end do
      end do
    case(1)
      !$acc parallel     loop collapse(3) default(present) async(1)
      !$omp target teams loop collapse(3)
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
  subroutine copy(a,b)
    real(rp), intent(in ), dimension(:,:,:) :: a
    real(rp), intent(out), dimension(:,:,:) :: b
    integer :: i,j,k,n1,n2,n3
    n1 = min(size(a,1),size(b,1))
    n2 = min(size(a,2),size(b,2))
    n3 = min(size(a,3),size(b,3))
    !$acc parallel     loop collapse(3) default(present) async(1)
    !$omp target teams loop collapse(3)
    do k=1,n3
      do j=1,n2
        do i=1,n1
          b(i,j,k) = a(i,j,k)
        end do
      end do
    end do
  end subroutine copy
#endif
end module mod_fft
