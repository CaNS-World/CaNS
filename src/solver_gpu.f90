! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2025 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_solver_gpu
#if defined(_OPENACC)
  use, intrinsic :: iso_c_binding, only: C_PTR
#if !defined(_USE_DIEZDECOMP)
  use cudecomp
#else
  use diezdecomp
#endif
  use mod_common_cudecomp, only: dtype_rp => cudecomp_real_rp, &
                                 cudecomp_is_t_in_place, &
                                 solver_buf_0,solver_buf_1, &
                                 pz_aux_1,work, &
                                 ap_x   => ap_x_poi, &
                                 ap_y   => ap_y_poi, &
                                 ap_z   => ap_z_poi, &
                                 ap_x_0 => ap_x    , &
                                 ap_z_0 => ap_z    , &
                                 ch => handle,gd => gd_poi, gd_io => gd_poi_io, &
                                 istream => istream_acc_queue_1_comm_lib
  use mod_fft            , only: signal_processing,fftf_gpu,fftb_gpu
  use mod_param          , only: ipencil_axis,is_poisson_pcr_tdma, &
                                 is_use_diezdecomp,is_diezdecomp_x2z_z2x_transposes
  use mod_types
  implicit none
  private
  public solver_gpu,solver_gaussel_z_gpu
  contains
  subroutine solver_gpu(n,ng,arrplan,normfft,lambdaxy,a,b,c,bc,c_or_f,p,is_ptdma_update,aa_z,cc_z)
    implicit none
    integer , intent(in), dimension(3) :: n,ng
#if !defined(_USE_HIP)
    integer    , intent(in), dimension(2,2) :: arrplan
#else
    type(C_PTR), intent(in), dimension(2,2) :: arrplan
#endif
    real(rp), intent(in) :: normfft
    real(rp), intent(in), dimension(:,:) :: lambdaxy
    real(rp), intent(in), dimension(:) :: a,b,c
    character(len=1), dimension(0:1,3), intent(in) :: bc
    character(len=1), intent(in), dimension(3) :: c_or_f
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    logical , intent(inout), target, optional :: is_ptdma_update
    real(rp), intent(inout), dimension(:,:,:), optional :: aa_z,cc_z
    real(rp), pointer, contiguous, dimension(:,:,:) :: px,py,pz
    integer :: i,j,k,q
    logical :: is_periodic_z
    integer, dimension(3) :: n_x,n_y,n_z,n_z_0,lo_z_0,hi_z_0,pad_io
    type(cudecompPencilInfo) :: ap_io
    integer :: istat
    logical :: is_ptdma_update_
    real(rp) :: norm
    !
    norm = normfft
    !
    is_ptdma_update_ = .true.
    if(present(is_ptdma_update)) is_ptdma_update_ = is_ptdma_update
    !
    n_z_0(:) = ap_z_0%shape(:)
    lo_z_0(:) = ap_z_0%lo(:)
    hi_z_0(:) = ap_z_0%hi(:)
    n_x(:) = ap_x%shape(:)
    n_y(:) = ap_y%shape(:)
    n_z(:) = ap_z%shape(:)
    if(is_poisson_pcr_tdma) then
      n_z(:) = n_z_0(:) ! equal to the (unpadded) ap_y%shape(:) under initmpi.f90
    end if
    px(1:n_x(1),1:n_x(2),1:n_x(3)) => solver_buf_0(1:product(n_x(:)))
    if(cudecomp_is_t_in_place) then
      py(1:n_y(1),1:n_y(2),1:n_y(3)) => solver_buf_0(1:product(n_y(:)))
    else
      py(1:n_y(1),1:n_y(2),1:n_y(3)) => solver_buf_1(1:product(n_y(:)))
    end if
    pz(1:n_z(1),1:n_z(2),1:n_z(3)) => solver_buf_0(1:product(n_z(:)))
    !
    select case(ipencil_axis)
    case(1)
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP parallel do   collapse(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            px(i,j,k) = p(i,j,k)
          end do
        end do
      end do
    case(2)
#if !defined(_USE_DIEZDECOMP)
      ap_io = ap_x_0
      pad_io(:) = ap_x%shape(:) - ap_io%shape(:)
#if !defined(_USE_DIEZDECOMP)
      !$acc host_data use_device(p,px,work)
#endif
      istat = cudecompTransposeYtoX(ch,gd_io,p ,px,work,dtype_rp,input_halo_extents  = [1,1,1], &
                                                                 output_halo_extents = [0,0,0], &
                                                                 input_padding       = [0,0,0], &
                                                                 output_padding      = pad_io, &
                                                                 stream=istream)
#if !defined(_USE_DIEZDECOMP)
      !$acc end host_data
#endif
    case(3)
      if(.not.is_diezdecomp_x2z_z2x_transposes) then
        istat = cudecompGetPencilInfo(ch,gd_io,ap_io,2) ! mem_order = [2,3,1]
        pad_io(:) = ap_y%shape(:) - ap_io%shape(:)
#if !defined(_USE_DIEZDECOMP)
        pad_io(ap_y%order(:)) = pad_io(:)
        !$acc host_data use_device(p,py,px,work)
#endif
        istat = cudecompTransposeZtoY(ch,gd_io,p ,py,work,dtype_rp,input_halo_extents  = [1,1,1], &
                                                                   output_halo_extents = [0,0,0], &
                                                                   input_padding       = [0,0,0], &
                                                                   output_padding      = pad_io, &
                                                                   stream=istream)
        istat = cudecompTransposeYtoX(ch,gd   ,py,px,work,dtype_rp,stream=istream)
#if !defined(_USE_DIEZDECOMP)
      !$acc end host_data
#endif
      else
#if defined(_USE_DIEZDECOMP)
        !$acc parallel loop collapse(3) default(present) async(1)
        !$OMP parallel do   collapse(3) DEFAULT(shared)
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              pz(i,j,k) = p(i,j,k)
            end do
          end do
        end do
        istat = diezdecompTransposeZtoX(ch,gd,pz,px,work,dtype_rp,stream=istream)
#endif
      end if
    end select
    !
    call signal_processing(0,'F',bc(0,1)//bc(1,1),c_or_f(1),ng(1),n_x,1,px)
    call fftf_gpu(arrplan(1,1),px)
    call signal_processing(1,'F',bc(0,1)//bc(1,1),c_or_f(1),ng(1),n_x,1,px)
    !
#if !defined(_USE_DIEZDECOMP)
    !$acc host_data use_device(px,py,work)
#endif
    istat = cudecompTransposeXtoY(ch,gd,px,py,work,dtype_rp,stream=istream)
#if !defined(_USE_DIEZDECOMP)
    !$acc end host_data
#endif
    !
    call signal_processing(0,'F',bc(0,2)//bc(1,2),c_or_f(2),ng(2),n_y,1,py)
    call fftf_gpu(arrplan(1,2),py)
    call signal_processing(1,'F',bc(0,2)//bc(1,2),c_or_f(2),ng(2),n_y,1,py)
    !
    q = merge(1,0,c_or_f(3) == 'f'.and.bc(1,3) == 'D'.and.hi_z_0(3) == ng(3))
    is_periodic_z = bc(0,3)//bc(1,3) == 'PP'
    if(.not.is_poisson_pcr_tdma) then
#if !defined(_USE_DIEZDECOMP)
      !$acc host_data use_device(py,pz,work)
#endif
      istat = cudecompTransposeYtoZ(ch,gd,py,pz,work,dtype_rp,stream=istream)
#if !defined(_USE_DIEZDECOMP)
      !$acc end host_data
#endif
      !
      call gaussel_gpu(n_z_0(1),n_z_0(2),n_z_0(3)-q,0,a,b,c,is_periodic_z,norm,pz,work,pz_aux_1,lambdaxy)
      !
#if !defined(_USE_DIEZDECOMP)
      !$acc host_data use_device(pz,py,work)
#endif
      istat = cudecompTransposeZtoY(ch,gd,pz,py,work,dtype_rp,stream=istream)
#if !defined(_USE_DIEZDECOMP)
      !$acc end host_data
#endif
    else
      block
        use mod_common_cudecomp, only: ap_y
        integer :: n_y_1,n_y_2,n_y_3
        !
        ! transpose py -> pz with non-axis-contiguous layout
        !
        n_y_3 = ap_y%shape(3)
        n_y_2 = ap_y%shape(2)
        n_y_1 = ap_y%shape(1)
        !$acc parallel loop collapse(3) default(present) async(1)
        !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
        do k=1,n_y_3
          do j=1,n_y_2
            do i=1,n_y_1
              pz(i,j,k) = py(j,k,i)
            end do
          end do
        end do
      end block
      !
      call gaussel_ptdma_gpu(n_z_0(1),n_z_0(2),n_z_0(3)-q,lo_z_0(3),0,a,b,c,is_periodic_z,norm,pz,work,pz_aux_1,is_ptdma_update_,lambdaxy,aa_z,cc_z)
      if(present(is_ptdma_update)) is_ptdma_update = is_ptdma_update_
      !
      block
        use mod_common_cudecomp, only: ap_y
        integer :: n_y_1,n_y_2,n_y_3
        !
        ! transpose pz -> py with axis-contiguous layout
        !
        n_y_3 = ap_y%shape(3)
        n_y_2 = ap_y%shape(2)
        n_y_1 = ap_y%shape(1)
        !$acc parallel loop collapse(3) default(present) async(1)
        !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
        do k=1,n_y_3
          do j=1,n_y_2
            do i=1,n_y_1
              py(j,k,i) = pz(i,j,k)
            end do
          end do
        end do
      end block
    end if
    !
    call signal_processing(0,'B',bc(0,2)//bc(1,2),c_or_f(2),ng(2),n_y,1,py)
    call fftb_gpu(arrplan(2,2),py)
    call signal_processing(1,'B',bc(0,2)//bc(1,2),c_or_f(2),ng(2),n_y,1,py)
    !
#if !defined(_USE_DIEZDECOMP)
    !$acc host_data use_device(py,px,work)
#endif
    istat = cudecompTransposeYtoX(ch,gd,py,px,work,dtype_rp,stream=istream)
#if !defined(_USE_DIEZDECOMP)
    !$acc end host_data
#endif
    !
    call signal_processing(0,'B',bc(0,1)//bc(1,1),c_or_f(1),ng(1),n_x,1,px)
    call fftb_gpu(arrplan(2,1),px)
    call signal_processing(1,'B',bc(0,1)//bc(1,1),c_or_f(1),ng(1),n_x,1,px)
    !
    select case(ipencil_axis)
    case(1)
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP parallel do   collapse(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            p(i,j,k) = px(i,j,k)
          end do
        end do
      end do
    case(2)
#if !defined(_USE_DIEZDECOMP)
      !$acc host_data use_device(px,p,work)
#endif
      istat = cudecompTransposeXtoY(ch,gd_io,px,p ,work,dtype_rp,input_halo_extents  = [0,0,0], &
                                                                 output_halo_extents = [1,1,1], &
                                                                 input_padding       = pad_io, &
                                                                 output_padding      = [0,0,0], &
                                                                 stream=istream)
#if !defined(_USE_DIEZDECOMP)
      !$acc end host_data
#endif
    case(3)
      if(.not.is_diezdecomp_x2z_z2x_transposes) then
#if !defined(_USE_DIEZDECOMP)
        !$acc host_data use_device(px,py,p,work)
#endif
        istat = cudecompTransposeXtoY(ch,gd   ,px,py,work,dtype_rp,stream=istream)
        istat = cudecompTransposeYtoZ(ch,gd_io,py,p ,work,dtype_rp,input_halo_extents  = [0,0,0], &
                                                                   output_halo_extents = [1,1,1], &
                                                                   input_padding       = pad_io, &
                                                                   output_padding      = [0,0,0], &
                                                                   stream=istream)
#if !defined(_USE_DIEZDECOMP)
      !$acc end host_data
#endif
      else
#if defined(_USE_DIEZDECOMP)
        istat = diezdecompTransposeXtoZ(ch,gd,px,pz,work,dtype_rp,stream=istream)
        !$acc parallel loop collapse(3) default(present) async(1)
        !$OMP parallel do   collapse(3) DEFAULT(shared)
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              p(i,j,k) = pz(i,j,k)
            end do
          end do
        end do
#endif
      end if
    end select
  end subroutine solver_gpu
  !
  subroutine gaussel_gpu(nx,ny,n,nh,a,b,c,is_periodic,norm,p,d,p2,lambdaxy)
    use mod_param, only: eps
    implicit none
    integer , intent(in) :: nx,ny,n,nh
    real(rp), intent(in), dimension(:) :: a,b,c
    logical , intent(in) :: is_periodic
    real(rp), intent(in) :: norm
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    real(rp),                dimension(nx,ny,n) :: d,p2
    real(rp), intent(in), dimension(:,:), optional :: lambdaxy
    real(rp) :: z,lxy
    integer :: i,j,k,nn
    real(rp), allocatable, save, dimension(:) :: dd,pp2
    !
    !solve tridiagonal system
    !
    nn = n
    if(is_periodic) nn = nn-1
    if(present(lambdaxy)) then
      !$acc parallel loop gang vector collapse(2) default(present) private(lxy,z) async(1)
      do j=1,ny
        do i=1,nx
          lxy = lambdaxy(i,j)
          !
          z = 1./(b(1)+lxy+eps)
          d(i,j,1) = c(1)*z
          p(i,j,1) = p(i,j,1)*norm*z
          !$acc loop seq
          do k=2,nn
            z        = 1./(b(k)+lxy-a(k)*d(i,j,k-1)+eps)
            p(i,j,k) = (p(i,j,k)*norm-a(k)*p(i,j,k-1))*z
            d(i,j,k) = c(k)*z
          end do
          !
          !$acc loop seq
          do k=nn-1,1,-1
            p(i,j,k) = p(i,j,k) - d(i,j,k)*p(i,j,k+1)
          end do
        end do
      end do
      if(is_periodic) then
        !$acc parallel loop gang vector collapse(2) default(present) private(lxy,z) async(1)
        do j=1,ny
          do i=1,nx
            lxy = lambdaxy(i,j)
            !
            !$acc loop seq
            do k=1,nn
              p2(i,j,k) = 0.
            end do
            p2(i,j,1 ) = -a(1 )
            p2(i,j,nn) = -c(nn)
            !
            z = 1./(b(1)+lxy+eps)
            d( i,j,1) = c(1)*z
            p2(i,j,1) = p2(i,j,1)*z
            !$acc loop seq
            do k=2,nn
              z         = 1./(b(k)+lxy-a(k)*d(i,j,k-1)+eps)
              p2(i,j,k) = (p2(i,j,k)-a(k)*p2(i,j,k-1))*z
              d(i,j,k)  = c(k)*z
            end do
            !
            !$acc loop seq
            do k=nn-1,1,-1
              p2(i,j,k) = p2(i,j,k) - d(i,j,k)*p2(i,j,k+1)
            end do
            !
            p(i,j,nn+1) = (p(i,j,nn+1)*norm       - c(nn+1)*p( i,j,1) - a(nn+1)*p( i,j,nn)) / &
                          (b(    nn+1)      + lxy + c(nn+1)*p2(i,j,1) + a(nn+1)*p2(i,j,nn)+eps)
            !$acc loop seq
            do k=1,nn
              p(i,j,k) = p(i,j,k) + p2(i,j,k)*p(i,j,nn+1)
            end do
          end do
        end do
      end if
    else
      if(.not.allocated(dd)) then
        allocate(dd(n+1)) ! needs to accomodate both face-centered and cell-centered variables
        !$acc enter data create(dd) async(1)
      end if
      !$acc parallel loop gang vector collapse(2) default(present) private(z) async(1)
      do j=1,ny
        do i=1,nx
          z = 1./(b(1)+eps)
          dd(1) = c(1)*z
          p(i,j,1) = p(i,j,1)*norm*z
          !$acc loop seq
          do k=2,nn
            z        = 1./(b(k)-a(k)*dd(k-1)+eps)
            p(i,j,k) = (p(i,j,k)*norm-a(k)*p(i,j,k-1))*z
            dd(k)    = c(k)*z
          end do
          !
          !$acc loop seq
          do k=nn-1,1,-1
            p(i,j,k) = p(i,j,k) - dd(k)*p(i,j,k+1)
          end do
        end do
      end do
      if(is_periodic) then
        if(.not.allocated(pp2)) then
          allocate(pp2(n+1)) ! needs to accomodate both face-centered and cell-centered variables
          !$acc enter data create(pp2) async(1)
        end if
        !$acc parallel loop gang vector collapse(2) default(present) private(z) async(1)
        do j=1,ny
          do i=1,nx
            !$acc loop seq
            do k=1,nn
              pp2(k) = 0.
            end do
            pp2(1 ) = -a(1 )
            pp2(nn) = -c(nn)
            !
            z = 1./(b(1)+eps)
            dd( 1) = c(1)*z
            pp2(1) = pp2(1)*z
            !$acc loop seq
            do k=2,nn
              z      = 1./(b(k)+lxy-a(k)*dd(k-1)+eps)
              pp2(k) = (pp2(k)-a(k)*pp2(k-1))*z
              dd(k)  = c(k)*z
            end do
            !
            !$acc loop seq
            do k=nn-1,1,-1
              pp2(k) = pp2(k) - dd(k)*pp2(k+1)
            end do
            !
            p(i,j,nn+1) = (p(i,j,nn+1)*norm       - c(nn+1)*p( i,j,1) - a(nn+1)*p( i,j,nn)) / &
                          (b(    nn+1)      + lxy + c(nn+1)*pp2(   1) + a(nn+1)*pp2(   nn)+eps)
            !$acc loop seq
            do k=1,nn-1
              p(i,j,k) = p(i,j,k) + pp2(k)*p(i,j,nn+1)
            end do
          end do
        end do
      end if
    end if
  end subroutine gaussel_gpu
  !
  subroutine gaussel_ptdma_gpu(nx,ny,n,lo,nh,a,b,c,is_periodic,norm,p,aa,cc,is_update,lambdaxy,aa_z_save,cc_z_save)
    !
    ! distributed TDMA solver
    !
    use mod_common_cudecomp, only: gd_ptdma,ap_z_ptdma,work => work_ptdma
    use mod_param          , only: eps
    !
    implicit none
    integer , intent(in) :: nx,ny,n,lo,nh
    real(rp), intent(in), dimension(:) :: a,b,c
    logical , intent(in) :: is_periodic
    real(rp), intent(in) :: norm
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    real(rp),                dimension(nx,ny,n) :: aa,cc
    logical , intent(inout), optional :: is_update
    real(rp), intent(in   ), dimension(:,:), optional :: lambdaxy
    real(rp), intent(inout), dimension(:,:,:), optional :: aa_z_save,cc_z_save
    real(rp), allocatable, dimension(:,:,:), save :: aa_y,cc_y,pp_y,aa_z,cc_z,pp_z
    real(rp), allocatable, dimension(:,:,:), save :: pp_z_2,cc_z_0
    real(rp) :: z,z1,z2,lxy
    integer :: i,j,k
    integer , dimension(3) :: nr_z
    integer :: nx_r,ny_r,nn,dk_g
    integer :: istat
    !
    nr_z(:) = ap_z_ptdma%shape(:)
    if(.not.allocated(pp_y)) then
      allocate(aa_y(nx,ny,2), &
               cc_y(nx,ny,2), &
               pp_y(nx,ny,2), &
               aa_z(nr_z(1),nr_z(2),nr_z(3)), &
               cc_z(nr_z(1),nr_z(2),nr_z(3)), &
               pp_z(nr_z(1),nr_z(2),nr_z(3)))
      if(is_periodic) then
        allocate(cc_z_0(nr_z(1),nr_z(2),nr_z(3)), &
                 pp_z_2(nr_z(1),nr_z(2),nr_z(3)))
      end if
      !$acc enter data create(aa_y,cc_y,pp_y,aa_z,cc_z,pp_z,cc_z_0,pp_z_2) async(1)
    end if
    !
    dk_g = lo-1
    if(present(lambdaxy)) then
      !
      ! factor inner rows of z-distributed systems so that they are only coupled to the boundaries:
      !
      !$acc parallel loop gang vector collapse(2) default(present) private(lxy,z,z1,z2) async(1)
      do j=1,ny
        do i=1,nx
          lxy = lambdaxy(i,j)
          !
          z1 = 1./(b(1+dk_g)+lxy+eps)
          z2 = 1./(b(2+dk_g)+lxy+eps)
          aa(i,j,1) = a(1+dk_g)*z1
          aa(i,j,2) = a(2+dk_g)*z2
          cc(i,j,1) = c(1+dk_g)*z1
          cc(i,j,2) = c(2+dk_g)*z2
          p(i,j,1) = p(i,j,1)*norm*z1
          p(i,j,2) = p(i,j,2)*norm*z2
          !
          ! elimination of lower diagonals
          !
          !$acc loop seq
          do k=3,n
            z = 1./(b(k+dk_g)+lxy-a(k+dk_g)*cc(i,j,k-1)+eps)
            p(i,j,k) = (p(i,j,k)*norm-a(k+dk_g)*p(i,j,k-1))*z
            aa(i,j,k) = -a(k+dk_g)*aa(i,j,k-1)*z
            cc(i,j,k) = c(k+dk_g)*z
          end do
          !
          ! elimination of upper diagonals
          !
          !$acc loop seq
          do k=n-2,2,-1
            p(i,j,k)  = p( i,j,k)-cc(i,j,k)*p( i,j,k+1)
            aa(i,j,k) = aa(i,j,k)-cc(i,j,k)*aa(i,j,k+1)
            cc(i,j,k) =          -cc(i,j,k)*cc(i,j,k+1)
          end do
          z = 1./(1.-aa(i,j,2)*cc(i,j,1)+eps)
          p(i,j,1) = (p(i,j,1)-cc(i,j,1)*p(i,j,2))*z
          aa(i,j,1) = aa(i,j,1)*z
          cc(i,j,1) = -cc(i,j,1)*cc(i,j,2)*z
          !
          ! gather reduced systems
          !
          aa_y(i,j,1) = aa(i,j,1); aa_y(i,j,2) = aa(i,j,n)
          cc_y(i,j,1) = cc(i,j,1); cc_y(i,j,2) = cc(i,j,n)
          pp_y(i,j,1) = p(i,j,1) ; pp_y(i,j,2) = p(i,j,n)
        end do
      end do
    else
      !
      ! factor inner rows of z-distributed systems so that they are only coupled to the boundaries:
      !
      !$acc parallel loop gang vector collapse(2) default(present) private(z,z1,z2) async(1)
      do j=1,ny
        do i=1,nx
          z1 = 1./(b(1+dk_g)+eps)
          z2 = 1./(b(2+dk_g)+eps)
          aa(i,j,1) = a(1+dk_g)*z1
          aa(i,j,2) = a(2+dk_g)*z2
          cc(i,j,1) = c(1+dk_g)*z1
          cc(i,j,2) = c(2+dk_g)*z2
          p(i,j,1) = p(i,j,1)*norm*z1
          p(i,j,2) = p(i,j,2)*norm*z2
          !
          ! elimination of lower diagonals
          !
          !$acc loop seq
          do k=3,n
            z = 1./(b(k+dk_g)-a(k+dk_g)*cc(i,j,k-1)+eps)
            p(i,j,k) = (p(i,j,k)*norm-a(k+dk_g)*p(i,j,k-1))*z
            aa(i,j,k) = -a(k+dk_g)*aa(i,j,k-1)*z
            cc(i,j,k) = c(k+dk_g)*z
          end do
          !
          ! elimination of upper diagonals
          !
          !$acc loop seq
          do k=n-2,2,-1
            p(i,j,k)  = p( i,j,k)-cc(i,j,k)*p( i,j,k+1)
            aa(i,j,k) = aa(i,j,k)-cc(i,j,k)*aa(i,j,k+1)
            cc(i,j,k) =          -cc(i,j,k)*cc(i,j,k+1)
          end do
          z = 1./(1.-aa(i,j,2)*cc(i,j,1)+eps)
          p(i,j,1) = (p(i,j,1)-cc(i,j,1)*p(i,j,2))*z
          aa(i,j,1) = aa(i,j,1)*z
          cc(i,j,1) = -cc(i,j,1)*cc(i,j,2)*z
          !
          ! gather reduced systems
          !
          aa_y(i,j,1) = aa(i,j,1); aa_y(i,j,2) = aa(i,j,n)
          cc_y(i,j,1) = cc(i,j,1); cc_y(i,j,2) = cc(i,j,n)
          pp_y(i,j,1) = p(i,j,1) ; pp_y(i,j,2) = p(i,j,n)
        end do
      end do
    end if
    !
    ! transpose to gather reduced subdomain boundary systems along z
    !
    nn   = nr_z(3)
    ny_r = nr_z(2)
    nx_r = nr_z(1)
    if(present(is_update) .and. present(aa_z_save) .and. present(cc_z_save)) then
      if(is_update) then
        is_update = .false.
#if !defined(_USE_DIEZDECOMP)
        !$acc host_data use_device(aa_y,cc_y,aa_z_save,cc_z_save,work)
#endif
        istat = cudecompTransposeYtoZ(ch,gd_ptdma,aa_y,aa_z_save,work,dtype_rp,stream=istream)
        istat = cudecompTransposeYtoZ(ch,gd_ptdma,cc_y,cc_z_save,work,dtype_rp,stream=istream)
#if !defined(_USE_DIEZDECOMP)
        !$acc end host_data
#endif
      end if
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1,nn
        do j=1,ny_r
          do i=1,nx_r
            aa_z(i,j,k) = aa_z_save(i,j,k)
            cc_z(i,j,k) = cc_z_save(i,j,k)
          end do
        end do
      end do
    else
#if !defined(_USE_DIEZDECOMP)
      !$acc host_data use_device(aa_y,cc_y,aa_z,cc_z,work)
#endif
      istat = cudecompTransposeYtoZ(ch,gd_ptdma,aa_y,aa_z,work,dtype_rp,stream=istream)
      istat = cudecompTransposeYtoZ(ch,gd_ptdma,cc_y,cc_z,work,dtype_rp,stream=istream)
#if !defined(_USE_DIEZDECOMP)
      !$acc end host_data
#endif
    end if
#if !defined(_USE_DIEZDECOMP)
    !$acc host_data use_device(pp_y,pp_z,work)
#endif
    istat = cudecompTransposeYtoZ(ch,gd_ptdma,pp_y,pp_z,work,dtype_rp,stream=istream)
#if !defined(_USE_DIEZDECOMP)
    !$acc end host_data
#endif
    !
    ! solve reduced systems
    !
    if(is_periodic) then
      nn = nn-1
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1,nn+1
        do j=1,ny_r
          do i=1,nx_r
            cc_z_0(i,j,k) = cc_z(i,j,k)
          end do
        end do
      end do
    end if
    !$acc parallel loop gang vector collapse(2) default(present) private(z) async(1)
    do j=1,ny_r
      do i=1,nx_r
        !$acc loop seq
        do k=2,nn
          z = 1./(1.-aa_z(i,j,k)*cc_z(i,j,k-1)+eps)
          pp_z(i,j,k) = (pp_z(i,j,k)-aa_z(i,j,k)*pp_z(i,j,k-1))*z
          cc_z(i,j,k) = cc_z(i,j,k)*z
        end do
        !$acc loop seq
        do k=nn-1,1,-1
          pp_z(i,j,k) = pp_z(i,j,k) - cc_z(i,j,k)*pp_z(i,j,k+1)
        end do
      end do
    end do
    if(is_periodic) then
      associate(cc_z => cc_z_0)
      !$acc parallel loop gang vector collapse(2) default(present) private(z) async(1)
      do j=1,ny_r
        do i=1,nx_r
          !$acc loop seq
          do k=1,nn
            pp_z_2(i,j,k) = 0.
          end do
          pp_z_2(i,j,1 ) = -aa_z(i,j,1 )
          pp_z_2(i,j,nn) = pp_z_2(i,j,nn) - cc_z(i,j,nn)
          !
          !$acc loop seq
          do k=2,nn
            z = 1./(1.-aa_z(i,j,k)*cc_z(i,j,k-1)+eps)
            pp_z_2(i,j,k) = (pp_z_2(i,j,k)-aa_z(i,j,k)*pp_z_2(i,j,k-1))*z
            cc_z(i,j,k) = cc_z(i,j,k)*z
          end do
          !
          !$acc loop seq
          do k=nn-1,1,-1
            pp_z_2(i,j,k) = pp_z_2(i,j,k) - cc_z(i,j,k)*pp_z_2(i,j,k+1)
          end do
          pp_z(i,j,nn+1) = (pp_z(i,j,nn+1) - cc_z(i,j,nn+1)*pp_z(  i,j,1) - aa_z(i,j,nn+1)*pp_z(  i,j,nn)) / &
                           (1.             + cc_z(i,j,nn+1)*pp_z_2(i,j,1) + aa_z(i,j,nn+1)*pp_z_2(i,j,nn)+eps)
          !$acc loop seq
          do k=1,nn
            pp_z(i,j,k) = pp_z(i,j,k) + pp_z_2(i,j,k)*pp_z(i,j,nn+1)
          end do
        end do
      end do
      end associate
    end if
    !
    ! transpose solution to the original z-distributed form
    !
#if !defined(_USE_DIEZDECOMP)
    !$acc host_data use_device(pp_z,pp_y,work)
#endif
    istat = cudecompTransposeZtoY(ch,gd_ptdma,pp_z,pp_y,work,dtype_rp,stream=istream)
#if !defined(_USE_DIEZDECOMP)
    !$acc end host_data
#endif
    !
    ! obtain final solution on the inner points
    !
    !$acc parallel loop gang vector collapse(2) default(present) private(z) async(1)
    do j=1,ny
      do i=1,nx
        p(i,j,1) = pp_y(i,j,1)
        p(i,j,n) = pp_y(i,j,2)
        !$acc loop seq
        do k=2,n-1
          p(i,j,k) = p(i,j,k) - aa(i,j,k)*p(i,j,1) - cc(i,j,k)*p(i,j,n)
        end do
      end do
    end do
  end subroutine gaussel_ptdma_gpu
  !
  subroutine gaussel_ptdma_gpu_fast_1d(nx,ny,n,lo,nh,a_g,b_g,c_g,is_periodic,norm,p)
    !
    ! distributed TDMA solver for many 1D systems on GPUs
    !
    ! original author - Rafael Diez (TU Delft)
    !
    use mod_common_cudecomp, only: gd_ptdma,ap_z_ptdma,ap_y_ptdma,work => work_ptdma
    use mod_common_cudecomp, only: buf => work
    use mod_common_mpi     , only: myid
    use mod_param, only: dims
    use mod_param, only: eps
    !
    implicit none
    integer , intent(in) :: nx,ny,n,lo,nh
    real(rp), intent(in), dimension(:) :: a_g,b_g,c_g
    logical , intent(in) :: is_periodic
    real(rp), intent(in) :: norm
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    real(rp), pointer, contiguous, dimension(:,:,:) :: pp_x,pp_y
    real(rp), allocatable, dimension(:,:,:), save :: pp_z
    real(rp), allocatable, dimension(:    ), save :: aa,bb,cc,aa_z,bb_z,cc_z,pp_z_2
    real(rp), pointer, contiguous, dimension(:,:) :: aa_all,bb_all,cc_all
    integer :: i,j,k,dk_g,nn
    integer :: islab,myslab,nranks_z,kg,llo
    integer , dimension(3) :: nr_z,nr_y
    integer :: nx_r,ny_r,nng
    integer :: istat
    !
    nr_y(:) = ap_y_ptdma%shape(:)
    nr_z(:) = ap_z_ptdma%shape(:)
    if(.not.allocated(pp_z)) then
      allocate(aa(n+1), & ! n+1 just in case one solves for a boundary normal variable in the first call
               bb(n+1), &
               cc(n+1), &
               aa_z(nr_z(3)), &
               bb_z(nr_z(3)), &
               cc_z(nr_z(3)), &
               pp_z(nr_z(1),nr_z(2),nr_z(3)))
      if(is_periodic) then
        allocate(pp_z_2(nr_z(3)))
      end if
      !$acc enter data create(aa,bb,cc,aa_z,bb_z,cc_z,pp_x,pp_y,pp_z,pp_z_2)
    end if
    !
    ! p_x <-> p_y transposes performed in-place, so that it is a no-op for (z-parallel) slab decomposition
    !
    pp_x(1:nx     ,1:ny     ,1:2      ) => buf(1:nx*ny*2         ) ! p_x <-> p_y transposes performed in place, so t
    pp_y(1:nr_y(1),1:nr_y(2),1:nr_y(3)) => buf(1:product(nr_y(:)))
    nng      = size(b_g)
    nranks_z = dims(2)
    myslab   = mod(myid,nranks_z)
    aa_all(1:n+1,0:nranks_z-1) => work(0*(n+1)*nranks_z+1:1*(n+1)*nranks_z)
    bb_all(1:n+1,0:nranks_z-1) => work(1*(n+1)*nranks_z+1:2*(n+1)*nranks_z)
    cc_all(1:n+1,0:nranks_z-1) => work(2*(n+1)*nranks_z+1:3*(n+1)*nranks_z)
    !$acc parallel loop gang default(present) private(nn,llo,k) async(1)
    do islab=0,nranks_z-1
      !
      ! manual partitioning to avoid an MPI_ALLGATHER
      ! initialization step and storing the result
      !
      ! for future reference, one can do, only once, in an initialization step:
      !  ```
      !  allocate(lo_y_all(3,0:nrank-1),hi_y_all(3,0:nrank-1))
      !  call MPI_ALLGATHER(lo_y,3,MPI_INTEGER,lo_all,3,MPI_INTEGER,MPI_COMM_WORLD)
      !  call MPI_ALLGATHER(hi_y,3,MPI_INTEGER,hi_all,3,MPI_INTEGER,MPI_COMM_WORLD)
      !
      !  lo = lo_y_all(3,islab)
      !  hi = hi_y_all(3,islab)
      !  ```
      !
      nn  = nng/nranks_z
      if(islab+1 <= mod(nng,nranks_z)) nn = nn+1
      llo = 1 + islab*nn
      if(islab+1 >  mod(nng,nranks_z)) then
        llo = llo + mod(nng,nranks_z)
      end if
      !
      !$acc loop private(kg)
      do k=1,nn
        kg = k + llo-1
        aa_all(k,islab) = a_g(kg)
        bb_all(k,islab) = b_g(kg)
        cc_all(k,islab) = c_g(kg)
      end do
      !$acc loop seq
      do k=3,nn
        bb_all(k,islab) = bb_all(k,islab) - aa_all(k,islab)/bb_all(k-1,islab)*cc_all(k-1,islab)
        aa_all(k,islab) =                 - aa_all(k,islab)/bb_all(k-1,islab)*aa_all(k-1,islab)
      end do
      !$acc loop seq
      do k=nn-2,2,-1
        aa_all(k,islab) = aa_all(k,islab) - cc_all(k,islab)/bb_all(k+1,islab)*aa_all(k+1,islab)
        cc_all(k,islab) =                 - cc_all(k,islab)/bb_all(k+1,islab)*cc_all(k+1,islab)
      end do
      if(nn > 1) then
        bb_all(1,islab) = bb_all(1,islab) - cc_all(1,islab)/bb_all(2,islab)*aa_all(2,islab)
        cc_all(1,islab) =                 - cc_all(1,islab)/bb_all(2,islab)*cc_all(2,islab)
      end if
      !
      k = 2*islab + 1
      aa_z(k  ) = aa_all(1 ,islab)
      aa_z(k+1) = aa_all(nn,islab)
      bb_z(k  ) = bb_all(1 ,islab)
      bb_z(k+1) = bb_all(nn,islab)
      cc_z(k  ) = cc_all(1 ,islab)
      cc_z(k+1) = cc_all(nn,islab)
    end do
    !
    !$acc parallel loop default(present) async(1)
    do k=1,n+1
      aa(k) = aa_all(k,myslab)
      bb(k) = bb_all(k,myslab)
      cc(k) = cc_all(k,myslab)
    end do
    !
    nn   = nr_z(3)
    ny_r = nr_z(2)
    nx_r = nr_z(1)
    dk_g = lo-1
    if(is_periodic) then
      nn = nn-1
      !$acc parallel loop default(present) async(1)
      do k=1,nn+1
        pp_z_2(k) = 0._rp
      end do
      !$acc parallel default(present) async(1)
      pp_z_2(1 ) = -aa_z(1 )
      pp_z_2(nn) = pp_z_2(nn) - cc_z(nn)
      pp_z_2(1 ) = pp_z_2(1)/bb_z(1)
      !$acc end parallel
    end if
    !
    !$acc parallel default(present) async(1)
    cc_z(1) = cc_z(1)/bb_z(1)
    !$acc loop seq
    do k=2,nn
      cc_z(k) = cc_z(k)/(bb_z(k) - aa_z(k)*cc_z(k-1))
    end do
    !$acc end parallel
    !$acc parallel loop default(present) async(1)
    do k=1,n
      bb(k) = bb(k)**(-1)
    end do
    !
    if(is_periodic) then
      !$acc parallel loop seq default(present) async(1)
      do k=2,nn
        pp_z_2(k) = (pp_z_2(k) - aa_z(k)*pp_z_2(k-1))/ &
                    (bb_z(  k) - aa_z(k)*cc_z(  k-1))
      end do
      !$acc parallel loop seq default(present) async(1)
      do k=nn-1,1,-1
        pp_z_2(k) = pp_z_2(k) - pp_z_2(k+1)*cc_z(k)
      end do
    end if
    !
    ! solve distributed TDMA problem
    !
    !$acc parallel loop gang vector collapse(2) default(present) async(1)
    do j=1,ny
      do i=1,nx
        p(i,j,1) = p(i,j,1)*norm
        p(i,j,2) = p(i,j,2)*norm
        !$acc loop seq
        do k=3,n
          p(i,j,k) = p(i,j,k)*norm - a_g(k+dk_g)*bb(k-1)*p(i,j,k-1)
        end do
        !$acc loop seq
        do k=n-2,1,-1
          p(i,j,k) = p(i,j,k)      - c_g(k+dk_g)*bb(k+1)*p(i,j,k+1)
        end do
      end do
    end do
    !
    !$acc parallel loop gang vector collapse(2) default(present) async(1)
    do j=1,ny
      do i=1,nx
        pp_x(i,j,1) = p(i,j,1)
        pp_x(i,j,2) = p(i,j,n)
      end do
    end do
    !
#if !defined(_USE_DIEZDECOMP)
    !$acc host_data use_device(pp_x,pp_y,pp_z,work)
#endif
    if(.not.is_diezdecomp_x2z_z2x_transposes) then
      istat = cudecompTransposeXtoY(ch,gd_ptdma,pp_x,pp_y,work,dtype_rp,stream=istream)
      istat = cudecompTransposeYtoZ(ch,gd_ptdma,pp_y,pp_z,work,dtype_rp,stream=istream)
    else
#if defined(_USE_DIEZDECOMP)
      istat = diezdecompTransposeXtoZ(ch,gd_ptdma,pp_x,pp_z,work,dtype_rp,stream=istream)
#endif
    end if
#if !defined(_USE_DIEZDECOMP)
    !$acc end host_data
#endif
    !
    !$acc parallel loop gang vector collapse(2) default(present) async(1)
    do j=1,ny_r
      do i=1,nx_r
        pp_z(i,j,1) = pp_z(i,j,1)/bb_z(1)
        !$acc loop seq
        do k=2,nn
          pp_z(i,j,k) = (pp_z(i,j,k) - aa_z(k)*pp_z(i,j,k-1))/ &
                        (bb_z(    k) - aa_z(k)*cc_z(    k-1))
        end do
        !$acc loop seq
        do k=nn-1,1,-1
          pp_z(i,j,k) = pp_z(i,j,k) - pp_z(i,j,k+1)*cc_z(k)
        end do
      end do
    end do
    !
    if(is_periodic) then
      !$acc parallel loop gang vector collapse(2) default(present) async(1)
      do j=1,ny_r
        do i=1,nx_r
          pp_z(i,j,nn+1) = (pp_z(i,j,nn+1) - cc_z(nn+1)*pp_z(  i,j,1) - aa_z(nn+1)*pp_z(  i,j,nn))/ &
                           (bb_z(    nn+1) + cc_z(nn+1)*pp_z_2(    1) + aa_z(nn+1)*pp_z_2(    nn))
          !$acc loop seq
          do k=1,nn
            pp_z(i,j,k) = pp_z(i,j,k) + pp_z_2(k)*pp_z(i,j,nn+1)
          end do
        end do
      end do
    end if
    !
#if !defined(_USE_DIEZDECOMP)
    !$acc host_data use_device(pp_z,pp_y,pp_x,work)
#endif
    if(.not.is_diezdecomp_x2z_z2x_transposes) then
      istat = cudecompTransposeZtoY(ch,gd_ptdma,pp_z,pp_y,work,dtype_rp,stream=istream)
      istat = cudecompTransposeYtoX(ch,gd_ptdma,pp_y,pp_x,work,dtype_rp,stream=istream)
    else
#if defined(_USE_DIEZDECOMP)
      istat = diezdecompTransposeZtoX(ch,gd_ptdma,pp_z,pp_x,work,dtype_rp,stream=istream)
#endif
    end if
#if !defined(_USE_DIEZDECOMP)
    !$acc end host_data
#endif
    !
    !$acc parallel loop gang vector collapse(2) default(present) async(1)
    do j=1,ny
      do i=1,nx
        p(i,j,1) = pp_x(i,j,1)
        p(i,j,n) = pp_x(i,j,2)
      end do
    end do
    !$acc parallel loop gang vector collapse(3) default(present) async(1)
    do k=2,n-1
      do j=1,ny
        do i=1,nx
          p(i,j,k) = (p(i,j,k) - aa(k)*p(i,j,1) &
                               - cc(k)*p(i,j,n))*bb(k)
        end do
      end do
    end do
  end subroutine gaussel_ptdma_gpu_fast_1d
  !
  subroutine solver_gaussel_z_gpu(n,ng,hi,a,b,c,bcz,c_or_f,norm,p)
    use mod_param, only: eps
    implicit none
    integer , intent(in), dimension(3) :: ng,n,hi
    real(rp), intent(in), dimension(:) :: a,b,c
    character(len=1), dimension(0:1), intent(in) :: bcz
    character(len=1), intent(in), dimension(3) :: c_or_f
    real(rp), intent(in) :: norm
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    real(rp), pointer, contiguous, dimension(:,:,:) :: px,py,pz
    integer :: i,j,k,q
    logical :: is_periodic_z
    integer, dimension(3) :: n_x,n_y,n_z,n_z_0
    integer               :: lo_z,hi_z
    logical :: is_no_decomp_z
    integer :: istat
    !
    n_z_0(:) = ap_z_0%shape(:)
    hi_z = hi(3)
    lo_z = hi(3)-n(3)+1
    if(.not.is_poisson_pcr_tdma) then
      n_x(:) = ap_x%shape(:)
      n_y(:) = ap_y%shape(:)
      n_z(:) = ap_z%shape(:)
      is_no_decomp_z = n_x(3) == n_z(3).or.ipencil_axis == 3 ! not decomposed along z: xsize(3) == ysize(3) == ng(3) when dims(2) = 1
      if(.not.is_no_decomp_z) then
        px(1:n_x(1),1:n_x(2),1:n_x(3)) => solver_buf_0(1:product(n_x(:)))
        if(cudecomp_is_t_in_place) then
          py(1:n_y(1),1:n_y(2),1:n_y(3)) => solver_buf_0(1:product(n_y(:)))
        else
          py(1:n_y(1),1:n_y(2),1:n_y(3)) => solver_buf_1(1:product(n_y(:)))
        end if
        pz(1:n_z(1),1:n_z(2),1:n_z(3)) => solver_buf_0(1:product(n_z(:)))
      end if
      !
      if(.not.is_no_decomp_z) then
        select case(ipencil_axis)
        !
        ! n.b.: since x- and y-aligned partitions are anyway sub-optimal with z-implicit diffusion,
        !       the first transposes below are not optimized to skip the copy of the
        !       input array `p` to `px`/`py`, as done above in `solver_gpu()`;
        !       the same holds for the reciprocate operations after the solution is obtained
        !
        case(1)
          !$acc parallel loop collapse(3) default(present) async(1)
          !$OMP parallel do   collapse(3) DEFAULT(shared)
          do k=1,n(3)
            do j=1,n(2)
              do i=1,n(1)
                px(i,j,k) = p(i,j,k)
              end do
            end do
          end do
#if !defined(_USE_DIEZDECOMP)
          !$acc host_data use_device(px,py,pz,work)
#endif
          if(.not.is_diezdecomp_x2z_z2x_transposes) then
            istat = cudecompTransposeXtoY(ch,gd,px,py,work,dtype_rp,stream=istream)
            istat = cudecompTransposeYtoZ(ch,gd,py,pz,work,dtype_rp,stream=istream)
          else
#if defined(_USE_DIEZDECOMP)
            istat = diezdecompTransposeXtoZ(ch,gd,px,pz,work,dtype_rp,stream=istream)
#endif
          end if
#if !defined(_USE_DIEZDECOMP)
          !$acc end host_data
#endif
        case(2)
          !
          ! transpose p -> py to axis-contiguous layout
          !
          !$acc parallel loop collapse(3) default(present) async(1)
          !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
          do k=1,n(3)
            do j=1,n(2)
              do i=1,n(1)
                py(j,k,i) = p(i,j,k)
              end do
            end do
          end do
#if !defined(_USE_DIEZDECOMP)
          !$acc host_data use_device(py,pz,work)
#endif
          istat = cudecompTransposeYtoZ(ch,gd,py,pz,work,dtype_rp,stream=istream)
#if !defined(_USE_DIEZDECOMP)
          !$acc end host_data
#endif
        case(3)
        end select
      end if
    end if
    !
    q = merge(1,0,c_or_f(3) == 'f'.and.bcz(1) == 'D'.and.hi_z == ng(3))
    is_periodic_z = bcz(0)//bcz(1) == 'PP'
    if(.not.is_no_decomp_z) then
      if(.not.is_poisson_pcr_tdma) then
        call gaussel_gpu(n_z_0(1),n_z_0(2),n_z_0(3)-q,0,a,b,c,is_periodic_z,norm,pz,work,pz_aux_1)
      else
        call gaussel_ptdma_gpu_fast_1d(n(1),n(2),n(3)-q,lo_z,1,a,b,c,is_periodic_z,norm,p)
      end if
    else
      call gaussel_gpu(n(1),n(2),n(3)-q,1,a,b,c,is_periodic_z,norm,p,work,pz_aux_1)
    end if
    !
    if(.not.is_poisson_pcr_tdma .and. .not.is_no_decomp_z) then
      select case(ipencil_axis)
      case(1)
#if !defined(_USE_DIEZDECOMP)
        !$acc host_data use_device(pz,py,px,work)
#endif
        if(.not.is_diezdecomp_x2z_z2x_transposes) then
          istat = cudecompTransposeZtoY(ch,gd,pz,py,work,dtype_rp,stream=istream)
          istat = cudecompTransposeYtoX(ch,gd,py,px,work,dtype_rp,stream=istream)
        else
#if defined(_USE_DIEZDECOMP)
          istat = diezdecompTransposeZtoX(ch,gd,pz,px,work,dtype_rp,stream=istream)
#endif
        end if
#if !defined(_USE_DIEZDECOMP)
        !$acc end host_data
#endif
        !$acc parallel loop collapse(3) default(present) async(1)
        !$OMP parallel do   collapse(3) DEFAULT(shared)
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              p(i,j,k) = px(i,j,k)
            end do
          end do
        end do
      case(2)
#if !defined(_USE_DIEZDECOMP)
        !$acc host_data use_device(pz,py,work)
#endif
        istat = cudecompTransposeZtoY(ch,gd,pz,py,work,dtype_rp,stream=istream)
#if !defined(_USE_DIEZDECOMP)
        !$acc end host_data
#endif
        !
        ! transpose py -> p to default layout
        !
        !$acc parallel loop collapse(3) default(present) async(1)
        !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              p(i,j,k) = py(j,k,i)
            end do
          end do
        end do
      case(3)
      end select
    end if
  end subroutine solver_gaussel_z_gpu
#if 0
  subroutine gaussel_ptdma_gpu_fast(nx,ny,n,lo,nh,a_g,b_g,c_g,is_periodic,norm,p,is_update,lambdaxy,aa,bb,cc,aa_z,bb_z,cc_z,pp_z_2)
    !
    ! distributed TDMA solver using pre-computed coefficients
    !
    ! original author - Rafael Diez (TU Delft)
    !
    use mod_common_cudecomp, only: gd_ptdma,ap_z_ptdma,work => work_ptdma
    use mod_param          , only: eps
    !
    implicit none
    integer , intent(in) :: nx,ny,n,lo,nh
    real(rp), intent(in), dimension(:) :: a_g,b_g,c_g
    logical , intent(in) :: is_periodic
    real(rp), intent(in) :: norm
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    logical , intent(inout) :: is_update
    real(rp), intent(in), dimension(:,:), optional :: lambdaxy
    real(rp), intent(inout), dimension(:,:,:), optional :: aa,bb,cc,aa_z,bb_z,cc_z,pp_z_2
    real(rp), allocatable, dimension(:,:,:), save :: aa_y,bb_y,cc_y,pp_y,pp_z
    real(rp) :: z,z1,z2,lxy
    integer :: i,j,k,dk_g
    integer , dimension(3) :: nr_z
    integer :: nx_r,ny_r,nn
    integer :: istat
    !
    nr_z(:) = ap_z_ptdma%shape(:)
    if(.not.allocated(pp_y)) then
      allocate(aa_y(nx,ny,2), &
               bb_y(nx,ny,2), &
               cc_y(nx,ny,2), &
               pp_y(nx,ny,2), &
               pp_z(nr_z(1),nr_z(2),nr_z(3)))
      !$acc enter data create(aa_y,bb_y,cc_y,pp_y,pp_z) async(1)
    end if
    nn   = nr_z(3)
    ny_r = nr_z(2)
    nx_r = nr_z(1)
    if(is_periodic) nn = nn-1
    dk_g = lo-1
    !
    if(is_update) then
      is_update = .false.
      !
      ! update pre-computed coefficients
      !
      !$acc parallel loop gang vector collapse(2) default(present) async(1)
      do j=1,ny
        do i=1,nx
          do k=1,n
            aa(i,j,k) = a_g(k+dk_g)
            bb(i,j,k) = b_g(k+dk_g) + lambdaxy(i,j)
            cc(i,j,k) = c_g(k+dk_g)
          end do
        end do
      end do
      !
      !$acc parallel loop gang vector collapse(2) default(present) async(1)
      do j=1,ny
        do i=1,nx
          !$acc loop seq
          do k=3,n
            bb(i,j,k) = bb(i,j,k) - aa(i,j,k)/bb(i,j,k-1)*cc(i,j,k-1)
            aa(i,j,k) =           - aa(i,j,k)/bb(i,j,k-1)*aa(i,j,k-1)
          end do
        end do
      end do
      !$acc parallel loop collapse(2) default(present) async(1)
      do j=1,ny
        do i=1,nx
          !$acc loop seq
          do k=n-2,2,-1
            aa(i,j,k) = aa(i,j,k) - cc(i,j,k)/bb(i,j,k+1)*aa(i,j,k+1)
            cc(i,j,k) =           - cc(i,j,k)/bb(i,j,k+1)*cc(i,j,k+1)
          end do
        end do
      end do
      if(n > 1) then
        !$acc parallel loop gang vector collapse(2) default(present) async(1)
        do j=1,ny
          do i=1,nx
            bb(i,j,1) = bb(i,j,1) - cc(i,j,1)/bb(i,j,2)*aa(i,j,2)
            cc(i,j,1) =           - cc(i,j,1)/bb(i,j,2)*cc(i,j,2)
          end do
        end do
      end if
      !
      !$acc parallel loop gang vector collapse(2) default(present) async(1)
      do j=1,ny
        do i=1,nx
          aa_y(i,j,1) = aa(i,j,1)
          aa_y(i,j,2) = aa(i,j,n)
          bb_y(i,j,1) = bb(i,j,1)
          bb_y(i,j,2) = bb(i,j,n)
          cc_y(i,j,1) = cc(i,j,1)
          cc_y(i,j,2) = cc(i,j,n)
        end do
      end do
      !
#if !defined(_USE_DIEZDECOMP)
      !$acc host_data use_device(aa_y,bb_y,cc_y,aa_z,bb_z,cc_z,work)
#endif
      istat = cudecompTransposeYtoZ(ch,gd_ptdma,aa_y,aa_z,work,dtype_rp,stream=istream)
      istat = cudecompTransposeYtoZ(ch,gd_ptdma,bb_y,bb_z,work,dtype_rp,stream=istream)
      istat = cudecompTransposeYtoZ(ch,gd_ptdma,cc_y,cc_z,work,dtype_rp,stream=istream)
#if !defined(_USE_DIEZDECOMP)
      !$acc end host_data
#endif
      !
      if(is_periodic) then
        !$acc parallel loop collapse(3) default(present) async(1)
        do k=1,nn+1
          do j=1,ny_r
            do i=1,nx_r
              pp_z_2(i,j,k) = 0._rp
            end do
          end do
        end do
        !$acc parallel loop gang vector collapse(2) default(present) async(1)
        do j=1,ny_r
          do i=1,nx_r
            pp_z_2(i,j,1 ) = -aa_z(i,j,1 )
            pp_z_2(i,j,nn) = pp_z_2(i,j,nn) - cc_z(i,j,nn)
            pp_z_2(i,j,1 ) = pp_z_2(i,j,1)/bb_z(i,j,1)
          end do
        end do
      end if
      !
      !$acc parallel loop gang vector collapse(2) default(present) async(1)
      do j=1,ny_r
        do i=1,nx_r
          cc_z(i,j,1) = cc_z(i,j,1)/bb_z(i,j,1)
          !$acc loop seq
          do k=2,nn
            cc_z(i,j,k) = cc_z(i,j,k)/(bb_z(i,j,k) - aa_z(i,j,k)*cc_z(i,j,k-1))
          end do
        end do
      end do
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1,n
        do j=1,ny
          do i=1,nx
            bb(i,j,k) = bb(i,j,k)**(-1)
          end do
        end do
      end do
      !
      if(is_periodic) then
        !$acc parallel loop gang vector collapse(2) default(present) async(1)
        do j=1,ny_r
          do i=1,nx_r
            !$acc loop seq
            do k=2,nn
              pp_z_2(i,j,k) = (pp_z_2(i,j,k) - aa_z(i,j,k)*pp_z_2(i,j,k-1))/ &
                              (bb_z(  i,j,k) - aa_z(i,j,k)*cc_z(  i,j,k-1))
            end do
            !$acc loop seq
            do k=nn-1,1,-1
              pp_z_2(i,j,k) = pp_z_2(i,j,k) - pp_z_2(i,j,k+1)*cc_z(i,j,k)
            end do
          end do
        end do
      end if
    end if
    !
    ! solve distributed TDMA problem
    !
    !$acc parallel loop gang vector collapse(2) default(present) async(1)
    do j=1,ny
      do i=1,nx
        p(i,j,1) = p(i,j,1)*norm
        p(i,j,2) = p(i,j,2)*norm
        !$acc loop seq
        do k=3,n
          p(i,j,k) = p(i,j,k)*norm - a_g(k+dk_g)*bb(i,j,k-1)*p(i,j,k-1)
        end do
        !$acc loop seq
        do k=n-2,1,-1
          p(i,j,k) = p(i,j,k)      - c_g(k+dk_g)*bb(i,j,k+1)*p(i,j,k+1)
        end do
      end do
    end do
    !
    !$acc parallel loop gang vector collapse(2) default(present) async(1)
    do j=1,ny
      do i=1,nx
        pp_y(i,j,1) = p(i,j,1)
        pp_y(i,j,2) = p(i,j,n)
      end do
    end do
    !
#if !defined(_USE_DIEZDECOMP)
    !$acc host_data use_device(pp_y,pp_z,work)
#endif
    istat = cudecompTransposeYtoZ(ch,gd_ptdma,pp_y,pp_z,work,dtype_rp,stream=istream)
#if !defined(_USE_DIEZDECOMP)
    !$acc end host_data
#endif
    !
    !$acc parallel loop gang vector collapse(2) default(present) async(1)
    do j=1,ny_r
      do i=1,nx_r
        pp_z(i,j,1) = pp_z(i,j,1)/bb_z(i,j,1)
        !$acc loop seq
        do k=2,nn
          pp_z(i,j,k) = (pp_z(i,j,k) - aa_z(i,j,k)*pp_z(i,j,k-1))/ &
                        (bb_z(i,j,k) - aa_z(i,j,k)*cc_z(i,j,k-1))
        end do
        !$acc loop seq
        do k=nn-1,1,-1
          pp_z(i,j,k) = pp_z(i,j,k) - pp_z(i,j,k+1)*cc_z(i,j,k)
        end do
      end do
    end do
    !
    if(is_periodic) then
      !$acc parallel loop gang vector collapse(2) default(present) async(1)
      do j=1,ny_r
        do i=1,nx_r
          pp_z(i,j,nn+1) = (pp_z(i,j,nn+1)*norm - cc_z(i,j,nn+1)*pp_z(  i,j,1) - aa_z(i,j,nn+1)*pp_z(  i,j,nn))/ &
                           (bb_z(i,j,nn+1)      + cc_z(i,j,nn+1)*pp_z_2(i,j,1) + aa_z(i,j,nn+1)*pp_z_2(i,j,nn))
          !$acc loop seq
          do k=1,nn
            pp_z(i,j,k) = pp_z(i,j,k) + pp_z_2(i,j,k)*pp_z(i,j,nn+1)
          end do
        end do
      end do
    end if
    !
#if !defined(_USE_DIEZDECOMP)
    !$acc host_data use_device(pp_z,pp_y,work)
#endif
    istat = cudecompTransposeZtoY(ch,gd_ptdma,pp_z,pp_y,work,dtype_rp,stream=istream)
#if !defined(_USE_DIEZDECOMP)
    !$acc end host_data
#endif
    !
    !$acc parallel loop gang vector collapse(2) default(present) async(1)
    do j=1,ny
      do i=1,nx
        p(i,j,1) = pp_y(i,j,1)
        p(i,j,n) = pp_y(i,j,2)
      end do
    end do
    !$acc parallel loop gang vector collapse(3) default(present) async(1)
    do k=2,n-1
      do j=1,ny
        do i=1,nx
          p(i,j,k) = (p(i,j,k) - aa(i,j,k)*p(i,j,1) &
                               - cc(i,j,k)*p(i,j,n))*bb(i,j,k)
        end do
      end do
    end do
  end subroutine gaussel_ptdma_gpu_fast
#endif
#endif
end module mod_solver_gpu
