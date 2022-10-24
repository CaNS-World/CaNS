! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_solver_gpu
#if defined(_OPENACC)
  use, intrinsic :: iso_c_binding, only: C_PTR
  use cudecomp
  use mod_fft  ,  only: signal_processing,fftf_gpu,fftb_gpu
  use mod_types
  use mod_common_mpi     , only: ipencil_axis
  use mod_common_cudecomp, only: dtype_rp => cudecomp_real_rp, &
                                 cudecomp_is_t_in_place, &
                                 solver_buf_0,solver_buf_1, &
                                 pz_aux_1,pz_aux_2,work, &
                                 ap_x   => ap_x_poi, &
                                 ap_y   => ap_y_poi, &
                                 ap_z   => ap_z_poi, &
                                 ap_z_0 => ap_z    , &
                                 ch => handle,gd => gd_poi, &
                                 istream => istream_acc_queue_1
  implicit none
  private
  public solver_gpu
#if defined(_IMPDIFF_1D)
  public solver_gaussel_z_gpu
#endif
  contains
  subroutine solver_gpu(n,ng,arrplan,normfft,lambdaxy,a,b,c,bc,c_or_f,p)
    implicit none
    integer , intent(in), dimension(3) :: n,ng
    integer , intent(in), dimension(2,2) :: arrplan
    real(rp), intent(in) :: normfft
    real(rp), intent(in), dimension(:,:) :: lambdaxy
    real(rp), intent(in), dimension(:) :: a,b,c
    character(len=1), dimension(0:1,3), intent(in) :: bc
    character(len=1), intent(in), dimension(3) :: c_or_f
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    real(rp), pointer, contiguous, dimension(:,:,:) :: px,py,pz
    integer :: q
    integer, dimension(3) :: n_x,n_y,n_z,n_z_0
    integer :: istat
    !
    n_z_0(:) = ap_z_0%shape(:)
    n_x(:) = ap_x%shape(:)
    n_y(:) = ap_y%shape(:)
    n_z(:) = ap_z%shape(:)
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
      !$acc kernels default(present) async(1)
      !$OMP PARALLEL WORKSHARE
      px(1:n(1),1:n(2),1:n(3)) = p(1:n(1),1:n(2),1:n(3))
      !$OMP END PARALLEL WORKSHARE
      !$acc end kernels
    case(2)
      block
        integer :: i,j,k
        !
        ! transpose p -> py to axis-contiguous layout
        !
        !$acc parallel loop collapse(3) default(present) async(1)
        !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(NONE) &
        !$OMP SHARED(n,p,py)
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              py(j,k,i) = p(i,j,k)
            end do
          end do
        end do
      end block
      !$acc host_data use_device(py,px,work)
      istat = cudecompTransposeYtoX(ch,gd,py,px,work,dtype_rp,stream=istream)
      !$acc end host_data
    case(3)
      !$acc kernels default(present) async(1)
      !$OMP PARALLEL WORKSHARE
      pz(1:n(1),1:n(2),1:n(3)) = p(1:n(1),1:n(2),1:n(3))
      !$OMP END PARALLEL WORKSHARE
      !$acc end kernels
      !$acc host_data use_device(pz,py,px,work)
      istat = cudecompTransposeZtoY(ch,gd,pz,py,work,dtype_rp,stream=istream)
      istat = cudecompTransposeYtoX(ch,gd,py,px,work,dtype_rp,stream=istream)
      !$acc end host_data
    end select
    !
    call signal_processing(0,'F',bc(0,1)//bc(1,1),c_or_f(1),ng(1),n_x,1,px)
    call fftf_gpu(arrplan(1,1),px)
    call signal_processing(1,'F',bc(0,1)//bc(1,1),c_or_f(1),ng(1),n_x,1,px)
    !
    !$acc host_data use_device(px,py,work)
    istat = cudecompTransposeXtoY(ch,gd,px,py,work,dtype_rp,stream=istream)
    !$acc end host_data
    !
    call signal_processing(0,'F',bc(0,2)//bc(1,2),c_or_f(2),ng(2),n_y,1,py)
    call fftf_gpu(arrplan(1,2),py)
    call signal_processing(1,'F',bc(0,2)//bc(1,2),c_or_f(2),ng(2),n_y,1,py)
    !
    !$acc host_data use_device(py,pz,work)
    istat = cudecompTransposeYtoZ(ch,gd,py,pz,work,dtype_rp,stream=istream)
    !$acc end host_data
    !
    q = 0
    if(c_or_f(3) == 'f'.and.bc(1,3) == 'D') q = 1
    if(bc(0,3)//bc(1,3) == 'PP') then
      call gaussel_periodic_gpu(n_z_0(1),n_z_0(2),n_z_0(3)-q,0,a,b,c,pz,work,pz_aux_1,pz_aux_2,lambdaxy)
    else
      call gaussel_gpu(         n_z_0(1),n_z_0(2),n_z_0(3)-q,0,a,b,c,pz,work,lambdaxy)
    end if
    !
    !$acc host_data use_device(pz,py,work)
    istat = cudecompTransposeZtoY(ch,gd,pz,py,work,dtype_rp,stream=istream)
    !$acc end host_data
    !
    call signal_processing(0,'B',bc(0,2)//bc(1,2),c_or_f(2),ng(2),n_y,1,py)
    call fftb_gpu(arrplan(2,2),py)
    call signal_processing(1,'B',bc(0,2)//bc(1,2),c_or_f(2),ng(2),n_y,1,py)
    !
    !$acc host_data use_device(py,px,work)
    istat = cudecompTransposeYtoX(ch,gd,py,px,work,dtype_rp,stream=istream)
    !$acc end host_data
    !
    call signal_processing(0,'B',bc(0,1)//bc(1,1),c_or_f(1),ng(1),n_x,1,px)
    call fftb_gpu(arrplan(2,1),px)
    call signal_processing(1,'B',bc(0,1)//bc(1,1),c_or_f(1),ng(1),n_x,1,px)
    !
    select case(ipencil_axis)
    case(1)
      !$acc kernels default(present) async(1)
      !$OMP PARALLEL WORKSHARE
      p(1:n(1),1:n(2),1:n(3)) = px(1:n(1),1:n(2),1:n(3))*normfft
      !$OMP END PARALLEL WORKSHARE
      !$acc end kernels
    case(2)
      !$acc host_data use_device(px,py,work)
      istat = cudecompTransposeXtoY(ch,gd,px,py,work,dtype_rp,stream=istream)
      !$acc end host_data
      block
        integer :: i,j,k
        !
        ! transpose py -> p to default layout
        !
        !$acc parallel loop collapse(3) default(present) async(1)
        !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(NONE) &
        !$OMP SHARED(n,p,py)
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              p(i,j,k) = py(j,k,i)*normfft
            end do
          end do
        end do
      end block
    case(3)
      !$acc host_data use_device(px,py,pz,work)
      istat = cudecompTransposeXtoY(ch,gd,px,py,work,dtype_rp,stream=istream)
      istat = cudecompTransposeYtoZ(ch,gd,py,pz,work,dtype_rp,stream=istream)
      !$acc end host_data
      !$acc kernels default(present) async(1)
      !$OMP PARALLEL WORKSHARE
      p(1:n(1),1:n(2),1:n(3)) = pz(1:n(1),1:n(2),1:n(3))*normfft
      !$OMP END PARALLEL WORKSHARE
      !$acc end kernels
    end select
  end subroutine solver_gpu
  !
  subroutine gaussel_gpu(nx,ny,n,nh,a,b,c,p,d,lambdaxy)
    use mod_param, only: eps
    implicit none
    integer , intent(in) :: nx,ny,n,nh
    real(rp), intent(in), dimension(:) :: a,b,c
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    real(rp),                dimension(nx,ny,n) :: d
    real(rp), intent(in), dimension(:,:), optional :: lambdaxy
    real(rp) :: z,lxy
    integer :: i,j,k
    !
    !solve tridiagonal system
    !
    if(present(lambdaxy)) then
      !$acc parallel loop gang vector collapse(2) default(present) private(lxy,z) async(1)
      do j=1,ny
        do i=1,nx
          lxy = lambdaxy(i,j)
          !
          ! call dgtsv_homebrewed(n,a,bb,c,p)
          !
          z = 1./(b(1)+lxy+eps)
          d(i,j,1) = c(1)*z
          p(i,j,1) = p(i,j,1)*z
          !$acc loop seq
          do k=2,n
            z = 1./(b(k)+lxy-a(k)*d(i,j,k-1)+eps)
            d(i,j,k) = c(k)*z
            p(i,j,k) = (p(i,j,k)-a(k)*p(i,j,k-1))*z
          end do
          !
          ! backward substitution
          !
          !$acc loop seq
          do k=n-1,1,-1
            p(i,j,k) = p(i,j,k) - d(i,j,k)*p(i,j,k+1)
          end do
        end do
      end do
    else
      !$acc parallel loop gang vector collapse(2) default(present) private(z) async(1)
      do j=1,ny
        do i=1,nx
          !
          ! call dgtsv_homebrewed(n,a,b,c,p)
          !
          z = 1./(b(1)+eps)
          d(i,j,1) = c(1)*z
          p(i,j,1) = p(i,j,1)*z
          !$acc loop seq
          do k=2,n
            z = 1./(b(k)-a(k)*d(i,j,k-1)+eps)
            d(i,j,k) = c(k)*z
            p(i,j,k) = (p(i,j,k)-a(k)*p(i,j,k-1))*z
          end do
          !
          ! backward substitution
          !
          !$acc loop seq
          do k=n-1,1,-1
            p(i,j,k) = p(i,j,k) - d(i,j,k)*p(i,j,k+1)
          end do
        end do
      end do
    end if
  end subroutine gaussel_gpu
  !
  subroutine gaussel_periodic_gpu(nx,ny,n,nh,a,b,c,p,d,p1,p2,lambdaxy)
    use mod_param, only: eps
    implicit none
    integer , intent(in) :: nx,ny,n,nh
    real(rp), intent(in), dimension(:) :: a,b,c
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    real(rp),                dimension(nx,ny,n) :: d,p1,p2
    real(rp), intent(in), dimension(:,:), optional :: lambdaxy
    real(rp) :: z,lxy
    integer :: i,j,k
    !
    ! solve tridiagonal system
    !
    if(present(lambdaxy)) then
      !$acc parallel loop gang vector collapse(2) default(present) private(lxy,z) async(1)
      do j=1,ny
        do i=1,nx
          lxy = lambdaxy(i,j)
          !$acc loop seq
          do k=1,n-1
            p1(i,j,k) = p(i,j,k)
          end do
          !
          ! call dgtsv_homebrewed(n-1,a,bb,c,p1)
          !
          z = 1./(b(1)+lxy+eps)
          d(i,j,1)  = c(1)*z
          p1(i,j,1) = p1(i,j,1)*z
          !$acc loop seq
          do k=2,n-1
            z         = 1./(b(k)+lxy-a(k)*d(i,j,k-1)+eps)
            d(i,j,k)  = c(k)*z
            p1(i,j,k) = (p1(i,j,k)-a(k)*p1(i,j,k-1))*z
          end do
          !
          ! backward substitution
          !
          !$acc loop seq
          do k=n-2,1,-1
            p1(i,j,k) = p1(i,j,k) - d(i,j,k)*p1(i,j,k+1)
          end do
          !
          !$acc loop seq
          do k=1,n
            p2(i,j,k) = 0.
          end do
          p2(i,j,1  ) = -a(1  )
          p2(i,j,n-1) = -c(n-1)
          !
          ! call dgtsv_homebrewed(n-1,a,bb,c,p2)
          !
          z = 1./(b(1)+lxy+eps)
          d( i,j,1) = c(1)*z
          p2(i,j,1) = p2(i,j,1)*z
          !$acc loop seq
          do k=2,n-1
            z        = 1./(b(k)+lxy-a(k)*d(i,j,k-1)+eps)
            d(i,j,k) = c(k)*z
            p2(i,j,k) = (p2(i,j,k)-a(k)*p2(i,j,k-1))*z
          end do
          !
          ! backward substitution
          !
          !$acc loop seq
          do k=n-2,1,-1
            p2(i,j,k) = p2(i,j,k) - d(i,j,k)*p2(i,j,k+1)
          end do
          p(i,j,n) = (p(i,j,n)       - c(n)*p1(i,j,1) - a(n)*p1(i,j,n-1)) / &
                     (b(    n) + lxy + c(n)*p2(i,j,1) + a(n)*p2(i,j,n-1)+eps)
          !$acc loop seq
          do k=1,n-1
            p(i,j,k) = p1(i,j,k) + p2(i,j,k)*p(i,j,n)
          end do
        end do
      end do
    else
      !$acc parallel loop gang vector collapse(2) default(present) private(z) async(1)
      do j=1,ny
        do i=1,nx
          !$acc loop seq
          do k=1,n-1
            p1(i,j,k) = p(i,j,k)
          end do
          !
          ! call dgtsv_homebrewed(n-1,a,b,c,p1)
          !
          z = 1./(b(1)+eps)
          d(i,j,1)  = c(1)*z
          p1(i,j,1) = p1(i,j,1)*z
          !$acc loop seq
          do k=2,n-1
            z         = 1./(b(k)-a(k)*d(i,j,k-1)+eps)
            d(i,j,k)  = c(k)*z
            p1(i,j,k) = (p1(i,j,k)-a(k)*p1(i,j,k-1))*z
          end do
          !
          ! backward substitution
          !
          !$acc loop seq
          do k=n-2,1,-1
            p1(i,j,k) = p1(i,j,k) - d(i,j,k)*p1(i,j,k+1)
          end do
          !
          !!$acc loop seq
          !do k=1,n
          !  p2(i,j,k) = 0.
          !end do
          p2(i,j,1  ) = -a(1  )
          p2(i,j,n-1) = -c(n-1)
          !
          ! call dgtsv_homebrewed(n-1,a,b,c,p2)
          !
          z = 1./(b(1)+eps)
          d( i,j,1) = c(1)*z
          p2(i,j,1) = p2(i,j,1)*z
          !$acc loop seq
          do k=2,n-1
            z        = 1./(b(k)-a(k)*d(i,j,k-1)+eps)
            d(i,j,k) = c(k)*z
           !p2(i,j,k) = (p2(i,j,k)-a(k)*p2(i,j,k-1))*z
            p2(i,j,k) = (         -a(k)*p2(i,j,k-1))*z
          end do
          !
          ! backward substitution
          !
          !$acc loop seq
          do k=n-2,1,-1
            p2(i,j,k) = p2(i,j,k) - d(i,j,k)*p2(i,j,k+1)
          end do
          p(i,j,n) = (p(i,j,n) - c(n)*p1(i,j,1) - a(n)*p1(i,j,n-1)) / &
                     (b(    n) + c(n)*p2(i,j,1) + a(n)*p2(i,j,n-1)+eps)
          !$acc loop seq
          do k=1,n-1
            p(i,j,k) = p1(i,j,k) + p2(i,j,k)*p(i,j,n)
          end do
        end do
      end do
    end if
  end subroutine gaussel_periodic_gpu
  !
#if defined(_IMPDIFF_1D)
  subroutine solver_gaussel_z_gpu(n,a,b,c,bcz,c_or_f,p)
    use mod_param, only: eps
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(:) :: a,b,c
    character(len=1), dimension(0:1), intent(in) :: bcz
    character(len=1), intent(in), dimension(3) :: c_or_f
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    real(rp), pointer, contiguous, dimension(:,:,:) :: px,py,pz
    integer :: q
    integer, dimension(3) :: n_x,n_y,n_z,n_z_0
    integer :: istat
    !
    n_z_0(:) = ap_z_0%shape(:)
#if !defined(_DECOMP_Z)
    n_x(:) = ap_x%shape(:)
    n_y(:) = ap_y%shape(:)
    n_z(:) = ap_z%shape(:)
    px(1:n_x(1),1:n_x(2),1:n_x(3)) => solver_buf_0(1:product(n_x(:)))
    if(cudecomp_is_t_in_place) then
      py(1:n_y(1),1:n_y(2),1:n_y(3)) => solver_buf_0(1:product(n_y(:)))
    else
      py(1:n_y(1),1:n_y(2),1:n_y(3)) => solver_buf_1(1:product(n_y(:)))
    end if
    pz(1:n_z(1),1:n_z(2),1:n_z(3)) => solver_buf_0(1:product(n_z(:)))
#endif
    select case(ipencil_axis)
    case(1)
      !$acc kernels default(present) async(1)
      !$OMP PARALLEL WORKSHARE
      px(:,:,:) = p(1:n(1),1:n(2),1:n(3))
      !$OMP END PARALLEL WORKSHARE
      !$acc end kernels
      !$acc host_data use_device(px,py,pz,work)
      istat = cudecompTransposeXtoY(ch,gd,px,py,work,dtype_rp,stream=istream)
      istat = cudecompTransposeYtoZ(ch,gd,py,pz,work,dtype_rp,stream=istream)
      !$acc end host_data
    case(2)
      block
        integer :: i,j,k
        !
        ! transpose p -> py to axis-contiguous layout
        !
        !$acc parallel loop collapse(3) default(present) async(1)
        !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(NONE) &
        !$OMP SHARED(n,p,py)
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              py(j,k,i) = p(i,j,k)
            end do
          end do
        end do
      end block
      !$acc host_data use_device(py,pz,work)
      istat = cudecompTransposeYtoZ(ch,gd,py,pz,work,dtype_rp,stream=istream)
      !$acc end host_data
    case(3)
    end select
    !
    if(ipencil_axis /= 3) then
      q = 0
      if(c_or_f(3) == 'f'.and.bcz(1) == 'D') q = 1
      if(bcz(0)//bcz(1) == 'PP') then
        call gaussel_periodic_gpu(n_z_0(1),n_z_0(2),n_z_0(3)-q,0,a,b,c,pz,work,pz_aux_1,pz_aux_2)
      else
        call gaussel_gpu(         n_z_0(1),n_z_0(2),n_z_0(3)-q,0,a,b,c,pz,work)
      end if
    else
      q = 0
      if(c_or_f(3) == 'f'.and.bcz(1) == 'D') q = 1
      if(bcz(0)//bcz(1) == 'PP') then
        call gaussel_periodic_gpu(n_z_0(1),n_z_0(2),n_z_0(3)-q,1,a,b,c,p ,work,pz_aux_1,pz_aux_2)
      else
        call gaussel_gpu(         n_z_0(1),n_z_0(2),n_z_0(3)-q,1,a,b,c,p ,work)
      end if
    end if
    !
    select case(ipencil_axis)
    case(1)
      !$acc host_data use_device(pz,py,px,work)
      istat = cudecompTransposeZtoY(ch,gd,pz,py,work,dtype_rp,stream=istream)
      istat = cudecompTransposeYtoX(ch,gd,py,px,work,dtype_rp,stream=istream)
      !$acc end host_data
      !$acc kernels default(present) async(1)
      !$OMP PARALLEL WORKSHARE
      p(1:n(1),1:n(2),1:n(3)) = px(:,:,:)
      !$OMP END PARALLEL WORKSHARE
      !$acc end kernels
    case(2)
      !$acc host_data use_device(pz,py,work)
      istat = cudecompTransposeZtoY(ch,gd,pz,py,work,dtype_rp,stream=istream)
      !$acc end host_data
      block
        integer :: i,j,k
        !
        ! transpose py -> p to default layout
        !
        !$acc parallel loop collapse(3) default(present) async(1)
        !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(NONE) &
        !$OMP SHARED(n,p,py)
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              p(i,j,k) = py(j,k,i)
            end do
          end do
        end do
      end block
    case(3)
    end select
  end subroutine solver_gaussel_z_gpu
#endif
#endif
end module mod_solver_gpu
