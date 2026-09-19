! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_solver
  use, intrinsic :: iso_c_binding, only: C_PTR
  use decomp_2d
  use mod_fft       , only: fft,prep_dctviii,posp_dctviii
  use mod_param     , only: ipencil_axis,is_poisson_dtdma
  use mod_types
  implicit none
  private
  public solver,solver_gaussel_z
  contains
  subroutine solver(n,ng,arrplan,normfft,lambdaxy,a,b,c,bc,c_or_f,p,is_dtdma_update,aa_z,cc_z)
    !
    ! note: some of the transposes below are suboptimal in slab decompositions,
    ! as they would be a no-op if done in-place (e.g., `px = py` for xy slabs)
    !
    implicit none
    integer , intent(in), dimension(3) :: n,ng
    type(C_PTR), intent(in), dimension(2,2) :: arrplan
    real(rp), intent(in) :: normfft
    real(rp), intent(in), dimension(:,:) :: lambdaxy
    real(rp), intent(in), dimension(:) :: a,b,c
    character(len=1), dimension(0:1,3), intent(in) :: bc
    character(len=1), intent(in), dimension(3) :: c_or_f
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    logical , intent(inout), optional :: is_dtdma_update
    real(rp), intent(inout), dimension(:,:,:), optional :: aa_z,cc_z
    real(rp), allocatable, target, dimension(:,:,:) :: px,py,pz
    real(rp), pointer, contiguous :: pfft(:,:,:)
    integer :: q,idir
    logical :: is_periodic_z,is_mixed(2)
    integer, dimension(3) :: n_z,hi_z
    logical :: is_dtdma_update_
    real(rp) :: norm
    !
    norm = normfft
    do idir=1,2
      is_mixed(idir) = (c_or_f(idir) == 'f').and.(any(bc(0,idir)//bc(1,idir) == ['ND','DN']))
    end do
    !
    is_dtdma_update_ = .true.
    if(present(is_dtdma_update)) is_dtdma_update_ = is_dtdma_update
    n_z(:)  = zsize(:)
    hi_z(:) = zend(:)
    if(is_poisson_dtdma) then
      n_z(:)  = ysize(:)
      hi_z(:) = yend(:)
    end if
    allocate(px(xsize(1),xsize(2),xsize(3)))
    allocate(py(ysize(1),ysize(2),ysize(3)))
    allocate(pz(zsize(1),zsize(2),zsize(3)))
    select case(ipencil_axis)
    case(1)
      px(:,:,:) = p(1:n(1),1:n(2),1:n(3))
    case(2)
      py(:,:,:) = p(1:n(1),1:n(2),1:n(3))
      call transpose_y_to_x(py,px)
    case(3)
      pz(:,:,:) = p(1:n(1),1:n(2),1:n(3))
      !call transpose_z_to_x(pz,px)
      call transpose_z_to_y(pz,py)
      call transpose_y_to_x(py,px)
    end select
    !
    pfft => px
    if(is_mixed(1)) call prep_dctviii('F',bc(0,1)//bc(1,1),1,px,pfft)
    call fft(arrplan(1,1),pfft) ! fwd transform in x
    if(is_mixed(1)) call posp_dctviii('F',bc(0,1)//bc(1,1),1,pfft,px)
    !
    call transpose_x_to_y(px,py)
    pfft => py
    if(is_mixed(2)) call prep_dctviii('F',bc(0,2)//bc(1,2),2,py,pfft)
    call fft(arrplan(1,2),pfft) ! fwd transform in y
    if(is_mixed(2)) call posp_dctviii('F',bc(0,2)//bc(1,2),2,pfft,py)
    !
    q = merge(1,0,(c_or_f(3) == 'f').and.(bc(1,3) /= 'P').and.(hi_z(3) == ng(3)))
    is_periodic_z = bc(0,3)//bc(1,3) == 'PP'
    if(.not.is_poisson_dtdma) then
      call transpose_y_to_z(py,pz)
      !
      call gaussel(n_z(1),n_z(2),n_z(3)-q,0,a,b,c,is_periodic_z,norm,pz,lambdaxy)
      !
      call transpose_z_to_y(pz,py)
    else
      call gaussel_dtdma(n_z(1),n_z(2),n_z(3)-q,0,a,b,c,is_periodic_z,norm,py,lambdaxy,is_dtdma_update_,aa_z,cc_z)
      if(present(is_dtdma_update)) is_dtdma_update = is_dtdma_update_
    end if
    pfft => py
    if(is_mixed(2)) call prep_dctviii('B',bc(0,2)//bc(1,2),2,py,pfft)
    call fft(arrplan(2,2),pfft) ! bwd transform in y
    if(is_mixed(2)) call posp_dctviii('B',bc(0,2)//bc(1,2),2,pfft,py)
    if((c_or_f(2) == 'f').and.(bc(0,2)//bc(1,2) == 'NN')) then
      py(:,ng(2),:) = py(:,ng(2)-1,:)
    end if
    !
    call transpose_y_to_x(py,px)
    pfft => px
    if(is_mixed(1)) call prep_dctviii('B',bc(0,1)//bc(1,1),1,px,pfft)
    call fft(arrplan(2,1),pfft) ! bwd transform in x
    if(is_mixed(1)) call posp_dctviii('B',bc(0,1)//bc(1,1),1,pfft,px)
    if((c_or_f(1) == 'f').and.(bc(0,1)//bc(1,1) == 'NN')) then
      px(ng(1),:,:) = px(ng(1)-1,:,:)
    end if
    !
    select case(ipencil_axis)
    case(1)
      p(1:n(1),1:n(2),1:n(3)) = px(:,:,:)
    case(2)
      call transpose_x_to_y(px,py)
      p(1:n(1),1:n(2),1:n(3)) = py(:,:,:)
    case(3)
      !call transpose_x_to_z(px,pz)
      call transpose_x_to_y(px,py)
      call transpose_y_to_z(py,pz)
      p(1:n(1),1:n(2),1:n(3)) = pz(:,:,:)
    end select
  end subroutine solver
  !
  subroutine gaussel(nx,ny,n,nh,a,b,c,is_periodic,norm,p,lambdaxy)
    implicit none
    integer , intent(in) :: nx,ny,n,nh
    real(rp), intent(in), dimension(:) :: a,b,c
    logical , intent(in) :: is_periodic
    real(rp), intent(in) :: norm
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    real(rp), intent(in), dimension(nx,ny), optional :: lambdaxy
    real(rp), allocatable, dimension(:,:,:) :: d,p2
    real(rp) :: den,pivot_tol,z
    integer :: i,j,k,nn
    !
    ! a single periodic point has no Z Laplacian
    !
    if(is_periodic.and.(n == 1)) then
      if(present(lambdaxy)) then
        do j=1,ny
          do i=1,nx
            den = b(1) + lambdaxy(i,j)
            pivot_tol = epsilon(den)*max(abs(b(1)),abs(lambdaxy(i,j)))
            if(abs(den) <= pivot_tol) then
              p(i,j,1) = 0.
            else
              p(i,j,1) = p(i,j,1)*norm/den
            end if
          end do
        end do
      else
        z = norm/b(1)
        do j=1,ny
          do i=1,nx
            p(i,j,1) = p(i,j,1)*z
          end do
        end do
      end if
      return
    end if
    !
    ! solve tridiagonal system
    !
    nn = n
    if(is_periodic) nn = n-1
    !
    ! allocate work arrays
    !
    allocate(d(nx,ny,nn))
    !
    ! forward elimination
    !
    if(present(lambdaxy)) then
      do j=1,ny
        do i=1,nx
          z = 1._rp/(b(1) + lambdaxy(i,j))
          d(i,j,1) = c(1)*z
          p(i,j,1) = p(i,j,1)*norm*z
        end do
      end do
      !
      do k=2,nn
        do j=1,ny
          do i=1,nx
            den = b(k) + lambdaxy(i,j) - a(k)*d(i,j,k-1)
            !
            ! pin the constant pressure mode instead of regularizing its
            ! singular final equation
            !
            pivot_tol = epsilon(den)*max(abs(b(k)+lambdaxy(i,j)),abs(a(k)*d(i,j,k-1)))
            if((k == nn).and.(abs(den) <= pivot_tol)) then
              d(i,j,k) = 0._rp
              p(i,j,k) = 0._rp
            else
              z = 1._rp/den
              d(i,j,k) = c(k)*z
              p(i,j,k) = (p(i,j,k)*norm - a(k)*p(i,j,k-1))*z
            end if
          end do
        end do
      end do
    else
      z = 1._rp/b(1)
      do j=1,ny
        do i=1,nx
          d(i,j,1) = c(1)*z
          p(i,j,1) = p(i,j,1)*norm*z
        end do
      end do
      !
      do k=2,nn
        z = 1._rp/(b(k) - a(k)*d(1,1,k-1))
        do j=1,ny
          do i=1,nx
            d(i,j,k) = c(k)*z
            p(i,j,k) = (p(i,j,k)*norm - a(k)*p(i,j,k-1))*z
          end do
        end do
      end do
    end if
    !
    ! backward substitution
    !
    do k=nn-1,1,-1
      do j=1,ny
        do i=1,nx
          p(i,j,k) = p(i,j,k) - d(i,j,k)*p(i,j,k+1)
        end do
      end do
    end do
    !
    ! handle periodic closure with an auxiliary tridiagonal solve
    !
    if(is_periodic) then
      allocate(p2(nx,ny,nn))
      !
      ! initialize the auxiliary right-hand side for the periodic correction
      !
      do j=1,ny
        do i=1,nx
          p2(i,j,1:nn) = 0.
          p2(i,j,1)  = -a(1)
          p2(i,j,nn) = p2(i,j,nn) - c(nn)
        end do
      end do
      !
      ! forward elimination for the auxiliary system
      !
      if(present(lambdaxy)) then
        do j=1,ny
          do i=1,nx
            z = 1._rp/(b(1) + lambdaxy(i,j))
            d(i,j,1) = c(1)*z
            p2(i,j,1) = p2(i,j,1)*z
          end do
        end do
        do k=2,nn
          do j=1,ny
            do i=1,nx
              z = 1._rp/(b(k) + lambdaxy(i,j) - a(k)*d(i,j,k-1))
              d(i,j,k) = c(k)*z
              p2(i,j,k) = (p2(i,j,k) - a(k)*p2(i,j,k-1))*z
            end do
          end do
        end do
      else
        z = 1._rp/b(1)
        do j=1,ny
          do i=1,nx
            d(i,j,1) = c(1)*z
            p2(i,j,1) = p2(i,j,1)*z
          end do
        end do
        do k=2,nn
          z = 1._rp/(b(k) - a(k)*d(1,1,k-1))
          do j=1,ny
            do i=1,nx
              d(i,j,k) = c(k)*z
              p2(i,j,k) = (p2(i,j,k) - a(k)*p2(i,j,k-1))*z
            end do
          end do
        end do
      end if
      !
      ! backward substitution for the auxiliary system
      !
      do k=nn-1,1,-1
        do j=1,ny
          do i=1,nx
            p2(i,j,k) = p2(i,j,k) - d(i,j,k)*p2(i,j,k+1)
          end do
        end do
      end do
      !
      ! solve for the periodic closure value and correct the interior solution
      !
      if(present(lambdaxy)) then
        do j=1,ny
          do i=1,nx
            den = b(nn+1) + lambdaxy(i,j) + c(nn+1)*p2(i,j,1) + a(nn+1)*p2(i,j,nn)
            pivot_tol = epsilon(den)*max(abs(b(nn+1)+lambdaxy(i,j)), &
                                         abs(c(nn+1)*p2(i,j,1)+a(nn+1)*p2(i,j,nn)))
            if(abs(den) <= pivot_tol) then
              p(i,j,nn+1) = 0._rp
            else
              p(i,j,nn+1) = (p(i,j,nn+1)*norm - c(nn+1)*p(i,j,1) - a(nn+1)*p(i,j,nn))/den
            end if
          end do
        end do
      else
        do j=1,ny
          do i=1,nx
            p(i,j,nn+1) = (p(i,j,nn+1)*norm - c(nn+1)*p( i,j,1) - a(nn+1)*p( i,j,nn)) / &
                          (b(nn+1)          + c(nn+1)*p2(i,j,1) + a(nn+1)*p2(i,j,nn))
          end do
        end do
      end if
      !
      ! apply the (Sherman-Morrison) periodic correction to all interior points
      !
      do k=1,nn
        do j=1,ny
          do i=1,nx
            p(i,j,k) = p(i,j,k) + p2(i,j,k)*p(i,j,nn+1)
          end do
        end do
      end do
    end if
  end subroutine gaussel
  !
  subroutine gaussel_dtdma(nx,ny,n,nh,a,b,c,is_periodic,norm,p,lambdaxy,is_update,aa_z_save,cc_z_save)
    !
    ! distributed TDMA solver
    !
    use mod_common_mpi, only: dinfo_dtdma
    !
    implicit none
    integer , intent(in) :: nx,ny,n,nh
    real(rp), intent(in), dimension(:) :: a,b,c
    logical , intent(in) :: is_periodic
    real(rp), intent(in) :: norm
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    real(rp), intent(in), dimension(:,:), optional :: lambdaxy
    logical , intent(inout), optional :: is_update
    real(rp), intent(inout), dimension(:,:,:), optional :: aa_z_save,cc_z_save
    real(rp),              dimension(nx,ny,n) :: aa,cc
    real(rp), allocatable, dimension(: ,: ,:) :: aa_y,cc_y,pp_y,aa_z,cc_z,pp_z
    real(rp), allocatable, dimension(: ,: ,:) :: pp_z_2,cc_z_0
    real(rp) :: z,zz(2),bb(n),den,pivot_tol
    integer :: i,j,k
    integer , dimension(3) :: nr_z
    integer :: nx_r,ny_r,nn
    logical :: is_present_lambdaxy
    !
    is_present_lambdaxy = present(lambdaxy)
    nr_z(:) = dinfo_dtdma%zsz(:)
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
    !
    if(present(lambdaxy)) then
      !
      ! factor inner rows of z-distributed systems so that they are only coupled to the boundaries:
      !
      do j=1,ny
        do i=1,nx
          !
          bb(:) = b(1:n) + lambdaxy(i,j)
          zz(:) = 1._rp/bb(1:2)
          aa(i,j,1:2) = a(1:2)*zz(:)
          cc(i,j,1:2) = c(1:2)*zz(:)
          p( i,j,1:2) = p(i,j,1:2)*norm*zz(:)
          !
          ! elimination of lower diagonals
          !
          do k=3,n
            z = 1._rp/(bb(k) - a(k)*cc(i,j,k-1))
            p(i,j,k) = (p(i,j,k)*norm-a(k)*p(i,j,k-1))*z
            aa(i,j,k) = -a(k)*aa(i,j,k-1)*z
            cc(i,j,k) = c(k)*z
          end do
          !
          ! elimination of upper diagonals
          !
          do k=n-2,2,-1
            p(i,j,k)  = p(i,j,k) - cc(i,j,k)*p(i,j,k+1)
            aa(i,j,k) =  aa(i,j,k)-cc(i,j,k)*aa(i,j,k+1)
            cc(i,j,k) = -cc(i,j,k)*cc(i,j,k+1)
          end do
          !
          ! with two rows both points already belong to the reduced system
          !
          if(n > 2) then
            z = 1._rp/(1._rp - aa(i,j,2)*cc(i,j,1))
            p(i,j,1) = (p(i,j,1)-cc(i,j,1)*p(i,j,2))*z
            aa(i,j,1) = aa(i,j,1)*z
            cc(i,j,1) = -cc(i,j,1)*cc(i,j,2)*z
          end if
          !
          ! gather reduced systems
          !
          aa_y(i,j,1) = aa(i,j,1); aa_y(i,j,2) = aa(i,j,n)
          cc_y(i,j,1) = cc(i,j,1); cc_y(i,j,2) = cc(i,j,n)
          pp_y(i,j,1) = p(i,j,1) ; pp_y(i,j,2) = p(i,j,n)
        end do
      end do
    else
      do j=1,ny
        do i=1,nx
          zz(:) = 1._rp/b(1:2)
          aa(i,j,1:2) = a(1:2)*zz(:)
          cc(i,j,1:2) = c(1:2)*zz(:)
          p( i,j,1:2) = p(i,j,1:2)*norm*zz(:)
          !
          ! elimination of lower diagonals
          !
          do k=3,n
            z = 1._rp/(b(k) - a(k)*cc(i,j,k-1))
            p(i,j,k) = (p(i,j,k)*norm-a(k)*p(i,j,k-1))*z
            aa(i,j,k) = -a(k)*aa(i,j,k-1)*z
            cc(i,j,k) = c(k)*z
          end do
          !
          ! elimination of upper diagonals
          !
          do k=n-2,2,-1
            p(i,j,k)  = p(i,j,k) - cc(i,j,k)*p(i,j,k+1)
            aa(i,j,k) =  aa(i,j,k)-cc(i,j,k)*aa(i,j,k+1)
            cc(i,j,k) = -cc(i,j,k)*cc(i,j,k+1)
          end do
          !
          ! with two rows both points already belong to the reduced system
          !
          if(n > 2) then
            z = 1._rp/(1._rp - aa(i,j,2)*cc(i,j,1))
            p(i,j,1) = (p(i,j,1)-cc(i,j,1)*p(i,j,2))*z
            aa(i,j,1) = aa(i,j,1)*z
            cc(i,j,1) = -cc(i,j,1)*cc(i,j,2)*z
          end if
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
    if(present(is_update) .and. present(aa_z_save) .and. present(cc_z_save)) then
      if(is_update) then
        is_update = .false.
        call transpose_y_to_z(aa_y,aa_z_save,dinfo_dtdma)
        call transpose_y_to_z(cc_y,cc_z_save,dinfo_dtdma)
      end if
      aa_z(:,:,:) = aa_z_save(:,:,:)
      cc_z(:,:,:) = cc_z_save(:,:,:)
    else
      call transpose_y_to_z(aa_y,aa_z,dinfo_dtdma)
      call transpose_y_to_z(cc_y,cc_z,dinfo_dtdma)
    end if
    call transpose_y_to_z(pp_y,pp_z,dinfo_dtdma)
    !
    ! solve reduced systems
    !
    nn   = nr_z(3)
    ny_r = nr_z(2)
    nx_r = nr_z(1)
    if(is_periodic) then
      nn = nn-1
      cc_z_0(:,:,:) = cc_z(:,:,:)
    end if
    do j=1,ny_r
      do i=1,nx_r
        do k=2,nn
          den = 1._rp - aa_z(i,j,k)*cc_z(i,j,k-1)
          pivot_tol = epsilon(den)*max(1._rp,abs(aa_z(i,j,k)*cc_z(i,j,k-1)))
          if(is_present_lambdaxy.and.(k == nn).and.(abs(den) <= pivot_tol)) then ! pin the constant pressure mode
            pp_z(i,j,k) = 0._rp
            cc_z(i,j,k) = 0._rp
          else
            z = 1._rp/den
            pp_z(i,j,k) = (pp_z(i,j,k)-aa_z(i,j,k)*pp_z(i,j,k-1))*z
            cc_z(i,j,k) = cc_z(i,j,k)*z
          end if
        end do
        do k=nn-1,1,-1
          pp_z(i,j,k) = pp_z(i,j,k) - cc_z(i,j,k)*pp_z(i,j,k+1)
        end do
      end do
    end do
    if(is_periodic) then
      associate(cc_z => cc_z_0)
      do j=1,ny_r
        do i=1,nx_r
          pp_z_2(i,j,1:nn) = 0.
          pp_z_2(i,j,1 ) = -aa_z(i,j,1 )
          pp_z_2(i,j,nn) = pp_z_2(i,j,nn) - cc_z(i,j,nn)
          !
          do k=2,nn
            z = 1._rp/(1._rp - aa_z(i,j,k)*cc_z(i,j,k-1))
            pp_z_2(i,j,k) = (pp_z_2(i,j,k)-aa_z(i,j,k)*pp_z_2(i,j,k-1))*z
            cc_z(i,j,k) = cc_z(i,j,k)*z
          end do
          !
          do k=nn-1,1,-1
            pp_z_2(i,j,k) = pp_z_2(i,j,k) - cc_z(i,j,k)*pp_z_2(i,j,k+1)
          end do
          den = 1._rp + cc_z(i,j,nn+1)*pp_z_2(i,j,1) + aa_z(i,j,nn+1)*pp_z_2(i,j,nn)
          pivot_tol = epsilon(den)*max(1._rp,abs(cc_z(i,j,nn+1)*pp_z_2(i,j,1)+aa_z(i,j,nn+1)*pp_z_2(i,j,nn)))
          if(is_present_lambdaxy.and.(abs(den) <= pivot_tol)) then
            pp_z(i,j,nn+1) = 0._rp
          else
            pp_z(i,j,nn+1) = (pp_z(i,j,nn+1)-cc_z(i,j,nn+1)*pp_z(i,j,1)-aa_z(i,j,nn+1)*pp_z(i,j,nn))/den
          end if
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
    call transpose_z_to_y(pp_z,pp_y,dinfo_dtdma)
    !
    ! obtain final solution on the inner points
    !
    do j=1,ny
      do i=1,nx
        p(i,j,1) = pp_y(i,j,1)
        p(i,j,n) = pp_y(i,j,2)
        do k=2,n-1
          p(i,j,k) = p(i,j,k) - aa(i,j,k)*p(i,j,1) - cc(i,j,k)*p(i,j,n)
        end do
      end do
    end do
  end subroutine gaussel_dtdma
  !
  subroutine dgtsv_homebrewed(n,a,b,c,norm,p)
    implicit none
    integer , intent(in) :: n
    real(rp), intent(in   ), dimension(:) :: a,b,c
    real(rp), intent(in   )               :: norm
    real(rp), intent(inout), dimension(:) :: p
    real(rp), dimension(n) :: d
    real(rp) :: z
    integer :: l
    !
    ! Gauss elimination
    !
    z = 1._rp/b(1)
    d(1) = c(1)*z
    p(1) = p(1)*norm*z
    do l=2,n
      z = 1._rp/(b(l) - a(l)*d(l-1))
      d(l) = c(l)*z
      p(l) = (p(l)*norm-a(l)*p(l-1))*z
    end do
    !
    ! backward substitution
    !
    do l=n-1,1,-1
      p(l) = p(l) - d(l)*p(l+1)
    end do
  end subroutine dgtsv_homebrewed
  !
  subroutine solver_gaussel_z(n,ng,hi,a,b,c,bcz,c_or_f,norm,p)
    implicit none
    integer , intent(in), dimension(3) :: n,ng,hi
    real(rp), intent(in), dimension(:) :: a,b,c
    character(len=1), dimension(0:1), intent(in) :: bcz
    character(len=1), intent(in), dimension(3) :: c_or_f
    real(rp), intent(in) :: norm
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    real(rp), allocatable, dimension(:,:,:) :: px,py,pz
    integer :: q
    logical :: is_periodic_z
    integer, dimension(3) :: n_z,hi_z
    logical :: is_no_decomp_z
    !
    n_z(:)  = zsize(:)
    hi_z(:) = zend(:)
    if(is_poisson_dtdma) then
      n_z(:)  = ysize(:)
      hi_z(:) = yend(:)
    end if
    is_no_decomp_z = n(3) == ng(3)
    if(.not.is_no_decomp_z) then
      allocate(py(ysize(1),ysize(2),ysize(3)))
      if(.not.is_poisson_dtdma) allocate(pz(zsize(1),zsize(2),zsize(3)))
      select case(ipencil_axis)
      case(1)
        allocate(px(xsize(1),xsize(2),xsize(3)))
        px(:,:,:) = p(1:n(1),1:n(2),1:n(3))
        call transpose_x_to_y(px,py)
      case(2)
        py(:,:,:) = p(1:n(1),1:n(2),1:n(3))
      end select
      if(.not.is_poisson_dtdma) call transpose_y_to_z(py,pz)
    end if
    !
    q = merge(1,0,(c_or_f(3) == 'f').and.(bcz(1) /= 'P').and.(hi_z(3) == ng(3)))
    is_periodic_z = bcz(0)//bcz(1) == 'PP'
    if(is_no_decomp_z) then
      call gaussel(n(1),n(2),n(3)-q,1,a,b,c,is_periodic_z,norm,p)
    else if(is_poisson_dtdma) then
      !
      ! the reduced-system descriptor expects Y pencils without halos
      !
      call gaussel_dtdma(n_z(1),n_z(2),n_z(3)-q,0,a,b,c,is_periodic_z,norm,py)
    else
      call gaussel(n_z(1),n_z(2),n_z(3)-q,0,a,b,c,is_periodic_z,norm,pz)
    end if
    !
    if(.not.is_no_decomp_z) then
      if(.not.is_poisson_dtdma) call transpose_z_to_y(pz,py)
      select case(ipencil_axis)
      case(1)
        call transpose_y_to_x(py,px)
        p(1:n(1),1:n(2),1:n(3)) = px(:,:,:)
      case(2)
        p(1:n(1),1:n(2),1:n(3)) = py(:,:,:)
      end select
    end if
  end subroutine solver_gaussel_z
  !
#if 0
  subroutine gaussel_lapack(nx,ny,n,a,b,c,p)
    implicit none
#if !defined(_SINGLE_PRECISION)
    external :: dgttrf,dgttrs
    procedure(), pointer :: gttrf => dgttrf, gttrs => dgttrs
#else
    external :: sgttrf,sgttrs
    procedure(), pointer :: gttrf => sgttrf, gttrs => sgttrs
#endif
    integer , intent(in) :: nx,ny,n
    real(rp), intent(in), dimension(:) :: a,b,c
    real(rp), intent(inout), dimension(:,:,:) :: p
    real(rp), allocatable, dimension(:) :: aa,bb,cc,ccc
    integer , allocatable, dimension(:) :: ipiv
    integer :: i,j,info
    !real(rp), dimension(n,nx,ny) :: p_t
    !
    allocate(aa,source=a(2:n  ))
    allocate(bb,source=b(1:n  ))
    allocate(cc,source=c(1:n-1))
    allocate(ccc(n-2),ipiv(n))
    call gttrf(n,aa,bb,cc,ccc,ipiv,info)
    do j=1,ny
      do i=1,nx
        call gttrs('N',n,1,aa,bb,cc,ccc,ipiv,p(i,j,1:n),n,info)
      end do
    end do
    !p_t = reshape(p(1:nx,1:ny,1:n),shape(p_t),order=[2,3,1])
    !call gttrs('N',n,nx*ny,aa,bb,cc,ccc,ipiv,p_t(1:n,:,:),n,info)
    !p(1:nx,1:ny,1:n) = reshape(p_t,shape(p(1:nx,1:ny,1:n)),order=[3,1,2])
  end subroutine gaussel_lapack
#endif
end module mod_solver
