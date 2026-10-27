! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_bound
  use mpi
  use mod_common_mpi, only: ierr,halo
  use mod_param     , only: ipencil_axis
  use mod_types
  implicit none
  private
  integer, parameter :: HALO_LOWER = -1, HALO_NONE = 0, HALO_UPPER = 1
  public boundp,bounduvw,updt_rhs_b,set_bc,bc_dist
  contains
  subroutine bounduvw(cbc,n,nh,nh_g,bc,nb,is_bound,dl,dzc,dzf,u,v,w,keep_norm_values,halos_norm,halos_tang)
    !
    ! imposes velocity boundary conditions
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3,3) :: cbc
    integer , intent(in), dimension(3) :: n
    integer , intent(in)               :: nh,nh_g
    real(rp), intent(in), dimension(0:1,3,3) :: bc
    integer , intent(in), dimension(0:1,3  ) :: nb
    logical , intent(in), dimension(0:1,3  ) :: is_bound
    real(rp), intent(in), dimension(3 ) :: dl
    real(rp), intent(in), dimension(1-nh_g:) :: dzc,dzf
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: u,v,w
    logical , intent(in), optional :: keep_norm_values
    integer , intent(in), optional :: halos_norm,halos_tang
    logical, dimension(0:1,3) :: apply_norm_bc
    logical, dimension(0:1) :: updt_norm
    logical :: keep_norm,updt_tang
    integer :: idir,ib
    real(rp), allocatable, dimension(:,:,:), save :: drc,drf
    !
    ! distances between the ghost layers and their mirror points, for centered (drc)
    ! and non-centered (drf) components; these depend only on the grid, so they are
    ! computed only once, for the widest halo
    !
    if(.not.allocated(drc) .and. .not.allocated(drf)) then
      allocate(drc(0:nh_g,0:1,3),drf(0:nh_g,0:1,3))
      block
        real(rp), dimension(1-nh_g:n(1)+nh_g) :: dsx
        real(rp), dimension(1-nh_g:n(2)+nh_g) :: dsy
        dsx(:) = dl(1)
        dsy(:) = dl(2)
        do ib=0,1
          drc(:,ib,1) = bc_dist(n(1),nh_g,nh_g,ib,.true. ,dsx)
          drf(:,ib,1) = bc_dist(n(1),nh_g,nh_g,ib,.false.,dsx)
          drc(:,ib,2) = bc_dist(n(2),nh_g,nh_g,ib,.true. ,dsy)
          drf(:,ib,2) = bc_dist(n(2),nh_g,nh_g,ib,.false.,dsy)
          drc(:,ib,3) = bc_dist(n(3),nh_g,nh_g,ib,.true. ,dzc)
          drf(:,ib,3) = bc_dist(n(3),nh_g,nh_g,ib,.false.,dzf)
        end do
      end block
    end if
    keep_norm = .false.
    if(present(keep_norm_values)) keep_norm = keep_norm_values
    updt_norm(:) = .true.
    if(present(halos_norm)) then
      select case(halos_norm)
      case(HALO_LOWER)
        updt_norm(1) = .false.
      case(HALO_NONE)
        updt_norm(:) = .false.
      case(HALO_UPPER)
        updt_norm(0) = .false.
      end select
    end if
    updt_tang = .true.
    if(present(halos_tang)) then
      select case(halos_tang)
      case(HALO_LOWER,HALO_UPPER)
      case(HALO_NONE)
        updt_tang = .false.
      end select
    end if
    do idir = 1,3
      apply_norm_bc(0,idir) = .not.keep_norm .or. &
                              (cbc(0,idir,idir)//cbc(1,idir,idir) == 'PP'.and.updt_norm(0))
      apply_norm_bc(1,idir) = .not.keep_norm .or. &
                              (cbc(0,idir,idir)//cbc(1,idir,idir) == 'PP'.and.updt_norm(1))
    end do
    !
#if !defined(_OPENACC)
    call updthalo(nh,halo(1,nh),nb(:,1),1,u,halos_norm)
    call updthalo(nh,halo(1,nh),nb(:,1),1,v,halos_tang)
    call updthalo(nh,halo(1,nh),nb(:,1),1,w,halos_tang)
    call updthalo(nh,halo(2,nh),nb(:,2),2,u,halos_tang)
    call updthalo(nh,halo(2,nh),nb(:,2),2,v,halos_norm)
    call updthalo(nh,halo(2,nh),nb(:,2),2,w,halos_tang)
    call updthalo(nh,halo(3,nh),nb(:,3),3,u,halos_tang)
    call updthalo(nh,halo(3,nh),nb(:,3),3,v,halos_tang)
    call updthalo(nh,halo(3,nh),nb(:,3),3,w,halos_norm)
#else
    if(any(updt_norm(:))) then
      call updthalo_gpu(nh,cbc(0,:,1)//cbc(1,:,1)==['PP','PP','PP'],u,1)
      call updthalo_gpu(nh,cbc(0,:,2)//cbc(1,:,2)==['PP','PP','PP'],v,2)
      call updthalo_gpu(nh,cbc(0,:,3)//cbc(1,:,3)==['PP','PP','PP'],w,3)
    end if
    if(updt_tang) then
      call updthalo_gpu(nh,cbc(0,:,1)//cbc(1,:,1)==['PP','PP','PP'],u,2)
      call updthalo_gpu(nh,cbc(0,:,1)//cbc(1,:,1)==['PP','PP','PP'],u,3)
      call updthalo_gpu(nh,cbc(0,:,2)//cbc(1,:,2)==['PP','PP','PP'],v,1)
      call updthalo_gpu(nh,cbc(0,:,2)//cbc(1,:,2)==['PP','PP','PP'],v,3)
      call updthalo_gpu(nh,cbc(0,:,3)//cbc(1,:,3)==['PP','PP','PP'],w,1)
      call updthalo_gpu(nh,cbc(0,:,3)//cbc(1,:,3)==['PP','PP','PP'],w,2)
    end if
#endif
    !
    if(is_bound(0,1)) then
      if(apply_norm_bc(0,1)) call set_bc(cbc(0,1,1),0,1,nh,.false.,bc(0,1,1),drf(:,0,1),u)
                             call set_bc(cbc(0,1,2),0,1,nh,.true. ,bc(0,1,2),drc(:,0,1),v)
                             call set_bc(cbc(0,1,3),0,1,nh,.true. ,bc(0,1,3),drc(:,0,1),w)
    end if
    if(is_bound(1,1)) then
      if(apply_norm_bc(1,1)) call set_bc(cbc(1,1,1),1,1,nh,.false.,bc(1,1,1),drf(:,1,1),u)
                             call set_bc(cbc(1,1,2),1,1,nh,.true. ,bc(1,1,2),drc(:,1,1),v)
                             call set_bc(cbc(1,1,3),1,1,nh,.true. ,bc(1,1,3),drc(:,1,1),w)
    end if
    if(is_bound(0,2)) then
                             call set_bc(cbc(0,2,1),0,2,nh,.true. ,bc(0,2,1),drc(:,0,2),u)
      if(apply_norm_bc(0,2)) call set_bc(cbc(0,2,2),0,2,nh,.false.,bc(0,2,2),drf(:,0,2),v)
                             call set_bc(cbc(0,2,3),0,2,nh,.true. ,bc(0,2,3),drc(:,0,2),w)
    end if
    if(is_bound(1,2)) then
                             call set_bc(cbc(1,2,1),1,2,nh,.true. ,bc(1,2,1),drc(:,1,2),u)
      if(apply_norm_bc(1,2)) call set_bc(cbc(1,2,2),1,2,nh,.false.,bc(1,2,2),drf(:,1,2),v)
                             call set_bc(cbc(1,2,3),1,2,nh,.true. ,bc(1,2,3),drc(:,1,2),w)
    end if
    if(is_bound(0,3)) then
                             call set_bc(cbc(0,3,1),0,3,nh,.true. ,bc(0,3,1),drc(:,0,3),u)
                             call set_bc(cbc(0,3,2),0,3,nh,.true. ,bc(0,3,2),drc(:,0,3),v)
      if(apply_norm_bc(0,3)) call set_bc(cbc(0,3,3),0,3,nh,.false.,bc(0,3,3),drf(:,0,3),w)
    end if
    if(is_bound(1,3)) then
                             call set_bc(cbc(1,3,1),1,3,nh,.true. ,bc(1,3,1),drc(:,1,3),u)
                             call set_bc(cbc(1,3,2),1,3,nh,.true. ,bc(1,3,2),drc(:,1,3),v)
      if(apply_norm_bc(1,3)) call set_bc(cbc(1,3,3),1,3,nh,.false.,bc(1,3,3),drf(:,1,3),w)
    end if
    !
    ! refresh normal halos after imposing upper physical boundary values
    !
    if(.not.keep_norm .and. updt_norm(1)) then
      do idir=1,3
        if(cbc(1,idir,idir) == 'P') cycle
        select case(idir)
#if !defined(_OPENACC)
        case(1)
          call updthalo(nh,halo(1,nh),nb(:,1),1,u,HALO_UPPER)
        case(2)
          call updthalo(nh,halo(2,nh),nb(:,2),2,v,HALO_UPPER)
        case(3)
          call updthalo(nh,halo(3,nh),nb(:,3),3,w,HALO_UPPER)
#else
        case(1)
          call updthalo_gpu(nh,cbc(0,:,1)//cbc(1,:,1)==['PP','PP','PP'],u,1)
        case(2)
          call updthalo_gpu(nh,cbc(0,:,2)//cbc(1,:,2)==['PP','PP','PP'],v,2)
        case(3)
          call updthalo_gpu(nh,cbc(0,:,3)//cbc(1,:,3)==['PP','PP','PP'],w,3)
#endif
        end select
      end do
    end if
  end subroutine bounduvw
  !
  subroutine boundp(cbc,n,nh,nh_g,bc,nb,is_bound,dl,dzc,p,halo_side)
    !
    ! imposes pressure boundary conditions
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    integer , intent(in), dimension(3) :: n
    integer , intent(in)               :: nh,nh_g
    real(rp), intent(in), dimension(0:1,3) :: bc
    integer , intent(in), dimension(0:1,3) :: nb
    logical , intent(in), dimension(0:1,3) :: is_bound
    real(rp), intent(in), dimension(3 ) :: dl
    real(rp), intent(in), dimension(1-nh_g:) :: dzc
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    integer , intent(in), optional :: halo_side
    logical, dimension(0:1) :: update_halo
    integer :: idir,ib
    real(rp), allocatable, dimension(:,:,:), save :: dr
    !
    ! distances between the ghost layers and their mirror points; these depend only
    ! on the grid, so they are computed only once, for the widest halo, which serves
    ! the pressure and all scalars
    !
    if(.not.allocated(dr)) then
      allocate(dr(0:nh_g,0:1,3))
      block
        real(rp), dimension(1-nh_g:n(1)+nh_g) :: dsx
        real(rp), dimension(1-nh_g:n(2)+nh_g) :: dsy
        dsx(:) = dl(1)
        dsy(:) = dl(2)
        do ib=0,1
          dr(:,ib,1) = bc_dist(n(1),nh_g,nh_g,ib,.true.,dsx)
          dr(:,ib,2) = bc_dist(n(2),nh_g,nh_g,ib,.true.,dsy)
          dr(:,ib,3) = bc_dist(n(3),nh_g,nh_g,ib,.true.,dzc)
        end do
      end block
    end if
    update_halo(:) = .true.
    if(present(halo_side)) then
      select case(halo_side)
      case(HALO_LOWER)
        update_halo(1) = .false.
      case(HALO_NONE)
        update_halo(:) = .false.
      case(HALO_UPPER)
        update_halo(0) = .false.
      end select
    end if
    !
#if !defined(_OPENACC)
    do idir = 1,3
      call updthalo(nh,halo(idir,nh),nb(:,idir),idir,p,halo_side)
    end do
#else
    if(any(update_halo(:))) call updthalo_gpu(nh,cbc(0,:)//cbc(1,:)==['PP','PP','PP'],p)
#endif
    !
    if(is_bound(0,1).and.update_halo(0)) then
      call set_bc(cbc(0,1),0,1,nh,.true.,bc(0,1),dr(:,0,1),p)
    end if
    if(is_bound(1,1).and.update_halo(1)) then
      call set_bc(cbc(1,1),1,1,nh,.true.,bc(1,1),dr(:,1,1),p)
    end if
    if(is_bound(0,2).and.update_halo(0)) then
      call set_bc(cbc(0,2),0,2,nh,.true.,bc(0,2),dr(:,0,2),p)
    end if
    if(is_bound(1,2).and.update_halo(1)) then
      call set_bc(cbc(1,2),1,2,nh,.true.,bc(1,2),dr(:,1,2),p)
    end if
    if(is_bound(0,3).and.update_halo(0)) then
      call set_bc(cbc(0,3),0,3,nh,.true.,bc(0,3),dr(:,0,3),p)
    end if
    if(is_bound(1,3).and.update_halo(1)) then
      call set_bc(cbc(1,3),1,3,nh,.true.,bc(1,3),dr(:,1,3),p)
    end if
  end subroutine boundp
  !
  function bc_dist(n,nh,nh_g,ibound,centered,ds) result(dr)
    !
    ! distances between the `nh` ghost layers next to boundary `ibound` and their
    ! mirror points across it, from the grid spacing `ds` (cell-center spacing for a
    ! centered variable, face spacing otherwise); for a non-centered variable, `dr(0)`
    ! is the spacing of the one-sided Neumann condition at the boundary node
    !
    implicit none
    integer , intent(in) :: n,nh,nh_g,ibound
    logical , intent(in) :: centered
    real(rp), intent(in), dimension(1-nh_g:) :: ds
    real(rp), dimension(0:nh)                :: dr
    integer :: dh,m,m0,m1
    !
    dr(:) = 0.
    do dh=0,nh
      if(centered) then
        if(dh == nh) cycle ! a centered variable has `nh` ghost layers
        if(ibound == 0) then
          m0 = -dh  ; m1 = dh
        else
          m0 = n-dh ; m1 = n+dh
        end if
      else
        if(     dh == 0     ) then
          m0 = ibound*n; m1 = m0
        else if(ibound == 0) then
          m0 = 1-dh  ; m1 = dh
        else
          m0 = n+1-dh; m1 = n+dh
        end if
      end if
      do m=m0,m1
        dr(dh) = dr(dh) + ds(m)
      end do
    end do
  end function bc_dist
  !
  subroutine set_bc(ctype,ibound,idir,nh,centered,rvalue,dr,p)
    !
    ! imposes a boundary condition on the `nh` ghost layers of `p` next to boundary
    ! `ibound` of direction `idir`, and on the boundary node of a non-centered variable;
    ! the ghost values mirror the interior about the boundary:
    !
    ! - Dirichlet: p(ghost) = 2*rvalue - p(mirror); p(boundary node) = rvalue;
    ! - Neumann  : p(ghost) = p(mirror) -/+ dr*rvalue at the lower/upper boundary;
    !              p(boundary node) = p(adjacent node) -/+ dr(0)*rvalue;
    !
    ! dr(dh) -> distance between ghost layer `dh` and its mirror point (see `bc_dist`)
    !
    ! each layer takes one kernel, p(ghost) = factor + sgn*p(mirror) (sgn = 0 for a Dirichlet
    ! boundary node); the boundary node of a non-centered variable at the upper boundary is set
    ! together with the first ghost layer
    !
    implicit none
    character(len=1), intent(in) :: ctype
    integer , intent(in) :: ibound,idir,nh
    logical , intent(in) :: centered
    real(rp), intent(in) :: rvalue
    real(rp), intent(in), dimension(0:) :: dr
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    real(rp) :: factor,sgn,factor_b,sgn_b
    integer  :: n,dh,ig,im,l
    logical  :: is_with_bnode
    integer  :: i,j,k
    !
    n = size(p,idir) - 2*nh
    do dh=0,nh-1
      is_with_bnode = .false.
      !
      ! index `ig` of the ghost layer (or boundary node) and index `im` of its mirror (or source) layer
      !
      select case(ctype)
      case('P')
        !
        ! note: this periodic BC imposition assumes that the subroutine is only
        ! called for non-decomposed directions, for which `n` is the domain
        ! length in index space;
        ! note that the `is_bound(:,:)` mask above (set under `initmpi.f90`)
        ! is only true along the undecomposed pencil direction;
        ! along decomposed directions, periodicity is naturally set via the
        ! halo exchange
        !
        if(ibound == 0) then
          ig = 0-dh  ; im = n-dh
        else
          ig = n+1+dh; im = 1+dh
        end if
      case('D','N')
        l = dh ! layer of the distance `dr`
        if(centered) then
          if(ibound == 0) then
            ig = 0-dh  ; im = 1+dh
          else
            ig = n+1+dh; im = n-dh
          end if
        else
          if(ibound == 0) then
            ig = 0-dh  ; im = max(dh,1) ! `dh = 0` -> boundary node
          else
            ig = n+1+dh; im = n-1-dh; l = dh+1
            is_with_bnode = dh == 0
          end if
        end if
        if(ctype == 'D') then
          factor = 2.*rvalue
          sgn    = -1.
          if(.not.centered.and.ibound == 0.and.dh == 0) then
            factor = rvalue
            sgn    = 0.
          end if
          factor_b = rvalue
          sgn_b    = 0.
        else
          if(ibound == 0) factor = -dr(l)*rvalue
          if(ibound == 1) factor =  dr(l)*rvalue
          sgn    = 1.
          factor_b = dr(0)*rvalue
          sgn_b    = 1.
        end if
      case default
        return
      end select
      !
      if(ctype == 'P') then
        select case(idir)
        case(1)
          !$acc parallel loop collapse(2) default(present) async(1)
          !$OMP parallel do   collapse(2) DEFAULT(shared)
          do k=1-nh,size(p,3)-nh
            do j=1-nh,size(p,2)-nh
              p(ig,j,k) = p(im,j,k)
            end do
          end do
        case(2)
          !$acc parallel loop collapse(2) default(present) async(1)
          !$OMP parallel do   collapse(2) DEFAULT(shared)
          do k=1-nh,size(p,3)-nh
            do i=1-nh,size(p,1)-nh
              p(i,ig,k) = p(i,im,k)
            end do
          end do
        case(3)
          !$acc parallel loop collapse(2) default(present) async(1)
          !$OMP parallel do   collapse(2) DEFAULT(shared)
          do j=1-nh,size(p,2)-nh
            do i=1-nh,size(p,1)-nh
              p(i,j,ig) = p(i,j,im)
            end do
          end do
        end select
      else if(is_with_bnode) then
        select case(idir)
        case(1)
          !$acc parallel loop collapse(2) default(present) async(1)
          !$OMP parallel do   collapse(2) DEFAULT(shared)
          do k=1-nh,size(p,3)-nh
            do j=1-nh,size(p,2)-nh
              p(n ,j,k) = factor_b+sgn_b*p(n-1,j,k)
              p(ig,j,k) = factor  +sgn  *p(im ,j,k)
            end do
          end do
        case(2)
          !$acc parallel loop collapse(2) default(present) async(1)
          !$OMP parallel do   collapse(2) DEFAULT(shared)
          do k=1-nh,size(p,3)-nh
            do i=1-nh,size(p,1)-nh
              p(i,n ,k) = factor_b+sgn_b*p(i,n-1,k)
              p(i,ig,k) = factor  +sgn  *p(i,im ,k)
            end do
          end do
        case(3)
          !$acc parallel loop collapse(2) default(present) async(1)
          !$OMP parallel do   collapse(2) DEFAULT(shared)
          do j=1-nh,size(p,2)-nh
            do i=1-nh,size(p,1)-nh
              p(i,j,n ) = factor_b+sgn_b*p(i,j,n-1)
              p(i,j,ig) = factor  +sgn  *p(i,j,im )
            end do
          end do
        end select
      else
        select case(idir)
        case(1)
          !$acc parallel loop collapse(2) default(present) async(1)
          !$OMP parallel do   collapse(2) DEFAULT(shared)
          do k=1-nh,size(p,3)-nh
            do j=1-nh,size(p,2)-nh
              p(ig,j,k) = factor+sgn*p(im,j,k)
            end do
          end do
        case(2)
          !$acc parallel loop collapse(2) default(present) async(1)
          !$OMP parallel do   collapse(2) DEFAULT(shared)
          do k=1-nh,size(p,3)-nh
            do i=1-nh,size(p,1)-nh
              p(i,ig,k) = factor+sgn*p(i,im,k)
            end do
          end do
        case(3)
          !$acc parallel loop collapse(2) default(present) async(1)
          !$OMP parallel do   collapse(2) DEFAULT(shared)
          do j=1-nh,size(p,2)-nh
            do i=1-nh,size(p,1)-nh
              p(i,j,ig) = factor+sgn*p(i,j,im)
            end do
          end do
        end select
      end if
    end do
  end subroutine set_bc
  !
  subroutine inflow(idir,is_bound,nh,vel2d,u,v,w)
    implicit none
    integer , intent(in   )  :: idir
    logical , intent(in   ), dimension(0:1,3) :: is_bound
    integer , intent(in   )                   :: nh
    real(rp), intent(in   ), dimension(0:,0:   ) :: vel2d
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: u,v,w
    integer :: i,j,k
    integer, dimension(3) :: n
    !
    select case(idir)
      case(1) ! x direction
        if(is_bound(0,1)) then
          n(:) = shape(u) - 2*nh
          i = 0
          !$acc parallel loop collapse(2) default(present) async(1)
          !$OMP parallel do   collapse(2) DEFAULT(shared)
          do k=1,n(3)
            do j=1,n(2)
              u(i,j,k) = vel2d(j,k)
            end do
          end do
        end if
      case(2) ! y direction
        if(is_bound(0,2)) then
          n(:) = shape(v) - 2*nh
          j = 0
          !$acc parallel loop collapse(2) default(present) async(1)
          !$OMP parallel do   collapse(2) DEFAULT(shared)
          do k=1,n(3)
            do i=1,n(1)
              v(i,j,k) = vel2d(i,k)
            end do
          end do
        end if
      case(3) ! z direction
        if(is_bound(0,3)) then
          n(:) = shape(w) - 2*nh
          k = 0
          !$acc parallel loop collapse(2) default(present) async(1)
          !$OMP parallel do   collapse(2) DEFAULT(shared)
          do j=1,n(2)
            do i=1,n(1)
              w(i,j,k) = vel2d(i,j)
            end do
          end do
        end if
    end select
  end subroutine inflow
  !
  subroutine updt_rhs_b(c_or_f,cbc,n,nh,is_bound,rhsbx,rhsby,rhsbz,p,alpha)
    implicit none
    character(len=1), intent(in), dimension(3    ) :: c_or_f
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    integer , intent(in), dimension(3) :: n
    integer , intent(in)               :: nh
    logical , intent(in), dimension(0:1,3) :: is_bound
    real(rp), intent(in), dimension(:,:,0:), optional :: rhsbx,rhsby,rhsbz
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    real(rp), intent(in), optional :: alpha
    integer , dimension(3) :: q
    integer :: idir
    integer :: nn
    integer :: i,j,k
    real(rp) :: norm
    q(:) = 0
    do idir = 1,3
      if(c_or_f(idir) == 'f'.and.cbc(1,idir) /= 'P') q(idir) = 1
    end do
    norm = 1.
    if(present(alpha)) norm = alpha
    !
    if(present(rhsbx)) then
      if(is_bound(0,1)) then
        !$acc parallel loop collapse(2) default(present) async(1)
        !$OMP parallel do   collapse(2) DEFAULT(shared)
        do k=1,n(3)
          do j=1,n(2)
            p(1 ,j,k) = p(1 ,j,k) + rhsbx(j,k,0)*norm
          end do
        end do
      end if
      if(is_bound(1,1)) then
        nn = n(1)-q(1)
        !$acc parallel loop collapse(2) default(present) async(1)
        !$OMP parallel do   collapse(2) DEFAULT(shared)
        do k=1,n(3)
          do j=1,n(2)
            p(nn,j,k) = p(nn,j,k) + rhsbx(j,k,1)*norm
          end do
        end do
      end if
    end if
    if(present(rhsby)) then
      if(is_bound(0,2)) then
        !$acc parallel loop collapse(2) default(present) async(1)
        !$OMP parallel do   collapse(2) DEFAULT(shared)
        do k=1,n(3)
          do i=1,n(1)
            p(i,1 ,k) = p(i,1 ,k) + rhsby(i,k,0)*norm
          end do
        end do
      end if
      if(is_bound(1,2)) then
        nn = n(2)-q(2)
        !$acc parallel loop collapse(2) default(present) async(1)
        !$OMP parallel do   collapse(2) DEFAULT(shared)
        do k=1,n(3)
          do i=1,n(1)
            p(i,nn,k) = p(i,nn,k) + rhsby(i,k,1)*norm
          end do
        end do
      end if
    end if
    if(present(rhsbz)) then
      if(is_bound(0,3)) then
        !$acc parallel loop collapse(2) default(present) async(1)
        !$OMP parallel do   collapse(2) DEFAULT(shared)
        do j=1,n(2)
          do i=1,n(1)
            p(i,j,1 ) = p(i,j,1 ) + rhsbz(i,j,0)*norm
          end do
        end do
      end if
      if(is_bound(1,3)) then
        nn = n(3)-q(3)
        !$acc parallel loop collapse(2) default(present) async(1)
        !$OMP parallel do   collapse(2) DEFAULT(shared)
        do j=1,n(2)
          do i=1,n(1)
            p(i,j,nn) = p(i,j,nn) + rhsbz(i,j,1)*norm
          end do
        end do
      end if
    end if
  end subroutine updt_rhs_b
  !
  subroutine updthalo(nh,halo,nb,idir,p,halo_side)
    implicit none
    integer , intent(in) :: nh ! number of ghost points
    integer , intent(in) :: halo
    integer , intent(in), dimension(0:1) :: nb
    integer , intent(in) :: idir
    real(rp), dimension(1-nh:,1-nh:,1-nh:), intent(inout) :: p
    integer , intent(in), optional :: halo_side
    logical :: update_lower,update_upper
    integer , dimension(3) :: lo,hi
#if defined(_ASYNC_HALO)
    integer :: requests(4)
#endif
    !
    !  this subroutine updates the halo that store info
    !  from the neighboring computational sub-domain
    !
    update_lower = .true.
    update_upper = .true.
    if(present(halo_side)) then
      select case(halo_side)
      case(HALO_LOWER)
        update_upper = .false.
      case(HALO_NONE)
        update_lower = .false.
        update_upper = .false.
      case(HALO_UPPER)
        update_lower = .false.
      end select
    end if
    if(idir == ipencil_axis.or..not.(update_lower.or.update_upper)) return
    lo(:) = lbound(p)+nh
    hi(:) = ubound(p)-nh
    select case(idir)
    case(1) ! x direction
#if !defined(_ASYNC_HALO)
      if(update_upper) &
        call MPI_SENDRECV(p(lo(1)     ,lo(2)-nh,lo(3)-nh),1,halo,nb(0),0, &
                          p(hi(1)+1   ,lo(2)-nh,lo(3)-nh),1,halo,nb(1),0, &
                          MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      if(update_lower) &
        call MPI_SENDRECV(p(hi(1)-nh+1,lo(2)-nh,lo(3)-nh),1,halo,nb(1),0, &
                          p(lo(1)-nh  ,lo(2)-nh,lo(3)-nh),1,halo,nb(0),0, &
                          MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
#else
      requests(:) = MPI_REQUEST_NULL
      if(update_upper) then
        call MPI_IRECV( p(hi(1)+1,lo(2)-nh,lo(3)-nh),1,halo,nb(1),0, &
                        MPI_COMM_WORLD,requests(1),ierr)
        call MPI_ISEND(p(lo(1)   ,lo(2)-nh,lo(3)-nh),1,halo,nb(0),0, &
                        MPI_COMM_WORLD,requests(3),ierr)
      end if
      if(update_lower) then
        call MPI_IRECV( p(lo(1)-nh  ,lo(2)-nh,lo(3)-nh),1,halo,nb(0),1, &
                        MPI_COMM_WORLD,requests(2),ierr)
        call MPI_ISEND(p(hi(1)-nh+1,lo(2)-nh,lo(3)-nh),1,halo,nb(1),1, &
                        MPI_COMM_WORLD,requests(4),ierr)
      end if
      call MPI_WAITALL(4,requests,MPI_STATUSES_IGNORE,ierr)
#endif
    case(2) ! y direction
#if !defined(_ASYNC_HALO)
      if(update_upper) &
        call MPI_SENDRECV(p(lo(1)-nh,lo(2)     ,lo(3)-nh),1,halo,nb(0),0, &
                          p(lo(1)-nh,hi(2)+1   ,lo(3)-nh),1,halo,nb(1),0, &
                          MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      if(update_lower) &
        call MPI_SENDRECV(p(lo(1)-nh,hi(2)-nh+1,lo(3)-nh),1,halo,nb(1),0, &
                          p(lo(1)-nh,lo(2)-nh  ,lo(3)-nh),1,halo,nb(0),0, &
                          MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
#else
      requests(:) = MPI_REQUEST_NULL
      if(update_upper) then
        call MPI_IRECV(p(lo(1)-nh,hi(2)+1,lo(3)-nh),1,halo,nb(1),0, &
                       MPI_COMM_WORLD,requests(1),ierr)
        call MPI_ISEND(p(lo(1)-nh,lo(2)  ,lo(3)-nh),1,halo,nb(0),0, &
                       MPI_COMM_WORLD,requests(3),ierr)
      end if
      if(update_lower) then
        call MPI_IRECV(p(lo(1)-nh,lo(2)-nh  ,lo(3)-nh),1,halo,nb(0),1, &
                       MPI_COMM_WORLD,requests(2),ierr)
        call MPI_ISEND(p(lo(1)-nh,hi(2)-nh+1,lo(3)-nh),1,halo,nb(1),1, &
                       MPI_COMM_WORLD,requests(4),ierr)
      end if
      call MPI_WAITALL(4,requests,MPI_STATUSES_IGNORE,ierr)
#endif
    case(3) ! z direction
#if !defined(_ASYNC_HALO)
      if(update_upper) &
        call MPI_SENDRECV(p(lo(1)-nh,lo(2)-nh,lo(3)     ),1,halo,nb(0),0, &
                          p(lo(1)-nh,lo(2)-nh,hi(3)+1   ),1,halo,nb(1),0, &
                          MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      if(update_lower) &
        call MPI_SENDRECV(p(lo(1)-nh,lo(2)-nh,hi(3)-nh+1),1,halo,nb(1),0, &
                          p(lo(1)-nh,lo(2)-nh,lo(3)-nh  ),1,halo,nb(0),0, &
                          MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
#else
      requests(:) = MPI_REQUEST_NULL
      if(update_upper) then
        call MPI_IRECV(p(lo(1)-nh,lo(2)-nh,hi(3)+1),1,halo,nb(1),0, &
                       MPI_COMM_WORLD,requests(1),ierr)
        call MPI_ISEND(p(lo(1)-nh,lo(2)-nh,lo(3)  ),1,halo,nb(0),0, &
                       MPI_COMM_WORLD,requests(3),ierr)
      end if
      if(update_lower) then
        call MPI_IRECV(p(lo(1)-nh,lo(2)-nh,lo(3)-nh  ),1,halo,nb(0),1, &
                       MPI_COMM_WORLD,requests(2),ierr)
        call MPI_ISEND(p(lo(1)-nh,lo(2)-nh,hi(3)-nh+1),1,halo,nb(1),1, &
                       MPI_COMM_WORLD,requests(4),ierr)
      end if
      call MPI_WAITALL(4,requests,MPI_STATUSES_IGNORE,ierr)
#endif
    end select
  end subroutine updthalo
#if defined(_OPENACC)
  subroutine updthalo_gpu(nh,periods,p,idir_only)
    use mod_types
#if defined(_USE_DIEZDECOMP)
    use diezdecomp
#elif defined(_USE_HIP)
    use hipdecomp
#else
    use cudecomp
#endif
    use mod_common_cudecomp, only: work => work_halo, &
                                   ch => handle,gd => gd_halo, &
                                   dtype => cudecomp_real_rp, &
                                   istream => istream_acc_queue_1_comm_lib
    implicit none
    integer , intent(in) :: nh
    logical , intent(in) :: periods(3)
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    integer , intent(in), optional :: idir_only
    integer :: idir,istat
#if !defined(_USE_DIEZDECOMP)
    !$acc host_data use_device(p,work)
#endif
    do idir=1,3
      if(idir == ipencil_axis) cycle
      if(present(idir_only)) then
        if(idir /= idir_only) cycle
      end if
      select case(ipencil_axis)
      case(1)
        istat = cudecompUpdateHalosX(ch,gd,p,work,dtype,[nh,nh,nh],periods,idir,stream=istream)
      case(2)
        istat = cudecompUpdateHalosY(ch,gd,p,work,dtype,[nh,nh,nh],periods,idir,stream=istream)
      case(3)
        istat = cudecompUpdateHalosZ(ch,gd,p,work,dtype,[nh,nh,nh],periods,idir,stream=istream)
      end select
    end do
#if !defined(_USE_DIEZDECOMP)
    !$acc end host_data
#endif
  end subroutine updthalo_gpu
#endif
end module mod_bound
