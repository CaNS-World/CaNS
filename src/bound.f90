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
  public boundp,bounduvw,updt_rhs_b
  contains
  subroutine bounduvw(cbc,n,bc,nb,is_bound,dl,dzc,dzf,u,v,w,keep_norm_values,halos_norm,halos_tang)
    !
    ! imposes velocity boundary conditions
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3,3) :: cbc
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(0:1,3,3) :: bc
    integer , intent(in), dimension(0:1,3  ) :: nb
    logical , intent(in), dimension(0:1,3  ) :: is_bound
    real(rp), intent(in), dimension(3 ) :: dl
    real(rp), intent(in), dimension(0:) :: dzc,dzf
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    logical , intent(in), optional :: keep_norm_values
    integer , intent(in), optional :: halos_norm,halos_tang
    logical, dimension(0:1,3) :: apply_norm_bc
    logical, dimension(0:1) :: updt_norm
    logical :: keep_norm,updt_tang
    integer :: idir,nh
    !
    nh = 1
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
    call updthalo(nh,halo(1),nb(:,1),1,u,halos_norm)
    call updthalo(nh,halo(1),nb(:,1),1,v,halos_tang)
    call updthalo(nh,halo(1),nb(:,1),1,w,halos_tang)
    call updthalo(nh,halo(2),nb(:,2),2,u,halos_tang)
    call updthalo(nh,halo(2),nb(:,2),2,v,halos_norm)
    call updthalo(nh,halo(2),nb(:,2),2,w,halos_tang)
    call updthalo(nh,halo(3),nb(:,3),3,u,halos_tang)
    call updthalo(nh,halo(3),nb(:,3),3,v,halos_tang)
    call updthalo(nh,halo(3),nb(:,3),3,w,halos_norm)
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
      if(apply_norm_bc(0,1)) call set_bc(cbc(0,1,1),0,1,nh,.false.,bc(0,1,1),dl(1),u)
                             call set_bc(cbc(0,1,2),0,1,nh,.true. ,bc(0,1,2),dl(1),v)
                             call set_bc(cbc(0,1,3),0,1,nh,.true. ,bc(0,1,3),dl(1),w)
    end if
    if(is_bound(1,1)) then
      if(apply_norm_bc(1,1)) call set_bc(cbc(1,1,1),1,1,nh,.false.,bc(1,1,1),dl(1),u)
                             call set_bc(cbc(1,1,2),1,1,nh,.true. ,bc(1,1,2),dl(1),v)
                             call set_bc(cbc(1,1,3),1,1,nh,.true. ,bc(1,1,3),dl(1),w)
    end if
    if(is_bound(0,2)) then
                             call set_bc(cbc(0,2,1),0,2,nh,.true. ,bc(0,2,1),dl(2),u)
      if(apply_norm_bc(0,2)) call set_bc(cbc(0,2,2),0,2,nh,.false.,bc(0,2,2),dl(2),v)
                             call set_bc(cbc(0,2,3),0,2,nh,.true. ,bc(0,2,3),dl(2),w)
    end if
    if(is_bound(1,2)) then
                             call set_bc(cbc(1,2,1),1,2,nh,.true. ,bc(1,2,1),dl(2),u)
      if(apply_norm_bc(1,2)) call set_bc(cbc(1,2,2),1,2,nh,.false.,bc(1,2,2),dl(2),v)
                             call set_bc(cbc(1,2,3),1,2,nh,.true. ,bc(1,2,3),dl(2),w)
    end if
    if(is_bound(0,3)) then
                             call set_bc(cbc(0,3,1),0,3,nh,.true. ,bc(0,3,1),dzc(0)   ,u)
                             call set_bc(cbc(0,3,2),0,3,nh,.true. ,bc(0,3,2),dzc(0)   ,v)
      if(apply_norm_bc(0,3)) call set_bc(cbc(0,3,3),0,3,nh,.false.,bc(0,3,3),dzf(0)   ,w)
    end if
    if(is_bound(1,3)) then
                             call set_bc(cbc(1,3,1),1,3,nh,.true. ,bc(1,3,1),dzc(n(3)),u)
                             call set_bc(cbc(1,3,2),1,3,nh,.true. ,bc(1,3,2),dzc(n(3)),v)
      if(apply_norm_bc(1,3)) call set_bc(cbc(1,3,3),1,3,nh,.false.,bc(1,3,3),dzf(n(3)),w)
    end if
  end subroutine bounduvw
  !
  subroutine boundp(cbc,n,bc,nb,is_bound,dl,dzc,p,halo_side)
    !
    ! imposes pressure boundary conditions
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(0:1,3) :: bc
    integer , intent(in), dimension(0:1,3) :: nb
    logical , intent(in), dimension(0:1,3) :: is_bound
    real(rp), intent(in), dimension(3 ) :: dl
    real(rp), intent(in), dimension(0:) :: dzc
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    integer , intent(in), optional :: halo_side
    logical, dimension(0:1) :: update_halo
    integer :: idir,nh
    !
    nh = 1
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
      call updthalo(nh,halo(idir),nb(:,idir),idir,p,halo_side)
    end do
#else
    if(any(update_halo(:))) call updthalo_gpu(nh,cbc(0,:)//cbc(1,:)==['PP','PP','PP'],p)
#endif
    !
    if(is_bound(0,1).and.update_halo(0)) then
      call set_bc(cbc(0,1),0,1,nh,.true.,bc(0,1),dl(1),p)
    end if
    if(is_bound(1,1).and.update_halo(1)) then
      call set_bc(cbc(1,1),1,1,nh,.true.,bc(1,1),dl(1),p)
    end if
    if(is_bound(0,2).and.update_halo(0)) then
      call set_bc(cbc(0,2),0,2,nh,.true.,bc(0,2),dl(2),p)
    end if
    if(is_bound(1,2).and.update_halo(1)) then
      call set_bc(cbc(1,2),1,2,nh,.true.,bc(1,2),dl(2),p)
    end if
    if(is_bound(0,3).and.update_halo(0)) then
      call set_bc(cbc(0,3),0,3,nh,.true.,bc(0,3),dzc(0)   ,p)
    end if
    if(is_bound(1,3).and.update_halo(1)) then
      call set_bc(cbc(1,3),1,3,nh,.true.,bc(1,3),dzc(n(3)),p)
    end if
  end subroutine boundp
  !
  subroutine set_bc(ctype,ibound,idir,nh,centered,rvalue,dr,p)
    implicit none
    character(len=1), intent(in) :: ctype
    integer , intent(in) :: ibound,idir,nh
    logical , intent(in) :: centered
    real(rp), intent(in) :: rvalue,dr
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    real(rp) :: factor,sgn
    integer  :: n,dh
    integer  :: i,j,k
    !
    n = size(p,idir) - 2*nh
    factor = rvalue
    if(ctype == 'D'.and.centered) then
      factor = 2.*factor
      sgn    = -1.
    end if
    if(ctype == 'N') then
      if(     ibound == 0) then
        factor = -dr*factor ! n.b.: only valid for nh /= 1 or factor /= 0
      else if(ibound == 1) then
        factor =  dr*factor ! n.b.: only valid for nh /= 1 or factor /= 0
      end if
      sgn    = 1.
    end if
    !
    do dh=0,nh-1
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
        select case(idir)
        case(1)
          if     (ibound == 0) then
            !$acc parallel loop collapse(2) default(present) async(1)
            !$OMP parallel do   collapse(2) DEFAULT(shared)
            do k=1-nh,size(p,3)-nh
              do j=1-nh,size(p,2)-nh
                p(  0-dh,j,k) = p(n-dh,j,k)
              end do
            end do
          else if(ibound == 1) then
            !$acc parallel loop collapse(2) default(present) async(1)
            !$OMP parallel do   collapse(2) DEFAULT(shared)
            do k=1-nh,size(p,3)-nh
              do j=1-nh,size(p,2)-nh
                p(n+1+dh,j,k) = p(1+dh,j,k)
              end do
            end do
          end if
        case(2)
          if     (ibound == 0) then
            !$acc parallel loop collapse(2) default(present) async(1)
            !$OMP parallel do   collapse(2) DEFAULT(shared)
            do k=1-nh,size(p,3)-nh
              do i=1-nh,size(p,1)-nh
                p(i,  0-dh,k) = p(i,n-dh,k)
              end do
            end do
          else if(ibound == 1) then
            !$acc parallel loop collapse(2) default(present) async(1)
            !$OMP parallel do   collapse(2) DEFAULT(shared)
            do k=1-nh,size(p,3)-nh
              do i=1-nh,size(p,1)-nh
                p(i,n+1+dh,k) = p(i,1+dh,k)
              end do
            end do
          end if
        case(3)
          if     (ibound == 0) then
            !$acc parallel loop collapse(2) default(present) async(1)
            !$OMP parallel do   collapse(2) DEFAULT(shared)
            do j=1-nh,size(p,2)-nh
              do i=1-nh,size(p,1)-nh
                p(i,j,  0-dh) = p(i,j,n-dh)
              end do
            end do
          else if(ibound == 1) then
            !$acc parallel loop collapse(2) default(present) async(1)
            !$OMP parallel do   collapse(2) DEFAULT(shared)
            do j=1-nh,size(p,2)-nh
              do i=1-nh,size(p,1)-nh
                p(i,j,n+1+dh) = p(i,j,1+dh)
              end do
            end do
          end if
        end select
      case('D','N')
        if(centered) then
          select case(idir)
          case(1)
            if     (ibound == 0) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do k=1-nh,size(p,3)-nh
                do j=1-nh,size(p,2)-nh
                  p(  0-dh,j,k) = factor+sgn*p(1+dh,j,k)
                end do
              end do
            else if(ibound == 1) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do k=1-nh,size(p,3)-nh
                do j=1-nh,size(p,2)-nh
                  p(n+1+dh,j,k) = factor+sgn*p(n-dh,j,k)
                end do
              end do
            end if
          case(2)
            if     (ibound == 0) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do k=1-nh,size(p,3)-nh
                do i=1-nh,size(p,1)-nh
                  p(i,  0-dh,k) = factor+sgn*p(i,1+dh,k)
                end do
              end do
            else if(ibound == 1) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do k=1-nh,size(p,3)-nh
                do i=1-nh,size(p,1)-nh
                  p(i,n+1+dh,k) = factor+sgn*p(i,n-dh,k)
                end do
              end do
            end if
          case(3)
            if     (ibound == 0) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do j=1-nh,size(p,2)-nh
                do i=1-nh,size(p,1)-nh
                  p(i,j,  0-dh) = factor+sgn*p(i,j,1+dh)
                end do
              end do
            else if(ibound == 1) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do j=1-nh,size(p,2)-nh
                do i=1-nh,size(p,1)-nh
                  p(i,j,n+1+dh) = factor+sgn*p(i,j,n-dh)
                end do
              end do
            end if
          end select
        else if(.not.centered.and.ctype == 'D') then
          select case(idir)
          case(1)
            if     (ibound == 0) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do k=1-nh,size(p,3)-nh
                do j=1-nh,size(p,2)-nh
                  p(0-dh,j,k) = factor
                end do
              end do
            else if(ibound == 1) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do k=1-nh,size(p,3)-nh
                do j=1-nh,size(p,2)-nh
                  p(n+1,j,k) = p(n-1,j,k) ! unused
                  p(n+dh,j,k) = factor
                end do
              end do
            end if
          case(2)
            if     (ibound == 0) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do k=1-nh,size(p,3)-nh
                do i=1-nh,size(p,1)-nh
                  p(i,0-dh,k) = factor
                end do
              end do
            else if(ibound == 1) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do k=1-nh,size(p,3)-nh
                do i=1-nh,size(p,1)-nh
                  p(i,n+1,k) = p(i,n-1,k) ! unused
                  p(i,n+dh,k) = factor
                end do
              end do
            end if
          case(3)
            if     (ibound == 0) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do j=1-nh,size(p,2)-nh
                do i=1-nh,size(p,1)-nh
                  p(i,j,0-dh) = factor
                end do
              end do
            else if(ibound == 1) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do j=1-nh,size(p,2)-nh
                do i=1-nh,size(p,1)-nh
                  p(i,j,n+1) = p(i,j,n-1) ! unused
                  p(i,j,n+dh) = factor
                end do
              end do
            end if
          end select
        else if(.not.centered.and.ctype == 'N') then
          select case(idir)
          case(1)
            if     (ibound == 0) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do k=1-nh,size(p,3)-nh
                do j=1-nh,size(p,2)-nh
                  !p(0-dh,j,k) = 1./3.*(-2.*factor+4.*p(1+dh,j,k)-p(2+dh,j,k))
                  p(0-dh,j,k) = 1.*factor + p(  1+dh,j,k)
                end do
              end do
            else if(ibound == 1) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do k=1-nh,size(p,3)-nh
                do j=1-nh,size(p,2)-nh
                  !p(n+1,j,k) = 1./3.*(-2.*factor+4.*p(n-1,j,k)-p(n-2,j,k))
                  p(n+1,j,k) = p(n,j,k) ! unused
                  p(n+dh,j,k) = 1.*factor + p(n-1-dh,j,k)
                end do
              end do
            end if
          case(2)
            if     (ibound == 0) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do k=1-nh,size(p,3)-nh
                do i=1-nh,size(p,1)-nh
                  !p(i,0-dh,k) = 1./3.*(-2.*factor+4.*p(i,1+dh,k)-p(i,2+dh,k))
                  p(i,0-dh,k) = 1.*factor + p(i,  1+dh,k)
                end do
              end do
            else if(ibound == 1) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do k=1-nh,size(p,3)-nh
                do i=1-nh,size(p,1)-nh
                  !p(i,n+1,k) = 1./3.*(-2.*factor+4.*p(i,n-1,k)-p(i,n-2,k))
                  p(i,n+1,k) = p(i,n,k) ! unused
                  p(i,n+dh,k) = 1.*factor + p(i,n-1-dh,k)
                end do
              end do
            end if
          case(3)
            if     (ibound == 0) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do j=1-nh,size(p,2)-nh
                do i=1-nh,size(p,1)-nh
                  !p(i,j,0-dh) = 1./3.*(-2.*factor+4.*p(i,j,1+dh)-p(i,j,2+dh))
                  p(i,j,0-dh) = 1.*factor + p(i,j,  1+dh)
                end do
              end do
            else if(ibound == 1) then
              !$acc parallel loop collapse(2) default(present) async(1)
              !$OMP parallel do   collapse(2) DEFAULT(shared)
              do j=1-nh,size(p,2)-nh
                do i=1-nh,size(p,1)-nh
                  !p(i,j,n+1) = 1./3.*(-2.*factor+4.*p(i,j,n-1)-p(i,j,n-2))
                  p(i,j,n+1) = p(i,j,n) ! unused
                  p(i,j,n+dh) = 1.*factor + p(i,j,n-1-dh)
                end do
              end do
            end if
          end select
        end if
      end select
    end do
  end subroutine set_bc
  !
  subroutine inflow(idir,is_bound,vel2d,u,v,w)
    implicit none
    integer , intent(in   )  :: idir
    logical , intent(in   ), dimension(0:1,3) :: is_bound
    real(rp), intent(in   ), dimension(0:,0:   ) :: vel2d
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    integer :: i,j,k
    integer, dimension(3) :: n
    !
    select case(idir)
      case(1) ! x direction
        if(is_bound(0,1)) then
          n(:) = shape(u) - 2*1
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
          n(:) = shape(v) - 2*1
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
          n(:) = shape(w) - 2*1
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
  subroutine updt_rhs_b(c_or_f,cbc,n,is_bound,rhsbx,rhsby,rhsbz,p,alpha)
    implicit none
    character(len=1), intent(in), dimension(3    ) :: c_or_f
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    integer , intent(in), dimension(3) :: n
    logical , intent(in), dimension(0:1,3) :: is_bound
    real(rp), intent(in), dimension(:,:,0:), optional :: rhsbx,rhsby,rhsbz
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    real(rp), intent(in), optional :: alpha
    integer , dimension(3) :: q
    integer :: idir
    integer :: nn
    integer :: i,j,k
    real(rp) :: norm
    q(:) = 0
    do idir = 1,3
      if(c_or_f(idir) == 'f'.and.cbc(1,idir) == 'D') q(idir) = 1
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
#if !defined(_USE_DIEZDECOMP)
    use cudecomp
#else
    use diezdecomp
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
