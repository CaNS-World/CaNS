module mod_bound
  use mpi
  use mod_common_mpi, only: ierr,comm_cart,halo
  use mod_types
  implicit none
  private
  public boundp,bounduvw,updt_rhs_b
  contains
  subroutine bounduvw(cbc,n,bc,is_correc,dl,dzc,dzf,u,v,w)
    !
    ! imposes velocity boundary conditions
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3,3) :: cbc
    integer , intent(in), dimension(3) :: n 
    real(rp), intent(in), dimension(0:1,3,3) :: bc
    logical , intent(in)                   :: is_correc
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(0:) :: dzc,dzf
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    logical :: impose_norm_bc
    !
    do idir = 1,3
      call updthalo(1,halo(idir),nb(:,idir),idir,u)
      call updthalo(1,halo(idir),nb(:,idir),idir,v)
      call updthalo(1,halo(idir),nb(:,idir),idir,w)
    enddo
    !
    impose_norm_bc = (.not.is_correc).or.(cbc(0,1,1)//cbc(1,1,1).eq.'PP')
    if(is_bound(0,1)) then
      if(impose_norm_bc) call set_bc(cbc(0,1,1),0,n(1),1,.false.,bc(0,1,1),dl(1),u)
                         call set_bc(cbc(0,1,2),0,n(1),1,.true. ,bc(0,1,2),dl(1),v)
                         call set_bc(cbc(0,1,3),0,n(1),1,.true. ,bc(0,1,3),dl(1),w)
    endif
    if(is_bound(1,1)) then
      if(impose_norm_bc) call set_bc(cbc(1,1,1),1,n(1),1,.false.,bc(1,1,1),dl(1),u)
                         call set_bc(cbc(1,1,2),1,n(1),1,.true. ,bc(1,1,2),dl(1),v)
                         call set_bc(cbc(1,1,3),1,n(1),1,.true. ,bc(1,1,3),dl(1),w)
    endif
    impose_norm_bc = (.not.is_correc).or.(cbc(0,2,2)//cbc(1,2,2).eq.'PP')
    if(is_bound(0,2)) then
                         call set_bc(cbc(0,2,1),0,n(2),2,.true. ,bc(0,2,1),dl(2),u)
      if(impose_norm_bc) call set_bc(cbc(0,2,2),0,n(2),2,.false.,bc(0,2,2),dl(2),v)
                         call set_bc(cbc(0,2,3),0,n(2),2,.true. ,bc(0,2,3),dl(2),w)
     endif
    if(is_bound(1,2)) then
                         call set_bc(cbc(1,2,1),1,n(2),2,.true. ,bc(1,2,1),dl(2),u)
      if(impose_norm_bc) call set_bc(cbc(1,2,2),1,n(2),2,.false.,bc(1,2,2),dl(2),v)
                         call set_bc(cbc(1,2,3),1,n(2),2,.true. ,bc(1,2,3),dl(2),w)
    endif
    impose_norm_bc = (.not.is_correc).or.(cbc(0,3,3)//cbc(1,3,3).eq.'PP')
    if(is_bound(0,3)) then
                         call set_bc(cbc(0,3,1),0,n(3),3,.true. ,bc(0,3,1),dzc(0)   ,u)
                         call set_bc(cbc(0,3,2),0,n(3),3,.true. ,bc(0,3,2),dzc(0)   ,v)
      if(impose_norm_bc) call set_bc(cbc(0,3,3),0,n(3),3,.false.,bc(0,3,3),dzf(0)   ,w)
    endif
      if(is_bound(1,3)) then
                         call set_bc(cbc(1,3,1),1,n(3),3,.true. ,bc(1,3,1),dzc(n(3)),u)
                         call set_bc(cbc(1,3,2),1,n(3),3,.true. ,bc(1,3,2),dzc(n(3)),v)
      if(impose_norm_bc) call set_bc(cbc(1,3,3),1,n(3),3,.false.,bc(1,3,3),dzf(n(3)),w)
    endif
  end subroutine bounduvw
  !
  subroutine boundp(cbc,n,bc,dl,dzc,dzf,p)
    !
    ! imposes pressure boundary conditions
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    integer , intent(in), dimension(3) :: n 
    real(rp)         , intent(in), dimension(0:1,3) :: bc
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(0:) :: dzc,dzf
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    !
    call updthalo([n(1),n(2)],1,p)
    call updthalo([n(1),n(2)],2,p)
    !
    if(left .eq.MPI_PROC_NULL) then
      call set_bc(cbc(0,1),0,n(1),1,.true.,bc(0,1),dl(1),p)
    endif
    if(right.eq.MPI_PROC_NULL) then
      call set_bc(cbc(1,1),1,n(1),1,.true.,bc(1,1),dl(1),p)
    endif
    if(front.eq.MPI_PROC_NULL) then
      call set_bc(cbc(0,2),0,n(2),2,.true.,bc(0,2),dl(2),p)
     endif
    if(back .eq.MPI_PROC_NULL) then
      call set_bc(cbc(1,2),1,n(2),2,.true.,bc(1,2),dl(2),p)
    endif
    call set_bc(cbc(0,3),0,n(3),3,.true.,bc(0,3),dzc(0)   ,p)
    call set_bc(cbc(1,3),1,n(3),3,.true.,bc(1,3),dzc(n(3)),p)
  end subroutine boundp
  !
  subroutine set_bc(ctype,ibound,n,idir,centered,rvalue,dr,p)
    implicit none
    character(len=1), intent(in) :: ctype
    integer , intent(in) :: ibound,n,idir
    logical , intent(in) :: centered
    real(rp), intent(in) :: rvalue,dr
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    real(rp) :: factor,sgn
    !
    factor = rvalue
    if(ctype.eq.'D'.and.centered) then
      factor = 2.*factor
      sgn    = -1.
    endif
    if(ctype.eq.'N') then
      if(    ibound.eq.0) then
        factor = -dr*factor
      elseif(ibound.eq.1) then
        factor =  dr*factor
      endif
      sgn    = 1.
    endif
    !
    select case(ctype)
    case('P')
      select case(idir)
      case(1)
        !p(0  ,:,:) = p(n,:,:)
        !p(n+1,:,:) = p(1,:,:)
      case(2)
        !p(:,0  ,:) = p(:,n,:)
        !p(:,n+1,:) = p(:,1,:)
      case(3)
        !$OMP WORKSHARE
        p(:,:,0  ) = p(:,:,n)
        p(:,:,n+1) = p(:,:,1)
        !$OMP END WORKSHARE
      end select
    case('D','N')
      if(centered) then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(0  ,:,:) = factor+sgn*p(1,:,:)
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(n+1,:,:) = factor+sgn*p(n,:,:)
            !$OMP END WORKSHARE
          endif
        case(2)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(:,0  ,:) = factor+sgn*p(:,1,:)
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(:,n+1,:) = factor+sgn*p(:,n,:)
            !$OMP END WORKSHARE
          endif
        case(3)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(:,:,0  ) = factor+sgn*p(:,:,1)
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(:,:,n+1) = factor+sgn*p(:,:,n)
            !$OMP END WORKSHARE
          endif
        end select
      elseif(.not.centered.and.ctype.eq.'D') then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(0,:,:) = factor 
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(n  ,:,:) = factor
            p(n+1,:,:) = p(n-1,:,:)
            !$OMP END WORKSHARE
          endif
        case(2)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(:,0,:) = factor 
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(:,n  ,:) = factor
            p(:,n+1,:) = p(:,n-1,:)
            !$OMP END WORKSHARE
          endif
        case(3)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(:,:,0) = factor 
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(:,:,n  ) = factor
            p(:,:,n+1) = p(:,:,n-1)
            !$OMP END WORKSHARE
          endif
        end select
      elseif(.not.centered.and.ctype.eq.'N') then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            !p(0,:,:) = 1./3.*(-2.*factor+4.*p(1  ,:,:)-p(2  ,:,:))
            p(0,:,:) = 1.*factor + p(1  ,:,:)
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            !p(n,:,:) = 1./3.*(-2.*factor+4.*p(n-1,:,:)-p(n-2,:,:))
            p(n,:,:) = 1.*factor + p(n-1,:,:)
            p(n+1,:,:) = p(n,:,:) ! not needed
            !$OMP END WORKSHARE
          endif
        case(2)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            !p(:,0  ,:) = 1./3.*(-2.*factor+4.*p(:,1,:)-p(:,2  ,:))
            p(:,0,:) = 1.*factor + p(:,1  ,:)
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            !p(:,n,:) = 1./3.*(-2.*factor+4.*p(:,n-1,:)-p(:,n-2,:))
            p(:,n,:) = 1.*factor + p(:,n-1,:)
            p(:,n+1,:) = p(:,n,:) ! not needed
            !$OMP END WORKSHARE
          endif
        case(3)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            !p(:,:,0) = 1./3.*(-2.*factor+4.*p(:,:,1  )-p(:,:,2  ))
            p(:,:,0) = 1.*factor + p(:,:,1  )
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            !p(:,:,n) = 1./3.*(-2.*factor+4.*p(:,:,n-1)-p(:,:,n-2))
            p(:,:,n) = 1.*factor + p(:,:,n-1)
            p(:,:,n+1) = p(:,:,n) ! not needed
            !$OMP END WORKSHARE
          endif
        end select
      endif
    end select
  end subroutine set_bc
  !
  subroutine inflow(n,idir,dl,dzf,vel2d,u,v,w)
    implicit none
    integer, intent(in), dimension(3) :: n
    integer, intent(in) :: idir
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(0:) :: dzf
    real(rp), dimension(0:,0:), intent(in) :: vel2d
    real(rp), dimension(0:,0:,0:), intent(inout) :: u,v,w
    real(rp) :: dx,dy,dxi,dyi
    real(rp), dimension(0:n(3)+1) :: dzfi
    integer :: i,j,k
    !
    dx   = dl(1)
    dxi  = dl(1)**(-1)
    dy   = dl(2)
    dyi  = dl(2)**(-1)
    dzfi = dzf**(-1)
    select case(idir)
      case(1) ! x direction
        if(left.eq.MPI_PROC_NULL) then
          i = 0
          do k=1,n(3)
            do j=1,n(2)
              u(i,j,k) = vel2d(j,k)
            enddo
          enddo 
        endif
      case(2) ! y direction
        j = 0
        if(front.eq.MPI_PROC_NULL) then
          do k=1,n(3)
            do i=1,n(1)
              v(i,j,k) = vel2d(i,k)
            enddo
          enddo 
        endif
      case(3) ! z direction
        k = 0
        do j=1,n(2)
          do i=1,n(1)
            w(i,j,k) = vel2d(i,j)
          enddo
        enddo 
    end select
  end subroutine inflow
  !
  subroutine updt_rhs_b(c_or_f,cbc,n,rhsbx,rhsby,rhsbz,p)
    implicit none
    character, intent(in), dimension(3) :: c_or_f
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(:,:,0:) :: rhsbx,rhsby,rhsbz
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    integer , dimension(3) :: q
    integer :: idir
    q(:) = 0
    do idir = 1,3
      if(c_or_f(idir).eq.'f'.and.cbc(1,idir).eq.'D') q(idir) = 1
    enddo
    if(left.eq.MPI_PROC_NULL) then
      !$OMP WORKSHARE
      p(1        ,1:n(2),1:n(3)) = p(1        ,1:n(2),1:n(3)) + rhsbx(:,:,0)
      !$OMP END WORKSHARE
    endif  
    if(right.eq.MPI_PROC_NULL) then
      !$OMP WORKSHARE
      p(n(1)-q(1),1:n(2),1:n(3)) = p(n(1)-q(1),1:n(2),1:n(3)) + rhsbx(:,:,1)
      !$OMP END WORKSHARE
    endif
    if(front.eq.MPI_PROC_NULL) then
      !$OMP WORKSHARE
      p(1:n(1),1        ,1:n(3)) = p(1:n(1),1        ,1:n(3)) + rhsby(:,:,0)
      !$OMP END WORKSHARE
    endif
    if(back.eq.MPI_PROC_NULL) then
      !$OMP WORKSHARE
      p(1:n(1),n(2)-q(2),1:n(3)) = p(1:n(1),n(2)-q(2),1:n(3)) + rhsby(:,:,1)
      !$OMP END WORKSHARE
    endif
    !$OMP WORKSHARE
    p(1:n(1),1:n(2),1        ) = p(1:n(1),1:n(2),1        ) + rhsbz(:,:,0)
    p(1:n(1),1:n(2),n(3)-q(3)) = p(1:n(1),1:n(2),n(3)-q(3)) + rhsbz(:,:,1)
    !$OMP END WORKSHARE
  end subroutine updt_rhs_b
  !
  subroutine updthalo(nh,halo,nb,idir,p)
    implicit none
    integer , dimension(3), intent(in) :: lo,hi
    integer , intent(in) :: nh ! number of ghost points
    integer , intent(in) :: halo
    integer , intent(in), dimension(0:1) :: nb
    integer , intent(in) :: idir
    real(rp), dimension(1-nh:,1-nh:,1-nh:), intent(inout) :: p
    integer , dimension(3) :: lo,hi
    !type(MPI_REQUEST) :: requests(4)
    !
    !  this subroutine updates the halo that store info
    !  from the neighboring computational sub-domain
    !
    lo(:) = lbound(p)+nh
    hi(:) = ubound(p)-nh
    select case(idir)
    case(1) ! x direction
      call MPI_SENDRECV(p(lo(1)     ,lo(2)-nh,lo(3)-nh),1,halo,nb(0),0, &
                        p(hi(1)+1   ,lo(2)-nh,lo(3)-nh),1,halo,nb(1),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE)
      call MPI_SENDRECV(p(hi(1)-nh+1,lo(2)-nh,lo(3)-nh),1,halo,nb(1),0, &
                        p(lo(1)-nh  ,lo(2)-nh,lo(3)-nh),1,halo,nb(0),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE)
      !call MPI_IRECV( p(hi(1)+1   ,lo(2)-nh,lo(3)-nh),1,halo,nb(1),0, &
      !                MPI_COMM_WORLD,requests(1))
      !call MPI_IRECV( p(lo(1)-nh  ,lo(2)-nh,lo(3)-nh),1,halo,nb(0),1, &
      !                MPI_COMM_WORLD,requests(2))
      !call MPI_ISSEND(p(lo(1)     ,lo(2)-nh,lo(3)-nh),1,halo,nb(0),0, &
      !                MPI_COMM_WORLD,requests(3))
      !call MPI_ISSEND(p(hi(1)-nh+1,lo(2)-nh,lo(3)-nh),1,halo,nb(1),1, &
      !                MPI_COMM_WORLD,requests(4))
      !call MPI_WAITALL(4,requests,MPI_STATUSES_IGNORE)
    case(2) ! y direction
      call MPI_SENDRECV(p(lo(1)-nh,lo(2)     ,lo(3)-nh),1,halo,nb(0),0, &
                        p(lo(1)-nh,hi(2)+1   ,lo(3)-nh),1,halo,nb(1),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE)
      call MPI_SENDRECV(p(lo(1)-nh,hi(2)-nh+1,lo(3)-nh),1,halo,nb(1),0, &
                        p(lo(1)-nh,lo(2)-nh  ,lo(3)-nh),1,halo,nb(0),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE)
      !call MPI_IRECV( p(lo(1)-nh,hi(2)+1   ,lo(3)-nh),1,halo,nb(1),0, &
      !                MPI_COMM_WORLD,requests(1))
      !call MPI_IRECV( p(lo(1)-nh,lo(2)-nh  ,lo(3)-nh),1,halo,nb(0),1, &
      !                MPI_COMM_WORLD,requests(2))
      !call MPI_ISSEND(p(lo(1)-nh,lo(2)     ,lo(3)-nh),1,halo,nb(0),0, &
      !                MPI_COMM_WORLD,requests(3))
      !call MPI_ISSEND(p(lo(1)-nh,hi(2)-nh+1,lo(3)-nh),1,halo,nb(1),1, &
      !                MPI_COMM_WORLD,requests(4))
      !call MPI_WAITALL(4,requests,MPI_STATUSES_IGNORE)
    case(3) ! z direction
      call MPI_SENDRECV(p(lo(1)-nh,lo(2)-nh,lo(3)     ),1,halo,nb(0),0, &
                        p(lo(1)-nh,lo(2)-nh,hi(3)+1   ),1,halo,nb(1),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE)
end module mod_bound
