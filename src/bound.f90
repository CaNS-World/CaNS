module mod_bound
  use mpi
  use mod_common_mpi, only: ierr,status,comm_cart,left,right,front,back,xhalo,yhalo
  implicit none
  private
  public boundp,bounduvw,updt_rhs_b
  contains
  subroutine bounduvw(cbc,n,bc,isoutflow,dl,dzc,dzf,u,v,w)
    !
    ! imposes velocity boundary conditions
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3,3) :: cbc
    integer, intent(in), dimension(3) :: n 
    real(8)         , intent(in), dimension(0:1,3,3) :: bc
    logical, intent(in), dimension(0:1,3) :: isoutflow
    real(8), intent(in), dimension(3) :: dl
    real(8), intent(in), dimension(0:) :: dzc,dzf
    real(8), intent(inout), dimension(0:,0:,0:) :: u,v,w
    integer :: q,idir,sgn,ioutflowdir
    !
    call updthalo((/n(1),n(2)/),1,u)
    call updthalo((/n(1),n(2)/),2,u)
    call updthalo((/n(1),n(2)/),1,v)
    call updthalo((/n(1),n(2)/),2,v)
    call updthalo((/n(1),n(2)/),1,w)
    call updthalo((/n(1),n(2)/),2,w)
    !
    if(left .eq.MPI_PROC_NULL) then
      call set_bc(cbc(0,1,1),0,n(1),1,.false.,bc(0,1,1),dl(1),u)
      call set_bc(cbc(0,1,2),0,n(1),1,.true. ,bc(0,1,2),dl(1),v)
      call set_bc(cbc(0,1,3),0,n(1),1,.true. ,bc(0,1,3),dl(1),w)
    endif
    if(right.eq.MPI_PROC_NULL) then
      call set_bc(cbc(1,1,1),1,n(1),1,.false.,bc(1,1,1),dl(1),u)
      call set_bc(cbc(1,1,2),1,n(1),1,.true. ,bc(1,1,2),dl(1),v)
      call set_bc(cbc(1,1,3),1,n(1),1,.true. ,bc(1,1,3),dl(1),w)
    endif
    if(front.eq.MPI_PROC_NULL) then
      call set_bc(cbc(0,2,1),0,n(2),2,.true. ,bc(0,2,1),dl(2),u)
      call set_bc(cbc(0,2,2),0,n(2),2,.false.,bc(0,2,2),dl(2),v)
      call set_bc(cbc(0,2,3),0,n(2),2,.true. ,bc(0,2,3),dl(2),w)
     endif
    if(back .eq.MPI_PROC_NULL) then
      call set_bc(cbc(1,2,1),1,n(2),2,.true. ,bc(1,2,1),dl(2),u)
      call set_bc(cbc(1,2,2),1,n(2),2,.false.,bc(1,2,2),dl(2),v)
      call set_bc(cbc(1,2,3),1,n(2),2,.true. ,bc(1,2,3),dl(2),w)
    endif
    call set_bc(cbc(0,3,1),0,n(3),3,.true. ,bc(0,3,1),dzc(0)   ,u)
    call set_bc(cbc(0,3,2),0,n(3),3,.true. ,bc(0,3,2),dzc(0)   ,v)
    call set_bc(cbc(0,3,3),0,n(3),3,.false.,bc(0,3,3),dzf(0)   ,w)
    call set_bc(cbc(1,3,1),1,n(3),3,.true. ,bc(1,3,1),dzc(n(3)),u)
    call set_bc(cbc(1,3,2),1,n(3),3,.true. ,bc(1,3,2),dzc(n(3)),v)
    call set_bc(cbc(1,3,3),1,n(3),3,.false.,bc(1,3,3),dzf(n(3)),w)
    !
    do q = 1,3
      do idir = 0,1
        if(isoutflow(idir,q)) then
          if(idir.eq.0) sgn = -1
          if(idir.eq.1) sgn = +1
          ioutflowdir = q*sgn
          call outflow(n,ioutflowdir,dl,dzf,u,v,w)
        endif
      enddo
    enddo
    return
  end subroutine bounduvw
  !
  subroutine boundp(cbc,n,bc,dl,dzc,dzf,p)
    !
    ! imposes pressure boundary conditions
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    integer, intent(in), dimension(3) :: n 
    real(8)         , intent(in), dimension(0:1,3) :: bc
    real(8), intent(in), dimension(3) :: dl
    real(8), intent(in), dimension(0:) :: dzc,dzf
    real(8), intent(inout), dimension(0:,0:,0:) :: p
    !
    call updthalo((/n(1),n(2)/),1,p)
    call updthalo((/n(1),n(2)/),2,p)
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
    return
  end subroutine boundp
  !
  subroutine set_bc(ctype,ibound,n,idir,stag,rvalue,dr,p)
    implicit none
    character(len=1), intent(in) :: ctype
    integer, intent(in) :: ibound,n,idir
    logical, intent(in) :: stag
    real(8), intent(in) :: rvalue,dr
    real(8), intent(inout), dimension(0:,0:,0:) :: p
    real(8) :: factor,sgn
    !
    factor = rvalue
    if(ctype.eq.'D'.and.stag) then
      factor = 2.d0*factor
      sgn    = -1.d0
    endif
    if(ctype.eq.'N'.and.stag) then
      if(    ibound.eq.0) then
        factor = -dr*factor
      elseif(ibound.eq.1) then
        factor =  dr*factor
      endif
      sgn    = 1.d0
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
        p(:,:,0  ) = p(:,:,n)
        p(:,:,n+1) = p(:,:,1)
      end select
    case('D','N')
      if(stag) then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            p(0  ,:,:) = factor+sgn*p(1,:,:)
          elseif(ibound.eq.1) then
            p(n+1,:,:) = factor+sgn*p(n,:,:)
          endif
        case(2)
          if    (ibound.eq.0) then
            p(:,0  ,:) = factor+sgn*p(:,1,:)
          elseif(ibound.eq.1) then
            p(:,n+1,:) = factor+sgn*p(:,n,:)
          endif
        case(3)
          if    (ibound.eq.0) then
            p(:,:,0  ) = factor+sgn*p(:,:,1)
          elseif(ibound.eq.1) then
            p(:,:,n+1) = factor+sgn*p(:,:,n)
          endif
        end select
      elseif(.not.stag.and.ctype.eq.'D') then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            p(0,:,:) = factor 
          elseif(ibound.eq.1) then
            p(n,:,:) = factor
          endif
        case(2)
          if    (ibound.eq.0) then
            p(:,0,:) = factor 
          elseif(ibound.eq.1) then
            p(:,n,:) = factor
          endif
        case(3)
          if    (ibound.eq.0) then
            p(:,:,0) = factor 
          elseif(ibound.eq.1) then
            p(:,:,n) = factor
          endif
        end select
      elseif(.not.stag.and.ctype.eq.'N') then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            p(0,  :,:) = 1.d0*factor + p(1  ,:,:)
          elseif(ibound.eq.1) then
            p(n+1,:,:) = 2.d0*factor + p(n-1,:,:)
          endif
        case(2)
          if    (ibound.eq.0) then
            p(:,0  ,:) = 1.d0*factor + p(:,1  ,:) 
          elseif(ibound.eq.1) then
            p(:,n+1,:) = 2.d0*factor + p(:,n-1,:)
          endif
        case(3)
          if    (ibound.eq.0) then
            p(:,:,0  ) = 1.d0*factor + p(:,:,1  )
          elseif(ibound.eq.1) then
            p(:,:,n+1) = 2.d0*factor + p(:,:,n-1)
          endif
        end select
      endif
    end select
    return
  end subroutine set_bc
  !
  subroutine outflow(n,idir,dl,dzf,u,v,w)
    implicit none
    integer, intent(in), dimension(3) :: n
    integer, intent(in) :: idir
    real(8), intent(in), dimension(3) :: dl
    real(8), intent(in), dimension(0:) :: dzf
    real(8), dimension(0:,0:,0:), intent(inout) :: u,v,w
    real(8) :: dx,dy,dxi,dyi
    real(8), dimension(0:n(3)+1) :: dzfi
    integer :: i,j,k
    !
    dx   = dl(1)     
    dxi  = dl(1)**(-1)
    dy   = dl(2)     
    dyi  = dl(2)**(-1)
    dzfi = dzf**(-1)
    !
    ! determine face velocity from zero divergence
    !
    select case(idir)
    case(1) ! x directio, right
      if(right.eq.MPI_PROC_NULL) then
        i = n(1) + 1
        do k=1,n(3)
          do j=1,n(2)
            u(i,j,k) = u(i-1,j,k) -dx*((v(i,j,k)-v(i,j-1,k))*dyi+(w(i,j,k)-w(i,j,k-1))*dzfi(k))
          enddo
        enddo 
      endif
    case(2) ! y direction, back
      j = n(2) + 1
      if(back.eq.MPI_PROC_NULL) then
        do k=1,n(3)
          do i=1,n(1)
            v(i,j,k) = v(i,j-1,k) -dy*((u(i,j,k)-u(i-1,j,k))*dxi+(w(i,j,k)-w(i,j,k-1))*dzfi(k))
          enddo
        enddo 
      endif
    case(3) ! z direction, top
      k = n(3) + 1
      do j=1,n(2)
        do i=1,n(1)
          w(i,j,k) = w(i,j,k-1) -dzf(k)*((u(i,j,k)-u(i-1,j,k))*dxi+(v(i,j,k)-v(i,j-1,k))*dyi)
        enddo
      enddo 
    case(-1) ! x direction, left
      if(left.eq.MPI_PROC_NULL) then
        i = 0
        do k=1,n(3)
          do j=1,n(2)
            u(i,j,k) = u(i+1,j,k) +dx*((v(i+1,j,k)-v(i+1,j-1,k))*dyi+(w(i+1,j,k)-w(i+1,j,k-1))*dzfi(k))
          enddo
        enddo 
      endif
    case(-2) ! y direction, front
      j = 0
      if(front.eq.MPI_PROC_NULL) then
        do k=1,n(3)
          do i=1,n(1)
            v(i,j,k) = v(i,j+1,k) +dy*((u(i,j+1,k)-u(i-1,j+1,k))*dxi+(w(i,j+1,k)-w(i,j+1,k-1))*dzfi(k))
          enddo
        enddo 
      endif
    case(-3) ! z direction, bottom
      k = 0
      do j=1,n(2)
        do i=1,n(1)
          w(i,j,k) = w(i,j,k+1) +dzf(k)*((u(i,j,k+1)-u(i-1,j,k+1))*dxi+(v(i,j,k+1)-v(i,j-1,k+1))*dyi)
        enddo
      enddo 
    end select
    return
  end subroutine outflow
  !
  subroutine inflow(n,idir,dl,dzf,vel2d,u,v,w)
    implicit none
    integer, intent(in), dimension(3) :: n
    integer, intent(in) :: idir
    real(8), intent(in), dimension(3) :: dl
    real(8), intent(in), dimension(0:) :: dzf
    real(8), dimension(0:,0:), intent(in) :: vel2d
    real(8), dimension(0:,0:,0:), intent(inout) :: u,v,w
    real(8) :: dx,dy,dxi,dyi
    real(8), dimension(0:n(3)+1) :: dzfi
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
    return
  end subroutine inflow
  !
  subroutine updt_rhs_b(c_or_f,cbc,n,rhsbx,rhsby,rhsbz,p)
    implicit none
    character, intent(in), dimension(3) :: c_or_f
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    integer, intent(in), dimension(3) :: n
    real(8), intent(in), dimension(:,:,0:) :: rhsbx,rhsby,rhsbz
    real(8), intent(inout), dimension(1:,1:,1:) :: p
    integer, dimension(3) :: q
    integer :: idir
    q(:) = 0
    do idir = 1,3
      if(c_or_f(idir).eq.'f'.and.cbc(1,idir).eq.'D') q(idir) = 1
    enddo
    if(left.eq.MPI_PROC_NULL) then
      p(1   ,:,:) = p(1   ,:,:) + rhsbx(:,:,0)
    endif  
    if(right.eq.MPI_PROC_NULL) then
      p(n(1)-q(1),:,:) = p(n(1)-q(1),:,:) + rhsbx(:,:,1)
    endif
    if(front.eq.MPI_PROC_NULL) then
      p(:,1   ,:) = p(:,1   ,:) + rhsby(:,:,0)
    endif
    if(back.eq.MPI_PROC_NULL) then
      p(:,n(2)-q(2),:) = p(:,n(2)-q(2),:) + rhsby(:,:,1)
    endif
    p(:,:,1   ) = p(:,:,1   ) + rhsbz(:,:,0)
    p(:,:,n(3)-q(3)) = p(:,:,n(3)-q(3)) + rhsbz(:,:,1)
    return
  end subroutine updt_rhs_b
  !
  subroutine updthalo(n,idir,p)
    implicit none
    integer, dimension(2), intent(in) :: n
    integer, intent(in) :: idir
    real(8), dimension(0:,0:,0:), intent(inout) :: p
    !integer :: requests(4), statuses(MPI_STATUS_SIZE,4)
    !
    !  This subroutine updates the halos that store info
    !  from the neighboring computational sub-domain
    !
    select case(idir)
    case(1) ! x direction
      call MPI_SENDRECV(p(1     ,0,0),1,xhalo,left ,0, &
                        p(n(1)+1,0,0),1,xhalo,right,0, &
                        comm_cart,status,ierr)
      call MPI_SENDRECV(p(n(1),0,0),1,xhalo,right,0, &
                        p(0   ,0,0),1,xhalo,left ,0, &
                        comm_cart,status,ierr)
         !call MPI_IRECV(p(0     ,0,0),1,xhalo,left ,1, &
         !               comm_cart,requests(2),error)
         !call MPI_IRECV(p(n(1)+1,0,0),1,xhalo,right,0, &
         !               comm_cart,requests(1),error)
         !call MPI_ISSEND(p(n(1),0,0),1,xhalo,right,1, &
         !               comm_cart,requests(4),error)
         !call MPI_ISSEND(p(1   ,0,0),1,xhalo,left ,0, &
         !               comm_cart,requests(3),error)
         !call MPI_WAITALL(4, requests, statuses, error)
    case(2) ! y direction
      call MPI_SENDRECV(p(0,1     ,0),1,yhalo,front,0, &
                        p(0,n(2)+1,0),1,yhalo,back ,0, &
                        comm_cart,status,ierr)
      call MPI_SENDRECV(p(0,n(2),0),1,yhalo,back ,0, &
                        p(0,0   ,0),1,yhalo,front,0, &
                        comm_cart,status,ierr)
         !call MPI_IRECV(p(0,n(2)+1,0),1,yhalo,back ,0, &
         !               comm_cart,requests(1),error)
         !call MPI_IRECV(p(0,0     ,0),1,yhalo,front,1, &
         !               comm_cart,requests(2),error)
         !call MPI_ISSEND(p(0,1   ,0),1,yhalo,front,0, &
         !               comm_cart,requests(3),error)
         !call MPI_ISSEND(p(0,n(2),0),1,yhalo,back ,1, &
         !               comm_cart,requests(4),error)
         !call MPI_WAITALL(4, requests, statuses, error)
    end select
    return
  end subroutine updthalo
end module mod_bound
