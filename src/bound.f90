module mod_bound
  use mpi
  use mod_common_mpi, only: ierr,status,comm_cart,left,right,front,back,xhalo,yhalo
  implicit none
  private
  public boundp,bounduvw
  contains
  subroutine bounduvw(cbc,n,bc,ioutflowdir,dl,dzc,dzf,u,v,w)
    implicit none
    character(len=1), intent(in), dimension(2,3,3) :: cbc
    integer, intent(in), dimension(3) :: n 
    real(8)         , intent(in), dimension(2,3,3) :: bc
    integer, intent(in) :: ioutflowdir
    real(8), intent(in), dimension(3) :: dl
    real(8), intent(in), dimension(0:) :: dzc,dzf
    real(8), intent(inout), dimension(0:,0:,0:) :: u,v,w
    !
    call updthalo((/n(1),n(2)/),1,u)
    call updthalo((/n(1),n(2)/),2,u)
    call updthalo((/n(1),n(2)/),1,v)
    call updthalo((/n(1),n(2)/),2,v)
    call updthalo((/n(1),n(2)/),1,w)
    call updthalo((/n(1),n(2)/),2,w)
    !
    if(left .eq.MPI_PROC_NULL) then
      call set_bc(cbc(1,1,1),0,n(1),1,.true. ,bc(1,1,1),dl(1),u)
      call set_bc(cbc(1,1,2),0,n(1),1,.false.,bc(1,1,2),dl(1),v)
      call set_bc(cbc(1,1,3),0,n(1),1,.false.,bc(1,1,3),dl(1),w)
    endif
    if(right.eq.MPI_PROC_NULL) then
      call set_bc(cbc(2,1,1),1,n(1),1,.true. ,bc(2,1,1),dl(1),u)
      call set_bc(cbc(2,1,2),1,n(1),1,.false.,bc(2,1,2),dl(1),v)
      call set_bc(cbc(2,1,3),1,n(1),1,.false.,bc(2,1,3),dl(1),w)
    endif
    if(front.eq.MPI_PROC_NULL) then
      call set_bc(cbc(1,2,1),0,n(2),2,.false.,bc(1,2,1),dl(2),u)
      call set_bc(cbc(1,2,2),0,n(2),2,.true. ,bc(1,2,2),dl(2),v)
      call set_bc(cbc(1,2,3),0,n(2),2,.false.,bc(1,2,3),dl(2),w)
     endif
    if(back .eq.MPI_PROC_NULL) then
      call set_bc(cbc(2,2,1),1,n(2),2,.false.,bc(2,2,1),dl(2),u)
      call set_bc(cbc(2,2,2),1,n(2),2,.true. ,bc(2,2,2),dl(2),v)
      call set_bc(cbc(2,2,3),1,n(2),2,.false.,bc(2,2,3),dl(2),w)
    endif
    call set_bc(cbc(1,3,1),0,n(3),3,.false.,bc(1,3,1),dzc(0)   ,u)
    call set_bc(cbc(1,3,2),0,n(3),3,.false.,bc(1,3,2),dzc(0)   ,v)
    call set_bc(cbc(1,3,3),0,n(3),3,.true. ,bc(1,3,3),dzf(0)   ,w)
    call set_bc(cbc(2,3,1),1,n(3),3,.false.,bc(2,3,1),dzc(n(3)),u) ! check
    call set_bc(cbc(2,3,2),1,n(3),3,.false.,bc(2,3,2),dzc(n(3)),v) ! check
    call set_bc(cbc(2,3,3),1,n(3),3,.true. ,bc(2,3,3),dzf(n(3)),w) ! check
    call outflow(n,ioutflowdir,dl,dzf,u,v,w)
    return
  end subroutine bounduvw
  !
  subroutine boundp(n,cbc,p)
    implicit none
    integer, intent(in), dimension(3) :: n
    character(len=1), intent(in), dimension(2,3) :: cbc
    real(8), intent(inout), dimension(0:,0:,0:) :: p
    !
    call updthalo((/n(1),n(2)/),1,p)
    call updthalo((/n(1),n(2)/),2,p)
    if(left .eq.MPI_PROC_NULL) then
      select case(cbc(1,1))
        !case('P')
        !  p(0     ,:,:) =  p(n(1),:,:)
        case('N')
          p(0     ,:,:) =  p(1   ,:,:)
        case('D')
          p(0     ,:,:) = -p(1   ,:,:)
      end select
    endif
    if(right.eq.MPI_PROC_NULL) then
      select case(cbc(2,1))
        !case('PP')
        !  p(n(1)+1,:,:) =  p(1   ,:,:)
        case('N')
          p(n(1)+1,:,:) =  p(n(1),:,:)
        case('D')
          p(n(1)+1,:,:) = -p(n(1),:,:)
      end select
    endif
    if(front.eq.MPI_PROC_NULL) then
      select case(cbc(1,2))
        !case('P')
        !   p(:,0     ,:) =  p(:,n(2),:)
        case('N')
          p(:,0     ,:) =  p(:,1   ,:)
        case('D')
          p(:,0     ,:) = -p(:,1   ,:)
      end select
    endif
    if(back .eq.MPI_PROC_NULL) then
      select case(cbc(2,2))
        !case('P')
        !  p(:,n(2)+1,:) =  p(:,1   ,:)
        case('N')
          p(:,n(2)+1,:) =  p(:,n(2),:)
        case('D')
          p(:,n(2)+1,:) = -p(:,n(2),:)
      end select
    endif
    select case(cbc(1,3))
      case('P')
        p(:,:,0) = p(:,:,n(3))
      case('N')
        p(:,:,0) =  p(:,:,1)
      case('D')
        p(:,:,0) = -p(:,:,1)
    end select
    select case(cbc(2,3))
      case('P')
        p(:,:,n(3)+1) = p(:,:,1   )
      case('N')
        p(:,:,n(3)+1) =  p(:,:,n(3))
      case('D')
        p(:,:,n(3)+1) = -p(:,:,n(3))
    end select
    return
  end subroutine boundp
  !
  subroutine set_bc(ctype,ibound,n,idir,iface,rvalue,dr,p)
    implicit none
    character(len=1), intent(in) :: ctype
    integer, intent(in) :: ibound,n,idir
    logical, intent(in) :: iface
    real(8), intent(in) :: rvalue,dr
    real(8), intent(inout), dimension(0:,0:,0:) :: p
    real(8) :: factor,sgn
    !
    factor = rvalue
    if(ctype.eq.'D'.and..not.iface) then
      factor = 2.*factor
      sgn    = -1.
    endif
    if(ctype.eq.'N'.and..not.iface) then
      factor = dr*factor
      sgn    = 1.
    endif
    !
    select case(ctype)
      case('P')
        select case(idir)
          case(1)
            p(0  ,:,:) = p(n,:,:)
            p(n+1,:,:) = p(1,:,:)
          case(2)
            p(:,0  ,:) = p(:,n,:)
            p(:,n+1,:) = p(:,1,:)
          case(3)
            p(:,:,0  ) = p(:,:,n)
            p(:,:,n+1) = p(:,:,1)
        end select
      ! end case('P')
      case('D','N')
        if(iface.and.ctype.eq.'D') then
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
        elseif(iface.and.ctype.eq.'N') then
          select case(idir)
            case(1)
              if    (ibound.eq.0) then
                p(0,  :,:) = factor + p(2  ,:,:) 
              elseif(ibound.eq.1) then
                p(n+1,:,:) = factor + p(n-1,:,:)
              endif
            case(2)
              if    (ibound.eq.0) then
                p(:,0  ,:) = factor + p(:,2  ,:) 
              elseif(ibound.eq.1) then
                p(:,n+1,:) = factor + p(:,n-1,:)
              endif
            case(3)
              if    (ibound.eq.0) then
                p(:,:,0  ) = factor + p(:,:,2  )
              elseif(ibound.eq.1) then
                p(:,:,n+1) = factor + p(:,:,n-1)
              endif
          end select
        elseif(.not.iface) then
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
        endif
      ! end case
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
    select case(idir)
      case(1) ! x direction
        if(right.eq.MPI_PROC_NULL) then
          i = n(1) + 1
          do k=1,n(3)
            do j=1,n(2)
              u(i,j,k) = u(i-1,j,k) -dx*((v(i,j,k)-v(i,j-1,k))*dyi+(w(i,j,k)-w(i,j,k-1))*dzfi(k))
            enddo
          enddo 
        endif
      case(2) ! y direction
        j = n(2) + 1
        if(back.eq.MPI_PROC_NULL) then
          do k=1,n(3)
            do i=1,n(1)
              v(i,j,k) = v(i,j-1,k) -dy*((u(i,j,k)-u(i-1,j,k))*dxi+(w(i,j,k)-w(i,j,k-1))*dzfi(k))
            enddo
          enddo 
        endif
      case(3) ! z direction
        k = n(3) + 1
        do j=1,n(2)
          do i=1,n(1)
            w(i,j,k) = w(i,j,k-1) -dzf(k)*((u(i,j,k)-v(i-1,j,k))*dxi+(v(i,j,k)-v(i,j-1,k))*dyi)
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
