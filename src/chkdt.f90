module mod_chkdt
  use mpi
  use mod_common_mpi, only:ierr
  use mod_types
  implicit none
  private
  public chkdt
  contains
  subroutine chkdt(n,dl,dzci,dzfi,visc,u,v,w,dtmax)
    !
    ! computes maximum allowed timestep
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), intent(in) :: visc
    real(rp), intent(in), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out) :: dtmax
    real(rp) :: dxi,dyi,dzi
    real(rp) :: ux,uy,uz,vx,vy,vz,wx,wy,wz
    real(rp) :: dtix,dtiy,dtiz,dti,dlmin
    integer :: i,j,k,im,ip,jm,jp,km,kp
    !
    dti = 0.
    dxi = 1./dl(1)
    dyi = 1./dl(2)
    dzi = 1./dl(3)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,u,v,w,dxi,dyi,dzi,dzci,dzfi) &
    !$OMP PRIVATE(i,j,k,im,jm,km,ip,jp,kp) &
    !$OMP PRIVATE(ux,uy,uz,vx,vy,vz,wx,wy,wz,dtix,dtiy,dtiz) &
    !$OMP REDUCTION(max:dti)
    do k=1,n(3)
      kp = k + 1
      km = k - 1
      do j=1,n(2)
        jp = j + 1
        jm = j - 1
        do i=1,n(1)
          ip = i + 1
          im = i - 1
          ux = abs(u(i,j,k))
          vx = 0.25*abs( v(i,j,k)+v(i,jm,k)+v(ip,j,k)+v(ip,jm,k) )
          wx = 0.25*abs( w(i,j,k)+w(i,j,km)+w(ip,j,k)+w(ip,j,km) )
          dtix = ux*dxi+vx*dyi+wx*dzfi(k)
          uy = 0.25*abs( u(i,j,k)+u(i,jp,k)+u(im,jp,k)+u(im,j,k) )
          vy = abs(v(i,j,k))
          wy = 0.25*abs( w(i,j,k)+w(i,jp,k)+w(i,jp,km)+w(i,j,km) )
          dtiy = uy*dxi+vy*dyi+wy*dzfi(k)
          uz = 0.25*abs( u(i,j,k)+u(im,j,k)+u(im,j,kp)+u(i,j,kp) )
          vz = 0.25*abs( v(i,j,k)+v(i,jm,k)+v(i,jm,kp)+v(i,j,kp) )
          wz = abs(w(i,j,k))
          dtiz = uz*dxi+vz*dyi+wz*dzci(k)
          dti = max(dti,dtix,dtiy,dtiz)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    call mpi_allreduce(MPI_IN_PLACE,dti,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(dti.eq.0.) dti = 1.
    dlmin     = minval(dl)
    dlmin     = min(dlmin,minval(1./dzfi)) ! minimum of dzf is an estimate on the safe side
#ifdef IMPDIFF
    dtmax = sqrt(3.)/dti
#else
    dtmax = min(1.65/12./visc*dlmin**2,sqrt(3.)/dti)
#endif
    return
  end subroutine chkdt
end module mod_chkdt
