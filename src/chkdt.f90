module mod_chkdt
  use mpi
  use mod_common_mpi, only:ierr
  implicit none
  private
  public chkdt
  contains
  subroutine chkdt(n,dl,dzci,dzfi,visc,u,v,w,dtmax)
    implicit none
    integer, intent(in), dimension(3) :: n
    real(8), intent(in), dimension(3) :: dl
    real(8), intent(in), dimension(0:) :: dzci,dzfi
    real(8), intent(in) :: visc
    real(8), intent(in), dimension(0:,0:,0:) :: u,v,w
    real(8), intent(out) :: dtmax
    real(8) :: dxi,dyi,dzi
    real(8) :: ux,uy,uz,vx,vy,vz,wx,wy,wz
    real(8) :: dtix,dtiy,dtiz,dti,dlmin
    integer :: i,j,k
    !
    dti = 0.
    dxi = 1./dl(1)
    dyi = 1./dl(2)
    dzi = 1./dl(3)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP  SHARED(n,u,v,w,dxi,dyi,dzi,dzci,dzfi) &
    !$OMP  PRIVATE(i,j,k,ux,uy,uz,vx,vy,vz,wx,wy,wz,dtix,dtiy,dtiz) &
    !$OMP  REDUCTION(max:dti)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          ux = abs(u(i,j,k))
          vx = 0.25d0*abs( v(i,j,k)+v(i,j-1,k)+v(i+1,j,k)+v(i+1,j-1,k) )
          wx = 0.25d0*abs( w(i,j,k)+w(i,j,k-1)+w(i+1,j,k)+w(i+1,j,k-1) )
          dtix = ux*dxi+vx*dyi+wx*dzfi(k)
          uy = 0.25d0*abs( u(i,j,k)+u(i,j+1,k)+u(i-1,j+1,k)+u(i-1,j,k) )
          vy = abs(v(i,j,k))
          wy = 0.25d0*abs( w(i,j,k)+w(i,j+1,k)+w(i,j+1,k-1)+w(i,j,k-1) )
          dtiy = uy*dxi+vy*dyi+wy*dzfi(k)
          uz = 0.25d0*abs( u(i,j,k)+u(i-1,j,k)+u(i-1,j,k+1)+u(i,j,k+1) )
          vz = 0.25d0*abs( v(i,j,k)+v(i,j-1,k)+v(i,j-1,k+1)+v(i,j,k+1) )
          wz = abs(w(i,j,k))
          dtiz = uz*dxi+vz*dyi+wz*dzci(k)
          dti = max(dti,dtix,dtiy,dtiz)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    call mpi_allreduce(MPI_IN_PLACE,dti,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(dti.eq.0.d0) dti = 1.
    dlmin     = minval(dl)
#ifdef IMPDIFF
    dtmax = sqrt(3.d0)/dti
#else
    dtmax = min(1.65d0/12.d0/visc*dlmin**2,sqrt(3.d0)/dti)
#endif
    return
  end subroutine chkdt
end module mod_chkdt
