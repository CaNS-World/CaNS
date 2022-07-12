module mod_scal
  use mod_types
  implicit none
  private
  public scal
  contains
  subroutine scal(nx,ny,nz,dxi,dyi,dzi,dzci,dzfi,visc,u,v,w,s,dsdt)
    !
    !
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,dzi,visc
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in) :: u,v,w,s
    real(rp), dimension(:,:,:), intent(out) :: dsdt
    integer :: i,j,k
    real(rp) :: usip,usim,vsjp,vsjm,wskp,wskm
    real(rp) :: dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(usip,usim,vsjp,vsjm,wskp,wskm) &
    !$OMP PRIVATE(dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzi,visc,u,v,w,s,dsdt,dzci,dzfi)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          usim  = 0.5*( s(i-1,j,k)+s(i,j,k) )*u(i-1,j,k)
          usip  = 0.5*( s(i+1,j,k)+s(i,j,k) )*u(i  ,j,k)
          vsjm  = 0.5*( s(i,j-1,k)+s(i,j,k) )*v(i,j-1,k)
          vsjp  = 0.5*( s(i,j+1,k)+s(i,j,k) )*v(i,j  ,k)
          wskm  = 0.5*( s(i,j,k-1)+s(i,j,k) )*w(i,j,k-1)
          wskp  = 0.5*( s(i,j,k+1)+s(i,j,k) )*w(i,j,k  )
          dsdxp = (s(i+1,j,k)-s(i  ,j,k))*dxi
          dsdxm = (s(i  ,j,k)-s(i-1,j,k))*dxi
          dsdyp = (s(i,j+1,k)-s(i,j  ,k))*dyi
          dsdym = (s(i,j  ,k)-s(i,j-1,k))*dyi
          dsdzp = (s(i,j,k+1)-s(i,j,k  ))*dzci(k  )
          dsdzm = (s(i,j,k  )-s(i,j,k-1))*dzci(k-1)
          !
          dsdt(i,j,k) = dxi*(     -usip + usim ) + (dsdxp-dsdxm)*visc*dxi + &
                        dyi*(     -vsjp + vsjm ) + (dsdyp-dsdym)*visc*dyi + &
                        dzfi(k)*( -wskp + wskm ) + (dsdzp-dsdzm)*visc*dzfi(k)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine scal
end module mod_scal
