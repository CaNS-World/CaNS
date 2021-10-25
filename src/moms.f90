module mod_moms
  use mod_types
  implicit none
  private
  public momsad
  contains
  subroutine momsad(nx,ny,nz,dxi,dyi,dzi,dzci,dzfi,visc,u,v,w,s,dsdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,dzi,visc
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in) :: u,v,w,s
    real(rp), dimension(:,:,:), intent(out) :: dsdt
    integer :: im,ip,jm,jp,km,kp,i,j,k
    real(rp) :: usip,usim,vsjp,vsjm,wskp,wskm
    real(rp) :: dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,im,jm,km,ip,jp,kp) &
    !$OMP PRIVATE(usip,usim,vsjp,vsjm,wskp,wskm) &
    !$OMP PRIVATE(dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzi,visc,u,v,w,s,dsdt,dzci,dzfi)
    do k=1,nz
      kp = k + 1
      km = k - 1
      do j=1,ny
        jp = j + 1
        jm = j - 1
        do i=1,nx
          ip = i + 1
          im = i - 1
          usim  = 0.5*( s(im,j,k)+s(i,j,k) )*u(im,j,k)
          usip  = 0.5*( s(ip,j,k)+s(i,j,k) )*u(i ,j,k)
          vsjm  = 0.5*( s(i,jm,k)+s(i,j,k) )*v(i,jm,k)
          vsjp  = 0.5*( s(i,jp,k)+s(i,j,k) )*v(i,j ,k)
          wskm  = 0.5*( s(i,j,km)+s(i,j,k) )*w(i,j,km)
          wskp  = 0.5*( s(i,j,kp)+s(i,j,k) )*w(i,j,k )
          dsdxp = (s(ip,j,k)-s(i ,j,k))*dxi
          dsdxm = (s(i ,j,k)-s(im,j,k))*dxi
          dsdyp = (s(i,jp,k)-s(i,j ,k))*dyi
          dsdym = (s(i,j ,k)-s(i,jm,k))*dyi
          dsdzp = (s(i,j,kp)-s(i,j,k ))*dzci(k)
          dsdzm = (s(i,j,k )-s(i,j,km))*dzci(km)
          !
          dsdt(i,j,k) = dxi*(     -usip + usim ) + (dsdxp-dsdxm)*visc*dxi + &
                        dyi*(     -vsjp + vsjm ) + (dsdyp-dsdym)*visc*dyi + &
                        dzfi(k)*( -wskp + wskm ) + (dsdzp-dsdzm)*visc*dzfi(k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine momsad
end module mod_moms
