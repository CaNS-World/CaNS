module mod_moms
  implicit none
  private
  public momsad
  contains
  subroutine momsad(nx,ny,nz,dxi,dyi,dzi,dzci,dzfi,dzflzi,visc,u,v,w,s,dsdt)
    implicit none
    integer, intent(in) :: nx,ny,nz
    real(8), intent(in) :: dxi,dyi,dzi,visc
    real(8), intent(in), dimension(0:) :: dzci,dzfi,dzflzi
    real(8), dimension(0:,0:,0:), intent(in) :: u,v,w,s
    real(8), dimension(:,:,:), intent(out) :: dsdt
    integer :: im,ip,jm,jp,km,kp,i,j,k
    real(8) :: usip,usim,vsjp,vsjm,wskp,wskm
    real(8) :: dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm
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
          usim  = 0.5d0*( u(im,j,k)+u(i,j,k) )*s(i,j,k)
          usip  = 0.5d0*( u(ip,j,k)+u(i,j,k) )*s(i,j,k)
          vsjm  = 0.5d0*( v(i,jm,k)+v(i,j,k) )*s(i,j,k)
          vsjp  = 0.5d0*( v(i,jp,k)+v(i,j,k) )*s(i,j,k)
          wskm  = 0.5d0*( w(i,j,km)+w(i,j,k) )*s(i,j,k)
          wskp  = 0.5d0*( w(i,j,kp)+w(i,j,k) )*s(i,j,k)
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
    return
  end subroutine momsad
end module mod_moms
