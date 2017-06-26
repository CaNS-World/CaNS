module mod_mom
  use mpi
  use mod_param     , only: dims
  use mod_common_mpi, only: ierr
  implicit none
  private
  public momxad,momyad,momzad,momxp,momyp,momzp
  contains
  subroutine momxad(nx,ny,nz,dxi,dyi,dzi,dzci,dzfi,dzflzi,visc,u,v,w,dudt,taux)
    implicit none
    integer, intent(in) :: nx,ny,nz
    real(8), intent(in) :: dxi,dyi,dzi,visc
    real(8), intent(in), dimension(0:) :: dzci,dzfi,dzflzi
    real(8), dimension(0:,0:,0:), intent(in) :: u,v,w
    real(8), dimension(:,:,:), intent(out) :: dudt
    real(8), dimension(3)  , intent(out) :: taux
    integer :: im,ip,jm,jp,km,kp,i,j,k
    real(8) :: uuip,uuim,uvjp,uvjm,uwkp,uwkm
    real(8) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm
    integer :: nxg,nyg,nzg
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,im,jm,km,ip,jp,kp) &
    !$OMP PRIVATE(uuip,uuim,uvjp,uvjm,uwkp,uwkm) &
    !$OMP PRIVATE(dudxp,dudxm,dudyp,dudym,dudzp,dudzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzi,visc,u,v,w,dudt,dzci,dzfi)
    do k=1,nz
      kp = k + 1
      km = k - 1
      do j=1,ny
        jp = j + 1
        jm = j - 1
        do i=1,nx
          ip = i + 1
          im = i - 1
          uuip  = 0.25*( u(ip,j,k)+u(i,j,k) )*( u(ip,j ,k )+u(i,j ,k ) )
          uuim  = 0.25*( u(im,j,k)+u(i,j,k) )*( u(im,j ,k )+u(i,j ,k ) )
          uvjp  = 0.25*( u(i,jp,k)+u(i,j,k) )*( v(ip,j ,k )+v(i,j ,k ) )
          uvjm  = 0.25*( u(i,jm,k)+u(i,j,k) )*( v(ip,jm,k )+v(i,jm,k ) )
          uwkp  = 0.25*( u(i,j,kp)+u(i,j,k) )*( w(ip,j ,k )+w(i,j ,k ) )
          uwkm  = 0.25*( u(i,j,km)+u(i,j,k) )*( w(ip,j ,km)+w(i,j ,km) )
          dudxp = (u(ip,j,k)-u(i,j,k))*dxi
          dudxm = (u(i,j,k)-u(im,j,k))*dxi
          dudyp = (u(i,jp,k)-u(i,j,k))*dyi
          dudym = (u(i,j,k)-u(i,jm,k))*dyi
          dudzp = (u(i,j,kp)-u(i,j,k))*dzci(k)
          dudzm = (u(i,j,k)-u(i,j,km))*dzci(km)
          !
          ! Momentum balance
          !
          dudt(i,j,k) = dxi*(     -uuip + uuim ) + (dudxp-dudxm)*visc*dxi + &
                        dyi*(     -uvjp + uvjm ) + (dudyp-dudym)*visc*dyi + &
                        dzfi(k)*( -uwkp + uwkm ) + (dudzp-dudzm)*visc*dzfi(k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    taux(:) = 0.
    do k=1,nz
      do i=1,nx
        dudyp = (u(i,1 ,k)-u(i,0   ,k))*dyi*visc*dzflzi(k)
        dudym = (u(i,ny,k)-u(i,ny+1,k))*dyi*visc*dzflzi(k)
        taux(2) = taux(2) + (dudyp+dudym)
      enddo
    enddo
    do j=1,ny
      do i=1,nx
        dudzp = (u(i,j,1 )-u(i,j,0   ))*dzci(0)*visc
        dudzm = (u(i,j,nz)-u(i,j,nz+1))*dzci(nz)*visc
        taux(3) = taux(3) + (dudzp+dudzm)
      enddo
    enddo
    call mpi_allreduce(MPI_IN_PLACE,taux(1),3,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    nxg = nx*dims(1)
    nyg = ny*dims(2)
    nzg = nz
    taux(1) = taux(1)/(1.d0*nyg)
    taux(2) = taux(2)/(1.d0*nxg)
    taux(3) = taux(3)/(1.d0*nxg*nyg)
    return
  end subroutine momxad
  !
  subroutine momyad(nx,ny,nz,dxi,dyi,dzi,dzci,dzfi,dzflzi,visc,u,v,w,dvdt,tauy)
    implicit none
    integer, intent(in) :: nx,ny,nz
    real(8), intent(in) :: dxi,dyi,dzi,visc
    real(8), intent(in), dimension(0:) :: dzci,dzfi,dzflzi
    real(8), dimension(0:,0:,0:), intent(in) :: u,v,w
    real(8), dimension(:,:,:), intent(out) :: dvdt
    real(8), dimension(3), intent(out) :: tauy
    integer :: im,ip,jm,jp,km,kp,i,j,k
    real(8) :: uvip,uvim,vvjp,vvjm,wvkp,wvkm
    real(8) :: dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm
    integer :: nxg,nyg,nzg
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,im,jm,km,ip,jp,kp) &
    !$OMP PRIVATE(uvip,uvim,vvjp,vvjm,wvkp,wvkm) &
    !$OMP PRIVATE(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzi,visc,u,v,w,dvdt,dzci,dzfi)
    do k=1,nz
      kp = k + 1
      km = k - 1
      do j=1,ny
        jp = j + 1
        jm = j - 1
        do i=1,nx
          ip = i + 1
          im = i - 1
          uvip  = 0.25*( u(i ,j,k)+u(i ,jp,k) )*( v(i,j,k )+v(ip,j ,k) )
          uvim  = 0.25*( u(im,j,k)+u(im,jp,k) )*( v(i,j,k )+v(im,j ,k) )
          vvjp  = 0.25*( v(i,j,k )+v(i,jp,k)  )*( v(i,j,k )+v(i ,jp,k) )
          vvjm  = 0.25*( v(i,j,k )+v(i,jm,k)  )*( v(i,j,k )+v(i ,jm,k) )
          wvkp  = 0.25*( w(i,j,k )+w(i,jp,k)  )*( v(i,j,kp)+v(i ,j ,k) )
          wvkm  = 0.25*( w(i,j,km)+w(i,jp,km) )*( v(i,j,km)+v(i ,j ,k) )
          dvdxp = (v(ip,j,k)-v(i,j,k))*dxi
          dvdxm = (v(i,j,k)-v(im,j,k))*dxi
          dvdyp = (v(i,jp,k)-v(i,j,k))*dyi
          dvdym = (v(i,j,k)-v(i,jm,k))*dyi
          dvdzp = (v(i,j,kp)-v(i,j,k))*dzci(k)
          dvdzm = (v(i,j,k)-v(i,j,km))*dzci(km)
          !
          ! Momentum balance
          !
          dvdt(i,j,k) = dxi*(     -uvip + uvim ) + (dvdxp-dvdxm)*visc*dxi+ &
                        dyi*(     -vvjp + vvjm ) + (dvdyp-dvdym)*visc*dyi+ &
                        dzfi(k)*( -wvkp + wvkm ) + (dvdzp-dvdzm)*visc*dzfi(k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    tauy(:) = 0.
    do k=1,nz
      do j=1,ny
        dvdxp = (v(1 ,j,k)-v(0   ,j,k))*dxi*visc*dzflzi(k)
        dvdxm = (v(nx,j,k)-v(nx+1,j,k))*dxi*visc*dzflzi(k)
        tauy(1) = tauy(1) + (dvdxp+dvdxm)
      enddo
    enddo
    do j=1,ny
      do i=1,nx
        dvdzp = (v(i,j,1 )-v(i,j,0   ))*dzci(0)*visc
        dvdzm = (v(i,j,nz)-v(i,j,nz+1))*dzci(nz)*visc
        tauy(3) = tauy(3) + (dvdzp+dvdzm)
      enddo
    enddo
    call mpi_allreduce(MPI_IN_PLACE,tauy(1),3,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    nxg = nx*dims(1)
    nyg = ny*dims(2)
    nzg = nz
    tauy(1) = tauy(1)/(1.d0*nyg)
    tauy(2) = tauy(2)/(1.d0*nxg)
    tauy(3) = tauy(3)/(1.d0*nxg*nyg)
    return
  end subroutine momyad
  !
  subroutine momzad(nx,ny,nz,dxi,dyi,dzi,dzci,dzfi,dzflzi,visc,u,v,w,dwdt,tauz)
    implicit none
    integer, intent(in) :: nx,ny,nz
    real(8), intent(in) :: dxi,dyi,dzi,visc
    real(8), intent(in), dimension(0:) :: dzci,dzfi,dzflzi
    real(8), dimension(0:,0:,0:), intent(in) :: u,v,w
    real(8), dimension(:,:,:), intent(out) :: dwdt
    real(8), dimension(3), intent(out) :: tauz
    integer :: im,ip,jm,jp,km,kp,i,j,k
    real(8) :: uwip,uwim,vwjp,vwjm,wwkp,wwkm
    real(8) :: dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
    integer :: nxg,nyg,nzg
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,im,jm,km,ip,jp,kp) &
    !$OMP PRIVATE(uwip,uwim,vwjp,vwjm,wwkp,wwkm) &
    !$OMP PRIVATE(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzi,visc,u,v,w,dwdt,dzci,dzfi)
    do k=1,nz
      kp = k + 1
      km = k - 1
      do j=1,ny
        jp = j + 1
        jm = j - 1
        do i=1,nx
          ip = i + 1
          im = i - 1
          uwip  = 0.25*( w(i,j,k)+w(ip,j,k) )*( u(i ,j ,k)+u(i ,j ,kp) )
          uwim  = 0.25*( w(i,j,k)+w(im,j,k) )*( u(im,j ,k)+u(im,j ,kp) )
          vwjp  = 0.25*( w(i,j,k)+w(i,jp,k) )*( v(i ,j ,k)+v(i ,j ,kp) )
          vwjm  = 0.25*( w(i,j,k)+w(i,jm,k) )*( v(i ,jm,k)+v(i ,jm,kp) )
          wwkp  = 0.25*( w(i,j,k)+w(i,j,kp) )*( w(i ,j ,k)+w(i ,j ,kp) )
          wwkm  = 0.25*( w(i,j,k)+w(i,j,km) )*( w(i ,j ,k)+w(i ,j ,km) )
          dwdxp = (w(ip,j,k)-w(i,j,k))*dxi
          dwdxm = (w(i,j,k)-w(im,j,k))*dxi
          dwdyp = (w(i,jp,k)-w(i,j,k))*dyi
          dwdym = (w(i,j,k)-w(i,jm,k))*dyi
          dwdzp = (w(i,j,kp)-w(i,j,k))*dzfi(kp)
          dwdzm = (w(i,j,k)-w(i,j,km))*dzfi(k)
          !
          ! Momentum balance
          !
          dwdt(i,j,k) = dxi*(     -uwip + uwim ) + (dwdxp-dwdxm)*visc*dxi+ &
                        dyi*(     -vwjp + vwjm ) + (dwdyp-dwdym)*visc*dyi+ &
                        dzci(k)*( -wwkp + wwkm ) + (dwdzp-dwdzm)*visc*dzci(k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    tauz(:) = 0.
    do k=1,nz
      do j=1,ny
        dwdxp = (w(1 ,j,k)-w(0   ,j,k))*dxi*visc*dzflzi(k)
        dwdxm = (w(nx,j,k)-w(nx+1,j,k))*dxi*visc*dzflzi(k)
        tauz(1) = tauz(1) + (dwdxp+dwdxm)
      enddo
    enddo
    do k=1,nz
      do i=1,nx
        dwdyp = (w(i,1,k )-w(i,0   ,k))*dyi*visc*dzflzi(k)
        dwdym = (w(i,ny,k)-w(i,ny+1,k))*dyi*visc*dzflzi(k)
        tauz(2) = tauz(2) + (dwdyp+dwdym)
      enddo
    enddo
    call mpi_allreduce(MPI_IN_PLACE,tauz(1),3,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    nxg = nx*dims(1)
    nyg = ny*dims(2)
    nzg = nz
    tauz(1) = tauz(1)/(1.d0*nyg)
    tauz(2) = tauz(2)/(1.d0*nxg)
    tauz(3) = tauz(3)/(1.d0*nxg*nyg)
    return
  end subroutine momzad
  !
  subroutine momxp(nx,ny,nz,dxi,p,dudt)
    implicit none
    integer, intent(in) :: nx,ny,nz
    real(8), intent(in) :: dxi 
    real(8), dimension(0:,0:,0:), intent(in) :: p
    real(8), dimension(:,:,:), intent(out) :: dudt
    integer :: i,j,k
    integer :: ip
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,ip) &
    !$OMP SHARED(nx,ny,nz,dxi,p,dudt)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          ip = i + 1
          dudt(i,j,k) = - dxi*( p(ip,j,k)-p(i,j,k) )
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    return
  end subroutine momxp
  !
  subroutine momyp(nx,ny,nz,dyi,p,dvdt)
    implicit none
    integer, intent(in) :: nx,ny,nz
    real(8), intent(in) :: dyi 
    real(8), dimension(0:,0:,0:), intent(in) :: p
    real(8), dimension(:,:,:), intent(out) :: dvdt
    integer :: i,j,k
    integer :: jp
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,jp) &
    !$OMP SHARED(nx,ny,nz,dyi,p,dvdt)
    do k=1,nz
       do j=1,ny
          jp = j + 1
          do i=1,nx
             !
             ! Momentum balance
             !
             dvdt(i,j,k) = - dyi*( p(i,jp,k)-p(i,j,k) )
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    return
  end subroutine momyp
  !
  subroutine momzp(nx,ny,nz,dzci,p,dwdt)
    implicit none
    integer, intent(in) :: nx,ny,nz
    real(8), intent(in), dimension(0:) :: dzci
    real(8), dimension(0:,0:,0:), intent(in) :: p
    real(8), dimension(:,:,:), intent(out) :: dwdt
    integer :: kp
    integer :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,kp) &
    !$OMP SHARED(nx,ny,nz,p,dwdt,dzci)
    do k=1,nz
       kp = k + 1
       do j=1,ny
          do i=1,nx
             !
             ! Momentum balance
             !
             dwdt(i,j,k) = - dzci(k)*( p(i,j,kp)-p(i,j,k) )
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    return
  end subroutine momzp
end module mod_mom
