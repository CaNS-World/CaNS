module mod_mom
  use mpi
  use mod_param     , only: bforce
  use mod_common_mpi, only: ierr,dims,is_bound
  use mod_types
  implicit none
  private
  public momxad,momyad,momzad,momxp,momyp,momzp
  contains
  subroutine momxad(nx,ny,nz,dxi,dyi,dzi,dzci,dzfi,dzflzi,visc,u,v,w,dudt,taux)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,dzi,visc
    real(rp), intent(in), dimension(0:) :: dzci,dzfi,dzflzi
    real(rp), dimension(0:,0:,0:), intent(in) :: u,v,w
    real(rp), dimension(:,:,:), intent(out) :: dudt
    real(rp), dimension(3)  , intent(out) :: taux
    integer :: im,ip,jm,jp,km,kp,i,j,k
    real(rp) :: uuip,uuim,uvjp,uvjm,uwkp,uwkm
    real(rp) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm
    integer :: nxg,nyg,nzg
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,im,jm,km,ip,jp,kp) &
    !$OMP PRIVATE(uuip,uuim,uvjp,uvjm,uwkp,uwkm) &
    !$OMP PRIVATE(dudxp,dudxm,dudyp,dudym,dudzp,dudzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzi,visc,u,v,w,dudt,dzci,dzfi,bforce)
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
                        dzfi(k)*( -uwkp + uwkm ) + (dudzp-dudzm)*visc*dzfi(k) + &
                        bforce(1)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    taux(:) = 0.
    if( is_bound(0,2) ) then
      do k=1,nz
        do i=1,nx
          dudyp = (u(i,1 ,k)-u(i,0   ,k))*dyi*visc*dzflzi(k)
          taux(2) = taux(2) + dudyp
        enddo
      enddo
    endif
    if( is_bound(1,2) ) then
      do k=1,nz
        do i=1,nx
          dudym = (u(i,ny,k)-u(i,ny+1,k))*dyi*visc*dzflzi(k)
          taux(2) = taux(2) + dudyp
        enddo
      enddo
    endif
    if( is_bound(0,3) ) then
      do j=1,ny
        do i=1,nx
          dudzp = (u(i,j,1 )-u(i,j,0   ))*dzci(0)*visc
          taux(3) = taux(3) + dudzp
        enddo
      enddo
    endif
    if( is_bound(1,3) ) then
      do j=1,ny
        do i=1,nx
          dudzm = (u(i,j,nz)-u(i,j,nz+1))*dzci(nz)*visc
          taux(3) = taux(3) + dudzm
        enddo
      enddo
    endif
    call mpi_allreduce(MPI_IN_PLACE,taux(1),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    nxg = nx*dims(1)
    nyg = ny*dims(2)
    nzg = nz*dims(3)
    taux(1) = taux(1)/(1.*nyg)
    taux(2) = taux(2)/(1.*nxg)
    taux(3) = taux(3)/(1.*nxg*nyg)
    return
  end subroutine momxad
  !
  subroutine momyad(nx,ny,nz,dxi,dyi,dzi,dzci,dzfi,dzflzi,visc,u,v,w,dvdt,tauy)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,dzi,visc
    real(rp), intent(in), dimension(0:) :: dzci,dzfi,dzflzi
    real(rp), dimension(0:,0:,0:), intent(in) :: u,v,w
    real(rp), dimension(:,:,:), intent(out) :: dvdt
    real(rp), dimension(3), intent(out) :: tauy
    integer :: im,ip,jm,jp,km,kp,i,j,k
    real(rp) :: uvip,uvim,vvjp,vvjm,wvkp,wvkm
    real(rp) :: dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm
    integer :: nxg,nyg,nzg
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,im,jm,km,ip,jp,kp) &
    !$OMP PRIVATE(uvip,uvim,vvjp,vvjm,wvkp,wvkm) &
    !$OMP PRIVATE(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzi,visc,u,v,w,dvdt,dzci,dzfi,bforce)
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
                        dzfi(k)*( -wvkp + wvkm ) + (dvdzp-dvdzm)*visc*dzfi(k)+ &
                        bforce(2)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    tauy(:) = 0.
    if( is_bound(0,1) ) then
      do k=1,nz
        do j=1,ny
          dvdxp = (v(1 ,j,k)-v(0   ,j,k))*dxi*visc*dzflzi(k)
          tauy(1) = tauy(1) + dvdxp
        enddo
      enddo
    endif
    if( is_bound(1,1) ) then
      do k=1,nz
        do j=1,ny
          dvdxm = (v(nx,j,k)-v(nx+1,j,k))*dxi*visc*dzflzi(k)
          tauy(1) = tauy(1) + dvdxm
        enddo
      enddo
    endif
    if(is_bound(0,3) ) then
      do j=1,ny
        do i=1,nx
          dvdzp = (v(i,j,1 )-v(i,j,0   ))*dzci(0)*visc
          tauy(3) = tauy(3) + dvdzp
        enddo
      enddo
    endif
    if(is_bound(1,3) ) then
      do j=1,ny
        do i=1,nx
          dvdzm = (v(i,j,nz)-v(i,j,nz+1))*dzci(nz)*visc
          tauy(3) = tauy(3) + dvdzm
        enddo
      enddo
    endif
    call mpi_allreduce(MPI_IN_PLACE,tauy(1),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    nxg = nx*dims(1)
    nyg = ny*dims(2)
    nzg = nz*dims(3)
    tauy(1) = tauy(1)/(1.*nyg)
    tauy(2) = tauy(2)/(1.*nxg)
    tauy(3) = tauy(3)/(1.*nxg*nyg)
    return
  end subroutine momyad
  !
  subroutine momzad(nx,ny,nz,dxi,dyi,dzi,dzci,dzfi,dzflzi,visc,u,v,w,dwdt,tauz)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,dzi,visc
    real(rp), intent(in), dimension(0:) :: dzci,dzfi,dzflzi
    real(rp), dimension(0:,0:,0:), intent(in) :: u,v,w
    real(rp), dimension(:,:,:), intent(out) :: dwdt
    real(rp), dimension(3), intent(out) :: tauz
    integer :: im,ip,jm,jp,km,kp,i,j,k
    real(rp) :: uwip,uwim,vwjp,vwjm,wwkp,wwkm
    real(rp) :: dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
    integer :: nxg,nyg,nzg
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,im,jm,km,ip,jp,kp) &
    !$OMP PRIVATE(uwip,uwim,vwjp,vwjm,wwkp,wwkm) &
    !$OMP PRIVATE(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzi,visc,u,v,w,dwdt,dzci,dzfi,bforce)
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
                        dzci(k)*( -wwkp + wwkm ) + (dwdzp-dwdzm)*visc*dzci(k)+ &
                        bforce(3)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    tauz(:) = 0.
    if(is_bound(0,1) ) then
      do k=1,nz
        do j=1,ny
          dwdxp = (w(1 ,j,k)-w(0   ,j,k))*dxi*visc*dzflzi(k)
          tauz(1) = tauz(1) + dwdxp
        enddo
      enddo
    endif
    if(is_bound(1,1) ) then
      do k=1,nz
        do j=1,ny
          dwdxm = (w(nx,j,k)-w(nx+1,j,k))*dxi*visc*dzflzi(k)
          tauz(1) = tauz(1) + dwdxm
        enddo
      enddo
    endif
    if(is_bound(0,2) ) then
      do k=1,nz
        do i=1,nx
          dwdyp = (w(i,1,k )-w(i,0   ,k))*dyi*visc*dzflzi(k)
          tauz(2) = tauz(2) + dwdyp
        enddo
      enddo
    endif
    if(is_bound(1,2) ) then
      do k=1,nz
        do i=1,nx
          dwdym = (w(i,ny,k)-w(i,ny+1,k))*dyi*visc*dzflzi(k)
          tauz(2) = tauz(2) + dwdym
        enddo
      enddo
    endif
    call mpi_allreduce(MPI_IN_PLACE,tauz(1),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    nxg = nx*dims(1)
    nyg = ny*dims(2)
    nzg = nz*dims(3)
    tauz(1) = tauz(1)/(1.*nyg)
    tauz(2) = tauz(2)/(1.*nxg)
    tauz(3) = tauz(3)/(1.*nxg*nyg)
    return
  end subroutine momzad
  !
  subroutine momxp(nx,ny,nz,dxi,p,dudt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi 
    real(rp), dimension(0:,0:,0:), intent(in) :: p
    real(rp), dimension(:,:,:), intent(out) :: dudt
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
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dyi 
    real(rp), dimension(0:,0:,0:), intent(in) :: p
    real(rp), dimension(:,:,:), intent(out) :: dvdt
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
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in), dimension(0:) :: dzci
    real(rp), dimension(0:,0:,0:), intent(in) :: p
    real(rp), dimension(:,:,:), intent(out) :: dwdt
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
          dwdt(i,j,k) = - dzci(k)*( p(i,j,kp)-p(i,j,k) )
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    return
  end subroutine momzp
end module mod_mom
