module mod_mom
  use mpi
  use decomp_2d     , only: nx_global,ny_global,nz_global
  use mod_param     , only: bforce
  use mod_common_mpi, only: ierr
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
    integer :: i,j,k
    real(rp) :: uuip,uuim,uvjp,uvjm,uwkp,uwkm
    real(rp) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm
    integer :: nxg,nyg,nzg
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(uuip,uuim,uvjp,uvjm,uwkp,uwkm) &
    !$OMP PRIVATE(dudxp,dudxm,dudyp,dudym,dudzp,dudzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzi,visc,u,v,w,dudt,dzci,dzfi,bforce)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          uuip  = 0.25*( u(i+1,j,k)+u(i,j,k) )*( u(i+1,j  ,k  )+u(i,j  ,k  ) )
          uuim  = 0.25*( u(i-1,j,k)+u(i,j,k) )*( u(i-1,j  ,k  )+u(i,j  ,k  ) )
          uvjp  = 0.25*( u(i,j+1,k)+u(i,j,k) )*( v(i+1,j  ,k  )+v(i,j  ,k  ) )
          uvjm  = 0.25*( u(i,j-1,k)+u(i,j,k) )*( v(i+1,j-1,k  )+v(i,j-1,k  ) )
          uwkp  = 0.25*( u(i,j,k+1)+u(i,j,k) )*( w(i+1,j  ,k  )+w(i,j  ,k  ) )
          uwkm  = 0.25*( u(i,j,k-1)+u(i,j,k) )*( w(i+1,j  ,k-1)+w(i,j  ,k-1) )
          dudxp = (u(i+1,j,k)-u(i,j,k))*dxi
          dudxm = (u(i,j,k)-u(i-1,j,k))*dxi
          dudyp = (u(i,j+1,k)-u(i,j,k))*dyi
          dudym = (u(i,j,k)-u(i,j-1,k))*dyi
          dudzp = (u(i,j,k+1)-u(i,j,k))*dzci(k  )
          dudzm = (u(i,j,k)-u(i,j,k-1))*dzci(k-1)
          !
          ! Momentum balance
          !
          dudt(i,j,k) = dxi*(     -uuip + uuim ) + (dudxp-dudxm)*visc*dxi + &
                        dyi*(     -uvjp + uvjm ) + (dudyp-dudym)*visc*dyi + &
                        dzfi(k)*( -uwkp + uwkm ) + (dudzp-dudzm)*visc*dzfi(k) + &
                        bforce(1)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    taux(:) = 0.
    do k=1,nz
      do i=1,nx
        dudyp = (u(i,1 ,k)-u(i,0   ,k))*dyi*visc*dzflzi(k)
        dudym = (u(i,ny,k)-u(i,ny+1,k))*dyi*visc*dzflzi(k)
        taux(2) = taux(2) + (dudyp+dudym)
      end do
    end do
    do j=1,ny
      do i=1,nx
        dudzp = (u(i,j,1 )-u(i,j,0   ))*dzci(0)*visc
        dudzm = (u(i,j,nz)-u(i,j,nz+1))*dzci(nz)*visc
        taux(3) = taux(3) + (dudzp+dudzm)
      end do
    end do
    call mpi_allreduce(MPI_IN_PLACE,taux(1),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    taux(1) = taux(1)/(1.*nyg)
    taux(2) = taux(2)/(1.*nxg)
    taux(3) = taux(3)/(1.*nxg*nyg)
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
    integer :: i,j,k
    real(rp) :: uvip,uvim,vvjp,vvjm,wvkp,wvkm
    real(rp) :: dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm
    integer :: nxg,nyg,nzg
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(uvip,uvim,vvjp,vvjm,wvkp,wvkm) &
    !$OMP PRIVATE(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzi,visc,u,v,w,dvdt,dzci,dzfi,bforce)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          uvip  = 0.25*( u(i  ,j,k  )+u(i  ,j+1,k  ) )*( v(i,j,k)+v(i+1,j,k) )
          uvim  = 0.25*( u(i-1,j,k  )+u(i-1,j+1,k  ) )*( v(i,j,k)+v(i-1,j,k) )
          vvjp  = 0.25*( v(i  ,j,k  )+v(i  ,j+1,k  ) )*( v(i,j,k)+v(i,j+1,k) )
          vvjm  = 0.25*( v(i  ,j,k  )+v(i  ,j-1,k  ) )*( v(i,j,k)+v(i,j-1,k) )
          wvkp  = 0.25*( w(i  ,j,k  )+w(i  ,j+1,k  ) )*( v(i,j,k+1)+v(i,j,k) )
          wvkm  = 0.25*( w(i  ,j,k-1)+w(i  ,j+1,k-1) )*( v(i,j,k-1)+v(i,j,k) )
          dvdxp = (v(i+1,j,k)-v(i,j,k))*dxi
          dvdxm = (v(i,j,k)-v(i-1,j,k))*dxi
          dvdyp = (v(i,j+1,k)-v(i,j,k))*dyi
          dvdym = (v(i,j,k)-v(i,j-1,k))*dyi
          dvdzp = (v(i,j,k+1)-v(i,j,k))*dzci(k  )
          dvdzm = (v(i,j,k)-v(i,j,k-1))*dzci(k-1)
          !
          ! Momentum balance
          !
          dvdt(i,j,k) = dxi*(     -uvip + uvim ) + (dvdxp-dvdxm)*visc*dxi+ &
                        dyi*(     -vvjp + vvjm ) + (dvdyp-dvdym)*visc*dyi+ &
                        dzfi(k)*( -wvkp + wvkm ) + (dvdzp-dvdzm)*visc*dzfi(k)+ &
                        bforce(2)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    tauy(:) = 0.
    do k=1,nz
      do j=1,ny
        dvdxp = (v(1 ,j,k)-v(0   ,j,k))*dxi*visc*dzflzi(k)
        dvdxm = (v(nx,j,k)-v(nx+1,j,k))*dxi*visc*dzflzi(k)
        tauy(1) = tauy(1) + (dvdxp+dvdxm)
      end do
    end do
    do j=1,ny
      do i=1,nx
        dvdzp = (v(i,j,1 )-v(i,j,0   ))*dzci(0)*visc
        dvdzm = (v(i,j,nz)-v(i,j,nz+1))*dzci(nz)*visc
        tauy(3) = tauy(3) + (dvdzp+dvdzm)
      end do
    end do
    call mpi_allreduce(MPI_IN_PLACE,tauy(1),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    nxg = nx_global
    nyg = ny_global
    nzg = nz_global
    tauy(1) = tauy(1)/(1.*nyg)
    tauy(2) = tauy(2)/(1.*nxg)
    tauy(3) = tauy(3)/(1.*nxg*nyg)
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
    integer :: i,j,k
    real(rp) :: uwip,uwim,vwjp,vwjm,wwkp,wwkm
    real(rp) :: dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
    integer :: nxg,nyg,nzg
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(uwip,uwim,vwjp,vwjm,wwkp,wwkm) &
    !$OMP PRIVATE(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzi,visc,u,v,w,dwdt,dzci,dzfi,bforce)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          uwip  = 0.25*( w(i,j,k)+w(i+1,j,k) )*( u(i  ,j  ,k)+u(i  ,j  ,k+1) )
          uwim  = 0.25*( w(i,j,k)+w(i-1,j,k) )*( u(i-1,j  ,k)+u(i-1,j  ,k+1) )
          vwjp  = 0.25*( w(i,j,k)+w(i,j+1,k) )*( v(i  ,j  ,k)+v(i  ,j  ,k+1) )
          vwjm  = 0.25*( w(i,j,k)+w(i,j-1,k) )*( v(i  ,j-1,k)+v(i  ,j-1,k+1) )
          wwkp  = 0.25*( w(i,j,k)+w(i,j,k+1) )*( w(i  ,j  ,k)+w(i  ,j  ,k+1) )
          wwkm  = 0.25*( w(i,j,k)+w(i,j,k-1) )*( w(i  ,j  ,k)+w(i  ,j  ,k-1) )
          dwdxp = (w(i+1,j,k)-w(i,j,k))*dxi
          dwdxm = (w(i,j,k)-w(i-1,j,k))*dxi
          dwdyp = (w(i,j+1,k)-w(i,j,k))*dyi
          dwdym = (w(i,j,k)-w(i,j-1,k))*dyi
          dwdzp = (w(i,j,k+1)-w(i,j,k))*dzfi(k+1)
          dwdzm = (w(i,j,k)-w(i,j,k-1))*dzfi(k  )
          !
          ! Momentum balance
          !
          dwdt(i,j,k) = dxi*(     -uwip + uwim ) + (dwdxp-dwdxm)*visc*dxi+ &
                        dyi*(     -vwjp + vwjm ) + (dwdyp-dwdym)*visc*dyi+ &
                        dzci(k)*( -wwkp + wwkm ) + (dwdzp-dwdzm)*visc*dzci(k)+ &
                        bforce(3)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    tauz(:) = 0.
    do k=1,nz
      do j=1,ny
        dwdxp = (w(1 ,j,k)-w(0   ,j,k))*dxi*visc*dzflzi(k)
        dwdxm = (w(nx,j,k)-w(nx+1,j,k))*dxi*visc*dzflzi(k)
        tauz(1) = tauz(1) + (dwdxp+dwdxm)
      end do
    end do
    do k=1,nz
      do i=1,nx
        dwdyp = (w(i,1,k )-w(i,0   ,k))*dyi*visc*dzflzi(k)
        dwdym = (w(i,ny,k)-w(i,ny+1,k))*dyi*visc*dzflzi(k)
        tauz(2) = tauz(2) + (dwdyp+dwdym)
      end do
    end do
    call mpi_allreduce(MPI_IN_PLACE,tauz(1),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    nxg = nx_global
    nyg = ny_global
    nzg = nz_global
    tauz(1) = tauz(1)/(1.*nyg)
    tauz(2) = tauz(2)/(1.*nxg)
    tauz(3) = tauz(3)/(1.*nxg*nyg)
  end subroutine momzad
  !
  subroutine momxp(nx,ny,nz,dxi,p,dudt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi
    real(rp), dimension(0:,0:,0:), intent(in) :: p
    real(rp), dimension(:,:,:), intent(out) :: dudt
    integer :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(nx,ny,nz,dxi,p,dudt)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dudt(i,j,k) = - dxi*( p(i+1,j,k)-p(i,j,k) )
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine momxp
  !
  subroutine momyp(nx,ny,nz,dyi,p,dvdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dyi
    real(rp), dimension(0:,0:,0:), intent(in) :: p
    real(rp), dimension(:,:,:), intent(out) :: dvdt
    integer :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(nx,ny,nz,dyi,p,dvdt)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dvdt(i,j,k) = - dyi*( p(i,j+1,k)-p(i,j,k) )
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine momyp
  !
  subroutine momzp(nx,ny,nz,dzci,p,dwdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in), dimension(0:) :: dzci
    real(rp), dimension(0:,0:,0:), intent(in) :: p
    real(rp), dimension(:,:,:), intent(out) :: dwdt
    integer :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(nx,ny,nz,p,dwdt,dzci)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          dwdt(i,j,k) = - dzci(k)*( p(i,j,k+1)-p(i,j,k) )
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine momzp
  subroutine cmpt_wallshear(n,is_bound,l,dl,dzc,dzf,visc,u,v,w,taux,tauy,tauz)
    implicit none
    integer , intent(in ), dimension(3) :: n
    logical , intent(in ), dimension(0:1,3) :: is_bound
    real(rp), intent(in ), dimension(3)     :: l,dl
    real(rp), intent(in ), dimension(0:)    :: dzc,dzf
    real(rp), intent(in )                   :: visc
    real(rp), intent(in ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out), dimension(3) :: taux,tauy,tauz
    real(rp) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm, &
                dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm, &
                dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
    integer :: i,j,k,nx,ny,nz
    !
    nx = n(1)
    ny = n(2)
    nz = n(3)
    taux(:) = 0.
    if(is_bound(0,2)) then
      do k=1,nz
        do i=1,nx
          dudyp = (u(i,1 ,k)-u(i,0   ,k))/dl(2)*visc
          taux(2) = taux(2) + dudyp*(dl(1)*dzf(k))/(l(1)*l(3))
        enddo
      enddo
    endif
    if(is_bound(1,2)) then
      do k=1,nz
        do i=1,nx
          dudym = (u(i,ny,k)-u(i,ny+1,k))/dl(2)*visc
          taux(2) = taux(2) + dudyp*(dl(1)*dzf(k))/(l(1)*l(3))
        enddo
      enddo
    endif
    if(is_bound(0,3)) then
      do j=1,ny
        do i=1,nx
          dudzp = (u(i,j,1 )-u(i,j,0   ))/dzc(0)*visc
          taux(3) = taux(3) + dudzp*(dl(1)*dl(2))/(l(1)*l(2))
        enddo
      enddo
    endif
    if(is_bound(1,3)) then
      do j=1,ny
        do i=1,nx
          dudzm = (u(i,j,nz)-u(i,j,nz+1))/dzc(nz)*visc
          taux(3) = taux(3) + dudzm*(dl(1)*dl(2))/(l(1)*l(2))
        enddo
      enddo
    endif
    call mpi_allreduce(MPI_IN_PLACE,taux(1),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    !
    tauy(:) = 0.
    if(is_bound(0,1)) then
      do k=1,nz
        do j=1,ny
          dvdxp = (v(1 ,j,k)-v(0   ,j,k))/dl(1)*visc
          tauy(1) = tauy(1) + dvdxp*(dl(2)*dzf(1))/(l(2)*l(3))
        enddo
      enddo
    endif
    if(is_bound(1,1)) then
      do k=1,nz
        do j=1,ny
          dvdxm = (v(nx,j,k)-v(nx+1,j,k))/dl(1)*visc
          tauy(1) = tauy(1) + dvdxm*(dl(2)*dzf(1))/(l(2)*l(3))
        enddo
      enddo
    endif
    if(is_bound(0,3)) then
      do j=1,ny
        do i=1,nx
          dvdzp = (v(i,j,1 )-v(i,j,0   ))/dzc(0)*visc
          tauy(3) = tauy(3) + dvdzp*(dl(1)*dl(2))/(l(1)*l(2))
        enddo
      enddo
    endif
    if(is_bound(1,3)) then
      do j=1,ny
        do i=1,nx
          dvdzm = (v(i,j,nz)-v(i,j,nz+1))/dzc(nz)*visc
          tauy(3) = tauy(3) + dvdzm*(dl(1)*dl(2))/(l(1)*l(2))
        enddo
      enddo
    endif
    call mpi_allreduce(MPI_IN_PLACE,tauy(1),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    !
    tauz(:) = 0.
    if(is_bound(0,1)) then
      do k=1,nz
        do j=1,ny
          dwdxp = (w(1 ,j,k)-w(0   ,j,k))/dl(1)*visc
          tauz(1) = tauz(1) + dwdxp*(dl(2)*dzf(k))/(l(2)*l(3))
        enddo
      enddo
    endif
    if(is_bound(1,1)) then
      do k=1,nz
        do j=1,ny
          dwdxm = (w(nx,j,k)-w(nx+1,j,k))/dl(1)*visc
          tauz(1) = tauz(1) + dwdxm*(dl(2)*dzf(k))/(l(2)*l(3))
        enddo
      enddo
    endif
    if(is_bound(0,2)) then
      do k=1,nz
        do i=1,nx
          dwdyp = (w(i,1,k )-w(i,0   ,k))/dl(2)*visc
          tauz(2) = tauz(2) + dwdyp*(dl(1)*dzf(k))/(l(1)*l(3))
        enddo
      enddo
    endif
    if(is_bound(1,2)) then
      do k=1,nz
        do i=1,nx
          dwdym = (w(i,ny,k)-w(i,ny+1,k))/dl(2)*visc
          tauz(2) = tauz(2) + dwdym*(dl(1)*dzf(k))/(l(1)*l(3))
        enddo
      enddo
    endif
    call mpi_allreduce(MPI_IN_PLACE,tauz(1),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  end subroutine cmpt_wallshear
end module mod_mom
