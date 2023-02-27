! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_initflow
  use mpi
  use mod_common_mpi, only: ierr,myid
  use mod_param     , only: pi,dx,dy,dz,lx,ly,lz,is_wallturb,is_forced,velf,bforce,bcvel
  use mod_types
  implicit none
  private
  public initflow,add_noise
  contains
  subroutine initflow(inivel,ng,lo,zc,dzc,dzf,visc,u,v,w,p)
    !
    ! computes initial conditions for the velocity field
    !
    implicit none
    character(len=3), intent(in) :: inivel
    integer , intent(in), dimension(3) :: ng,lo
    real(rp), intent(in), dimension(0:) :: zc,dzc,dzf
    real(rp), intent(in) :: visc
    real(rp), dimension(0:,0:,0:), intent(out) :: u,v,w,p
    real(rp), allocatable, dimension(:) :: u1d
    !real(rp), allocatable, dimension(:,:) :: u2d
    integer :: i,j,k
    logical :: is_noise,is_mean,is_pair
    real(rp) :: xc,yc,zcc,xf,yf,zff
    real(rp), allocatable, dimension(:) :: zc2
    real(rp) :: uref,lref
    real(rp) :: ubulk,reb,retau
    integer, dimension(3) :: n
    !
    n(:) = shape(p) - 2*1
    allocate(u1d(n(3)))
    is_noise = .false.
    is_mean  = .false.
    is_pair  = .false.
    uref  = 1.
    ubulk = uref
    if(is_forced(1)) ubulk = velf(1)
    select case(trim(inivel))
    case('cou')
      uref = (bcvel(0,3,1)-bcvel(1,3,1)) ! uref = (ubot - utop)
      call couette(   n(3),zc/lz,uref ,u1d)
      uref = abs(uref)
    case('poi')
      call poiseuille(n(3),zc/lz,ubulk,u1d)
      is_mean = .true.
    case('tbl')
      call temporal_bl(n(3),zc,1._rp,visc,uref,u1d)
    case('iop') ! reversed 'poi'
      !
      ! convective reference frame moving with velocit `ubulk`;
      ! walls have negative velocity equal to `ubulk` in the laboratory frame
      !
      ubulk = 0.5*abs(bcvel(0,3,1)+bcvel(1,3,1))
      call poiseuille(n(3),zc/lz,ubulk,u1d)
      u1d(:) = u1d(:) - ubulk
      is_mean = .true.
    case('zer')
      u1d(:) = 0.
    case('uni')
      u1d(:) = uref
    case('log')
      reb = ubulk*lz/visc
      call log_profile(n(3),zc/lz,reb,u1d)
      is_noise = .true.
      is_mean = .true.
    case('hcl')
      deallocate(u1d)
      allocate(u1d(2*n(3)))
      allocate(zc2(0:2*n(3)+1))
      zc2(1     :  n(3)) =        zc(1   :n(3): 1)
      zc2(n(3)+1:2*n(3)) = 2*lz - zc(n(3):1   :-1)
      zc2(0)        = -zc(0)
      zc2(2*n(3)+1) = 2*lz + zc(0)
      reb = ubulk*(2*lz)/visc
      call log_profile(2*n(3),zc2/(2*lz),reb,u1d)
      is_noise = .true.
      is_mean = .true.
    case('hcp')
      deallocate(u1d)
      allocate(u1d(2*n(3)))
      allocate(zc2(0:2*n(3)+1))
      zc2(1     :  n(3)) =        zc(1   :n(3): 1)
      zc2(n(3)+1:2*n(3)) = 2*lz - zc(n(3):1   :-1)
      zc2(0)        = -zc(0)
      zc2(2*n(3)+1) = 2*lz + zc(0)
      call poiseuille(2*n(3),zc2/(2*lz),ubulk,u1d)
      is_mean = .true.
    case('tgv')
      do k=1,n(3)
        zcc = zc(k)/lz*2.*pi
        do j=1,n(2)
          yc = (j+lo(2)-1-.5)*dy/ly*2.*pi
          yf = (j+lo(2)-1-.0)*dy/ly*2.*pi
          do i=1,n(1)
            xc = (i+lo(1)-1-.5)*dx/lx*2.*pi
            xf = (i+lo(1)-1-.0)*dx/lx*2.*pi
            u(i,j,k) =  sin(xf)*cos(yc)*cos(zcc)*uref
            v(i,j,k) = -cos(xc)*sin(yf)*cos(zcc)*uref
            w(i,j,k) = 0.
            p(i,j,k) = 0.!(cos(2.*xc)+cos(2.*yc))*(cos(2.*zcc)+2.)/16.*uref**2
          end do
        end do
      end do
    case('tgw')
      do k=1,n(3)
        do j=1,n(2)
          yc = (j+lo(2)-1-.5)*dy
          yf = (j+lo(2)-1-.0)*dy
          do i=1,n(1)
            xc = (i+lo(1)-1-.5)*dx
            xf = (i+lo(1)-1-.0)*dx
            u(i,j,k) =  cos(xf)*sin(yc)*uref
            v(i,j,k) = -sin(xc)*cos(yf)*uref
            w(i,j,k) = 0.
            p(i,j,k) = -(cos(2.*xc)+cos(2.*yc))/4.*uref**2
          end do
        end do
      end do
    case('pdc','hdc')
      lref  = lz/2.
      if(trim(inivel) /= 'pdc') lref = 2.*lref
      if(is_wallturb) then ! turbulent flow
        uref  = (bforce(1)*lref)**(0.5) ! utau = sqrt(-dpdx*h)
        retau = uref*lref/visc
        reb   = (retau/.09)**(1./.88)
        ubulk = reb*visc/(2*lref)
      else                 ! laminar flow
        ubulk = (bforce(1)*lref**2/(3.*visc))
      end if
      if(trim(inivel) == 'pdc') then
        call poiseuille(n(3),zc/lz,ubulk,u1d)
      else
        deallocate(u1d)
        allocate(u1d(2*n(3)))
        allocate(zc2(0:2*n(3)+1))
        zc2(1     :  n(3)) =        zc(1   :n(3): 1)
        zc2(n(3)+1:2*n(3)) = 2*lz - zc(n(3):1   :-1)
        zc2(0)        = -zc(0)
        zc2(2*n(3)+1) = 2*lz + zc(0)
        call poiseuille(2*n(3),zc2/(2*lz),ubulk,u1d)
      end if
      is_mean = .true.
    case default
      if(myid == 0) print*, 'ERROR: invalid name for initial velocity field'
      if(myid == 0) print*, ''
      if(myid == 0) print*, '*** Simulation aborted due to errors in the case file ***'
      if(myid == 0) print*, '    check INFO_INPUT.md'
      call MPI_FINALIZE(ierr)
      error stop
    end select
    if(.not.any(inivel == ['tgv','tgw'])) then
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            u(i,j,k) = u1d(k)
            v(i,j,k) = 0.
            w(i,j,k) = 0.
            p(i,j,k) = 0.
          end do
        end do
      end do
    end if
    if(is_noise) then
      call add_noise(ng,lo,123,.05_rp,u(1:n(1),1:n(2),1:n(3)))
      call add_noise(ng,lo,456,.05_rp,v(1:n(1),1:n(2),1:n(3)))
      call add_noise(ng,lo,789,.05_rp,w(1:n(1),1:n(2),1:n(3)))
    end if
    if(is_mean) then
      if(trim(inivel) /= 'iop') then
        call set_mean(n,ubulk,dzf/lz*(dx/lx)*(dy/ly),u(1:n(1),1:n(2),1:n(3)))
      else
        call set_mean(n,0._rp,dzf/lz*(dx/lx)*(dy/ly),u(1:n(1),1:n(2),1:n(3)))
      end if
    end if
    if(is_wallturb) is_pair = .true.
    if(is_pair) then
      !
      ! initialize a streamwise vortex pair for a fast transition
      ! to turbulence in a pressure-driven channel:
      !        psi(x,y,z)  = f(z)*g(x,y), with
      !        f(z)        = (1-z**2)**2, and
      !        g(x,y)      = y*exp[-(16x**2-4y**2)]
      ! (x,y,z) --> (streamwise, spanwise, wall-normal) directions
      !
      ! see Henningson and Kim, JFM 1991
      !
      do k=1,n(3)
        zcc = 2.*zc(k)/lz - 1. ! z rescaled to be between -1 and +1
        zff = 2.*(zc(k)/lz + .5*dzf(k)/lz) - 1.
        do j=1,n(2)
          yc = ((lo(2)-1+j-0.5)*dy-.5*ly)*2./lz
          yf = ((lo(2)-1+j-0.0)*dy-.5*ly)*2./lz
          do i=1,n(1)
            xc = ((lo(1)-1+i-0.5)*dx-.5*lx)*2./lz
            xf = ((lo(1)-1+i-0.0)*dx-.5*lx)*2./lz
            !u(i,j,k) = u1d(k)
            v(i,j,k) = -1.*gxy(yf,xc)*dfz(zcc)*ubulk*1.5
            w(i,j,k) =  1.*fz(zff)*dgxy(yc,xc)*ubulk*1.5
            p(i,j,k) = 0.
          end do
        end do
      end do
      !
      ! alternatively, using a Taylor-Green vortex
      ! for the cross-stream velocity components
      ! (commented below)
      !
      !do k=1,n(3)
      !  zcc = (zc(k)/lz              )*2.*pi
      !  zff = (zc(k)/lz+0.5*dzc(k)/lz)*2.*pi
      !  do j=1,n(2)
      !    yc = (j+lo(2)-1-.5)*dy/ly*2.*pi
      !    yf = (j+lo(2)-1-.0)*dy/ly*2.*pi
      !    do i=1,n(1)
      !      xc = (i+lo(1)-1-.5)*dx/lx*2.*pi
      !      xf = (i+lo(1)-1-.0)*dx/lx*2.*pi
      !      !u(i,j,k) = u1d(k)
      !      v(i,j,k) =  sin(xc)*cos(yf)*cos(zcc)*ubulk
      !      w(i,j,k) = -cos(xc)*sin(yc)*cos(zff)*ubulk
      !      p(i,j,k) = 0.!(cos(2.*xc)+cos(2.*yc))*(cos(2.*zcc)+2.)/16.
      !    end do
      !  end do
      !end do
    end if
  end subroutine initflow
  !
  subroutine add_noise(ng,lo,iseed,norm,p)
    implicit none
    integer , intent(in), dimension(3) :: ng,lo
    integer , intent(in) :: iseed
    real(rp), intent(in) :: norm
    real(rp), intent(inout), dimension(:,:,:) :: p
    integer(4), allocatable, dimension(:) :: seed
    real(rp) :: rn
    integer, dimension(3) :: n
    integer :: i,j,k,ii,jj,kk
    !
    n(:) = shape(p)
    allocate(seed(64))
    seed(:) = iseed
    call random_seed( put = seed )
    do k=1,ng(3)
      kk = k-(lo(3)-1)
      do j=1,ng(2)
        jj = j-(lo(2)-1)
        do i=1,ng(1)
          ii = i-(lo(1)-1)
          call random_number(rn)
          if(ii >= 1.and.ii <= n(1) .and. &
             jj >= 1.and.jj <= n(2) .and. &
             kk >= 1.and.kk <= n(3) ) then
             p(ii,jj,kk) = p(ii,jj,kk) + 2.*(rn-.5)*norm
          end if
        end do
      end do
    end do
  end subroutine add_noise
  !
  subroutine set_mean(n,mean,grid_vol_ratio,p)
  implicit none
  integer , intent(in), dimension(3) :: n
  real(rp), intent(in), dimension(0:) :: grid_vol_ratio
  real(rp), intent(in) :: mean
  real(rp), intent(inout), dimension(:,:,:) :: p
  real(rp) :: meanold
  integer :: i,j,k
  meanold = 0.
  !$OMP PARALLEL DO DEFAULT(none) &
  !$OMP SHARED(n,p,grid_vol_ratio) &
  !$OMP REDUCTION(+:meanold)
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        meanold = meanold + p(i,j,k)*grid_vol_ratio(k)
      end do
    end do
  end do
  call MPI_ALLREDUCE(MPI_IN_PLACE,meanold,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  !
  if(meanold /= 0.) then
    !$OMP PARALLEL WORKSHARE
    p(:,:,:) = p(:,:,:)/meanold*mean
    !$OMP END PARALLEL WORKSHARE
  end if
  end subroutine set_mean
  !
  subroutine couette(n,zc,norm,p)
    !
    ! plane couette profile normalized by the wall velocity difference
    !
    implicit none
    integer , intent(in)   :: n
    real(rp), intent(in), dimension(0:) :: zc
    real(rp), intent(in)   :: norm
    real(rp), intent(out), dimension(n) :: p
    integer :: k
    real(rp) :: z
    do k=1,n
      z    = zc(k)
      p(k) = .5*(1.-2.*z)*norm
    end do
  end subroutine couette
  !
  subroutine poiseuille(n,zc,norm,p)
    implicit none
    integer , intent(in)   :: n
    real(rp), intent(in), dimension(0:) :: zc
    real(rp), intent(in)   :: norm
    real(rp), intent(out), dimension(n) :: p
    integer :: k
    real(rp) :: z
    !
    ! plane poiseuille profile normalized by the bulk velocity
    !
    do k=1,n
      z    = zc(k)
      p(k) = 6.*z*(1.-z)*norm
    end do
  end subroutine poiseuille
  !
  subroutine temporal_bl(n,zc,d,nu,norm,p)
    implicit none
    integer , intent(in )   :: n
    real(rp), intent(in ), dimension(0:) :: zc
    real(rp), intent(in )   :: d,nu,norm
    real(rp), intent(out), dimension(n) :: p
    integer  :: k
    real(rp) :: theta
    !
    ! temporal boudary layer profile
    ! with thickness d, viscosity nu, and wall velocity norm (at z=0)
    !
    theta = 54.*nu/norm
    do k=1,n
      p(k)=(0.5+(0.5)*tanh((d/(2.*theta))*(1.-zc(k)/d)))*norm
    end do
  end subroutine temporal_bl
  !
  subroutine log_profile(n,zc,reb,p)
    implicit none
    integer , intent(in)   :: n
    real(rp), intent(in), dimension(0:) :: zc
    real(rp), intent(in)   :: reb
    real(rp), intent(out), dimension(n) :: p
    integer :: k
    real(rp) :: z,retau ! z/lz and bulk Reynolds number
    retau = 0.09*reb**(0.88) ! from Pope's book
    do k=1,n
      z = zc(k)*2.*retau
      if(z >= retau) z = 2.*retau-z
      p(k) = 2.5*log(z) + 5.5
      if (z <= 11.6) p(k)=z
    end do
  end subroutine log_profile
  !
  ! functions to initialize the streamwise vortex pair
  ! (explained above)
  !
  function fz(zc)
  real(rp), intent(in) :: zc
  real(rp) :: fz
    fz = ((1.-zc**2)**2)
  end function
  !
  function dfz(zc)
  real(rp), intent(in) :: zc
  real(rp) :: dfz
    dfz = -4.*zc*(1.-zc**2)
  end function
  !
  function gxy(xc,yc)
  real(rp), intent(in) :: xc,yc
  real(rp) :: gxy
    gxy = yc*exp(-4.*(4.*xc**2+yc**2))
  end function
  !
  function dgxy(xc,yc)
  real(rp), intent(in) :: xc,yc
  real(rp) :: dgxy
    dgxy = exp(-4.*(4.*xc**2+yc**2))*(1.-8.*yc**2)
  end function
end module mod_initflow
