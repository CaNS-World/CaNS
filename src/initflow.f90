module mod_initflow
  use mpi
  use decomp_2d
  use mod_common_mpi, only: ierr,ijk_start,myid
  use mod_param     , only: dims,pi,dx,dy,dz,lx,ly,lz,uref,lref,is_wallturb,bforce
  use mod_types
  implicit none
  private
  public initflow,add_noise
  contains
  subroutine initflow(inivel,n,zclzi,dzclzi,dzflzi,visc,u,v,w,p)
    !
    ! computes initial conditions for the velocity field
    !
    implicit none
    character(len=3), intent(in) :: inivel
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(0:) :: zclzi,dzclzi,dzflzi
    real(rp), intent(in) :: visc
    real(rp), dimension(0:,0:,0:), intent(out) :: u,v,w,p
    real(rp), allocatable, dimension(:) :: u1d
    !real(rp), allocatable, dimension(:,:) :: u2d
    integer :: i,j,k,kk
    real(rp) :: q
    logical :: is_noise,is_mean,is_pair
    real(rp) :: xc,yc,zc,xf,yf,zf
    real(rp) :: reb,retau
    !
    allocate(u1d(n(3)))
    is_noise = .false.
    is_mean  = .false.
    is_pair  = .false.
    q = .5
    select case(trim(inivel))
    case('cou')
      call couette(   q,n(3),zclzi,uref,u1d)
    case('poi')
      call poiseuille(q,n(3),zclzi,uref,u1d)
      is_mean=.true.
    case('zer')
      u1d(:) = 0.
    case('log')
      call log_profile(q,n(3),zclzi,visc,u1d)
      is_noise = .true.
      is_mean = .true.
    case('hcl')
      deallocate(u1d)
      allocate(u1d(2*n(3)))
      call log_profile(q,2*n(3),zclzi,visc,u1d)
      is_noise = .true.
      is_mean=.true.
    case('hcp')
      deallocate(u1d)
      allocate(u1d(2*n(3)))
      call poiseuille(q,2*n(3),zclzi,uref,u1d)
      is_mean = .true.
    case('tgv')
      do k=1,n(3)
        zc = zclzi(k)*2.*pi
        do j=1,n(2)
          yc = (j+ijk_start(2)-.5)*dy/ly*2.*pi
          yf = (j+ijk_start(2)-.0)*dy/ly*2.*pi
          do i=1,n(1)
            xc = (i+ijk_start(1)-.5)*dx/lx*2.*pi
            xf = (i+ijk_start(1)-.0)*dx/lx*2.*pi
            u(i,j,k) =  sin(xf)*cos(yc)*cos(zc)
            v(i,j,k) = -cos(xc)*sin(yf)*cos(zc)
            w(i,j,k) = 0.
            p(i,j,k) = 0.!(cos(2.*xc)+cos(2.*yc))*(cos(2.*zc)+2.)/16.
          enddo
        enddo
      enddo
    case('pdc')
      if(is_wallturb) then ! turbulent flow
        retau  = (bforce(1)*lref)**.5*uref/visc
        reb    = (retau/.09)**(1./.88)
        uref   = (reb/2.)/retau
      else                 ! laminar flow
        uref = (bforce(1)*lref**2/(3.*visc))
      endif
      call poiseuille(q,n(3),zclzi,uref,u1d)
      is_mean=.true.
    case default
      if(myid.eq.0) print*, 'ERROR: invalid name for initial velocity field'
      if(myid.eq.0) print*, ''
      if(myid.eq.0) print*, '*** Simulation abortited due to errors in the case file ***'
      if(myid.eq.0) print*, '    check INFO_INPUT.md'
      call decomp_2d_finalize
      call MPI_FINALIZE(ierr)
      call exit
    end select
    if(inivel.ne.'tgv') then
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            u(i,j,k) = u1d(k)
            v(i,j,k) = 0.
            w(i,j,k) = 0.
            p(i,j,k) = 0.
          enddo
        enddo
      enddo
    endif
    if(is_noise) then
      call add_noise(n,123,.5_rp,u(1:n(1),1:n(2),1:n(3)))
      call add_noise(n,456,.5_rp,v(1:n(1),1:n(2),1:n(3)))
      call add_noise(n,789,.5_rp,w(1:n(1),1:n(2),1:n(3)))
    endif
    if(is_mean) then
      call set_mean(n,uref,dzflzi,u(1:n(1),1:n(2),1:n(3)))
    endif
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
        zc = 2.*zclzi(k) - 1. ! z rescaled to be between -1 and +1
        zf = 2.*(zclzi(k) + .5*dzflzi(k)) - 1.
        do j=1,n(2)
          yc = ((ijk_start(2)+j-0.5)*dy-.5*ly)*2./lz
          yf = ((ijk_start(2)+j-0.0)*dy-.5*ly)*2./lz
          do i=1,n(1)
            xc = ((ijk_start(1)+i-0.5)*dx-.5*lx)*2./lz
            xf = ((ijk_start(1)+i-0.0)*dx-.5*lx)*2./lz
            !u(i,j,k) = u1d(k)
            v(i,j,k) = -1. * gxy(yf,xc)*dfz(zc) * uref
            w(i,j,k) =  1. * fz(zf)*dgxy(yc,xc) * uref
            p(i,j,k) = 0.
          enddo
        enddo
      enddo
      !
      ! alternatively, using a Taylor-Green vortex 
      ! for the cross-stream velocity components
      ! (commented below)
      !
      !do k=1,n(3)
      !  zc = (zclzi(k)              )*2.*pi
      !  zf = (zclzi(k)+0.5*dzclzi(k))*2.*pi
      !  do j=1,n(2)
      !    yc = (j+ijk_start(2)-.5)*dy/ly*2.*pi
      !    yf = (j+ijk_start(2)-.0)*dy/ly*2.*pi
      !    do i=1,n(1)
      !      xc = (i+ijk_start(1)-.5)*dx/lx*2.*pi
      !      xf = (i+ijk_start(1)-.0)*dx/lx*2.*pi
      !      !u(i,j,k) = u1d(k)
      !      v(i,j,k) =  sin(xc)*cos(yf)*cos(zc)
      !      w(i,j,k) = -cos(xc)*sin(yc)*cos(zf)
      !      p(i,j,k) = 0.!(cos(2.*xc)+cos(2.*yc))*(cos(2.*zc)+2.)/16.
      !    enddo
      !  enddo
      !enddo
    endif
    deallocate(u1d)
    return
  end subroutine initflow
  !
  subroutine add_noise(n,iseed,norm,p)
    implicit none
    integer , intent(in), dimension(3) :: n
    integer , intent(in) :: iseed
    real(rp), intent(in) :: norm 
    real(rp), intent(inout), dimension(n(1),n(2),n(3)) :: p
    integer(4), allocatable, dimension(:) :: seed
    real(rp) :: rn
    integer, dimension(3) :: ng
    integer :: i,j,k,ii,jj,kk
    allocate(seed(64))
    seed(:) = iseed
    call random_seed( put = seed )
    ng(1:3) = n(1:3)*dims(1:3)
    do k=1,ng(3)
      kk = k - ijk_start(3)
      do j=1,ng(2)
        jj = j - ijk_start(2)
        do i=1,ng(1)
          ii = i - ijk_start(1)
          call random_number(rn)
          if(ii.ge.1.and.ii.le.n(1) .and. &
             jj.ge.1.and.jj.le.n(2) .and. &
             kk.ge.1.and.kk.le.n(3) ) then
             p(ii,jj,kk) = p(ii,jj,kk) + 2.*(rn-.5)*norm
          endif
        enddo
      enddo
    enddo
    return
  end subroutine add_noise
  !
  subroutine set_mean(n,mean,dzlzi,p)
  implicit none
  integer , intent(in), dimension(3) :: n
  real(rp), intent(in), dimension(0:) :: dzlzi 
  real(rp), intent(in) :: mean
  real(rp), intent(inout), dimension(n(1),n(2),n(3)) :: p
  real(rp) :: meanold
  integer :: i,j,k
    meanold = 0.
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,p,dzlzi) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP REDUCTION(+:meanold)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          meanold = meanold + p(i,j,k)*dzlzi(k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    call mpi_allreduce(MPI_IN_PLACE,meanold,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    meanold = meanold/(1.*n(1)*dims(1)*n(2)*dims(2))
    !
    if(meanold.ne.0.) then
      !$OMP WORKSHARE
      p(:,:,:) = p(:,:,:)/meanold*mean
      !$OMP END WORKSHARE
    endif
    return
  end subroutine set_mean
  !
  subroutine couette(q,n,zc,norm,p)
    !
    ! plane couette profile normalized by the wall velocity difference
    !
    implicit none
    real(rp), intent(in)   :: q
    integer , intent(in)   :: n
    real(rp), intent(in), dimension(0:) :: zc
    real(rp), intent(in)   :: norm
    real(rp), intent(out), dimension(n) :: p
    integer :: k
    real(rp) :: z
    do k=1,n
      z    = zc(k)!1.*((k-1)+q)/(1.*n)
      p(k) = .5*(1.-2.*z)*norm
    enddo
    return
  end subroutine couette
  !
  subroutine poiseuille(q,n,zc,norm,p)
    implicit none
    real(rp), intent(in)   :: q
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
      z    = zc(k)!1.*((k-1)+q)/(1.*n)
      p(k) = 6.*z*(1.-z)*norm
    enddo
    return
  end subroutine poiseuille
  !
  subroutine log_profile(q,n,zc,visc,p)
    implicit none
    real(rp), intent(in)   :: q
    integer , intent(in)   :: n
    real(rp), intent(in), dimension(0:) :: zc
    real(rp), intent(in)   :: visc
    real(rp), intent(out), dimension(n) :: p
    integer :: k
    real(rp) :: z,reb,retau ! z/lz and bulk Reynolds number
    reb = lref*uref/visc
    retau = 0.09*reb**(0.88) ! from Pope's book
    do k=1,n/2
      z    = zc(k)*2.*retau!1.*((k-1)+q)/(1.*n)*2.*retau
      p(k) = 2.5*log(z) + 5.5
      if (z.le.11.6) p(k)=z
      p(n+1-k) = p(k)
    enddo
    return
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
