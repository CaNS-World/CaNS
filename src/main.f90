!
!        CCCCCCCCCCCCC                    NNNNNNNN        NNNNNNNN    SSSSSSSSSSSSSSS
!     CCC::::::::::::C                    N:::::::N       N::::::N  SS:::::::::::::::S
!   CC:::::::::::::::C                    N::::::::N      N::::::N S:::::SSSSSS::::::S
!  C:::::CCCCCCCC::::C                    N:::::::::N     N::::::N S:::::S     SSSSSSS
! C:::::C       CCCCCC   aaaaaaaaaaaaa    N::::::::::N    N::::::N S:::::S
!C:::::C                 a::::::::::::a   N:::::::::::N   N::::::N S:::::S
!C:::::C                 aaaaaaaaa:::::a  N:::::::N::::N  N::::::N  S::::SSSS
!C:::::C                          a::::a  N::::::N N::::N N::::::N   SS::::::SSSSS
!C:::::C                   aaaaaaa:::::a  N::::::N  N::::N:::::::N     SSS::::::::SS
!C:::::C                 aa::::::::::::a  N::::::N   N:::::::::::N        SSSSSS::::S
!C:::::C                a::::aaaa::::::a  N::::::N    N::::::::::N             S:::::S
! C:::::C       CCCCCC a::::a    a:::::a  N::::::N     N:::::::::N             S:::::S
!  C:::::CCCCCCCC::::C a::::a    a:::::a  N::::::N      N::::::::N SSSSSSS     S:::::S
!   CC:::::::::::::::C a:::::aaaa::::::a  N::::::N       N:::::::N S::::::SSSSSS:::::S
!     CCC::::::::::::C  a::::::::::aa:::a N::::::N        N::::::N S:::::::::::::::SS
!        CCCCCCCCCCCCC   aaaaaaaaaa  aaaa NNNNNNNN         NNNNNNN  SSSSSSSSSSSSSSS
!-------------------------------------------------------------------------------------
! CaNS -- Canonical Navier-Stokes Solver
! Pedro Costa (p.simoes.costa@gmail.com)
!-------------------------------------------------------------------------------------
program fhs
  use iso_c_binding  , only: C_PTR
  use mpi
  use decomp_2d
  use mod_bound      , only: boundp,bounduvw
  use mod_chkdiv     , only: chkdiv
  use mod_chkdt      , only: chkdt
  use mod_common_mpi , only: myid,ierr
  use mod_correc     , only: correc
  use mod_debug      , only: chkmean,chkhelmholtz
  use mod_fft        , only: fftini,fftend
  use mod_fillps     , only: fillps
  use mod_initflow   , only: initflow
  use mod_initgrid   , only: initgrid
  use mod_initmpi    , only: initmpi
  use mod_initsolver , only: initsolver
  use mod_load       , only: load
  use mod_rk         , only: rk
  use mod_output     , only: out0d,out1d,out1d_2,out2d,out3d
  use mod_param      , only: itot,jtot,ktot,lx,ly,lz,dx,dy,dz,dxi,dyi,dzi,visc,pi,cbcvel,bcvel,cbcpre, &
                             icheck,iout0d,iout1d,iout2d,iout3d,isave, &
                             nstep,restart, &
                             rkcoeff, &
                             datadir, &
                             cfl,     &
                             inivel,  &
                             uref,lref, &
                             forceinx,forceiny,forceinz, &
                             imax,jmax,dims, &
                             nthreadsmax, &
                             gr, &
                             ioutflowdir
  use mod_solver     , only: solver
  !$ use omp_lib
  implicit none
  integer, parameter, dimension(3) :: ng = (/itot,jtot,ktot/)
  integer, parameter, dimension(3) :: n  = (/imax,jmax,ktot/)
  real(8), parameter, dimension(3) :: l   = (/lx,ly,lz/)
  real(8), parameter, dimension(3) :: dl  = (/dx,dy,dz/)
  real(8), parameter, dimension(3) :: dli = (/dxi,dyi,dzi/)
  logical, parameter, dimension(3) :: isforced = (/forceinx,forceiny,forceinz/)
  real(8), dimension(0:imax+1,0:jmax+1,0:ktot+1) :: u,v,w,p,up,vp,wp,pp
  real(8), dimension(imax,jmax,ktot)    :: dudtrko,dvdtrko,dwdtrko
  real(8), dimension(3) :: tauxo,tauyo,tauzo
  real(8), dimension(3) :: f
  type(C_PTR), dimension(2,2) :: arrplanp
  real(8), dimension(imax,jmax) :: lambdaxyp
  real(8), dimension(ktot) :: ap,bp,cp
  real(8) :: normfftp
#ifdef IMPDIFF
  type(C_PTR), dimension(2,2) :: arrplanu,arrplanv,arrplanw
  real(8), dimension(imax,jmax) :: lambdaxyu,lambdaxyv,lambdaxyw
  real(8), dimension(ktot) :: au,av,aw,bu,bv,bw,cu,cv,cw
  real(8), dimension(ktot) :: aa,bb,cc
  real(8) :: normfftu,normfftv,normfftw
  real(8) :: alpha,alphai
  integer :: i,j,k,im,ip,jm,jp,km,kp
  real(8) :: mean
#endif
  real(8) :: ristep
  real(8) :: dt,dti,dtmax,time,dtrk,dtrki
  integer :: irk,istep
  real(8), dimension(0:ktot+1) :: dzc,dzf,zc,zf,dzci,dzfi
  real(8) :: meanvel
  real(8), dimension(3) :: dpdl
  !real(8), allocatable, dimension(:) :: var
  real(8), dimension(10) :: var
#ifdef TIMING
  real(8) :: dt12,dt12av,dt12min,dt12max
#endif
  character(len=7) :: fldnum
  integer :: lenr
  !
  !$call omp_set_num_threads(nthreadsmax)
  call initmpi(ng,cbcpre)
  call initgrid(inivel,n(3),dz,gr,lz,dzc,dzf,zc,zf)
  if(myid.eq.0) then
    inquire (iolength=lenr) dzc(1)
    open(99,file=trim(datadir)//'grid.bin',access='direct',recl=4*n(3)*lenr)
    write(99,rec=1) dzc(1:n(3)),dzf(1:n(3)),zc(1:n(3)),zf(1:n(3))
    close(99)
  endif
  dzci = dzc**(-1)
  dzfi = dzf**(-1)
  if(.not.restart) then
    istep = 0
    time = 0.d0
    call initflow(inivel,n,zc/lz,dzc/lz,dzf/lz,visc,uref,u,v,w,p)
  else
    call load('r',trim(datadir)//'fld.bin',n,u(1:n(1),1:n(2),1:n(3)), &
                                             v(1:n(1),1:n(2),1:n(3)), &
                                             w(1:n(1),1:n(2),1:n(3)), &
                                             p(1:n(1),1:n(2),1:n(3)), &
                                             time,ristep)
    istep = nint(ristep)
  endif
  call bounduvw(cbcvel,n,bcvel,ioutflowdir,dl,dzc,dzf,u,v,w)
  call boundp(n,cbcpre,pp)
  call chkdt(n,dl,dzci,dzfi,visc,u,v,w,dtmax)
  dt = cfl*dtmax
  if(myid.eq.0) print*, 'dtmax = ', dtmax, 'dt = ',dt
  dti = 1.d0/dt
  !
  ! initialize Poisson solver
  !
  call initsolver(n,dli,dzci,dzfi,cbcpre       ,lambdaxyp,(/'c','c','c'/),ap,bp,cp,arrplanp,normfftp)
#ifdef IMPDIFF
  call initsolver(n,dli,dzci,dzfi,cbcvel(:,:,1),lambdaxyu,(/'f','c','c'/),au,bu,cu,arrplanu,normfftu)
  call initsolver(n,dli,dzci,dzfi,cbcvel(:,:,2),lambdaxyv,(/'c','f','c'/),av,bv,cv,arrplanv,normfftv)
  call initsolver(n,dli,dzci,dzfi,cbcvel(:,:,3),lambdaxyw,(/'c','c','f'/),aw,bw,cw,arrplanw,normfftw)
#endif
  !
  ! main loop
  !
  do while(istep.lt.nstep)
#ifdef TIMING
    dt12 = MPI_WTIME()
#endif
    istep = istep + 1
    time = time + dt
    if(myid.eq.0) print*, 'Timestep #', istep, 'Time = ', time
    dpdl(:) = 0.
    tauxo(:) = 0.
    tauyo(:) = 0.
    tauzo(:) = 0.
    do irk=1,3
      dtrk = sum(rkcoeff(:,irk))*dt
      dtrki = dtrk**(-1)
      call rk(rkcoeff(:,irk),n,dli,dzci,dzfi,dzf/lz,visc,dt,l, &
              u,v,w,p,dudtrko,dvdtrko,dwdtrko,tauxo,tauyo,tauzo,up,vp,wp,f)
#ifdef IMPDIFF
      alpha = -1.d0/(.5d0*visc*dtrk)
      bb(:) = bu(:) + alpha
      call solver(n,arrplanu,normfftu,lambdaxyu,au,bb,cu,cbcvel(:,3,1),(/'f','c','c'/),up(1:imax,1:jmax,1:ktot))
      bb(:) = bv(:) + alpha
      call solver(n,arrplanv,normfftv,lambdaxyv,av,bb,cv,cbcvel(:,3,2),(/'c','f','c'/),vp(1:imax,1:jmax,1:ktot))
      bb(:) = bw(:) + alpha
      call solver(n,arrplanw,normfftw,lambdaxyw,aw,bb,cw,cbcvel(:,3,3),(/'c','c','f'/),wp(1:imax,1:jmax,1:ktot))
#else
    if(isforced(1)) up(1:n(1),1:n(2),1:n(3)) = up(1:n(1),1:n(2),1:n(3)) + f(1)
    if(isforced(2)) vp(1:n(1),1:n(2),1:n(3)) = vp(1:n(1),1:n(2),1:n(3)) + f(2)
    if(isforced(3)) wp(1:n(1),1:n(2),1:n(3)) = wp(1:n(1),1:n(2),1:n(3)) + f(3)
#endif
    dpdl(:) = dpdl(:) + f(:)
#ifdef DEBUG
      if(isforced(1)) then
        call chkmean(n,dzf/lz,up,meanvel)
        if(myid.eq.0) print*,'Mean u = ', meanvel
      endif
      if(isforced(2)) then
        call chkmean(n,dzf/lz,v,meanvel)
        if(myid.eq.0) print*,'Mean v = ', meanvel
      endif
      if(isforced(3)) then
        call chkmean(n,dzf/lz,w,meanvel)
        if(myid.eq.0) print*,'Mean w = ', meanvel
      endif
#endif
      call bounduvw(cbcvel,n,bcvel,0,dl,dzc,dzf,up,vp,wp)
      call fillps(n,dli,dzfi,dtrki,up,vp,wp,pp)
      call solver(n,arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre(:,3),(/'c','c','c'/),pp(1:imax,1:jmax,1:ktot))
      call boundp(n,cbcpre,pp)
      call correc(n,dli,dzci,dtrk,pp,up,vp,wp,u,v,w)
      call bounduvw(cbcvel,n,bcvel,ioutflowdir,dl,dzc,dzf,u,v,w)
#ifdef IMPDIFF
      alphai = alpha**(-1)
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP PRIVATE(i,j,k,im,jm,km,ip,jp,kp) &
      !$OMP SHARED(p,pp,dzfi,dzci,alphai)
      do k=1,n(3)
        kp = k + 1
        km = k - 1
        do j=1,n(2)
          jp = j + 1
          jm = j - 1
          do i=1,n(1)
            ip = i + 1
            im = i - 1
            p(i,j,k) = p(i,j,k) + pp(i,j,k) + alphai*( &
                        (pp(ip,j,k)-2.*pp(i,j,k)+pp(im,j,k))*(dxi**2) + &
                        (pp(i,jp,k)-2.*pp(i,j,k)+pp(i,jm,k))*(dyi**2) + &
                        ((pp(i,j,kp)-pp(i,j,k))*dzci(k) - &
                         (pp(i,j,k)-pp(i,j,km))*dzci(km))*dzfi(k) )
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
#else
      !$OMP WORKSHARE
      p(:,:,:) = p(:,:,:) + pp(:,:,:)
      !$OMP END WORKSHARE
#endif
      call boundp(n,cbcpre,p)
    enddo
    dpdl(:) = -dpdl(:)*dti
    if(mod(istep,icheck).eq.0) then
      if(myid.eq.0) print*, 'Checking stability and divergence...'
      call chkdt(n,dl,dzci,dzfi,visc,u,v,w,dtmax)
      dt  = cfl*dtmax
      if(myid.eq.0) print*, 'dtmax = ', dtmax, 'dt = ',dt
      dti = 1.d0/dt
      call chkdiv(n,dli,dzfi,u,v,w)
    endif
    if(mod(istep,iout0d).eq.0) then
      !allocate(var(4))
      var(1) = 1.*istep
      var(2) = dt
      var(3) = time 
      call out0d(trim(datadir)//'time.out',3,var)
      if(any(isforced(:))) then
        var(1)   = time
        var(2:4) = dpdl(1:3)
        call out0d(trim(datadir)//'forcing.out',4,var)
      endif
      !deallocate(var)
    endif
    write(fldnum,'(i7.7)') istep
    if(mod(istep,iout1d).eq.0) then
      !call out1d(trim(datadir)//'umean_z_fld_'   //fldnum//'.out',n,3,zc/lz,dzf/lz,u)
      !call out1d(trim(datadir)//'vmean_z_fld_'   //fldnum//'.out',n,3,zc/lz,dzf/lz,v)
      !call out1d(trim(datadir)//'wmean_z_fld_'   //fldnum//'.out',n,3,zf/lz,dzc/lz,w)
      !call out1d(trim(datadir)//'umean_y_fld_'   //fldnum//'.out',n,2,zc/lz,dzf/lz,u)
      !call out1d(trim(datadir)//'vmean_y_fld_'   //fldnum//'.out',n,2,zf/lz,dzc/lz,v)
      !call out1d(trim(datadir)//'wmean_y_fld_'   //fldnum//'.out',n,2,zc/lz,dzf/lz,w)
      !call out1d_2(trim(datadir)//'velstats_fld_'//fldnum//'.out',n,3,zc/lz,u,v,w)
      call out1d(trim(datadir)//'umean_z.out',n,3,zc/lz,dzf/lz,u)
      call out1d(trim(datadir)//'vmean_z.out',n,3,zc/lz,dzf/lz,v)
      call out1d(trim(datadir)//'wmean_z.out',n,3,zf/lz,dzc/lz,w)
      call out1d(trim(datadir)//'umean_y.out',n,2,zc/lz,dzf/lz,u)
      call out1d(trim(datadir)//'vmean_y.out',n,2,zf/lz,dzc/lz,v)
      call out1d(trim(datadir)//'wmean_y.out',n,2,zc/lz,dzf/lz,w)
      call out1d_2(trim(datadir)//'velstats.out',n,3,zc/lz,u,v,w)
    endif
    if(mod(istep,iout2d).eq.0) then
      !call out2d(trim(datadir)//'fld_u_slice_fld_'//fldnum//'.bin',2,n(3)/2,u(1:n(1),1:n(2),1:n(3)))
      call out2d(trim(datadir)//'fld_u_slice.bin',2,n(2)/2,u(1:n(1),1:n(2),1:n(3)))
    endif
    if(mod(istep,iout3d).eq.0) then
      !call out3d(trim(datadir)//'fld_u.bin_fld_'//fldnum//'.bin',(/1,1,1/),u(1:n(1),1:n(2),1:n(3)))
      call out3d(trim(datadir)//'fld_u.bin',(/1,1,1/),u(1:n(1),1:n(2),1:n(3)))
    endif
    if(mod(istep,isave ).eq.0) then
      ristep = 1.d0*istep
      call load('w',trim(datadir)//'fld.bin',n,u(1:n(1),1:n(2),1:n(3)), &
                                               v(1:n(1),1:n(2),1:n(3)), &
                                               w(1:n(1),1:n(2),1:n(3)), &
                                               p(1:n(1),1:n(2),1:n(3)), &
                                               time,ristep)
    endif
#ifdef TIMING
      dt12 = MPI_WTIME()-dt12
      call MPI_ALLREDUCE(dt12,dt12av ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(dt12,dt12min,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(dt12,dt12max,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) print*, 'Avrg, min & max elapsed time: '
      if(myid.eq.0) print*, dt12av/(1.d0*product(dims)),dt12min,dt12max
#endif
  enddo
  !
  ! clear ffts
  !
  call fftend(arrplanp)
  if(myid.eq.0) print*, '*** Fim ***'
  call decomp_2d_finalize
  call MPI_FINALIZE(ierr)
  call exit
end program fhs
