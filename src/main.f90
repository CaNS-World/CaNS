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
program cans
  use iso_c_binding  , only: C_PTR
  use mpi
  use decomp_2d
  use mod_bound      , only: boundp,bounduvw,updt_rhs_b
  use mod_chkdiv     , only: chkdiv
  use mod_chkdt      , only: chkdt
  use mod_common_mpi , only: myid,ierr
  use mod_correc     , only: correc
  use mod_debug      , only: chkmean
  use mod_fft        , only: fftini,fftend
  use mod_fillps     , only: fillps
  use mod_initflow   , only: initflow
  use mod_initgrid   , only: initgrid
  use mod_initmpi    , only: initmpi
  use mod_initsolver , only: initsolver
  use mod_load       , only: load
  use mod_rk         , only: rk,rk_id
  use mod_output     , only: out0d,out1d,out1d_2,out2d,out3d
  use mod_param      , only: itot,jtot,ktot,lx,ly,lz,dx,dy,dz,dxi,dyi,dzi,visc,small, &
                             cbcvel,bcvel,cbcpre,bcpre, &
                             icheck,iout0d,iout1d,iout2d,iout3d,isave, &
                             nstep,restart, &
                             rkcoeff, &
                             datadir, &
                             cfl,     &
                             inivel,  &
                             uref,lref, &
                             imax,jmax,dims, &
                             nthreadsmax, &
                             gr, &
                             is_outflow,no_outflow,is_forced
  use mod_sanity     , only: test_sanity
  use mod_solver     , only: solver
  !$ use omp_lib
  implicit none
  integer, parameter, dimension(3) :: ng = (/itot,jtot,ktot/)
  integer, parameter, dimension(3) :: n  = (/imax,jmax,ktot/)
  real(8), parameter, dimension(3) :: l   = (/lx,ly,lz/)
  real(8), parameter, dimension(3) :: dl  = (/dx,dy,dz/)
  real(8), parameter, dimension(3) :: dli = (/dxi,dyi,dzi/)
  real(8), dimension(0:imax+1,0:jmax+1,0:ktot+1) :: u,v,w,p,up,vp,wp,pp
  real(8), dimension(imax,jmax,ktot)    :: dudtrko,dvdtrko,dwdtrko
  real(8), dimension(3) :: tauxo,tauyo,tauzo
  real(8), dimension(3) :: f
  type(C_PTR), dimension(2,2) :: arrplanp
  real(8), dimension(imax,jmax) :: lambdaxyp
  real(8), dimension(ktot) :: ap,bp,cp
  real(8) :: normfftp
  type rhs_bound
    real(8), dimension(n(2),n(3),0:1) :: x
    real(8), dimension(n(1),n(3),0:1) :: y
    real(8), dimension(n(1),n(2),0:1) :: z
  end type rhs_bound 
#ifdef IMPDIFF
  type(C_PTR), dimension(2,2) :: arrplanu,arrplanv,arrplanw
  real(8), dimension(imax,jmax) :: lambdaxyu,lambdaxyv,lambdaxyw
  real(8), dimension(ktot) :: au,av,aw,bu,bv,bw,bb,cu,cv,cw
  real(8) :: normfftu,normfftv,normfftw
  real(8) :: alpha,alphai
  integer :: i,j,k,im,ip,jm,jp,km,kp
  type(rhs_bound) :: rhsbu,rhsbv,rhsbw
#endif
  type(rhs_bound) :: rhsbp
  real(8) :: ristep
  real(8) :: dt,dti,dtmax,time,dtrk,dtrki,divtot,divmax
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
  integer :: lenr,kk
  logical :: kill
  !
  !$call omp_set_num_threads(nthreadsmax)
  call initmpi(ng,cbcpre)
  if(myid.eq.0) print*, '******************************'
  if(myid.eq.0) print*, '*** Beginning of simulation ***'
  if(myid.eq.0) print*, '******************************'
  if(myid.eq.0) print*, ''
  call initgrid(inivel,n(3),gr,lz,dzc,dzf,zc,zf)
  if(myid.eq.0) then
    inquire (iolength=lenr) dzc(1)
    open(99,file=trim(datadir)//'grid.bin',access='direct',recl=4*n(3)*lenr)
    write(99,rec=1) dzc(1:n(3)),dzf(1:n(3)),zc(1:n(3)),zf(1:n(3))
    close(99)
    open(99,file=trim(datadir)//'grid.out')
    do kk=0,ktot+1
      write(99,'(5E15.7)') 0.d0,zf(kk),zc(kk),dzf(kk),dzc(kk)
    enddo
    close(99)
  endif
  !
  ! test input files before proceeding with the calculation
  !
  call test_sanity(ng,n,dims,cbcvel,cbcpre,bcvel,bcpre,is_outflow,is_forced, &
                   dli,dzci,dzfi)
  !
  dzci = dzc**(-1)
  dzfi = dzf**(-1)
  if(.not.restart) then
    istep = 0
    time = 0.d0
    call initflow(inivel,n,zc/lz,dzc/lz,dzf/lz,visc,uref,u,v,w,p)
    if(myid.eq.0) print*, '*** Initial condition succesfully set ***'
  else
    call load('r',trim(datadir)//'fld.bin',n,u(1:n(1),1:n(2),1:n(3)), &
                                             v(1:n(1),1:n(2),1:n(3)), &
                                             w(1:n(1),1:n(2),1:n(3)), &
                                             p(1:n(1),1:n(2),1:n(3)), &
                                             time,ristep)
    istep = nint(ristep)
    if(myid.eq.0) print*, '*** Checkpoint loaded at time = ', time, 'time step = ', istep, '. ***'
  endif
  call bounduvw(cbcvel,n,bcvel,is_outflow,dl,dzc,dzf,u,v,w)
  call boundp(cbcpre,n,bcpre,dl,dzc,dzf,p)
  call chkdt(n,dl,dzci,dzfi,visc,u,v,w,dtmax)
  dt = cfl*dtmax
  if(myid.eq.0) print*, 'dtmax = ', dtmax, 'dt = ',dt
  dti = 1.d0/dt
  kill = .false.
  !
  ! initialize Poisson solver
  !
  call initsolver(n,dli,dzci,dzfi,cbcpre,bcpre(:,:),lambdaxyp,(/'c','c','c'/),ap,bp,cp,arrplanp,normfftp,rhsbp%x,rhsbp%y,rhsbp%z)
#ifdef IMPDIFF
  call initsolver(n,dli,dzci,dzfi,cbcvel(:,:,1),bcvel(:,:,1),lambdaxyu,(/'f','c','c'/),au,bu,cu,arrplanu,normfftu, &
                  rhsbu%x,rhsbu%y,rhsbu%z)
  call initsolver(n,dli,dzci,dzfi,cbcvel(:,:,2),bcvel(:,:,2),lambdaxyv,(/'c','f','c'/),av,bv,cv,arrplanv,normfftv, &
                  rhsbv%x,rhsbv%y,rhsbv%z)
  call initsolver(n,dli,dzci,dzfi,cbcvel(:,:,3),bcvel(:,:,3),lambdaxyw,(/'c','c','f'/),aw,bw,cw,arrplanw,normfftw, &
                  rhsbw%x,rhsbw%y,rhsbw%z)
#endif
  !
  ! main loop
  !
  if(myid.eq.0) print*, '*** Calculation loop starts now ***'
  do while(istep.lt.nstep)
#ifdef TIMING
    dt12 = MPI_WTIME()
#endif
    istep = istep + 1
    time = time + dt
    if(myid.eq.0) print*, 'Timestep #', istep, 'Time = ', time
    dpdl(:)  = 0.d0
    tauxo(:) = 0.d0
    tauyo(:) = 0.d0
    tauzo(:) = 0.d0
    do irk=1,3
      dtrk = sum(rkcoeff(:,irk))*dt
      dtrki = dtrk**(-1)
#ifndef IMPDIFF
      call rk(rkcoeff(:,irk),n,dli,dzci,dzfi,dzf/lz,dzc/lz,visc,dt,l, &
                 u,v,w,p,dudtrko,dvdtrko,dwdtrko,tauxo,tauyo,tauzo,up,vp,wp,f)
#else
      call rk_id(rkcoeff(:,irk),n,dli,dzci,dzfi,dzf/lz,dzc/lz,visc,dt,l, &
                 u,v,w,p,dudtrko,dvdtrko,dwdtrko,tauxo,tauyo,tauzo,up,vp,wp,f)
#endif
      if(is_forced(1)) up(1:n(1),1:n(2),1:n(3)) = up(1:n(1),1:n(2),1:n(3)) + f(1)
      if(is_forced(2)) vp(1:n(1),1:n(2),1:n(3)) = vp(1:n(1),1:n(2),1:n(3)) + f(2)
      if(is_forced(3)) wp(1:n(1),1:n(2),1:n(3)) = wp(1:n(1),1:n(2),1:n(3)) + f(3)
#ifdef IMPDIFF
      alpha = -1.d0/(.5d0*visc*dtrk)
      !$OMP WORKSHARE
      up(:,:,:) = up(:,:,:)*alpha
      !$OMP END WORKSHARE
      bb(:) = bu(:) + alpha
      call updt_rhs_b((/'f','c','c'/),cbcvel(:,:,1),n,rhsbu%x,rhsbu%y,rhsbu%z,up(1:imax,1:jmax,1:ktot))
      call solver(n,arrplanu,normfftu,lambdaxyu,au,bb,cu,cbcvel(:,3,1),(/'f','c','c'/),up(1:imax,1:jmax,1:ktot))
      !$OMP WORKSHARE
      vp(:,:,:) = vp(:,:,:)*alpha
      !$OMP END WORKSHARE
      bb(:) = bv(:) + alpha
      call updt_rhs_b((/'c','f','c'/),cbcvel(:,:,2),n,rhsbv%x,rhsbv%y,rhsbv%z,vp(1:imax,1:jmax,1:ktot))
      call solver(n,arrplanv,normfftv,lambdaxyv,av,bb,cv,cbcvel(:,3,2),(/'c','f','c'/),vp(1:imax,1:jmax,1:ktot))
      !$OMP WORKSHARE
      wp(:,:,:) = wp(:,:,:)*alpha
      !$OMP END WORKSHARE
      bb(:) = bw(:) + alpha
      call updt_rhs_b((/'c','c','f'/),cbcvel(:,:,3),n,rhsbw%x,rhsbw%y,rhsbw%z,wp(1:imax,1:jmax,1:ktot))
      call solver(n,arrplanw,normfftw,lambdaxyw,aw,bb,cw,cbcvel(:,3,3),(/'c','c','f'/),wp(1:imax,1:jmax,1:ktot))
#endif
      dpdl(:) = dpdl(:) + f(:)
#ifdef DEBUG
      if(is_forced(1)) then
        call chkmean(n,dzf/lz,up,meanvel)
        if(myid.eq.0) print*,'Mean u = ', meanvel
      endif
      if(is_forced(2)) then
        call chkmean(n,dzf/lz,vp,meanvel)
        if(myid.eq.0) print*,'Mean v = ', meanvel
      endif
      if(is_forced(3)) then
        call chkmean(n,dzc/lz,wp,meanvel)
        if(myid.eq.0) print*,'Mean w = ', meanvel
      endif
#endif
      call bounduvw(cbcvel,n,bcvel,no_outflow,dl,dzc,dzf,up,vp,wp) ! outflow BC only at final velocity
      call fillps(n,dli,dzfi,dtrki,up,vp,wp,pp)
      call updt_rhs_b((/'c','c','c'/),cbcpre,n,rhsbp%x,rhsbp%y,rhsbp%z,pp(1:imax,1:jmax,1:ktot))
      call solver(n,arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre(:,3),(/'c','c','c'/),pp(1:imax,1:jmax,1:ktot))
      call boundp(cbcpre,n,bcpre,dl,dzc,dzf,pp)
      call correc(n,dli,dzci,dtrk,pp,up,vp,wp,u,v,w)
      call bounduvw(cbcvel,n,bcvel,is_outflow,dl,dzc,dzf,u,v,w)
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
                        (pp(ip,j,k)-2.d0*pp(i,j,k)+pp(im,j,k))*(dxi**2) + &
                        (pp(i,jp,k)-2.d0*pp(i,j,k)+pp(i,jm,k))*(dyi**2) + &
                        ((pp(i,j,kp)-pp(i,j,k ))*dzci(k ) - &
                         (pp(i,j,k )-pp(i,j,km))*dzci(km))*dzfi(k) )
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
#else
      !$OMP WORKSHARE
      p(:,:,:) = p(:,:,:) + pp(:,:,:)
      !$OMP END WORKSHARE
#endif
      call boundp(cbcpre,n,bcpre,dl,dzc,dzf,p)
    enddo
    dpdl(:) = -dpdl(:)*dti
    if(mod(istep,icheck).eq.0) then
      if(myid.eq.0) print*, 'Checking stability and divergence...'
      call chkdt(n,dl,dzci,dzfi,visc,u,v,w,dtmax)
      dt  = cfl*dtmax
      if(myid.eq.0) print*, 'dtmax = ', dtmax, 'dt = ',dt
      if(dtmax.lt.small) then
        if(myid.eq.0) print*, 'ERROR: timestep is too small.'
        if(myid.eq.0) print*, 'Aborting ...'
        istep = nstep + 1 ! i.e. exit main loop
        kill = .true.
      endif
      dti = 1.d0/dt
      call chkdiv(n,dli,dzfi,u,v,w,divtot,divmax)
      if(divmax.gt.small.or.divtot.ne.divtot) then
        if(myid.eq.0) print*, 'ERROR: maximum divergence is too large.'
        if(myid.eq.0) print*, 'Aborting ...'
        istep = nstep + 1 ! i.e. exit main loop
        kill = .true.
      endif
    endif
    !
    ! output routines below
    !
    if(mod(istep,iout0d).eq.0) then
      !allocate(var(4))
      var(1) = 1.d0*istep
      var(2) = dt
      var(3) = time
      call out0d(trim(datadir)//'time.out',3,var)
      if(any(is_forced(:))) then
        var(1)   = time
        var(2:4) = dpdl(1:3)
        call out0d(trim(datadir)//'forcing.out',4,var)
      endif
      !deallocate(var)
    endif
    write(fldnum,'(i7.7)') istep
    if(mod(istep,iout1d).eq.0) then
      include 'out1d.h90'
    endif
    if(mod(istep,iout2d).eq.0) then
      include 'out2d.h90'
    endif
    if(mod(istep,iout3d).eq.0) then
      include 'out3d.h90'
    endif
    if(mod(istep,isave ).eq.0) then
      ristep = 1.d0*istep
      call load('w',trim(datadir)//'fld.bin',n,u(1:n(1),1:n(2),1:n(3)), &
                                               v(1:n(1),1:n(2),1:n(3)), &
                                               w(1:n(1),1:n(2),1:n(3)), &
                                               p(1:n(1),1:n(2),1:n(3)), &
                                               time,ristep)
    if(myid.eq.0) print*, '*** Checkpoint saved at time = ', time, 'time step = ', istep, '. ***'
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
#ifdef IMPDIFF
  call fftend(arrplanu)
  call fftend(arrplanv)
  call fftend(arrplanw)
#endif
  if(myid.eq.0.and.(.not.kill)) print*, '*** Fim ***'
  call decomp_2d_finalize
  call MPI_FINALIZE(ierr)
  call exit
end program cans
