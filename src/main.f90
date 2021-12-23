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
  use mod_debug      , only: chk_mean
  use mod_fft        , only: fftini,fftend
  use mod_fillps     , only: fillps
  use mod_initflow   , only: initflow
  use mod_initgrid   , only: initgrid
  use mod_initmpi    , only: initmpi
  use mod_initsolver , only: initsolver
  use mod_load       , only: load
  use mod_rk         , only: rk
  use mod_output     , only: out0d,out1d,out1d_2,out2d,out3d,write_log_output,write_visu_2d,write_visu_3d
  use mod_param      , only: lx,ly,lz,dx,dy,dz,dxi,dyi,dzi,uref,lref,rey,visc,small, &
                             nb,is_bound,cbcvel,bcvel,cbcpre,bcpre, &
                             icheck,iout0d,iout1d,iout2d,iout3d,isave, &
                             nstep,time_max,tw_max,stop_type,restart,is_overwrite_save, &
                             rkcoeff,   &
                             datadir,   &
                             cfl,dtmin, &
                             inivel,    &
                             dims, &
                             nthreadsmax, &
                             gr, &
                             is_forced,velf,bforce, &
                             n,n_z,ng,lo,hi,l,dl,dli, &
                             read_input
  use mod_sanity     , only: test_sanity
  use mod_solver     , only: solver
#if defined(_IMPDIFF) && defined(_IMPDIFF_1D)
  use mod_solver     , only: solver_gaussel_z
#endif
  use mod_types
  !$ use omp_lib
  implicit none
  real(rp), allocatable, dimension(:,:,:) :: u,v,w,p,up,vp,wp,pp
  real(rp), allocatable, dimension(:,:,:)    :: dudtrko,dvdtrko,dwdtrko
  real(rp), dimension(3) :: tauxo,tauyo,tauzo
  real(rp), dimension(3) :: f
  type(C_PTR), dimension(2,2) :: arrplanp
  real(rp), allocatable, dimension(:,:) :: lambdaxyp
  real(rp), allocatable, dimension(:) :: ap,bp,cp
  real(rp) :: normfftp
  type rhs_bound
    real(rp), allocatable, dimension(:,:,:) :: x
    real(rp), allocatable, dimension(:,:,:) :: y
    real(rp), allocatable, dimension(:,:,:) :: z
  end type rhs_bound
  type(rhs_bound) :: rhsbp
#if defined(_IMPDIFF)
  type(C_PTR), dimension(2,2) :: arrplanu,arrplanv,arrplanw
  real(rp), allocatable, dimension(:,:) :: lambdaxyu,lambdaxyv,lambdaxyw
  real(rp), allocatable, dimension(:) :: au,av,aw,bu,bv,bw,bb,cu,cv,cw
  real(rp) :: normfftu,normfftv,normfftw
  real(rp) :: alpha,alphai
  integer :: i,j
  type(rhs_bound) :: rhsbu,rhsbv,rhsbw
#endif
  real(rp) :: dt,dti,dtmax,time,dtrk,dtrki,divtot,divmax
  integer :: irk,istep
  real(rp), allocatable, dimension(:) :: dzc  ,dzf  ,zc  ,zf  ,dzci  ,dzfi, &
                                         dzc_g,dzf_g,zc_g,zf_g,dzci_g,dzfi_g
  real(rp) :: meanvelu,meanvelv,meanvelw
  real(rp), dimension(3) :: dpdl
  !real(rp), allocatable, dimension(:) :: var
  real(rp), dimension(10) :: var
#if defined(_TIMING)
  real(rp) :: dt12,dt12av,dt12min,dt12max
#endif
  real(rp) :: twi,tw
  character(len=7  ) :: fldnum
  character(len=100) :: filename
  integer :: k,kk
  logical :: is_done,kill
  integer :: rlen
  !
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  !
  ! read parameter file
  !
  call read_input(myid)
  !
  ! initialize MPI/OpenMP
  !
  !$call omp_set_num_threads(nthreadsmax)
  call initmpi(ng,dims,cbcpre,n_z,lo,hi,n,nb,is_bound)
  twi = MPI_WTIME()
  !
  ! allocate variables
  !
  allocate(u( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           v( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           w( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           p( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           up(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           vp(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           wp(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           pp(0:n(1)+1,0:n(2)+1,0:n(3)+1))
  allocate(dudtrko(n(1),n(2),n(3)), &
           dvdtrko(n(1),n(2),n(3)), &
           dwdtrko(n(1),n(2),n(3)))
  allocate(lambdaxyp(n_z(1),n_z(2)))
  allocate(ap(n_z(3)),bp(n_z(3)),cp(n_z(3)))
  allocate(dzc( 0:n(3)+1), &
           dzf( 0:n(3)+1), &
           zc(  0:n(3)+1), &
           zf(  0:n(3)+1), &
           dzci(0:n(3)+1), &
           dzfi(0:n(3)+1))
  allocate(dzc_g( 0:ng(3)+1), &
           dzf_g( 0:ng(3)+1), &
           zc_g(  0:ng(3)+1), &
           zf_g(  0:ng(3)+1), &
           dzci_g(0:ng(3)+1), &
           dzfi_g(0:ng(3)+1))
  allocate(rhsbp%x(n(2),n(3),0:1), &
           rhsbp%y(n(1),n(3),0:1), &
           rhsbp%z(n(1),n(2),0:1))
#if defined(_IMPDIFF)
  allocate(lambdaxyu(n_z(1),n_z(2)), &
           lambdaxyv(n_z(1),n_z(2)), &
           lambdaxyw(n_z(1),n_z(2)))
  allocate(au(n_z(3)),bu(n_z(3)),cu(n_z(3)), &
           av(n_z(3)),bv(n_z(3)),cv(n_z(3)), &
           aw(n_z(3)),bw(n_z(3)),cw(n_z(3)), &
           bb(n_z(3)))
  allocate(rhsbu%x(n(2),n(3),0:1), &
           rhsbu%y(n(1),n(3),0:1), &
           rhsbu%z(n(1),n(2),0:1), &
           rhsbv%x(n(2),n(3),0:1), &
           rhsbv%y(n(1),n(3),0:1), &
           rhsbv%z(n(1),n(2),0:1), &
           rhsbw%x(n(2),n(3),0:1), &
           rhsbw%y(n(1),n(3),0:1), &
           rhsbw%z(n(1),n(2),0:1))
#endif
  if(myid == 0) print*, '*******************************'
  if(myid == 0) print*, '*** Beginning of simulation ***'
  if(myid == 0) print*, '*******************************'
  if(myid == 0) print*, ''
  call initgrid(inivel,ng(3),gr,lz,dzc_g,dzf_g,zc_g,zf_g)
  if(myid == 0) then
    inquire(iolength=rlen) 1._rp
    open(99,file=trim(datadir)//'grid.bin',access='direct',recl=4*ng(3)*rlen)
    write(99,rec=1) dzc_g(1:ng(3)),dzf_g(1:ng(3)),zc_g(1:ng(3)),zf_g(1:ng(3))
    close(99)
    open(99,file=trim(datadir)//'grid.out')
    do kk=0,ng(3)+1
      write(99,'(5E15.7)') 0.,zf_g(kk),zc_g(kk),dzf_g(kk),dzc_g(kk)
    end do
    close(99)
    open(99,file=trim(datadir)//'geometry.out')
      write(99,*) ng(1),ng(2),ng(3)
      write(99,*) l(1),l(2),l(3)
    close(99)
  end if
  do kk=lo(3)-1,hi(3)+1
    k = kk-(lo(3)-1)
    zc( k) = zc_g(kk)
    zf( k) = zf_g(kk)
    dzc(k) = dzc_g(kk)
    dzf(k) = dzf_g(kk)
  end do
  dzci(:) = dzc(:)**(-1)
  dzfi(:) = dzf(:)**(-1)
  dzci_g(:) = dzc_g(:)**(-1)
  dzfi_g(:) = dzf_g(:)**(-1)
  !
  ! test input files before proceeding with the calculation
  !
  call test_sanity(ng,dims,n,n_z,lo,hi,stop_type,cbcvel,cbcpre,bcvel,bcpre,is_forced, &
                   nb,is_bound,dli,dzci_g,dzfi_g,dzci,dzfi)
  !
  if(.not.restart) then
    istep = 0
    time = 0.
    call initflow(inivel,ng,lo,zc/lz,dzc/lz,dzf/lz,visc,u,v,w,p)
    if(myid == 0) print*, '*** Initial condition succesfully set ***'
  else
    call load('r',trim(datadir)//'fld.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,u,v,w,p,time,istep)
    if(myid == 0) print*, '*** Checkpoint loaded at time = ', time, 'time step = ', istep, '. ***'
  end if
  call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,u,v,w)
  call boundp(cbcpre,n,bcpre,nb,is_bound,dl,dzc,p)
  !
  ! post-process and write initial condition
  !
  write(fldnum,'(i7.7)') istep
  include 'out1d.h90'
  include 'out2d.h90'
  include 'out3d.h90'
  !
  dudtrko(:,:,:) = 0.
  dvdtrko(:,:,:) = 0.
  dwdtrko(:,:,:) = 0.
  call chkdt(n,dl,dzci,dzfi,visc,u,v,w,dtmax)
  dt = min(cfl*dtmax,dtmin)
  if(myid == 0) print*, 'dtmax = ', dtmax, 'dt = ',dt
  dti = 1./dt
  kill = .false.
  !
  ! initialize Poisson solver
  !
  call initsolver(ng,zstart,zend,dli,dzci_g,dzfi_g,cbcpre,bcpre(:,:),lambdaxyp,['c','c','c'],ap,bp,cp,arrplanp,normfftp, &
                  rhsbp%x,rhsbp%y,rhsbp%z)
#if defined(_IMPDIFF)
  call initsolver(ng,zstart,zend,dli,dzci_g,dzfi_g,cbcvel(:,:,1),bcvel(:,:,1),lambdaxyu,['f','c','c'],au,bu,cu,arrplanu,normfftu, &
                  rhsbu%x,rhsbu%y,rhsbu%z)
  call initsolver(ng,zstart,zend,dli,dzci_g,dzfi_g,cbcvel(:,:,2),bcvel(:,:,2),lambdaxyv,['c','f','c'],av,bv,cv,arrplanv,normfftv, &
                  rhsbv%x,rhsbv%y,rhsbv%z)
  call initsolver(ng,zstart,zend,dli,dzci_g,dzfi_g,cbcvel(:,:,3),bcvel(:,:,3),lambdaxyw,['c','c','f'],aw,bw,cw,arrplanw,normfftw, &
                  rhsbw%x,rhsbw%y,rhsbw%z)
#if defined(_IMPDIFF_1D)
  deallocate(lambdaxyu,lambdaxyv,lambdaxyw)
  call fftend(arrplanu)
  call fftend(arrplanv)
  call fftend(arrplanw)
#endif
#endif
  !
  ! main loop
  !
  if(myid == 0) print*, '*** Calculation loop starts now ***'
  is_done = .false.
  do while(.not.is_done)
#if defined(_TIMING)
    dt12 = MPI_WTIME()
#endif
    istep = istep + 1
    time = time + dt
    if(myid == 0) print*, 'Timestep #', istep, 'Time = ', time
    dpdl(:)  = 0.
    tauxo(:) = 0.
    tauyo(:) = 0.
    tauzo(:) = 0.
    do irk=1,3
      dtrk = sum(rkcoeff(:,irk))*dt
      dtrki = dtrk**(-1)
      call rk(rkcoeff(:,irk),n,dli,l,dzci,dzfi,visc,dt,u,v,w,p,is_bound,is_forced,velf,bforce, &
              dudtrko,dvdtrko,dwdtrko,tauxo,tauyo,tauzo,up,vp,wp,f)
      if(is_forced(1)) up(1:n(1),1:n(2),1:n(3)) = up(1:n(1),1:n(2),1:n(3)) + f(1)
      if(is_forced(2)) vp(1:n(1),1:n(2),1:n(3)) = vp(1:n(1),1:n(2),1:n(3)) + f(2)
      if(is_forced(3)) wp(1:n(1),1:n(2),1:n(3)) = wp(1:n(1),1:n(2),1:n(3)) + f(3)
#if defined(_IMPDIFF)
      alpha = -1./(.5*visc*dtrk)
      !$OMP WORKSHARE
      up(1:n(1),1:n(2),1:n(3)) = up(1:n(1),1:n(2),1:n(3))*alpha
      !$OMP END WORKSHARE
      bb(:) = bu(:) + alpha
      call updt_rhs_b(['f','c','c'],cbcvel(:,:,1),n,is_bound,rhsbu%x,rhsbu%y,rhsbu%z,up)
#if !defined(_IMPDIFF_1D)
      call solver(n,arrplanu,normfftu,lambdaxyu,au,bb,cu,cbcvel(:,3,1),['f','c','c'],up)
#else
      call solver_gaussel_z(n                  ,au,bb,cu,cbcvel(:,3,1),['f','c','c'],up)
#endif
      !$OMP WORKSHARE
      vp(1:n(1),1:n(2),1:n(3)) = vp(1:n(1),1:n(2),1:n(3))*alpha
      !$OMP END WORKSHARE
      bb(:) = bv(:) + alpha
      call updt_rhs_b(['c','f','c'],cbcvel(:,:,2),n,is_bound,rhsbv%x,rhsbv%y,rhsbv%z,vp)
#if !defined(_IMPDIFF_1D)
      call solver(n,arrplanv,normfftv,lambdaxyv,av,bb,cv,cbcvel(:,3,2),['c','f','c'],vp)
#else
      call solver_gaussel_z(n                  ,av,bb,cv,cbcvel(:,3,2),['c','f','c'],vp)
#endif
      !$OMP WORKSHARE
      wp(1:n(1),1:n(2),1:n(3)) = wp(1:n(1),1:n(2),1:n(3))*alpha
      !$OMP END WORKSHARE
      bb(:) = bw(:) + alpha
      call updt_rhs_b(['c','c','f'],cbcvel(:,:,3),n,is_bound,rhsbw%x,rhsbw%y,rhsbw%z,wp)
#if !defined(_IMPDIFF_1D)
      call solver(n,arrplanw,normfftw,lambdaxyw,aw,bb,cw,cbcvel(:,3,3),['c','c','f'],wp)
#else
      call solver_gaussel_z(n                  ,aw,bb,cw,cbcvel(:,3,3),['c','c','f'],wp)
#endif
#endif
      dpdl(:) = dpdl(:) + f(:)
      call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,up,vp,wp)
#if !defined(_IMPDIFF)
#if defined(_ONE_PRESS_CORR)
      dtrk  = dt
      dtrki = dt**(-1)
      if(irk < 3) then ! pressure correction only at the last RK step
        !$OMP WORKSHARE
        u(:,:,:) = up(:,:,:)
        v(:,:,:) = vp(:,:,:)
        w(:,:,:) = wp(:,:,:)
        !$OMP END WORKSHARE
        cycle
      end if
#endif
#endif
      call fillps(n,dli,dzfi,dtrki,up,vp,wp,pp)
      call updt_rhs_b(['c','c','c'],cbcpre,n,is_bound,rhsbp%x,rhsbp%y,rhsbp%z,pp)
      call solver(n,arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre(:,3),['c','c','c'],pp)
      call boundp(cbcpre,n,bcpre,nb,is_bound,dl,dzc,pp)
      call correc(n,dli,dzci,dtrk,pp,up,vp,wp,u,v,w)
      call bounduvw(cbcvel,n,bcvel,nb,is_bound,.true.,dl,dzc,dzf,u,v,w)
#if defined(_IMPDIFF)
      alphai = alpha**(-1)
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP PRIVATE(i,j,k) &
      !$OMP SHARED(n,p,pp,dxi,dyi,dzfi,dzci,alphai)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            p(i,j,k) = p(i,j,k) + pp(i,j,k) + alphai*( &
#if !defined(_IMPDIFF_1D)
                        (pp(i+1,j,k)-2.*pp(i,j,k)+pp(i-1,j,k))*(dxi**2) + &
                        (pp(i,j+1,k)-2.*pp(i,j,k)+pp(i,j-1,k))*(dyi**2) + &
#endif
                        ((pp(i,j,k+1)-pp(i,j,k  ))*dzci(k  ) - &
                         (pp(i,j,k  )-pp(i,j,k-1))*dzci(k-1))*dzfi(k) )
          end do
        end do
      end do
      !$OMP END PARALLEL DO
#else
      !$OMP WORKSHARE
      p(1:n(1),1:n(2),1:n(3)) = p(1:n(1),1:n(2),1:n(3)) + pp(1:n(1),1:n(2),1:n(3))
      !$OMP END WORKSHARE
#endif
      call boundp(cbcpre,n,bcpre,nb,is_bound,dl,dzc,p)
    end do
    dpdl(:) = -dpdl(:)*dti
    !
    ! check simulation stopping criteria
    !
    if(stop_type(1)) then ! maximum number of time steps reached
      if(istep >= nstep   ) is_done = is_done.or..true.
    end if
    if(stop_type(2)) then ! maximum simulation time reached
      if(time  >= time_max) is_done = is_done.or..true.
    end if
    if(stop_type(3)) then ! maximum wall-clock time reached
      tw = (MPI_WTIME()-twi)/3600.
      if(tw    >= tw_max  ) is_done = is_done.or..true.
    end if
    if(mod(istep,icheck) == 0) then
      if(myid == 0) print*, 'Checking stability and divergence...'
      call chkdt(n,dl,dzci,dzfi,visc,u,v,w,dtmax)
      dt  = min(cfl*dtmax,dtmin)
      if(myid == 0) print*, 'dtmax = ', dtmax, 'dt = ',dt
      if(dtmax < small) then
        if(myid == 0) print*, 'ERROR: timestep is too small.'
        if(myid == 0) print*, 'Aborting...'
        is_done = .true.
        kill = .true.
      end if
      dti = 1./dt
      call chkdiv(lo,hi,dli,dzfi,u,v,w,divtot,divmax)
      if(myid == 0) print*, 'Total divergence = ', divtot, '| Maximum divergence = ', divmax
      if(divmax > small.or.divtot.ne.divtot) then
        if(myid == 0) print*, 'ERROR: maximum divergence is too large.'
        if(myid == 0) print*, 'Aborting...'
        is_done = .true.
        kill = .true.
      end if
    end if
    !
    ! output routines below
    !
    if(mod(istep,iout0d) == 0) then
      !allocate(var(4))
      var(1) = 1.*istep
      var(2) = dt
      var(3) = time
      call out0d(trim(datadir)//'time.out',3,var)
      !
      if(any(is_forced(:)).or.any(abs(bforce(:)) > 0.)) then
        meanvelu = 0.
        meanvelv = 0.
        meanvelw = 0.
        if(is_forced(1).or.abs(bforce(1)) > 0.) then
          call chk_mean(n,dl(1)*dl(2)*dzf/(l(1)*l(2)*l(3)),up,meanvelu)
        end if
        if(is_forced(2).or.abs(bforce(2)) > 0.) then
          call chk_mean(n,dl(1)*dl(2)*dzf/(l(1)*l(2)*l(3)),vp,meanvelv)
        end if
        if(is_forced(3).or.abs(bforce(3)) > 0.) then
          call chk_mean(n,dl(1)*dl(2)*dzf/(l(1)*l(2)*l(3)),wp,meanvelw)
        end if
        if(.not.any(is_forced(:))) dpdl(:) = -bforce(:) ! constant pressure gradient
        var(1)   = time
        var(2:4) = dpdl(1:3)
        var(5:7) = [meanvelu,meanvelv,meanvelw]
        call out0d(trim(datadir)//'forcing.out',7,var)
      end if
      !deallocate(var)
    end if
    write(fldnum,'(i7.7)') istep
    if(mod(istep,iout1d) == 0) then
      include 'out1d.h90'
    end if
    if(mod(istep,iout2d) == 0) then
      include 'out2d.h90'
    end if
    if(mod(istep,iout3d) == 0) then
      include 'out3d.h90'
    end if
    if(mod(istep,isave ) == 0.or.(is_done.and..not.kill)) then
      if(is_overwrite_save) then
        filename = 'fld.bin'
      else
        filename = 'fld_'//fldnum//'.bin'
      end if
      call load('w',trim(datadir)//trim(filename),MPI_COMM_WORLD,ng,[1,1,1],lo,hi,u,v,w,p,time,istep)
      if(.not.is_overwrite_save) then
        !
        ! fld.bin -> last checkpoint file (symbolic link)
        !
        if(myid == 0) call execute_command_line('ln -sf '//trim(filename)//' '//trim(datadir)//'fld.bin')
      end if
      if(myid == 0) print*, '*** Checkpoint saved at time = ', time, 'time step = ', istep, '. ***'
    end if
#if defined(_TIMING)
      dt12 = MPI_WTIME()-dt12
      call MPI_ALLREDUCE(dt12,dt12av ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(dt12,dt12min,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(dt12,dt12max,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
      if(myid == 0) print*, 'Avrg, min & max elapsed time: '
      if(myid == 0) print*, dt12av/(1.*product(dims)),dt12min,dt12max
#endif
  end do
  !
  ! clear ffts
  !
  call fftend(arrplanp)
#if defined(_IMPDIFF) && !defined(_IMPDIFF_1D)
  call fftend(arrplanu)
  call fftend(arrplanv)
  call fftend(arrplanw)
#endif
  !
  ! deallocate variables
  !
  deallocate(u,v,w,p,up,vp,wp,pp)
  deallocate(dudtrko,dvdtrko,dwdtrko)
  deallocate(lambdaxyp)
  deallocate(ap,bp,cp)
  deallocate(dzc,dzf,zc,zf,dzci,dzfi)
  deallocate(rhsbp%x,rhsbp%y,rhsbp%z)
#if defined(_IMPDIFF)
#if !defined(_IMPDIFF_1D)
  deallocate(lambdaxyu,lambdaxyv,lambdaxyw)
#endif
  deallocate(au,bu,cu,av,bv,cv,aw,bw,cw,bb)
  deallocate(rhsbu%x,rhsbu%y,rhsbu%z, &
             rhsbv%x,rhsbv%y,rhsbv%z, &
             rhsbw%x,rhsbw%y,rhsbw%z)
#endif
  if(myid == 0.and.(.not.kill)) print*, '*** Fim ***'
  call decomp_2d_finalize
  call MPI_FINALIZE(ierr)
  stop
end program cans
