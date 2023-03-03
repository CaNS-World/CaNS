! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
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
!-------------------------------------------------------------------------------------
program cans
#if defined(_DEBUG)
  use, intrinsic :: iso_fortran_env, only: compiler_version,compiler_options
#endif
  use, intrinsic :: iso_c_binding  , only: C_PTR
  use, intrinsic :: ieee_arithmetic, only: is_nan => ieee_is_nan
  use mpi
  use decomp_2d
  use mod_bound          , only: boundp,bounduvw,updt_rhs_b
  use mod_chkdiv         , only: chkdiv
  use mod_chkdt          , only: chkdt
  use mod_common_mpi     , only: myid,ierr
  use mod_correc         , only: correc
  use mod_fft            , only: fftini,fftend
  use mod_fillps         , only: fillps
  use mod_initflow       , only: initflow
  use mod_initgrid       , only: initgrid
  use mod_initmpi        , only: initmpi
  use mod_initsolver     , only: initsolver
  use mod_load           , only: load
  use mod_mom            , only: bulk_forcing
  use mod_rk             , only: rk,rk_cmpt_bulk_forcing
  use mod_output         , only: out0d,gen_alias,out1d,out1d_chan,out2d,out3d,write_log_output,write_visu_2d,write_visu_3d
  use mod_param          , only: lz,visc,small, &
                                 nb,is_bound,cbcvel,bcvel,cbcpre,bcpre, &
                                 icheck,iout0d,iout1d,iout2d,iout3d,isave, &
                                 nstep,time_max,tw_max,stop_type,restart,is_overwrite_save,nsaves_max, &
                                 rkcoeff,   &
                                 datadir,   &
                                 cfl,dtmin, &
                                 inivel,    &
                                 dims, &
                                 gr, &
                                 is_forced,velf,bforce, &
                                 ng,l,dl,dli, &
                                 read_input
  use mod_sanity         , only: test_sanity_input,test_sanity_solver
#if !defined(_OPENACC)
  use mod_solver         , only: solver
#if defined(_IMPDIFF_1D)
  use mod_solver         , only: solver_gaussel_z
#endif
#else
  use mod_solver_gpu     , only: solver => solver_gpu
#if defined(_IMPDIFF_1D)
  use mod_solver_gpu     , only: solver_gaussel_z => solver_gaussel_z_gpu
#endif
  use mod_workspaces     , only: init_wspace_arrays,set_cufft_wspace
  use mod_common_cudecomp, only: istream_acc_queue_1
#endif
  use mod_timer          , only: timer_tic,timer_toc,timer_print
  use mod_updatep        , only: updatep
  use mod_utils          , only: bulk_mean
  !@acc use mod_utils    , only: device_memory_footprint
  use mod_types
  use omp_lib
  implicit none
  integer , dimension(3) :: lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z
  real(rp), allocatable, dimension(:,:,:) :: u,v,w,p,pp
  real(rp), dimension(3) :: tauxo,tauyo,tauzo
  real(rp), dimension(3) :: f
#if !defined(_OPENACC)
  type(C_PTR), dimension(2,2) :: arrplanp
#else
  integer    , dimension(2,2) :: arrplanp
#endif
  real(rp), allocatable, dimension(:,:) :: lambdaxyp
  real(rp), allocatable, dimension(:) :: ap,bp,cp
  real(rp) :: normfftp
  type rhs_bound
    real(rp), allocatable, dimension(:,:,:) :: x
    real(rp), allocatable, dimension(:,:,:) :: y
    real(rp), allocatable, dimension(:,:,:) :: z
  end type rhs_bound
  type(rhs_bound) :: rhsbp
  real(rp) :: alpha
#if defined(_IMPDIFF)
#if !defined(_OPENACC)
  type(C_PTR), dimension(2,2) :: arrplanu,arrplanv,arrplanw
#else
  integer    , dimension(2,2) :: arrplanu,arrplanv,arrplanw
#endif
  real(rp), allocatable, dimension(:,:) :: lambdaxyu,lambdaxyv,lambdaxyw,lambdaxy
  real(rp), allocatable, dimension(:) :: au,av,aw,bu,bv,bw,cu,cv,cw,aa,bb,cc
  real(rp) :: normfftu,normfftv,normfftw
  type(rhs_bound) :: rhsbu,rhsbv,rhsbw
  real(rp), allocatable, dimension(:,:,:) :: rhsbx,rhsby,rhsbz
#endif
  real(rp) :: dt,dti,dtmax,time,dtrk,dtrki,divtot,divmax
  integer :: irk,istep
  real(rp), allocatable, dimension(:) :: dzc  ,dzf  ,zc  ,zf  ,dzci  ,dzfi, &
                                         dzc_g,dzf_g,zc_g,zf_g,dzci_g,dzfi_g, &
                                         grid_vol_ratio_c,grid_vol_ratio_f
  real(rp) :: meanvelu,meanvelv,meanvelw
  real(rp), dimension(3) :: dpdl
  !real(rp), allocatable, dimension(:) :: var
  real(rp), dimension(42) :: var
#if defined(_TIMING)
  real(rp) :: dt12,dt12av,dt12min,dt12max
#endif
  real(rp) :: twi,tw
  integer  :: savecounter
  character(len=7  ) :: fldnum
  character(len=4  ) :: chkptnum
  character(len=100) :: filename
  integer :: k,kk
  logical :: is_done,kill
  integer :: rlen
  !
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  !
  ! read parameter file
  !
  call read_input(myid)
  !
  ! initialize MPI/OpenMP
  !
  !$ call omp_set_num_threads(omp_get_max_threads())
  call initmpi(ng,dims,cbcpre,lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z,nb,is_bound)
  twi = MPI_WTIME()
  savecounter = 0
  !
  ! allocate variables
  !
  allocate(u( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           v( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           w( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           p( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           pp(0:n(1)+1,0:n(2)+1,0:n(3)+1))
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
  allocate(grid_vol_ratio_c,mold=dzc)
  allocate(grid_vol_ratio_f,mold=dzf)
  allocate(rhsbp%x(n(2),n(3),0:1), &
           rhsbp%y(n(1),n(3),0:1), &
           rhsbp%z(n(1),n(2),0:1))
#if defined(_IMPDIFF)
  allocate(lambdaxyu(n_z(1),n_z(2)), &
           lambdaxyv(n_z(1),n_z(2)), &
           lambdaxyw(n_z(1),n_z(2)), &
           lambdaxy( n_z(1),n_z(2)))
  allocate(au(n_z(3)),bu(n_z(3)),cu(n_z(3)), &
           av(n_z(3)),bv(n_z(3)),cv(n_z(3)), &
           aw(n_z(3)),bw(n_z(3)),cw(n_z(3)), &
           aa(n_z(3)),bb(n_z(3)),cc(n_z(3)))
  allocate(rhsbu%x(n(2),n(3),0:1), &
           rhsbu%y(n(1),n(3),0:1), &
           rhsbu%z(n(1),n(2),0:1), &
           rhsbv%x(n(2),n(3),0:1), &
           rhsbv%y(n(1),n(3),0:1), &
           rhsbv%z(n(1),n(2),0:1), &
           rhsbw%x(n(2),n(3),0:1), &
           rhsbw%y(n(1),n(3),0:1), &
           rhsbw%z(n(1),n(2),0:1), &
           rhsbx(  n(2),n(3),0:1), &
           rhsby(  n(1),n(3),0:1), &
           rhsbz(  n(1),n(2),0:1))
#endif
#if defined(_DEBUG)
  if(myid == 0) print*, 'This executable of CaNS was built with compiler: ', compiler_version()
  if(myid == 0) print*, 'Using the options: ', compiler_options()
  block
    character(len=MPI_MAX_LIBRARY_VERSION_STRING) :: mpi_version
    integer :: ilen
    call MPI_GET_LIBRARY_VERSION(mpi_version,ilen,ierr)
    if(myid == 0) print*, 'MPI Version: ', trim(mpi_version)
  end block
  if(myid == 0) print*, ''
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
      write(99,'(5E16.7e3)') 0.,zf_g(kk),zc_g(kk),dzf_g(kk),dzc_g(kk)
    end do
    close(99)
    open(99,file=trim(datadir)//'geometry.out')
      write(99,*) ng(1),ng(2),ng(3)
      write(99,*) l(1),l(2),l(3)
    close(99)
  end if
  !$acc enter data copyin(lo,hi,n) async
  !$acc enter data copyin(bforce,dl,dli,l) async
  !$acc enter data copyin(zc_g,zf_g,dzc_g,dzf_g) async
  !$acc enter data create(zc,zf,dzc,dzf,dzci,dzfi,dzci_g,dzfi_g) async
  !
  !$acc parallel loop default(present) private(k) async
  do kk=lo(3)-1,hi(3)+1
    k = kk-(lo(3)-1)
    zc( k) = zc_g(kk)
    zf( k) = zf_g(kk)
    dzc(k) = dzc_g(kk)
    dzf(k) = dzf_g(kk)
    dzci(k) = dzc(k)**(-1)
    dzfi(k) = dzf(k)**(-1)
  end do
  !$acc kernels default(present) async
  dzci_g(:) = dzc_g(:)**(-1)
  dzfi_g(:) = dzf_g(:)**(-1)
  !$acc end kernels
  !$acc enter data create(grid_vol_ratio_c,grid_vol_ratio_f) async
  !$acc kernels default(present) async
  grid_vol_ratio_c(:) = dl(1)*dl(2)*dzc(:)/(l(1)*l(2)*l(3))
  grid_vol_ratio_f(:) = dl(1)*dl(2)*dzf(:)/(l(1)*l(2)*l(3))
  !$acc end kernels
  !$acc update self(dzci,dzfi) async
  !$acc exit data copyout(zc_g,zf_g,dzc_g,dzf_g,dzci_g,dzfi_g) async ! not needed on the device
  !$acc wait
  !
  ! test input files before proceeding with the calculation
  !
  call test_sanity_input(ng,dims,stop_type,cbcvel,cbcpre,bcvel,bcpre,is_forced)
  !
  ! initialize Poisson solver
  !
  call initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcpre,bcpre(:,:), &
                  lambdaxyp,['c','c','c'],ap,bp,cp,arrplanp,normfftp,rhsbp%x,rhsbp%y,rhsbp%z)
  !$acc enter data copyin(lambdaxyp,ap,bp,cp) async
  !$acc enter data copyin(rhsbp,rhsbp%x,rhsbp%y,rhsbp%z) async
  !$acc wait
#if defined(_IMPDIFF)
  call initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcvel(:,:,1),bcvel(:,:,1), &
                  lambdaxyu,['f','c','c'],au,bu,cu,arrplanu,normfftu,rhsbu%x,rhsbu%y,rhsbu%z)
  call initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcvel(:,:,2),bcvel(:,:,2), &
                  lambdaxyv,['c','f','c'],av,bv,cv,arrplanv,normfftv,rhsbv%x,rhsbv%y,rhsbv%z)
  call initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcvel(:,:,3),bcvel(:,:,3), &
                  lambdaxyw,['c','c','f'],aw,bw,cw,arrplanw,normfftw,rhsbw%x,rhsbw%y,rhsbw%z)
#if defined(_IMPDIFF_1D)
  deallocate(lambdaxyu,lambdaxyv,lambdaxyw,lambdaxy)
  call fftend(arrplanu)
  call fftend(arrplanv)
  call fftend(arrplanw)
  deallocate(rhsbu%x,rhsbu%y,rhsbv%x,rhsbv%y,rhsbw%x,rhsbw%y,rhsbx,rhsby)
#endif
  !$acc enter data copyin(lambdaxyu,au,bu,cu,lambdaxyv,av,bv,cv,lambdaxyw,aw,bw,cw) async
  !$acc enter data copyin(rhsbu,rhsbu%x,rhsbu%y,rhsbu%z) async
  !$acc enter data copyin(rhsbv,rhsbv%x,rhsbv%y,rhsbv%z) async
  !$acc enter data copyin(rhsbw,rhsbw%x,rhsbw%y,rhsbw%z) async
  !$acc enter data create(lambdaxy,aa,bb,cc) async
  !$acc enter data create(rhsbx,rhsby,rhsbz) async
  !$acc wait
#endif
#if defined(_OPENACC)
  !
  ! determine workspace sizes and allocate the memory
  !
  call init_wspace_arrays()
  call set_cufft_wspace(pack(arrplanp,.true.),istream_acc_queue_1)
#if defined(_IMPDIFF) && !defined(_IMPDIFF_1D)
  call set_cufft_wspace(pack(arrplanu,.true.),istream_acc_queue_1)
  call set_cufft_wspace(pack(arrplanv,.true.),istream_acc_queue_1)
  call set_cufft_wspace(pack(arrplanw,.true.),istream_acc_queue_1)
#endif
  if(myid == 0) print*,'*** Device memory footprint (Gb): ', &
                  device_memory_footprint(n,n_z)/(1._sp*1024**3), ' ***'
#endif
#if defined(_DEBUG_SOLVER)
  call test_sanity_solver(ng,lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z,dli,dzc,dzf,dzci,dzfi,dzci_g,dzfi_g, &
                          nb,is_bound,cbcvel,cbcpre,bcvel,bcpre)
#endif
  !
  if(.not.restart) then
    istep = 0
    time = 0.
    !$acc update self(zc,dzc,dzf)
    call initflow(inivel,ng,lo,zc,dzc,dzf,visc,u,v,w,p)
    if(myid == 0) print*, '*** Initial condition succesfully set ***'
  else
    call load('r',trim(datadir)//'fld.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,u,v,w,p,time,istep)
    if(myid == 0) print*, '*** Checkpoint loaded at time = ', time, 'time step = ', istep, '. ***'
  end if
  !$acc enter data copyin(u,v,w,p) create(pp)
  call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,u,v,w)
  call boundp(cbcpre,n,bcpre,nb,is_bound,dl,dzc,p)
  !
  ! post-process and write initial condition
  !
  write(fldnum,'(i7.7)') istep
  !$acc update self(u,v,w,p)
  include 'out1d.h90'
  include 'out2d.h90'
  include 'out3d.h90'
  !
  call chkdt(n,dl,dzci,dzfi,visc,u,v,w,dtmax)
  dt = min(cfl*dtmax,dtmin)
  if(myid == 0) print*, 'dtmax = ', dtmax, 'dt = ',dt
  dti = 1./dt
  kill = .false.
  !
  ! main loop
  !
  if(myid == 0) print*, '*** Calculation loop starts now ***'
  is_done = .false.
  do while(.not.is_done)
#if defined(_TIMING)
    !$acc wait(1)
    dt12 = MPI_WTIME()
#endif
    istep = istep + 1
    time = time + dt
    if(myid == 0) print*, 'Time step #', istep, 'Time = ', time
    tauxo(:) = 0.
    tauyo(:) = 0.
    tauzo(:) = 0.
    dpdl(:)  = 0.
    do irk=1,3
      dtrk = sum(rkcoeff(:,irk))*dt
      dtrki = dtrk**(-1)
      call rk_cmpt_bulk_forcing(rkcoeff(:,irk),n,dli,l,dzci,dzfi,visc,dt,is_bound,is_forced,u,v,w,tauxo,tauyo,tauzo,f,.true.)
      call rk(rkcoeff(:,irk),n,dli,dzci,dzfi,visc,dt,p,bforce,u,v,w)
#if defined(_IMPDIFF)
      alpha = -.5*visc*dtrk
      !$OMP PARALLEL WORKSHARE
      !$acc kernels present(rhsbx,rhsby,rhsbz,rhsbu) async(1)
#if !defined(_IMPDIFF_1D)
      rhsbx(:,:,0:1) = rhsbu%x(:,:,0:1)*alpha
      rhsby(:,:,0:1) = rhsbu%y(:,:,0:1)*alpha
#endif
      rhsbz(:,:,0:1) = rhsbu%z(:,:,0:1)*alpha
      !$acc end kernels
      !$OMP END PARALLEL WORKSHARE
      call updt_rhs_b(['f','c','c'],cbcvel(:,:,1),n,is_bound,rhsbx,rhsby,rhsbz,u)
      !$acc kernels default(present) async(1)
      !$OMP PARALLEL WORKSHARE
      aa(:) = au(:)*alpha
      bb(:) = bu(:)*alpha + 1.
      cc(:) = cu(:)*alpha
#if !defined(_IMPDIFF_1D)
      lambdaxy(:,:) = lambdaxyu(:,:)*alpha
#endif
      !$OMP END PARALLEL WORKSHARE
      !$acc end kernels
#if !defined(_IMPDIFF_1D)
      call solver(n,ng,arrplanu,normfftu,lambdaxy,aa,bb,cc,cbcvel(:,:,1),['f','c','c'],u)
#else
      call solver_gaussel_z(n                    ,aa,bb,cc,cbcvel(:,3,1),['f','c','c'],u)
#endif
      !$OMP PARALLEL WORKSHARE
      !$acc kernels present(rhsbx,rhsby,rhsbz,rhsbv) async(1)
#if !defined(_IMPDIFF_1D)
      rhsbx(:,:,0:1) = rhsbv%x(:,:,0:1)*alpha
      rhsby(:,:,0:1) = rhsbv%y(:,:,0:1)*alpha
#endif
      rhsbz(:,:,0:1) = rhsbv%z(:,:,0:1)*alpha
      !$acc end kernels
      !$OMP END PARALLEL WORKSHARE
      call updt_rhs_b(['c','f','c'],cbcvel(:,:,2),n,is_bound,rhsbx,rhsby,rhsbz,v)
      !$acc kernels default(present) async(1)
      !$OMP PARALLEL WORKSHARE
      aa(:) = av(:)*alpha
      bb(:) = bv(:)*alpha + 1.
      cc(:) = cv(:)*alpha
#if !defined(_IMPDIFF_1D)
      lambdaxy(:,:) = lambdaxyv(:,:)*alpha
#endif
      !$OMP END PARALLEL WORKSHARE
      !$acc end kernels
#if !defined(_IMPDIFF_1D)
      call solver(n,ng,arrplanv,normfftv,lambdaxy,aa,bb,cc,cbcvel(:,:,2),['c','f','c'],v)
#else
      call solver_gaussel_z(n                    ,aa,bb,cc,cbcvel(:,3,2),['c','f','c'],v)
#endif
      !$OMP PARALLEL WORKSHARE
      !$acc kernels present(rhsbx,rhsby,rhsbz,rhsbw) async(1)
#if !defined(_IMPDIFF_1D)
      rhsbx(:,:,0:1) = rhsbw%x(:,:,0:1)*alpha
      rhsby(:,:,0:1) = rhsbw%y(:,:,0:1)*alpha
#endif
      rhsbz(:,:,0:1) = rhsbw%z(:,:,0:1)*alpha
      !$acc end kernels
      !$OMP END PARALLEL WORKSHARE
      call updt_rhs_b(['c','c','f'],cbcvel(:,:,3),n,is_bound,rhsbx,rhsby,rhsbz,w)
      !$acc kernels default(present) async(1)
      !$OMP PARALLEL WORKSHARE
      aa(:) = aw(:)*alpha
      bb(:) = bw(:)*alpha + 1.
      cc(:) = cw(:)*alpha
#if !defined(_IMPDIFF_1D)
      lambdaxy(:,:) = lambdaxyw(:,:)*alpha
#endif
      !$OMP END PARALLEL WORKSHARE
      !$acc end kernels
#if !defined(_IMPDIFF_1D)
      call solver(n,ng,arrplanw,normfftw,lambdaxy,aa,bb,cc,cbcvel(:,:,3),['c','c','f'],w)
#else
      call solver_gaussel_z(n                    ,aa,bb,cc,cbcvel(:,3,3),['c','c','f'],w)
#endif
      call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,u,v,w)
      call rk_cmpt_bulk_forcing(rkcoeff(:,irk),n,dli,l,dzci,dzfi,visc,dt,is_bound,is_forced,u,v,w,tauxo,tauyo,tauzo,f,.false.)
#endif
      call bulk_forcing(n,is_forced,f,u,v,w)
      dpdl(:) = dpdl(:) + f(:)
      call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,u,v,w)
      call fillps(n,dli,dzfi,dtrki,u,v,w,pp)
      call updt_rhs_b(['c','c','c'],cbcpre,n,is_bound,rhsbp%x,rhsbp%y,rhsbp%z,pp)
      call solver(n,ng,arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre,['c','c','c'],pp)
      call boundp(cbcpre,n,bcpre,nb,is_bound,dl,dzc,pp)
      call correc(n,dli,dzci,dtrk,pp,u,v,w)
      call bounduvw(cbcvel,n,bcvel,nb,is_bound,.true.,dl,dzc,dzf,u,v,w)
      call updatep(n,dli,dzci,dzfi,alpha,pp,p)
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
        if(myid == 0) print*, 'ERROR: time step is too small.'
        if(myid == 0) print*, 'Aborting...'
        is_done = .true.
        kill = .true.
      end if
      dti = 1./dt
      call chkdiv(lo,hi,dli,dzfi,u,v,w,divtot,divmax)
      if(myid == 0) print*, 'Total divergence = ', divtot, '| Maximum divergence = ', divmax
#if !defined(_MASK_DIVERGENCE_CHECK)
      if(divmax > small.or.is_nan(divtot)) then
        if(myid == 0) print*, 'ERROR: maximum divergence is too large.'
        if(myid == 0) print*, 'Aborting...'
        is_done = .true.
        kill = .true.
      end if
#endif
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
          call bulk_mean(n,grid_vol_ratio_f,u,meanvelu)
        end if
        if(is_forced(2).or.abs(bforce(2)) > 0.) then
          call bulk_mean(n,grid_vol_ratio_f,v,meanvelv)
        end if
        if(is_forced(3).or.abs(bforce(3)) > 0.) then
          call bulk_mean(n,grid_vol_ratio_c,w,meanvelw)
        end if
        if(.not.any(is_forced(:))) dpdl(:) = -bforce(:) ! constant pressure gradient
        var(1)   = time
        var(2:4) = dpdl(1:3)
        var(5:7) = [meanvelu,meanvelv,meanvelw]
        call out0d(trim(datadir)//'forcing.out',7,var)
      end if
    end if
    write(fldnum,'(i7.7)') istep
    if(mod(istep,iout1d) == 0) then
      !$acc update self(u,v,w,p)
      include 'out1d.h90'
    end if
    if(mod(istep,iout2d) == 0) then
      !$acc update self(u,v,w,p)
      include 'out2d.h90'
    end if
    if(mod(istep,iout3d) == 0) then
      !$acc update self(u,v,w,p)
      include 'out3d.h90'
    end if
    if(mod(istep,isave ) == 0.or.(is_done.and..not.kill)) then
      if(is_overwrite_save) then
        filename = 'fld.bin'
      else
        filename = 'fld_'//fldnum//'.bin'
        if(nsaves_max > 0) then
          if(savecounter >= nsaves_max) savecounter = 0
          savecounter = savecounter + 1
          write(chkptnum,'(i4.4)') savecounter
          filename = 'fld_'//chkptnum//'.bin'
          var(1) = 1.*istep
          var(2) = time
          var(3) = 1.*savecounter
          call out0d(trim(datadir)//'log_checkpoints.out',3,var)
        end if
      end if
      !$acc update self(u,v,w,p)
      call load('w',trim(datadir)//trim(filename),MPI_COMM_WORLD,ng,[1,1,1],lo,hi,u,v,w,p,time,istep)
      if(.not.is_overwrite_save) then
        !
        ! fld.bin -> last checkpoint file (symbolic link)
        !
        call gen_alias(myid,trim(datadir),trim(filename),'fld.bin')
      end if
      if(myid == 0) print*, '*** Checkpoint saved at time = ', time, 'time step = ', istep, '. ***'
    end if
#if defined(_TIMING)
      !$acc wait(1)
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
  if(myid == 0.and.(.not.kill)) print*, '*** Fim ***'
  call decomp_2d_finalize
  call MPI_FINALIZE(ierr)
end program cans
